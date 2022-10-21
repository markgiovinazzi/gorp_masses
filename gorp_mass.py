# import libraries
import ssl, numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
from scipy import interpolate
from scipy.special import erf
from astroquery.gaia import Gaia

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context

# small function to simplify magnitude conversions
def absmag(appmag, plx):

    return appmag + 5 - 5 * np.log10(1000. / plx)

# function useed to estimate masses from Gaia
def gaia_posterior(ids, N = 100000):

    '''
    Estimate a posterior on mass given GAIA Source ID(s)
    Relation derived in Giovinazzi, M. R. & Blake, C. H. (2022)
    
    If you use this code, please cite the relevant paper:
    https://ui.adsabs.harvard.edu/abs/2022AJ....164..164G
    
    PARAMETERS
        ids: str, array
            GAIA Source ID(s)
        N: int
            number of random samples used in estimating uncertainty in mass; default = 100,000
            
    RETURNS
        source_id: str, array
            input source ids that are within the eligible photometric range
        mass: float, array
            point estimate on the mass(es) corresponding to the provided GAIA Source IDs
        mass_sig: float, array
            standard error on the mass(es)
        is_MS: bool, array
            flag determining whether a star is main sequence
            True if the star is main sequence, False otherwise
            this relation is NOT applicable to stars that return `is_MS = False`
        is_ss: bool, array
            flag relaying whether Gaia registers the star as a good fit to a single-star astrometric solution
            True if the star is likley a single star with a good astrometric solution from Gaia, False otherwise
            if False, user should be skeptical of returned mass estimate
        phot_par: bool, array
            flag for photometric and astrometric parameters that may bias our mass estimates
            True if apparent `g_mag > 5`, `rp_flux_over_err > 10`, `rel_plx_err < 0.2`, `plx_err < 2` mas, `plx > 1` mas, False otherwise
            if False, photometry may be poor or saturated, or absolute magnitude (due to plx) is uncertain; use mass estimate with caution
        ext_flag: bool, array
            flag determining whether or not the star may be subjected to significant extinction
            True if galactic latitude `|b| > 10` degrees, False otherwise
            if False, be mindful that extinction may make star appear fainter, causing potential underestimate in mass
    '''

    # reformat and determine th enumber of ids provided based on the input type
    if isinstance(ids, list): N_ids = len(ids); ids = str(ids)[1:-1]
    elif isinstance(ids, str): N_ids = 1
    else: print('Input "ids" must be either \'str\' or list/tuple of \'str\''); exit(1)

    # set up our query to the Gaia DR3 catalog
    query = """
    SELECT g.source_id, g.phot_rp_mean_mag, g.parallax, g.parallax_error,
           g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.astrometric_excess_noise,
           g.b, g.phot_g_mean_mag, g.bp_rp, g.phot_rp_mean_flux_over_error
    FROM gaiadr3.gaia_source as g
    WHERE g.source_id IN (""" + ids + """)
    """

    # excecute the query and store the results
    print('Querying the GAIA Archive...')
    job = Gaia.launch_job(query)
    r = job.get_results()
    
    # first, check which stars are not eligible for our relation
    source_ids  = np.array(r['source_id'])
    rp_mags     = np.array(r['phot_rp_mean_mag'])
    plxs        = np.array(r['parallax'])
    abs_rp_mags = absmag(rp_mags, plxs)
    ineligible_mask = (abs_rp_mags < 4) | (abs_rp_mags > 14.5)

    # the program cannot estimate masses for stars that are outside the recommended range
    if np.sum(ineligible_mask) == N_ids:
        print('Sorry, all ids provided are outside the eligible range (4.0 < MGRP < 14.5).'); exit(1)
    elif np.sum(ineligible_mask) > 0:
        print('The following ' + str(int(np.sum(ineligible_mask))) + ' ids are outside the eligible range (4.0 < MGRP < 14.5).')
        for i, abs_rp in zip(source_ids[ineligible_mask], abs_rp_mags[ineligible_mask]):
            print(str(i) + ' (abs_rp_mag = ' + str(round(abs_rp, 2)) + ')')
        N_ids -= np.sum(ineligible_mask)
        # adjust queried results by remeoving ineligble stars
        r = r[~ineligible_mask]
    else:
        print('All ids provided are within the eligible photometric range (4.0 < MGRP < 14.5).')

    # re-itemize downloaded parameters for gaia source id(s)
    source_ids        = np.array(r['source_id'])
    rp_mags           = np.array(r['phot_rp_mean_mag'])
    plxs              = np.array(r['parallax'])
    plx_errs          = np.array(r['parallax_error'])
    ruwes             = np.array(r['ruwe'])
    frac_multi_peaks  = np.array(r['ipd_frac_multi_peak'])
    harmonic_amps     = np.array(r['ipd_gof_harmonic_amplitude'])
    astr_noise        = np.array(r['astrometric_excess_noise'])
    bs                = np.array(r['b'])
    g_mags            = np.array(r['phot_g_mean_mag'])
    bp_rps            = np.array(r['bp_rp'])
    rp_flux_over_errs = np.array(r['phot_rp_mean_flux_over_error'])

    # create additional quantities based on those downloaded above
    abs_rp_mags  = absmag(rp_mags, plxs)
    abs_g_mags   = absmag(g_mags, plxs)
    rel_plx_errs = plx_errs / plxs
    
    # tests for main sequence -- source id(s) that do not pass cannot have their masses photometrically estimated
    is_MS = (abs_g_mags > 3) & (abs_g_mags < 3.1 * bp_rps + 5) & (abs_g_mags > 4)
    
    # stars in Gaia that do not meet these criteria are not good fits to single-star astrometric solutions
    is_ss = (ruwes < 1.4) & (frac_multi_peaks <= 10) & (harmonic_amps <= 0.1) & (astr_noise <= 1)
    
    # additional photometric and parallactic cuts; mass estimates for stars not passing these should be used cautiously
    phot_par = (g_mags > 5) & (rp_flux_over_errs > 10) & (rel_plx_errs < 0.2) & (plx_errs < 2) & (plxs > 1)
    
    # stars in the galactic plane may be subjected to extinction
    ext_flag = np.abs(bs) > 10

    # load in uncertainties from file
    with fits.open('resources/massmag_errs_new.fits') as hdul:

        points = hdul[1].data['X'][0]
        massmag_errs = hdul[1].data['ERF_DATA'][0]
    
    massmag_errs = np.std(massmag_errs, axis = 1)
    f = interpolate.interp1d(points, massmag_errs)

    # define the relation from the paper
    A, B, C = 0.444687, -0.0968913, 0.0746491
    log10masses = A + B * abs_rp_mags - C * (1. + erf(abs_rp_mags - 9.5))
    masses = 10**log10masses
    try: mass_errs = np.array([np.std(10**np.random.normal(log10masses[i], f(abs_rp_mags[i]), N)) for i in range(N_ids)])
    except: mass_errs = np.std(10**np.random.normal(log10masses, f(abs_rp_mags), N))
    
    # columnize the results and write to an ascii table
    t = Table()
    t['source_id'] = source_ids
    t['mass'] = masses
    t['mass_err'] = mass_errs
    t['is_MS'] = is_MS
    t['is_ss'] = is_ss
    t['phot_par'] = phot_par
    t['ext_flag'] = ext_flag
    
    print('Writing to file...')
    ascii.write(t, 'gorp_id_results.dat', overwrite = True)
    print('Done!')
    
# function useed to estimate masses from RP magnitudes
def rp_posterior(abs_rp_mags, N = 100000):

    '''
    Estimate a posterior on mass given absolute RP magnitudes
    Relation derived in Giovinazzi, M. R. & Blake, C. H. (2022)
    
    If you use this code, please cite the relevant paper:
    https://ui.adsabs.harvard.edu/abs/2022AJ....164..164G
    
    PARAMETERS
        abs_rp_mags: float, array
            GAIA Source ID(s)
        N: int
            number of random samples used in estimating uncertainty in mass; default = 100,000
            
    RETURNS
        mass: float, array
            point estimate on the mass(es) corresponding to the provided GAIA Source IDs
        mass_sig: float, array
            standard error on the mass(es)
    '''

    if isinstance(abs_rp_mags, float): abs_rp_mags = np.array(abs_rp_mags)
    ineligible_mask = (abs_rp_mags < 4) | (abs_rp_mags > 14.5)

    # the program cannot estimate masses for stars that are outside the recommended range
    if np.sum(ineligible_mask) == len(abs_rp_mags):
        print('Sorry, all RP mags provided are outside the eligible range (4.0 < MGRP < 14.5). Make sure you are using absolute mags!'); exit(1)
    elif np.sum(ineligible_mask) > 0:
        print('Some RP mags are outside the eligible range (4.0 < MGRP < 14.5). Ignoring those.')
        # adjust queried results by remeoving ineligble stars
        abs_rp_mags = abs_rp_mags[~ineligible_mask]
    else:
        print('All RP mags provided are within the eligible photometric range (4.0 < MGRP < 14.5).')

    # load in uncertainties from file
    with fits.open('resources/massmag_errs_new.fits') as hdul:

        points = hdul[1].data['X'][0]
        massmag_errs = hdul[1].data['ERF_DATA'][0]
    
    massmag_errs = np.std(massmag_errs, axis = 1)
    f = interpolate.interp1d(points, massmag_errs)

    # define the relation from the paper
    A, B, C = 0.444687, -0.0968913, 0.0746491
    log10masses = A + B * abs_rp_mags - C * (1. + erf(abs_rp_mags - 9.5))
    masses = 10**log10masses
    try: mass_errs = np.array([np.std(10**np.random.normal(log10masses[i], f(abs_rp_mags[i]), N)) for i in range(N_ids)])
    except: mass_errs = np.std(10**np.random.normal(log10masses, f(abs_rp_mags), N))
    
    # columnize the results and write to an ascii table
    t = Table()
    t['abs_rp_mag'] = abs_rp_mags
    t['mass'] = masses
    t['mass_err'] = mass_errs
    
    print('Writing to file...')
    ascii.write(t, 'gorp_rp_results.dat', overwrite = True)
    print('Done!')
