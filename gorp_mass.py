# import libraries
import os, ssl, pkg_resources, numpy as np, matplotlib.pyplot as plt
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

# function used to estimate masses from Gaia
def gaia_posterior(ids, N = 100000, plot_1d = False, plot_2d = False, plot_path = 'gorp_plots', plot_ext = 'pdf', lit_vals = None, lit_errs = None, lit_labs = None,
    colors = ['darkorange', 'darkgreen', 'darkred', 'purple', 'black']):

    '''
    Estimate a posterior on mass given GAIA Source ID(s)
    Relation derived in Giovinazzi, M. R. & Blake, C. H. (2022)
    
    If you use this code, please cite the relevant paper:
    https://ui.adsabs.harvard.edu/abs/2022AJ....164..164G
    
    PARAMETERS
        ids: str, array
            GAIA Source ID(s)
            if `ids` format is array, array values must be strings
        N: int
            number of random samples used in estimating uncertainty in mass; default = 100,000
        plot_1d: bool
            if True, plot of 1d distribution of masses will be saved
        plot_2d: bool
            if True, plot of 2d distribution of masses and absolute magnitudes will be saved
        plot_path: str
            path to directory where plots will be saved
            if path does not already exist, a new directory `plot_path` will be created
            by default, the program will create a new directory called 'gorp_plots'
        plot_ext: str
            file extension for plots
            default is pdf to be used in publications
        lit_vals: float, array
            optional entry to overlay vertical lines over 1d mass posterior
            currently only used for 1d plotting
            if `lit_vals` format is array, array values must be floats
        lit_errs: float, array
            optional entry to overlay vertical 1-sigma error regions over 1d mass posterior
            must provide input for `lit_vals`
            currently only used for 1d plotting
            if `lit_errs` format is array, array values must be floats
        lit_labs: str, array
            optional entry to add legend to plot with text labels corresponding to each literature value
            must provide input for `lit_vals`
            currently only used for 1d plotting
            if `lit_labs` format is array, array values must be strings
        colors: str, array
            colors of the vertical lines + error regions used to visualize supplied literature values
            currently only used for 1d plotting
            `colors` is pre-defined to offer well-contrasted colors for up to the first six literature values/errors/labels
            if `colors` format is array, array values must be strings
            
    RETURNS
        ASCII data table with the following columns
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

    # reformat and determine the number of ids provided based on the input type
    if isinstance(ids, list): N_ids = len(ids); ids = str(ids)[1:-1]
    elif isinstance(ids, str): N_ids = 1
    else: print('Input "ids" must be either \'str\' or list/tuple of \'str\''); os._exit(1)

    # set up our query to the Gaia DR3 catalog
    query = """
    SELECT g.source_id, g.phot_rp_mean_mag, g.parallax, g.parallax_error, g.phot_rp_mean_flux, g.phot_rp_mean_flux_error,
           g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.astrometric_excess_noise,
           g.b, g.phot_g_mean_mag, g.bp_rp, g.phot_rp_mean_flux_over_error
    FROM gaiadr3.gaia_source as g
    WHERE g.source_id IN (""" + ids + """)
    """

    # excecute the query and store the results
    print('Querying the GAIA Archive...')
    job = Gaia.launch_job(query)
    r = job.get_results()
    
    # first, check that stars have valid parallax and rp mag measurements
    valid_mask = ~(r['parallax'] > 0) | np.isnan(r['phot_rp_mean_mag'])

    # the program cannot estimate masses for stars without vaild parallax or rp mag measurements; remove them
    if np.sum(valid_mask) == N_ids:
        print('Sorry, none of the ids provided have valid parallax or RP mag measurements.'); os._exit(1)
    elif np.sum(valid_mask) > 0:
        print('The following ' + str(int(np.sum(valid_mask))) + ' ids have invalid parallax or RP mag measurements.')
        for i, plx, app_rp in zip(r['source_id'][valid_mask], r['parallax'][valid_mask], r['phot_rp_mean_mag'][valid_mask]):
            try:
                print(str(i) + ' (plx = ' + str(plx) + ', rp_mag = ' + str(round(app_rp, 2)) + ')')
            except:
                try:
                    print(str(i) + ' (plx = ' + str(round(plx, 2)) + ', rp_mag = ' + str(app_rp) + ')')
                except:
                    print(str(i) + ' (plx = ' + str(plx) + ', rp_mag = ' + str(app_rp) + ')')
        # adjust queried results by removing invalid stars
        r = r[~valid_mask]
        N_ids = len(r)
    else:
        print('All ids provided have valid parallax and RP mag measurements.')

    # next, check which stars are not eligible for our relation
    source_ids  = np.array(r['source_id'])
    rp_mags     = np.array(r['phot_rp_mean_mag'])
    plxs        = np.array(r['parallax'])
    abs_rp_mags = absmag(rp_mags, plxs)
    ineligible_mask = (abs_rp_mags < 4) | (abs_rp_mags > 14.5)

    # the program cannot estimate masses for stars that are outside the recommended range; remove them
    if np.sum(ineligible_mask) == N_ids:
        print('Sorry, all ids provided are outside the eligible range (4.0 < MGRP < 14.5).'); os._exit(1)
    elif np.sum(ineligible_mask) > 0:
        print('The following ' + str(int(np.sum(ineligible_mask))) + ' ids are outside the eligible range (4.0 < MGRP < 14.5).')
        for i, abs_rp in zip(source_ids[ineligible_mask], abs_rp_mags[ineligible_mask]):
            print(str(i) + ' (abs_rp_mag = ' + str(round(abs_rp, 2)) + ')')
        # adjust queried results by removing ineligble stars
        r = r[~ineligible_mask]
        N_ids = len(r)
    else:
        print('All ids provided are within the eligible photometric range (4.0 < MGRP < 14.5).')

    # re-itemize downloaded parameters for gaia source id(s)
    source_ids        = np.array(r['source_id'])
    rp_mags           = np.array(r['phot_rp_mean_mag'])
    plxs              = np.array(r['parallax'])
    plx_errs          = np.array(r['parallax_error'])
    rp_fluxes         = np.array(r['phot_rp_mean_flux'])
    rp_flux_errs      = np.array(r['phot_rp_mean_flux_error'])
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
    with fits.open(os.path.join(pkg_resources.resource_filename('gorp_mass', 'resources'), 'massmag_errs_new.fits')) as hdul:

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
    
    # if True, 1d histograms will be saved for each input Gaia source id
    if plot_1d:
    
        # check to see if `plot_path` exists. This is helpful for default `gorp_plots` directory
        if not os.path.exists(plot_path): os.mkdir(plot_path)
    
        # provide the number of x-bins in the histogram
        nbins = 100
    
        # create and save 1d plots for each given Gaia source id
        for i in range(len(source_ids)):
        
            print('Plotting 1d histogram of ' + str(source_ids[i]) + '...')
        
            # plot a 1d probability density of masses for the given Gaia source
            plt.figure(figsize = (6.40, 3.95))
            data = np.random.normal(masses[i], mass_errs[i], N)
            bins, edges = plt.hist(data, histtype = 'step', color = 'darkblue', lw = 2, bins = nbins, density = True, zorder = 10)[:2]
            mask = (edges > masses[i] - mass_errs[i]) & (edges < masses[i] + mass_errs[i])
            fx = np.concatenate((np.array([masses[i] - mass_errs[i]]), sum([(j, j) for j in edges[mask]], ()), (np.array([masses[i] + mass_errs[i]]))))
            fy = np.concatenate((np.array(2 * [bins[np.argmax(edges[edges < masses[i] - mass_errs[i]])]]), sum([(j, j) for j in bins[mask[:-1]]], ())))
            plt.fill_between(fx, y1 = 0, y2 = fy, color = 'gray', alpha = 0.3)
            plt.axvline(x = masses[i], ls = '--', lw = 2, color = 'gray', ymax = bins[np.argmin(np.abs(edges - masses[i]))] / plt.axis()[-1], label = 'Giovinazzi & Blake 2022')
            plt.xlabel(r'$M_\mathrm{GORP}~\left[\mathrm{M_\odot}\right]$')
            plt.ylabel('Probability Density')
            plt.ylim(0, plt.axis()[-1])
            
            # ~automatic legend handling; legend placement is known to not be optimal, but for most cases this seeems to do a good job
            if lit_vals is not None:
                if isinstance(lit_vals, float): lit_vals = np.array([lit_vals])
                if lit_labs is not None:
                    if isinstance(lit_labs, str): lit_labs = np.array([lit_labs])
                    for j in range(len(lit_vals)):
                        plt.axvline(x = lit_vals[j], label = lit_labs[j], ls = '--', lw = 2, zorder = 1000, color = colors[j])
                    m = max(len(max(max(lit_labs, key = len))), 23)
                    if (np.abs(masses[i] - np.max(lit_vals)) >= np.abs(masses[i] - np.min(lit_vals))) or (np.max(lit_vals) <= (masses[i] + mass_errs[i])):
                        plt.xlim(xmax = max(max(lit_vals) + mass_errs[i], plt.xlim()[1], masses[i] + mass_errs[i] + m * masses[i]**0.75 / 100))
                        loc = 'upper right'
                    elif (np.abs(masses[i] - np.max(lit_vals)) < np.abs(masses[i] - np.min(lit_vals))) or (np.min(lit_vals) > (masses[i] - mass_errs[i])):
                        plt.xlim(xmin = min(min(lit_vals) - mass_errs[i], plt.xlim()[0], masses[i] - mass_errs[i] - m * masses[i]**0.75 / 100))
                        loc = 'upper left'
                    leg = plt.legend(loc = loc, edgecolor = None).set_zorder(10000)
                else:
                    for j in range(len(lit_vals)):
                        plt.axvline(x = lit_vals[j], ls = '--', lw = 2, zorder = 1000, color = colors[j])
                if lit_errs is not None:
                    if isinstance(lit_errs, float): lit_errs = np.array([lit_errs])
                    for j in range(len(lit_vals)):
                        plt.fill_between(x = np.array([lit_vals[j] - lit_errs[j], lit_vals[j] + lit_errs[j]]), y1 = 0, y2 = plt.axis()[-1], alpha = 0.3, zorder = 100 + j, color = colors[j])
                        
            plt.savefig(os.path.join(plot_path, str(source_ids[i]) + '_1d.' + plot_ext), bbox_inches = 'tight', dpi = 1000)
            plt.close('all')
            
    # if True, 2d histograms (mass vs magnitude) will be saved for each input Gaia source id
    if plot_2d:
    
        # check to see if `plot_path` exists. This is helpful for default `gorp_plots` directory
        if not os.path.exists(plot_path): os.mkdir(plot_path)
        
        # define the Vega mag system RP zeropoint from Gaia
        rp_zp = np.random.normal(24.7478955012, 0.0037793818, N)
    
        # provide the number of xy-bins in the histogram
        nbins = 100
    
        # create and save 2d plots for each given Gaia source id
        for i in range(len(source_ids)):
        
            print('Plotting 2d histogram of ' + str(source_ids[i]) + '. This may take a minute.')

            # construct empty mass-magnitude grid given input information from Gaia and standard errors from our relation
            app_rps = -2.5 * np.log10(np.random.normal(rp_fluxes[i], rp_flux_errs[i], N)) + rp_zp
            if np.sum(np.isnan(app_rps)) > N / 10.: print('2d plot of ' + str(source_ids[i]) + ' failed. RP_flux / error is too small.'); continue
            else: app_rps = app_rps[~np.isnan(app_rps)]
            abs_rps = absmag(app_rps, np.random.normal(plxs[i], plx_errs[i], N))
            if np.sum(np.isnan(abs_rps))  > N / 10.: print('2d plot of ' + str(source_ids[i]) + ' failed. plx / error is too small.'); continue
            else: abs_rps = abs_rps[~np.isnan(abs_rps)]
            abs_rps = abs_rps[(abs_rps > 4) & (abs_rps < 14.5)]
            x_bins, x_edges = np.histogram(abs_rps, bins = nbins)
            bin_centers = (x_edges[1:] + x_edges[:-1]) / 2.
            grid = np.zeros(len(bin_centers)**2).reshape(len(bin_centers), len(bin_centers))
            log10masses = A + B * abs_rps - C * (1. + erf(abs_rps - 9.5))
            log10mass_errs = f(abs_rps)
            mass_bins = np.linspace(10**(np.min(log10masses) - 3 * log10mass_errs[np.argmin(log10masses)]),
                                    10**(np.max(log10masses) + 3 * log10mass_errs[np.argmax(log10masses)]), nbins + 1)

            # populate the mass-magnitude grid with mass distributions for each sample of absolute magnitude
            for j in range(len(abs_rps)):
            
                masses = 10**np.random.normal(log10masses[j], log10mass_errs[j], int(np.sqrt(N)))
                y_bins = np.histogram(masses, bins = mass_bins)[0]
                grid_ind = np.argmin(np.abs(bin_centers - abs_rps[j]))
                grid[:, grid_ind] += y_bins
                
            # plot a 2d probability density of masses and magnitudes for the given Gaia source
            fig, axs = plt.subplots(2, 2, figsize = (6.40, 3.95), gridspec_kw = {'width_ratios': [1, 0.67], 'height_ratios': [0.67, 1]})
            fig.delaxes(axs[0, 1])
            fig.subplots_adjust(hspace = 0, wspace = 0)
            axs[1, 0].imshow(grid, cmap = 'Blues', extent = (np.min(abs_rps), np.max(abs_rps), np.max(mass_bins), np.min(mass_bins)), aspect = 'auto')
            axs[0, 0].stairs(np.sum(grid, axis = 0), x_edges, color = 'darkblue', lw = 2)
            axs[1, 1].stairs(np.sum(grid, axis = 1), mass_bins, color = 'darkblue', lw = 2, orientation = 'horizontal')
            axs[1, 0].set_xlim(np.max(abs_rps), np.min(abs_rps))
            axs[1, 0].set_ylim(np.min(mass_bins), np.max(mass_bins))
            axs[0, 0].set_xlim(np.max(abs_rps), np.min(abs_rps))
            axs[1, 1].set_ylim(np.min(mass_bins), np.max(mass_bins))
            axs[1, 0].set_xlabel(r'$M_{G_\mathrm{RP}}$')
            axs[1, 0].set_ylabel(r'$M_\mathrm{GORP}~\left[\mathrm{M_\odot}\right]$')
            for side in ['top', 'right', 'bottom', 'left']:
                for s in [0, 1]:
                    axs[s, s].spines[side].set_visible(False)
                    axs[s, s].tick_params(left = False, right = False, labelleft = False, labelbottom = False, bottom = False)
                    
            plt.savefig(os.path.join(plot_path, str(source_ids[i]) + '_2d.' + plot_ext), bbox_inches = 'tight', dpi = 1000)
            plt.close('all')
            
    # return a table of output values describing the mass posterior(s)
    return t
    
# function useed to estimate masses from RP magnitudes
def rp_posterior(abs_rp_mags, N = 100000):

    '''
    Estimate a posterior on mass given absolute magnitudes in the Gaia RP bandpass
    Relation derived in Giovinazzi, M. R. & Blake, C. H. (2022)
    
    If you use this code, please cite the relevant paper:
    https://ui.adsabs.harvard.edu/abs/2022AJ....164..164G
    
    PARAMETERS
        abs_rp_mags: float, array
            absolute magnitudes in the Gaia RP bandpass
            if `abs_rp_mags` format is array, array values must be floats
        N: int
            number of random samples used in estimating uncertainty in mass; default = 100,000
            
    RETURNS
        ASCII data table with the following columns
            mass: float, array
                point estimate on the mass(es) corresponding to the provided GAIA Source ID(s)
            mass_sig: float, array
                standard error on the mass(es)
    '''

    if isinstance(abs_rp_mags, float): abs_rp_mags = np.array(abs_rp_mags)
    N_rps = len(abs_rp_mags)
    ineligible_mask = (abs_rp_mags < 4) | (abs_rp_mags > 14.5)

    # the program cannot estimate masses for stars that are outside the recommended range
    if np.sum(ineligible_mask) == len(abs_rp_mags):
        print('Sorry, all RP mags provided are outside the eligible range (4.0 < MGRP < 14.5). Make sure you are using absolute mags!'); os._exit(1)
    elif np.sum(ineligible_mask) > 0:
        print('Some RP mags are outside the eligible range (4.0 < MGRP < 14.5). Ignoring those.')
        # adjust queried results by remeoving ineligble stars
        abs_rp_mags = abs_rp_mags[~ineligible_mask]
    else:
        print('All RP mags provided are within the eligible photometric range (4.0 < MGRP < 14.5).')

    # load in uncertainties from file
    with fits.open(os.path.join(pkg_resources.resource_filename('gorp_mass', 'resources'), 'massmag_errs_new.fits')) as hdul:

        points = hdul[1].data['X'][0]
        massmag_errs = hdul[1].data['ERF_DATA'][0]
    
    massmag_errs = np.std(massmag_errs, axis = 1)
    f = interpolate.interp1d(points, massmag_errs)

    # define the relation from the paper
    A, B, C = 0.444687, -0.0968913, 0.0746491
    log10masses = A + B * abs_rp_mags - C * (1. + erf(abs_rp_mags - 9.5))
    masses = 10**log10masses
    try: mass_errs = np.array([np.std(10**np.random.normal(log10masses[i], f(abs_rp_mags[i]), N)) for i in range(N_rps)])
    except: mass_errs = np.std(10**np.random.normal(log10masses, f(abs_rp_mags), N))
    
    # columnize the results and write to an ascii table
    t = Table()
    t['abs_rp_mag'] = abs_rp_mags
    t['mass'] = masses
    t['mass_err'] = mass_errs
    
    print('Writing to file...')
    ascii.write(t, 'gorp_rp_results.dat', overwrite = True)
    print('Done!')

    # return a table of output values describing the mass posterior(s)
    return t
