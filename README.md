README

This program estimates masses of stars in the Gaia catalog using the photometric mass relation derived in Giovinazzi, M. R. & Blake, C. H. (2022).

If you use this code, please cite the relevant paper:
[https://ui.adsabs.harvard.edu/abs/2022AJ....164..164G](url)

To get started, you should clone a copy of the code to your home machine.

```python
git clone https://github.com/markgiovinazzi/gorp_masses
cd gorp_masses
python3 setup.py sdist
pip3 install dist/gorp_mass-1.0.tar.gz
```

There are two functions available for use within this repository, `rp_posterior()` and `gaia_posterior()`. The mass-magnitude relation is built on top of Gaia photometry, specifically using the instrument's RP bandpass. If you have a single absolute RP magnitude, or an array of them, you can feed that into `rp_posterior()`, which takes a singular argument in the form of absolute RP magnitude(s) and returns an ASCII data table with the following columns: absolute RP magnitude, GORP mass estimate, GORP mass error.

The `gaia_posterior()` function, however, is designed to accept Gaia DR3 source ID(s), download the photometry, and return a GORP mass estimate. As an example, we can query the exoplanet archive for systems where host stars have published masses and registered Gaia source IDs. Suppose those IDs are stored in a variable called `gaia_ids`. We can execute the following code.

```python
from gorp_mass import gaia_posterior

posterior_table = gaia_posterior(gaia_ids)
gorp_masses = posterior_table['mass']
gorp_mass_errs = posterior_table['mass_err']
```

Then, our GORP masses and corresponding errors can be plotted against those published in the exoplanet archive, as a comparison. Here is just that.

![test](https://user-images.githubusercontent.com/14206224/211969877-5eae4ee3-f63f-4e9f-8e0a-a9a73e1515c2.jpeg)

Not bad! The mass-magnitude relation is applicable to stars with `4 < absolute_RP_mag < 14.5`, which corresponds roughly to 0-1 solar mass. The code will alert you of stars with ineligible absolute RP magnitudes, invalid parallaxes, or other indicators of poor measurement, and will skip over them. The data table also returns quality flags to help identify the trustworthiness of a given GORP mass estimate. The flags are described below.

```
is_MS: bool, array
  flag determining whether a star is main sequence
  True if the star is main sequence, False otherwise
  this relation is NOT applicable to stars that return `is_MS = False`
is_ss: bool, array
  flag relaying whether Gaia registers the star as a good fit to a single-star astrometric solution
  True if the star is likely a single star with a good astrometric solution from Gaia, False otherwise
  if False, user should be skeptical of returned mass estimate
phot_par: bool, array
  flag for photometric and astrometric parameters that may bias our mass estimates
  True if apparent `g_mag > 5`, `rp_flux_over_err > 10`, `rel_plx_err < 0.2`, `plx_err < 2` mas, `plx > 1` mas, False otherwise
  if False, photometry may be poor or saturated, or absolute magnitude (due to plx) is uncertain; use mass estimate with caution
ext_flag: bool, array
  flag determining whether or not the star may be subjected to significant extinction
  True if galactic latitude `|b| > 10` degrees, False otherwise
  if False, be mindful that extinction may make star appear fainter, causing potential underestimate in mass
```

Additionally, the program has some bonus plotting features included. For example, what if we want to make the 1d GORP mass posterior for TRAPPIST-1 (Gaia DR3 2635476908753563008)? The next two lines of code will save the following plot to a user-specified directory that defaults to `gorp_plots`.

```
from gorp_mass import gaia_posterior
gaia_posterior('2635476908753563008', plot_1d = True)
```

![2635476908753563008_1d](https://user-images.githubusercontent.com/14206224/211975098-2bb889a7-9732-45f4-b77e-70eadadfe611.jpeg)

Ducrot et al. 2020 derived a mass estimate for TRAPPIST-1 of `m = 0.0898 +/- 0.0023` solar masses. We can include it in our plot by passing information through the arguments `lit_vals`, `lit_errs`, and `lit_labs` like so.

```
from gorp_mass import gaia_posterior
gaia_posterior('2635476908753563008', plot_1d = True, lit_vals = 0.0898, lit_errs = 0.0023, lit_labs = 'Ducrot et al. 2020')
```

![2635476908753563008_1d](https://user-images.githubusercontent.com/14206224/211978009-e78a445e-186e-43cd-bc7a-87aa4b6abc71.jpeg)

The arguments `lit_vals`, `lit_errs`, and `lit_labs` can accept either singular values or arrays of multiple values to include additional literature values on the plot for direct comparison. If `lit_vals` is given, `lit_errs` and `lit_labs` can be used optionally. If `lit_labs` is given, a legend will be added to the plot (automatically including this paper's reference as well). The program will attempt to find an optimal placement for the legend given the range of `lit_vals`, but legend placement is currently known to not always be perfect. Future updates to this code should correct for this.

And, while Gaia photometry is fantastic, some stars can still be poorly measured. So, how trustworthy is a star's photometric mass estimate if its flux, parallax, etc has large errors? To understand this, we can visualize the GORP mass estimate in a 2d sense. Let's first explore this with TRAPPIST-1 again, this time by turning on `plot_2d`.

```
from gorp_mass import gaia_posterior
gaia_posterior('2635476908753563008', plot_2d = True)
```

![2635476908753563008_2d](https://user-images.githubusercontent.com/14206224/211978840-0e95dd09-bfdb-453e-906f-2d82c89e09c1.jpeg)

Beautiful! TRAPPIST-1 is pristinely measured in Gaia's RP bandpass. But, not all stars are. For example, Gaia DR3 1061855729642662144's absolute RP magnitude has a width that spans nearly a full magnitude, so we should be skeptical of the GORP mass derived from its absolute RP magnitude point estimate.

```
from gorp_mass import gaia_posterior
gaia_posterior('1061855729642662144', plot_2d = True)
```

![1061855729642662144_2d](https://user-images.githubusercontent.com/14206224/211983210-101de1d6-0538-4bd4-bf8e-af0437086df3.jpeg)

If any bugs are found, or if there are suggestions on how this code can be improved, please feel free to contact me at mgiovinazzi@amherst.edu.
