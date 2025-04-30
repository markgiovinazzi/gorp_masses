from setuptools import setup, find_packages

setup(
    name='gorp_mass',
    version='1.0',
    packages=find_packages(),
    include_package_data=True,
    package_data={'gorp_mass': ['resources/*']},
    zip_safe=False,
)
