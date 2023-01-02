from setuptools import setup, find_packages

setup(name = 'gorp_mass',
      version = '1.0',
      py_modules = ['gorp_mass'],
      packages = find_packages() + ['resources'],
      include_package_data = True
      )
