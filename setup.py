from setuptools import setup, find_packages

setup(name='gorp',
      version='1.0',
      py_modules=['gorp'],
      packages=find_packages() + ['err_file'],
      include_package_data=True
      )
