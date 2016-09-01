from setuptools import setup
import numpy as np


setup(name='qdune',
      version='0.1.6',
      description='dune point cloud processing tools',
      author='Thomas Ashley',
      author_email='tashley22@gmail.com',
      url='https://github.com/tashley/qDune',
      license='MIT',
      packages=['qdune'],
      include_package_data=True,
      include_dirs = [np.get_include()],
      package_data={'qdune': ['*.txt','*.csv',]},
      install_requires = ['numpy', 'pandas', 'matplotlib']
)
