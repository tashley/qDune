from setuptools import setup


setup(name='qdune',
      version='0.1.1',
      description='dune point cloud processing tools',
      author='Thomas Ashley',
      author_email='tashley22@gmail.com',
      url='https://github.com/tashley/qDune',
      license='MIT',
      packages=['qdune'],
      package_data={'qdune': ['*.xyz', '*.csv']})
