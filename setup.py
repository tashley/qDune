

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def setupPackage():
#    setup(name='qDune',
#          description='dune point cloud processing tools',
#          classifiers=[
#              'Intended Audience :: Science/Research',
#              'Intended Audience :: Developers',
#              'License :: MIT)',
#              'Programming Language :: Python',
#              'Programming Language :: Python :: 3.5',
#              'Topic :: Scientific/Engineering'],
#          author='Thomas Ashley',
#          url = 'URL to get it at.',
#          author_email='tashley22@gmail.com',
#          version='1.0',
#          license='MIT',
#          packages=['qDune'],
#          keywords=['dunes', 'bedload'],
#          platforms='Windows',
#          package_data={'qDune': ['*.xyz', '*.csv']})

    setup(name='qdune',
          version='0.1.0',
          description='dune point cloud processing tools',
          author='Thomas Ashley',
          author_email='tashley22@gmail.com',
          license='MIT',
          packages=['qdune'],
          zip_safe=False,
          package_data={'qdune': ['*.xyz', '*.csv']})


if __name__ == '__main__':
    setupPackage()
