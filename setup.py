

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def setupPackage():

    setup(name='qdune',
          version='0.1.0',
          description='dune point cloud processing tools',
          url='https://github.com/tashley/qDune'
          author='Thomas Ashley',
          author_email='tashley22@gmail.com',
          license='MIT',
          packages=['qdune'],
          zip_safe=False,
          package_data={'qdune': ['*.xyz', '*.csv']})


if __name__ == '__main__':
    setupPackage()
