from distutils.core import setup

setup(
    name='biopy',
    version='0.1.1',
    author='B.A. van den Berg',
    author_email='b.a.vandenberg@gmail.com',
    packages=['biopy'],
    package_dir={'biopy': 'biopy'},
    package_data={'biopy': ['data/*/*']},
    url='http://pypi.python.org/pypi/biopy/',
    license='LICENSE.txt',
    description='A little Swiss Army knife for bioinformaticians',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.7.1",
        "scipy >= 0.12.0",
        "matplotlib >= 1.2.2"
    ]
)
