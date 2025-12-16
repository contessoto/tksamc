from setuptools import setup, find_packages
import os

# Read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='tksamc',
    version='1.0.0',
    description='TKSA - Electrostatic Free Energy calculation for each ionizable residue',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Vinícius G. Contessoto',
    url='https://github.com/contessoto/tksamc',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'mdtraj',
        'numba',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    entry_points={
        'console_scripts': [
            'tksamc=tksamc.cli:main',
        ],
    },
    python_requires='>=3.6',
)
