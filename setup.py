from setuptools import setup, find_packages

setup(
    name='tksamc',
    version='1.0.0',
    description='TKSA - Electrostatic Free Energy calculation for each ionizable residue',
    author='Vinícius G. Contessoto',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'mdtraj',
        'numba',
    ],
    entry_points={
        'console_scripts': [
            'tksamc=tksamc.cli:main',
        ],
    },
)
