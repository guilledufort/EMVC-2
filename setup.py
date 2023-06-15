from setuptools import setup, find_packages

setup(
    name='EMVC-2',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'cython>=0.29.17',
        'numpy>=1.16.6,<=1.20.3',
        'argparse>=1.1',
        'pysam==0.20.0',
        'scipy>=1.1.0,<1.5.4',
        'tqdm>=4.46.0',
        'scikit-learn>=1.2.1'
    ],
)
# En lnano funciona 3.8 pero 3.9 da conflictos con samtools
# En condor funciona 3.9 y no da conflicto con nada ,pero las nuevas versiones de los paquetes dan conflicto en 3.8. Espec'ificamente numpy.
