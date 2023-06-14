from setuptools import setup, find_packages

setup(
    name='EMVC-2',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'cython>=3.0.0',
        'numpy>=1.25.0',
        'argparse>=1.4.0',
        'pysam>=0.21.0',
        'scipy>=1.11.0',
        'tqdm>=4.65.0',
        'scikit-learn>=1.2.2'
    ],
)