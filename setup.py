from setuptools import setup, find_packages

setup(
    name='EMVC-2',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'argparse==1.1',
        'pysam==0.20.0',
        'scipy==1.4.1',
        'tqdm==4.46.0',
    ],
)