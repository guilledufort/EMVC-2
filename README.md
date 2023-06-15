# EMVC-2
## An efficient SNV variant caller based on the expectation maximization algorithm. EMVC-2 is implemented in C and uses a python wrapper.
### Authors: Guillermo Dufort y Álvarez, Martí Xargay, Idoia Ochoa, and Alba Pages-Zamora
### Contact: gdufort@fing.edu.uy

## Install with Conda
To install directly from source, follow the instructions in the next section.

EMVC-2 is available on conda via the bioconda channel. See [this](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) page for installation instructions for conda. Once conda is installed, do the following to install emvc-2.
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install emvc-2
```
Note that if emvc-2 is installed this way, it should be invoked with the command `emvc-2` rather than `./emvc-2`. The bioconda [help page](https://bioconda.github.io/tutorials/gcb2020.html#conda-environments) shows the commands if you wish to install emvc-2 in an environment.

## Install from source code

### Download repository
```bash
git clone https://github.com/guilledufort/EMVC-2.git
```
### Requirements

Software requirements
1. python ( >= 3.8.1, <= 3.8.5 )
2. samtools ( == 1.9 )

Compiler requirement
1. gcc ( >= 4.8.1 )

Python libraries requirement
1. cython ( >=0.29.17 ),
2. numpy ( >=1.16.6,<=1.20.3 ),
3. argparse ( >=1.1 ),
4. pysam ( >=0.15.4,<=0.16.0.1 ),
5. scipy ( >=1.1.0,<1.5.4 ),
6. tqdm ( >=4.46.0 ),
7. scikit-learn ( >=0.22.2,<=0.24.2 ),

### Compiling the *candidate_variants_finder* 

The following instructions will create the *candidate_variants_finder* executable in the root directory, which is needed to run EMVC-2.
To compile *candidate_variants_finder* you need to have the gcc compiler. 

On Linux (Ubuntu or CentOS) gcc usually comes installed by default, but if not run the following:
```bash
sudo apt update
sudo apt-get install gcc
```

On macOS, install GCC compiler:
- Install HomeBrew (https://brew.sh/)
- Install GCC (this step will be faster if Xcode command line tools are already installed using ```xcode-select --install```):
```bash
brew update
brew install gcc@9
```

To check if the gcc compiler is properly installed in your system run:

On Linux
```bash
gcc --version
```
On MacOS:
```bash
gcc-9 --version
```
The output should be the description of the installed software.

To compile *candidate_variants_finder* run:
```bash
cd EMVC-2/
make
```
### Install python dependencies

To install the *python dependencies* run:
```bash
cd EMVC-2/
python setup.py install
```

### Install samtools

To install *samtools*, you can use *conda*:
```bash
conda install -c bioconda samtools==1.9
```
or follow the instructions in [the github repository](https://github.com/samtools/samtools).

## Usage

```console 

python EMVC-2.py [-h] -i BAM_FILE -r REF_FILE [-p THREADS] [-t ITERATIONS] [-m LEARNERS] [-v VERBOSE] -o OUT_FILE

optional arguments:
  -h, --help            show this help message and exit
  -i BAM_FILE, --bam_file BAM_FILE
                        The bam file
  -r REF_FILE, --ref_file REF_FILE
                        The reference fasta file
  -p THREADS, --threads THREADS
                        The number of parallel threads (default 8)
  -t ITERATIONS, --iterations ITERATIONS
                        The number of EM iterations (default 5)
  -m LEARNERS, --learners LEARNERS
                        The number of learners (default 7)
  -v VERBOSE, --verbose VERBOSE
                        Make output verbose (default 0)
  -o OUT_FILE, --out_file OUT_FILE
                        The output file name

```


## Usage example
We add an *example* folder with a test file to run a simple example of the tool. The hs37d5 reference file must be downloaded following the instructions detailed in the previous section for the example to work.
If installed using conda, use the command `emvc-2` instead of `python EMVC-2.py`.

To run the variant caller with 8 threads on the example file *example.bam*:
```bash
cd EMVC-2
python EMVC-2.py -i example/example.bam -r reference/hs37d5/hs37d5.fa.gz -p 8 -o example/example.vcf
```

## Original paper datasets information

To test the performance of the EMVC-2 SNV variant caller we ran experiments on the following datasets.

| Dataset        | Reference | Size (GB) | Coverage | Sequencing Method     | Download link                                                                                                                                                       |
|----------------|-----------|-----------|----------|-----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ERR262997      | HG001     | 104       | 30       | Illumina HiSeq 2000   | [link](https://www.ebi.ac.uk/ena/browser/view/ERA207860?show=reads)                                                                                            |
| NovaSeq        | HG001     | 49        | 25       | Illumina NovaSeq 6000 | is not available for download                                                                                                                                                                   |
| Ashkenazim son | HG002     | 48        | 25       | Illumina HiSeq 2500   | [link](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/10XGenomics)                                                          |
| pangenomics2   | HG002     | 61        | 30       | Illumina HiSeq 2500   | [link](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/)               |
| pangenomics3   | HG003     | 66        | 30       | Illumina HiSeq 2500   | [link](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG003/) |
| pangenomics4   | HG004     | 61        | 30       | Illumina HiSeq 2500   | [link](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG004/) |
| Chinese Son    | HG005     | 34        | 15       | Illumina HiSeq 2500   | [link](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/NIST_Stanford_Illumina_6kb_matepair/fastqs/)         


### Downloading the datasets and the reference genome

To download a dataset you have to run the *download_script.sh* with the specific dataset name as a parameter.
For example, to download *ERR262997* run:
```bash
cd EMVC-2/datasets
./download_script.sh ERR262997
```

To download the human reference genome version hs37d5 run:
```bash
cd EMVC-2/reference
./download_script.sh hs37d5
```

The scripts use the command *curl* to perform the download. 
To install *curl* on macOS run:
 ```bash
brew install curl
```
To install *curl* on Ubuntu or CentOS run:
 ```bash
sudo apt-get install curl
```

## Alignment information

To obtain alignment information in [BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) for each pair of FASTQ files we recommend using the tool [BWA](https://github.com/lh3/bwa).

To install bwa with conda run:
 ```bash
conda install bwa
```

To align a pair of FASTQ files against a reference genome using BWA run:
 ```bash
bwa mem -t  [-@THREADS] [REF] [FASTQ_R1] [FASTQ_R2] \
    | samtools sort [-@THREADS] -o [BAM_FILE] 
```
