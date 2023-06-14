# MIT License

# Copyright (c) 2023 Guillermo Dufort y Álvarez

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# Import the necessary modules
import argparse # 1.1
import pysam # 0.20.0
import scipy.stats as stats # 1.4.1
import tqdm # 4.46.0
import os
import multiprocessing
import shutil
import time
import datetime
import pickle

# Define a function to sort a bam file if it is not sorted and return the sorted file name
def sort_bam(bam_file):
    # Check if the bam file is sorted by checking the header
    header = pysam.view("-H", bam_file)
    # Check if the header contains both @HD and SO:coordinate in the same line
    if "@HD" in header and "SO:coordinate" in header.split("@HD")[1]:
        # print("The bam file is already sorted.")
        # If sorted, return the original file name
        sorted_file = bam_file
    else:
        # If not sorted, sort it using pysam.sort and save it as a new file with _sorted suffix
        print("Sorting the bam file...")
        sorted_file = bam_file.replace(".bam", "_sorted.bam")
        pysam.sort("-o", sorted_file, bam_file)
    # Return the sorted file name
    return sorted_file

# Define a function to call variants for a given contig using the system command
def variant_call(args):
    
    contig, temp_folder, ref_file, bam_file, niterations, learners, verbose = args
    # Use the command to generate a vcf file for the contig
    # print(f"Calling variants for contig {contig}...")
    out_name = f"{temp_folder}/{contig}.vcf"
    # Get the path of the current script file
    script_path = os.path.abspath(__file__)
    # Get the path of the binary file
    binary_path = os.path.join(os.path.dirname(script_path), 'candidate_variants_finder')
    # Replace the file names and parameters in the command with the corresponding variables
    command = f"samtools mpileup -B -C 0 -d 0 -f {ref_file} -Q 0 --ff 4,8,10,20,40,80,100,200,400,800 {bam_file} -r {contig} 2> /dev/null | {binary_path} - {out_name} {learners} {niterations} - - - - - - {verbose}"
    if verbose:
        print(command)
    # Execute the command using os.system
    os.system(command)
    # Return the vcf file name
    return out_name

def features_from_line(line, features):
    line_s = line.split("\t")
    # print (line_s)
    POS = int(line_s[1])
    QUAL = int(line_s[5])
    line_s_f = line_s[9].split(":")
    # print (line_s_f)
    GT = int(line_s_f[0].split('/')[0] + line_s_f[0].split('/')[1])
    DP = int(line_s_f[1])
    FR = int(line_s_f[2])
    FA = int(line_s_f[3])
    RR = int(line_s_f[4])
    RA = int(line_s_f[5])
    PROB_AA = float(line_s_f[6])
    PROB_CC = float(line_s_f[7])
    PROB_TT = float(line_s_f[8])
    PROB_GG = float(line_s_f[9])
    PROB_AC = float(line_s_f[10])
    PROB_AT = float(line_s_f[11])
    PROB_AG = float(line_s_f[12])
    PROB_CT = float(line_s_f[13])
    PROB_CG = float(line_s_f[14])
    PROB_TG = float(line_s_f[15][:-1])
    x = []
    for f in features:
        if f == 'QUAL':
            x.append(QUAL)
        elif f == 'GT':
            x.append(GT)
        elif f == 'DP':
            x.append(DP)
        elif f == 'FR':
            x.append(FR)
        elif f == 'FA':
            x.append(FA)
        elif f == 'RR':
            x.append(RR)
        elif f == 'RA':
            x.append(RA)
        elif f == 'PROB_AA':
            x.append(PROB_AA)
        elif f == 'PROB_CC':
            x.append(PROB_CC)
        elif f == 'PROB_TT':
            x.append(PROB_TT)
        elif f == 'PROB_GG':
            x.append(PROB_GG)
        elif f == 'PROB_AC':
            x.append(PROB_AC)
        elif f == 'PROB_AT':
            x.append(PROB_AT)
        elif f == 'PROB_AG':
            x.append(PROB_AG)
        elif f == 'PROB_CT':
            x.append(PROB_CT)
        elif f == 'PROB_CG':
            x.append(PROB_CG)
        elif f == 'PROB_TG':
            x.append(PROB_TG)
        elif f == 'ALT_%':
            x.append(float(FA + RA) * 100 / float(DP))
        elif f == 'CLASS_Entropy':
            probs = [PROB_AA, PROB_CC,  PROB_TT, PROB_GG, PROB_AC, PROB_AT,  PROB_AG, PROB_CT, PROB_CG, PROB_TG]
            # Calculate entropy of probs with log2
            x.append(stats.entropy(probs, base=2, axis=0))

    return x

def filter_variants(args):
    original_vcf, clf, features, filtered_name = args
    o_vcf = open(original_vcf, 'r')
    out_f = open(filtered_name, 'w')
    line_o = ""
    while "#CHROM" not in line_o:
        line_o = o_vcf.readline()
        out_f.write(line_o)

    line_o = o_vcf.readline()
    x = []
    while line_o:
        x.append(features_from_line(line_o, features))
        line_o = o_vcf.readline()
   
    o_vcf.close()

    # exit if there are no variants
    if len(x) == 0:
        return filtered_name, 0, 0
    
    o_vcf = open(original_vcf, 'r')

    line_o = ""
    while "#CHROM" not in line_o:
        line_o = o_vcf.readline()
    
    preds = clf.predict(x)

    i = 0

    variants_kept = 0
    variants_filtered = 0

    line_o = o_vcf.readline()
    
    while line_o:    
        if preds[i] == 'TP':
            out_f.write(line_o)
            variants_kept += 1
        else:
            variants_filtered += 1
            
        i += 1
        line_o = o_vcf.readline()

    o_vcf.close()
    out_f.close()

    return filtered_name, variants_kept, variants_filtered

# Define a function that takes a database name and a version name as parameters
def join_vcf_files(files, output_file):
    # Set the first variable to 1 outside the for loop
    first = 1
    for file_name in files:
        # Open the file using the open function
        with open(file_name) as f:
            # If it is the first chromosome, copy the whole file to the output file
            if first == 1:
                # Change the first variable to 0 after processing the first file
                first = 0
                with open(output_file, "w") as o:
                    for line in f:
                        o.write(line)
            # Otherwise, skip the header lines and append the variant lines to the output file
            else:
                with open(output_file, "a") as o:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        else:
                            o.write(line)

def main():
    try:
        # Create an argument parser
        parser = argparse.ArgumentParser(description="EMVC-2 v1.0")

        # Add the arguments
        parser.add_argument("-i", "--bam_file", help="The bam file", required=True)
        parser.add_argument("-r", "--ref_file", help="The reference fasta file", required=True)
        parser.add_argument("-p", "--threads", help="The number of parallel threads (default 8)", default=8, type=int)
        parser.add_argument("-t", "--iterations", help="The number of EM iterations (default 5)", default=5, type=int)
        parser.add_argument("-m", "--learners", help="The number of learners (default 7)", default=7, type=int)
        parser.add_argument("-v", "--verbose", help="Make output verbose (default 0)", default=0, type=int)
        parser.add_argument("-o", "--out_file", help="The output file name", required=True)

        # Parse the arguments
        args = parser.parse_args()

        # Check if the bam file exists
        if not os.path.isfile(args.bam_file):
            raise FileNotFoundError("The bam file does not exist.")

        # Check if the reference file exists
        if not os.path.isfile(args.ref_file):
            raise FileNotFoundError("The reference file does not exist.")

        # Check if samtools is installed
        if not shutil.which("samtools"):
            raise FileNotFoundError("Samtools is not installed. Try running 'conda install -c bioconda samtools'.")

        # # Check if the output file already exists
        # if os.path.isfile(args.out_file):
        #     raise FileExistsError("The output file already exists.")

        # Start timer
        start = time.time()

        # Print the program name and version
        print("/*  EMVC-2 v1.0 */")
        print("Authors: Guillermo Dufort y Álvarez, Martí Xargay, Idoia Ochoa, and Alba Pages-Zamora")
        print("Contact: gdufort@fing.edu.uy")
        print("/*****************/ \n")

        # Print arguments
        print("Arguments:")
        print("BAM file: " + args.bam_file)
        print("Reference file: " + args.ref_file)
        print("Threads: " + str(args.threads))
        print("Output file: " + args.out_file + "\n")

        # Check if the bam file is indexed
        if not os.path.isfile(args.bam_file + ".bai"):
            # Index the bam file
            print("Indexing the bam file...")
            pysam.index(args.bam_file)

        # Create a temporary folder in the bam_file directory to store the intermediate files
        temp_folder = os.path.join(os.path.dirname(args.bam_file), "tmp")
        os.makedirs(temp_folder, exist_ok=True)

        # Sort the bam file if it is not sorted and get the sorted file name
        bam_file = sort_bam(args.bam_file)

        # Get the contig names from the header of the bam file using pysam.AlignmentFile.header.references and store them in a list
        bam = pysam.AlignmentFile(bam_file, "rb")
        contigs = list(bam.header.references)
        # Print all contig names in one line
        print("Step 1- Finding candidate variants with the EM algorithm for {} contigs:".format(len(contigs)), ", ".join(contigs))

        # Create a multiprocessing pool with the number of threads
        pool = multiprocessing.Pool(args.threads)

        # Create a list of tuples with the arguments for each variant_call function call
        args_list = [(contig, temp_folder, args.ref_file, bam_file, args.iterations, args.learners, args.verbose) for contig in contigs]

        # Run the variant_call function for each tuple of arguments in parallel and get the results as a list of bcf files
        # Wrap the result of imap_unordered with the tqdm object pool.imap_unordered
        vcf_files = list(tqdm.tqdm(pool.imap_unordered(variant_call, args_list), total=len(args_list), bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}"))

        # Close the pool and wait for the processes to finish
        pool.close()
        pool.join()

        model_name = "RF_depth5_test0.01_GT_DP_ALT_%_CLASS_Entropy_ERR174324_1_22_41"
        clf = pickle.load(open("dt_model/{}.pkl".format(model_name), 'rb'))
        features = ["GT", "DP", "ALT_%", "CLASS_Entropy"]

        args_list = [(vcf_file, clf, features, vcf_file[:-4] + "_RF.vcf") for vcf_file in vcf_files]

        # Create a multiprocessing pool with the number of threads
        pool = multiprocessing.Pool(args.threads)

        # Step 2: Filter variants using the trained model
        print("\nStep 2- Filtering variants with the trained decision tree model...")
        # Run the filter_variants function for each tuple of arguments in parallel and get the results as a list of tuples
        filtered_data = list(tqdm.tqdm(pool.imap_unordered(filter_variants, args_list), total=len(args_list), bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}"))
    
        # Close the pool and wait for the processes to finish
        pool.close()
        pool.join()
        # Total number of filtered variants
        total_filtered = 0
        total_kept = 0
        filtered_vcfs_files = []
        for (filtered_vcf, kept, filtered) in filtered_data:
            total_filtered += filtered
            total_kept += kept
            filtered_vcfs_files.append(filtered_vcf)

        new_order_vcfs = []
        # Sort the filterted vcf files in the same order as the original contigs
        for i in range(len(contigs)):
            for vcf_file in filtered_vcfs_files:
                if "/{}_RF.vcf".format(contigs[i]) in vcf_file:
                    new_order_vcfs.append(vcf_file)
                    break
        
        filtered_vcfs_files = new_order_vcfs
        
        print("Total variants: {}, filtered: {}".format(total_kept, total_filtered))

        # Join the vcf files into one vcf file
        print("\nJoining the filtered vcf files...")
        join_vcf_files(filtered_vcfs_files, args.out_file)

        # Delete the temporary folder and all its files
        shutil.rmtree(temp_folder, ignore_errors=True)

        # End timer
        end = time.time()

        # Print output file name
        print("\nOutput file: {}".format(args.out_file))
        # Print the total time in hh:mm:ss format
        print("Total time: {}".format(datetime.timedelta(seconds=end-start)))
    except Exception as e:
        # Print the error message
        print("Error: {}".format(str(e)))
        # Remove the temporary folder and all its files
        shutil.rmtree(temp_folder, ignore_errors=True)

if __name__ == "__main__":
    main()