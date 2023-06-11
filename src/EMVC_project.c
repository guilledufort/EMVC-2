#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "EMVC_functions_project.h"

int main(int argc, char* argv[]){

    // print_init();

    //The arguments introduced after the flags can be extracted from argv[]
    char *gt_file_name = argv[1];
    char *file_name = argv[2];
    // char *learners = argv[3];
    char *niterations = argv[4];
    char *results_file_name = argv[5];
    char *s_file_name = argv[6];
    char *nobv_file_name = argv[7];
    char *obv_var_file_name = argv[8];
    char *obv_invar_file_name = argv[9];
    char *truth_bed_file_name = argv[10];
    char *verbose = argv[11];

    // omp_set_num_threads(8);

    int T = atoi(niterations);
    int verbose_bool = atoi(verbose);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // ***PART 1 - MPILEUP READ*****************************************************

    /*
    VARIABLE DECLARATION AND INITIALIZATION

    This is where all variables for the 1st PART are declared and initialized if needed. Here is also where memory is initially allocated for the 
    dynamically updated memory units (with the calloc function). Also here is where all .txt files will be openned, if requested, for later analysis.
    
    Variable explanation:

        -   char *chr       Pointer to a chunk of memory storing the name of the chromosome we are currently reading from. When reading
                            from the pileup, this information is written in the first column. It will be necessary to generate the vcf file, where 
                            the chromosome of origin appears in the first column.

        -   int pos         Integer in charge of saving the position corresponding to each pileup read. It is written in the second column in pileup format
                            (outputing from the mpileup from samtools).

        -   char ref        Saves the reference nucleobase from each pileup line. It is written in the third column in pileup format

        -   int len         Saves the length of the sequence (and its corresponding Q-scores), that is, the number of reads that cover a specific 
                            position. This len is also described as read depth, and it is outputed in the final .vcf file. The number of mapping reads shall be one
                            of the criteria to evaluate the reliability of decided variants.

        -   int fwd_ref     Saves the number of reads mapping REF in "forward".

        -   int fwd_alt     Saves the number of reads mapping ALT in "forward".

        -   int rev_ref     Saves the number of reads mapping REF in "reverse".

        -   int rev_alt     Saves the number of reads mapping ALT in "reverse".

                *** Note that fwd_x shall always be smaller than len. 
                        fwd_ref + fwd_alt + rev_ref + rev_alt = len
                    Ideally, approximately half of all reads mapping REF should be "forward" strand, and the other half, "reverse" strand.
                    On the other hand, reads mapping alternative versions should have similar parity. This parity shall be one of the criteria to evaluate
                    the reliability of decided variants.

        -   char *seq       Pointer to a chunk of memory big enough to store the whole nucleobase sequence, therefore acting like a len-position-long string.

        -   char *qs        Pointer to a chunk of memory big enough to store the whole Q-Score sequence, therefore acting like a len-position-long string.

        -   char *line      Pointer that will point to the memory holding the line of the pileup file that is being read at the moment.

        -   char *token     Pointer that will point to the memory holding the token of the line that is being processed at the moment. A token is defined as every '\t'-separated
                            element of a pileup line.

        -   char *str_pos   Pointer that will point to the memory holding the position of the current pileup line, in a string format. It will later be converted to integer and
                            saved into 'pos'.

        -   char *str_ref   Pointer that will point to the memory holding the reference nucleobase of the current pileup line, in a string format. The corresponding char will
                            later be saved into 'ref'.

        -   char *str_len   Pointer that will point to the memory holding the length of the nucleobase sequence of the current pileup line, in a string format. It will later be
                            converted to integer and saved into 'pos'.

        -   char *ptr1  Pointer used to turn 'str_pos' into integer 'pos'.

        -   char *ptr2  Pointer used to turn 'str_len' into integer 'len'.

        -   FILE *obv_var_indx_file     File pointer used to store all relevant information regarding OBVIOUS VARIANT positions. It is saved in the
                                        file named: "obv_var_indx_file.txt".
                    It contains (in columns, tab-separated):
                        - 1st column: The position/index of such obvious variant line.
                        - 2nd column: The reference base, in upper or lower case letters.
                        - 3rd column: The length of the base sequence, that is, read depth.
                        - 4th column: The number of reads mapping REF that have been mapped forward.
                        - 5th column: The number of reads mapping ALT that have been mapped forward.
                        - 6th column: The number of reads mapping REF that have been mapped reverse.
                        - 7th column: The number of reads mapping ALT that have been mapped reverse.
                        - 8th column: The nucleobase sequence, in upper or lower case letters.

        -   FILE *obv_invar_indx_file   File pointer used to store all relevant information regarding OBVIOUS VARIANT positions. It is saved in the
                                        file named: "obv_invar_indx_file.txt". Note that this will be by far the longest (and "least interesting") 
                                        file created, so information is very limited. 
                    It contains (in columns, tab-separated):
                        - 1st column: The position/index of such obvious invariant line.
                        - 2nd column: The reference base, in upper or lower case letters.

        -   FILE *nobv_file             File pointer used to store all relevant information regarding NOT OBVIOUS positions. It is specially useful
                                        when checking all those positions entering the EM algorithm. It is saved in the file named: "nobv_indx_file.txt"
                    It contains (in columns, tab-separated): 
                        - 1st column: The position/index of such not obvious line.
                        - 2nd column: The reference base, in upper or lower case letters.
                        - 3rd column: The length of the base sequence, that is, read depth.
                        - 4th column: The number of reads mapping REF that have been mapped forward.
                        - 5th column: The number of reads mapping ALT that have been mapped forward.
                        - 6th column: The number of reads mapping REF that have been mapped reverse.
                        - 7th column: The number of reads mapping ALT that have been mapped reverse.
                        - 8th column: The nucleobase sequence, in upper case letters {A,C,T,G}.
                        - 9th column: The QS sequence, comma-separated, in Phred+33 scale.
            
            ***Note that even if these files are created, this information should be saved in memory, to avoid having to read from the file if needed for
               later use. To do so, we may use "*obv_var_indxs", "*obv_invar_indxs" and "*nobv_indxs".

        -   int *obv_var_indxs      Dynamically updated (that is, memory is extended when it gets full) 3-column matrix, storing information from OBVIOUS VARIANT positions.
                    It contains:
                        - 1st column: the positions/locus of such obvious variant.
                        - 2nd column: the corresponding reference nucleobase in integers: {A,C,T,G}={0,1,2,3}.
                        - 3rd column: the sequence length, that is, the read depth. Will come in handy for generating the .vcf file.
                        - 4th column: the number of reads mapping REF "forward". Will come in handy for generating the .vcf file.
                        - 5th column: the number of reads mapping ALT "forward". Will come in handy for generating the .vcf file.
                        - 6th column: the number of reads mapping REF "reverse". Will come in handy for generating the .vcf file.
                        - 7th column: the number of reads mapping ALT "reverse". Will come in handy for generating the .vcf file.                        
                        - 8th column: the decided class {0,1,2,3,4,5,6,7,8,9}. In this case, no algorithm shall run, since there is no conflict on the decided class.
                        - 9th column: the QS of the decided variant. In this case, since the algorithm does not run, an arbitrary QS=1000 is assigned to these positions
                                      (the probability of this decision being correct would be virtually alpha=1, and -10*log10(1-alpha) = infinite) 
    
        -   int *obv_invar_indxs    Dynamically updated (that is, memory is extended when it gets full) vector, storing the OBVIOUS INVARIANT positions. It is accessed with
                                    the offset_2d() function.      

        -   int *nobv_indxs     Dynamically updated (that is, memory is extended when it gets full) 3-column matrix, storing information from NOT OBVIOUS positions.
                                It is accessed with the offset_2d() function.
                    It contains:
                        - 1st column: the positions/locus of such obvious variant.
                        - 2nd column: the corresponding reference nucleobase in integers: {A,C,T,G}={0,1,2,3}.
                        - 3rd column: the sequence length, that is, the read depth. Will come in handy for generating the .vcf file.
                        - 4th column: the number of reads mapping REF "forward". Will come in handy for generating the .vcf file.
                        - 5th column: the number of reads mapping ALT "forward". Will come in handy for generating the .vcf file.
                        - 6th column: the number of reads mapping REF "reverse". Will come in handy for generating the .vcf file.
                        - 7th column: the number of reads mapping ALT "reverse". Will come in handy for generating the .vcf file. 
                        - 8th column: the decided class {0,1,2,3,4,5,6,7,8,9}. Will be filled once the EM algorithm has run
                        - 9th column: the QS of the decided variant, corresponding to the a posteriori probability (alpha) of the decided class. Will be filled once the EM 
                                      algorithm has run. It is calculated as -10*log10(1-alpha). Note that if alpha is so close to 1 that it is impossible for the compiler
                                      to distinguish it from an actual 1, the value 1000 is introduced in this column, similarly to how it is done in the obvious variants.

        -   int *s  Dynamically updated (that is, memory is extended when it gets full) 3D tensor saving the number of times a particular position
                    has been targeted with a determined nucelobase {A,C,T,G} by a concrete learner. Therefore, if the m'th learner provides the l'th base
                    {A,C,T,G}={0,1,2,3} for the n'th position in 'x' different reads, one may find the number 'x' in the m'th x-axis position, l'th y-axis 
                    position, n'th z-axis position: s[m,l,n]=x.

                    It is accessed with the offset_3d() function.

                    Note that the n'th position refers to the index number 'n' from the NOT OBVIOUS list (that is, every 'n' refers to a not obvious position),
                    not to the locus of this line in the whole chromosome. A not obvious position might have n=4, and locus/position=60237.

        -   int obv_var_cnt     Counts the number of detected OBVIOUS VARIANT positions, plus also serving as index to continue filling the
                                'obv_var_indcs' matrix.

        -   int obv_invar_cnt   Counts the number of detected OBVIOUS INVARIANT positions, plus also serving as index to continue filling the
                                'obv_invar_indxs' vector.

        -   int nobv_cnt        Counts the number of detected NOT OBVIOUS positions, plus also serving as index to continue filling the 'nobv_indxs' 
                                matrix.

        -   int disc_cnt        Counts the number of discarded positions. This is going to be either because the reference nucleobase of a particular position
                                had an 'N' or 'n', or the nucleobase sequence of a particular line of the pileup consisted only of 'N', 'n' or '*'.
                                It's going to be useful to check that the number of accounted positions corresponds to the number of read positions.

        -   int token_cnt       When reading each line of the pileup, this counter is used to keep track on the column a particular token is extracted from.

        -   int ret     Saves the value returned by the getline() function, which reads from stdin. If the end of file is reached, ret==-1.

        -   int i       Counts the number of lines from the pileup file that are read. Every DISPLAY_PERIOD positions, it is displayed on streen the current line 
                        that's being read. This "i" counter is specially useful to keep count of the current pileup line. 

        -   size_t length   Saves the length of the pileup line being read (including all tokens).

    */

    // ***VARIABLE DECLARATION AND INITIALIZATION***

    int pos=0, len=0, fwd_ref=0, fwd_alt=0, rev_ref=0, rev_alt=0, obv_invar_cnt=0, obv_var_cnt=0, nobv_cnt=0, disc_cnt=0, token_cnt=0, ret=1, i=0, *s=NULL, *obv_var_indxs=NULL, *obv_invar_indxs=NULL, *nobv_indxs=NULL;
    char ref=0, *chr=NULL, *seq=NULL, *qs=NULL, *line=NULL, *token=NULL, *str_pos=NULL, *str_ref=NULL, *str_len=NULL, *ptr1=NULL, *ptr2=NULL;
    FILE *obv_var_indx_file=NULL, *obv_invar_indx_file=NULL, *nobv_file=NULL, *truth_bed_file = NULL;
    size_t length = 0;
 

    /*
    Here, memory is initially allocated for the dynamically managed memory units: 's', 'nobv_indxs', 'obv_var_indxs' and 'obv_invar_indxs'.
    Note that in every memory unit, a multiple of P positions are allocated. That is because P is the arbitrarily assigned initial n-axis size
    (more information in the "int* update_matrix()" function description).
    */

    s = (int*) calloc(P*M*L, sizeof(int));
    nobv_indxs = (int*) calloc(P*9, sizeof(int));   //four columns. One for the index, another for the corresponding reference base, the next one for the read depth, one for the 
                                                    //number of forward mapped reads, one for the decision and the last for the QS of the decision
    obv_var_indxs = (int*) calloc(P*9, sizeof(int));    //four columns. One for the index, another for the corresponding reference base, the next one for the read depth, one for 
                                                        //the number of forward mapped reads, one for the obvious variant and the last one for QS of decision (which is always 1)
    obv_invar_indxs = (int*) calloc(P, sizeof(int)); //one column, saving the indexes

    //The files described above are created (to be written in) and assigned to the corresponding pointer, if so has been specified with the corresponding command flags
    if(strcmp(obv_var_file_name,"-")!=0){ //If the obv_var_file_name was specified in the command, create the corresponding file with the expressed name
        obv_var_indx_file = fopen(obv_var_file_name,"w");
    }
    if(strcmp(obv_invar_file_name,"-")!=0){ //If the obv_invar_file_name was specified in the command, create the corresponding file with the expressed name
        obv_invar_indx_file = fopen(obv_invar_file_name,"w");
    }
    if(strcmp(nobv_file_name,"-")!=0){ //If the nobv_file_name was specified in the command, create the corresponding file with the expressed name
        nobv_file = fopen(nobv_file_name,"w");
    }
    unsigned long * truth_bed = NULL;
    if(strcmp(truth_bed_file_name,"-")!=0){ //If the truth_bed_file_name was specified in the command, create the corresponding file with the expressed name
        // Count the number of lines in the truth bed file for the size of the array
        truth_bed_file = fopen(truth_bed_file_name,"r"); 
        // Counting TRUTH BED FILE lines
        unsigned long bed_regions = count_lines(truth_bed_file);
        fclose(truth_bed_file);

        truth_bed_file = fopen(truth_bed_file_name,"r"); 
        truth_bed = read_truth_bed_file(truth_bed_file, bed_regions);
        fclose(truth_bed_file);
    }

    // ***End of variable declaration and initialization***

    /*
    MPILEUP READ

    Since we already have the mpileup function from samtools library that does the pileup for us, it is now time to extract all relevant data
    from the mpileup output format (http://www.htslib.org/doc/samtools-mpileup.html for more information). In every line, 
        - The 1st column (tab-separated) corresponds to the chromosome name (saved in 'chr')
        - The 2nd column corresponds to the line position/locus in the chromosome (saved in 'len')
        - The 3rd column corresponds to the line reference nucleobase (saved in 'ref')
        - The 4th column corresponds to the number of reads mapping this particular position, that is, the nucleobase and qs sequence length (saved 
          in 'len')
        - The 5th column corresponds to the nucleobase sequence for this particular position (saved in 'seq')
        - The 6th column corresponds to the Q-Score sequence for this particular position (saved in 'qs')
    Note that after reading the chromosome name, the position, the reference base and the length, memory is allocated for the 'seq' and 'qs' strings
    according to the length specified before.

    The reading loop repeats indefinetely until the end of file (ret == -1) is reached, that is, the end of the pileup file is reached: as long as the
    return from the getline() is different from -1. The way how the the pileup file is read is very similar to how it is done in the read_ground_truth()
    function, but still worthy of explanation.
        A token_cnt is initialized to 0 in every new line. The strtok() function is used to return every piece of "line" that is '\t'-separated.
        This small strings, named tokens are individually analyzed. Until the end of "line"

    The dynamically updated memory units save information until full. Then, the memory is expanded to accomodate more information (for more information,
    check "int* update_matrix()"). Then, the matrix is trimmed from the unused allocated space (for more information, check "int* trim_matrix()").
    */

    // ***MPILEUP READ***

    unsigned int current_region = 0;

    while(ret!=-1){

        //READING THE FILE

        //Returns the next line, until the \n
        ret = getline(&line, &length, stdin);

        //If end of file is reached, continue. It will then check the looping condition, and since it will not fulfill it, the loop will stop
        if(ret==-1){
            continue;
        }

        token_cnt = 0;
        token = strtok(line, "\t");

        while (token != NULL){//As long as there are tokens to extract
            if(token_cnt>5){ //A token counter is used to check the column from which the token has been extracted.
                printf("\nError! Too many columns for a pileup...\n");
            }
            switch(token_cnt){ //Depending on the particular token we are extracting, we shall save (or not) the information in the corresponding variable
                case 0: chr=token; break;   //chromosome
                case 1: 
                    str_pos=token; 
                    pos = strtol(str_pos, &ptr1, 10); 
                    break;   //position
                case 2: 
                    str_ref=token; 
                    ref = str_ref[0];
                    break;   //reference nucleobase
                case 3: 
                    str_len=token; 
                    len = strtol(str_len, &ptr2, 10); 
                    break;   //read depth
                case 4: 
                    seq=token; 
                    break;   //nucleobase sequence
                case 5: qs=token; break;    //qs sequence
                default: printf("Error when reading. Should not be more than 5 tokens\n");
            }
            
            token = strtok(NULL, "\t"); //Extract the next token
            token_cnt++; //Update the token counter
        }

        if (truth_bed != NULL) {
            while (pos > truth_bed[current_region + 1] && current_region < BED_REGIONS - 1) {
                current_region += 2; 
            }
        }

        switch(decide_type(pos, ref, &len, &fwd_ref, &fwd_alt, &rev_ref, &rev_alt, seq, qs, truth_bed + current_region)){ //Deciding whether it is NOT OBVIOUS, OBVIOUS VARIANT or OBVIOUS INVARIANT
            case 1:
                add_to_nobv(pos, len, fwd_ref, fwd_alt, rev_ref, rev_alt, ref, seq, qs, nobv_file, s, nobv_indxs, &nobv_cnt, nobv_file_name);
                break;
            case 2:
                add_to_obv_var(pos, len, fwd_ref, fwd_alt, rev_ref, rev_alt, ref, seq, obv_var_indx_file, obv_var_indxs, &obv_var_cnt, obv_var_file_name);
                break;
            case 3:
                add_to_obv_invar(pos, ref, obv_invar_indx_file, obv_invar_indxs, &obv_invar_cnt, obv_invar_file_name);
                break;
            case 4:
                disc_cnt++;
                continue;
                break;
            default:
                exit(EXIT_FAILURE);
        }
        //Every DISPLAY_PERIOD positions, it is displayed on streen the current line that's being read. Useful to give some feedback on the user about the programs
        //stage of execution
        
        if(verbose_bool && i%DISPLAY_PERIOD==0){
            printf("Processing %d'th line of mpileup\n", i);
        }
        i++;

        s = update_matrix(s, nobv_cnt, M*L); //checking if the dynamically updated memory units should be expanded. More information in the "int* update_matrix()" description
        nobv_indxs = update_matrix(nobv_indxs, nobv_cnt, 9);
        obv_var_indxs = update_matrix(obv_var_indxs, obv_var_cnt, 9);
        obv_invar_indxs = update_matrix(obv_invar_indxs, obv_invar_cnt, 1);
    }
    
    if (verbose_bool)
        printf ("DISCARDED: %d \n" , disc_cnt);

    s = trim_matrix(s, nobv_cnt, M*L); //After all consecutive expantions throughout the reading of the pileup file, the dynamically updated memory might
                                       //have extra allocated space that should be properly trimmed. More information in the "int* trim_matrix()" description
    nobv_indxs = trim_matrix(nobv_indxs, nobv_cnt, 9);
    obv_var_indxs = trim_matrix(obv_var_indxs, obv_var_cnt, 9);
    obv_invar_indxs = trim_matrix(obv_invar_indxs, obv_invar_cnt, 1);

    print_s_file(s, nobv_cnt, s_file_name);

    if(strcmp(obv_var_file_name,"-")!=0){
        fclose(obv_var_indx_file);
    }
    if(strcmp(obv_invar_file_name,"-")!=0){
        fclose(obv_invar_indx_file);
    }
    if(strcmp(nobv_file_name,"-")!=0){
        fclose(nobv_file);
    }

    // ***End of mpileup reading and processing***

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // ***PART 2 - EM ALGORITHM RUN*************************************************

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    VARIABLE DECLARATION AND INITIALIZATION

    This is where all variables for the 2nd PART are declared and initialized if needed. Mainly, here is where key components of de described EM
    algorithm (see corresponding paper) are declared and initialized.

    Variable explanation:

        -   int N   Integer appearing several times throughout the code: it refers to the total number of NOT OBVIOUS positions, that is, the number
                    of positions actually entering the EM algorithm: that is why it is directly initialized with the 'nobv_cnt' value. Therefore, 
                    N is also the length of the 's' tensor and the 'nobv_indxs' matrix.

        -   double *alpha               pointer to a K x N matrix, storing the corresponding a posteriori probabilities, for each iteration. Note that, just like
                                        'pi' and 'gamma', it could have been designed with an extra dimension with T+1 positions, storing the different
                                        a posteriori probabilities calculated in each iteration. Its computation is part of the E-Step of the algorithm.
                                        For clarity, it could have been declared as double alpha[K][N]: however, since N integer is usually quite large, such
                                        significant memory allocations must be assigned with the malloc() family.

        -   double pi[K][T+1]            K x (T+1) matrix, which is better understood as a K positions long vector, storing its content over the T 
                                        iterations (T iterations + the initial state). This 'pi' vector stores the a priori probability for the K classes
                                        (for more information check the corresponding paper). Note that storing the evolution of the a priority probability
                                        over different iterations can be specially interesting when analysing the algorithm. Its computation is part of the
                                        M-Step of the algorithm.

        -   double gamma[L][K][M][T+1]   L x K x M x (T+1) tensor. It stores the confusion matrix (L x K) for each learner (x M) for each iteration (T iterations
                                        + the initial state). For more information over the meaning of these confusion matrixes, check the corresponding 
                                        paper). Note that storing the evolution over different iterations can be specially interesting when analysing the 
                                        algorithm. Its computation is part of the M-Step of the algorithm.
    */

    // ***VARIABLE DECLARATION AND INITIALIZATION***

    int N = nobv_cnt;
    double log_pi[K][T+1], log_gamma[L][K][M][T+1];
    double *alpha = (double*) calloc(K*N, sizeof(double));

    //Initialization heavily conditions the outcome of the algorithm: therefore, the initial values of pi and gamma have to be thoroughly considered
    initialize_pi(T, log_pi);
    initialize_gamma(T, log_gamma);

    //Just to print the initial state of the pi vector and the confusion matrixes of each learner
    if (verbose_bool){
        printf("Initial state of pi vector:\n");
        print_pi(T, log_pi, 0);
        printf("Initial state of gamma tensor:\n");
        print_gamma(T, log_gamma, 0);
    }

    // ***End of variable declaration and initialization***


    /*
    ALGORITHM RUN

    The function "void run_EMVC()" computes the EM algorithm, and the corresponding decisions {AA,CC,GG,TT,AC,AG,AT,CG,CT,GT} = {0,1,2,3,4,5,6,7,8,9}
    are stored in the third column of the 'nobv_indxs' vector.
    */

    // ***EM ALGORITHM RUN***

    run_EMVC(N, T, alpha, log_pi, log_gamma, s, nobv_indxs);

    // ***End of EM algorithm run***


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ***PART 3 - RESULTS ASSESSMENT***********************************************

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    VARIABLE DECLARATION AND INITIALIZATION

    Here is where variables are declared and initialized for the 3rd part. Also, the ground truth and reference files are opened in reading mode.
    Finally, memory is allocated to save the newly obtained NOT OBVIOUS VARIANT and NOT OBVIOUS INVARIANT information.    

    Variable explanation:

        -   int tp  Integer counting the number of "true positives" detected in our OBVIOUS VARIANT and NOT OBVIOUS VARIANT, when comparing the results
                    to the ground_truth: that is, correctly detected variants.

        -   int fp  Integer counting the number of "false positives" detected in our OBVIOUS VARIANT and NOT OBVIOUS VARIANT, when comparing the results
                    to the ground_truth: that is, incorrectly detected variants.

        -   int fn  Integer counting the number of "false negatives" detected in our OBVIOUS INVARIANT and NOT OBVIOUS INVARIANT, when comparing the 
                    results to the ground_truth: that is, not detected variants.

                ***Note that tn (which stands for "true negatives") is not defined. That's because it's a big number, which contributes little information,
                   and can easily be calculated knowing that:
                        tn = obv_var_cnt + obv_invar_cnt + nobv_var_cnt + nobv_invar_cnt - tp - fp - fn 


        -   int nobv_var_cnt    Integer counting number of NOT OBVIOUS VARIANTS decided through our algorithm, when comparing to the reference file.
                                It also serves as the index to keep filling the 'nobv_var_indxs' matrix.

        -   int nobv_invar_cnt  Integer counting number of NOT OBVIOUS INVARIANTS decided through our algorithm, when comparing to the reference file.
                                It also serves as the index to keep filling the 'nobv_invar_indxs' matrix.

                ***Note that nobv_var_cnt + nobv_invar_cnt shall be equal to nobv_cnt. Nontheless, separate adding is preferred to avoid miscounting


        -   int *nobv_var_indxs     2-column matrix of nobv_var_cnt lines. It stores the NOT OBVIOUS VARIANT positions/locus (first column) and their 
                                    corresponding decision {0,1,2,3,4,5,6,7,8,9} (second column).

        -   int *nobv_invar_indxs   2-column matrix of nobv_var_cnt lines. It stores the NOT OBVIOUS INVARIANT positions/locus (first column) and their 
                                    corresponding decision {0,1,2,3,4,5,6,7,8,9} (second column). 

                ***Note that: 
                    Unlike what was done with similar memory units (such as 'obv_var_indxs'), where, unknowing the final length of the matrix, a 
                    dynamic approach was followed, we can now initialize 'nobv_var_indxs' and 'nobv_invar_indxs' with N lines. Since positions
                    appearing in 'nobv_indxs' (which is N lines long) will go to either 'nobv_var_indxs' or 'nobv_invar_indxs', neither of them
                    will exceed the threshold of N lines. This is a much more comfortable technique, allowing to initially allocate the corresponding
                    memory and then trimming (with "int* trim_matrix()") the part that has been left unused. 

                    Also note that, as done previously with similar information, this tables are not printed onto a .txt file: that's
                    because it shall be represented in the .vcf file, with the corresponding format.


        -   double prec  double saving the calculation of the precision = tp/(tp+fp).

        -   double sens  double saving the calculation of the sensitivity = tp/(tp+fn).

        -   FILE *gt_file   File pointer used to open the ground truth (.vcf file) of the corresponding chromosome.

        -   FILE *results_file  File pointer to a prospective file containing all rellevant algorithm results. Note that this file will only be completed as long as
                                the corresponding flag has been activated (-v). In that case, the file will be named as the user has specified (the information it
                                will save is specified in the README.md)
    */

    // ***VARIABLE DECLARATION AND INITIALIZATION***

    int tp=0, fp=0, fn=0, nobv_var_cnt=0, nobv_invar_cnt=0, *nobv_var_indxs, *nobv_invar_indxs;
    double prec=0, sens=0;
    FILE *gt_file=NULL, *results_file=NULL;
    
    //Ground truth and reference files are assigned to the corresonding pointers
    if(strcmp(gt_file_name,"-")!=0){
        gt_file = fopen(gt_file_name,"r");   
    }
    if(strcmp(results_file_name,"-")!=0){
        results_file = fopen(results_file_name, "w");
    }

    //Memory for the NOT OBVIOUS VARIANT and NOT OBVIOUS INVARIANT positions is allocated
    nobv_var_indxs = (int*) calloc(N*10, sizeof(int));
    nobv_invar_indxs = (int*) calloc(N*2, sizeof(int));

    // ***End of variable declaration and initialization***


    /*
    RESULTS ASSESSMENT
    Compare the results with reference and with ground truth. 

    Firstly, the decisions made by the EM algorithm are compared to the reference file, to check whether the decisions correspond to actual variants
    or no (to know whether the NOT OBVIOUS are actually NOT OBVIOUS VARIANT or NOT OBVIOUS INVARIANT). Luckily, we had previously saved both the reference
    nucleobase and the corresponding decisions in the 'nobv_indxs' memory unit, together with the position/locus of such NOT OBVIOUS line. This is done 
    through the "void compare_to_ref()" function. This function also fills the 'nobv_var_indxs' and 'nobv_invar_indxs' memory units and updates their 
    corresponding counters. Note that unused memory space from this memory units is properly deleted through the "int* trim_matrix()" function.

    Secondly, the function "void compare_to_ground_truth()" compares the overall results to ground truth, to extract the exact number of true positives,
    false positives and false negatives obtained. For more information, check the "void compare_to_ground_truth()" description.
    */

    // ***RESULTS ASSESSMENT***

    compare_to_ref(N, nobv_var_indxs, nobv_invar_indxs, &nobv_var_cnt, &nobv_invar_cnt, nobv_indxs);
    nobv_var_indxs = trim_matrix(nobv_var_indxs, nobv_var_cnt, 10);
    nobv_invar_indxs = trim_matrix(nobv_invar_indxs, nobv_invar_cnt, 2);
    free(nobv_indxs);   

    if(strcmp(gt_file_name,"-")!=0){
        compare_to_ground_truth(&tp, &fp, &fn, obv_var_cnt, nobv_var_cnt, obv_var_indxs, nobv_var_indxs, gt_file);
        prec = compute_precision(tp, fp);
        sens = compute_sensitivity(tp, fn);
    }

    //Print relevant information
    if (verbose_bool)
        print_info(tp, fp, fn, obv_var_cnt, obv_invar_cnt, nobv_cnt, nobv_var_cnt, nobv_invar_cnt, disc_cnt, prec, sens, gt_file_name, results_file_name, results_file);
    
    if(strcmp(gt_file_name,"-")!=0){
        fclose(gt_file);   
    }
    if(strcmp(results_file_name,"-")!=0){
        fclose(results_file);
    }

    // ***End of results assessment***

    /*
    VCF FILE GENERATION: The .vcf file is generated accordingly, following the standard https://samtools.github.io/hts-specs/VCFv4.2.pdf
    */

    // ***VCF FILE GENERATION***
    generate_vcf(chr, obv_var_cnt, nobv_var_cnt, obv_var_indxs, nobv_var_indxs, alpha, file_name);
    // ***End of VCF file generation***

    return 0;
}
