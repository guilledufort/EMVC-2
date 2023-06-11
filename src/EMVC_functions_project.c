#define _GNU_SOURCE
#include "EMVC_functions_project.h"

unsigned int count_lines(FILE* file_ptr) {
    int line_count = 0;
    char* line = NULL;
    size_t length = 0;
    while (getline(&line, &length, file_ptr) != -1 ) {
        ++line_count;
    }

    return line_count;
}

unsigned long *read_truth_bed_file(FILE *truth_bed_file, unsigned long bed_regions) {
    unsigned long *res = (unsigned long *)calloc(bed_regions * 2, sizeof(unsigned long));
    unsigned int i = 0;
    int token_cnt, ret = 1;
    char *line = NULL, *token, *ptr1, *ptr2;
    ;
    size_t length = 0;

    while (ret != -1) {
        // READING THE FILE

        // Returns the next line, until the \n
        ret = getline(&line, &length, truth_bed_file);

        // If end of file is reached, continue. It will then check the looping condition, and since it will not fulfill it, the loop will stop
        if (ret == -1) {
            continue;
        }

        token_cnt = 0;
        token = strtok(line, "\t");

        while (token != NULL) {   // As long as there are tokens to extract
            if (token_cnt > 3) {  // A token counter is used to check the column from which the token has been extracted.
                printf("\nError! Too many columns for a truth bed file...\n");
            }
            switch (token_cnt) {  // Depending on the particular token we are extracting, we shall save (or not) the information in the corresponding variable
                case 0:
                    break;  // chromosome
                case 1:
                    res[i++] = strtol(token, &ptr1, 10);
                    break;  // initial region position
                case 2:
                    res[i++] = strtol(token, &ptr2, 10);
                    break;  // final region position
                default:
                    printf("Error when reading. Should not be more than 3 tokens\n");
            }

            token = strtok(NULL, "\t");  // Extract the next token
            token_cnt++;                 // Update the token counter
        }
    }
    return res;
}

int decide_type(int pos, char ref, int *len, int *fwd_ref, int *fwd_alt, int *rev_ref, int *rev_alt, char *seq, char *qs, unsigned long *truth_bed) {
    /*
    It is used in order to classify the mpileup line considering the contents of the seq and qs strings.
    Note that here is were mpileup lines are discarded (outputing 0 and 4) if they are not useable.

    According to the mpileup format, the '.' and ',' characters in the sequence string represent reads whose nucleobase in this position does not
    differ from the reference.

    The function:
        Returns 1 if one or more elements from the input 'seq' differ from the reference nucleobase. Therefore, it will be
        checked if more than one character from 'seq' are different to '.' or ','. Only at these loci will the EM algorithm apply.
        They will be accounted as NOT OBVIOUS. The algorithm shall apply to this position.

        Returns 2 if all of the elements from the input 'seq' differ from the reference nucleobase and are equal with one another.
        They will be accounted as OBVIOUS VARIANT. The algorithm shall not apply to this position.

        Returns 3 if no elements from the input 'seq' differ from the reference nucleobase. They will be accounted
        as OBVIOUS INVARIANT. The algorithm shall not apply to this position.

        Returns 4 if the mpileup line is unuseable and should be discarded. This may happen when the reference nucleobase is 'n' or 'N', or when
        the sequence only contained 'n', 'N' or '*', being therefore unusable. Therefor, the algorithm shall not apply to this position.
    */
    int count = 0, flag = 0;  // Counts the number of differing reads from the reference base in this particular mpileup line
    char first_base;
    if (ref == 'N' || ref == 'n') {  // The pileup offers bases that the reference did not have
        return 4;
    }

    clean_seq(seq, qs, len, fwd_ref, fwd_alt, rev_ref, rev_alt);  // In order to extract those characters {N, n, *, $} and segments (i.e. +5ACCTG, ^K) that should not appear
                                                                  // in the 'seq' and 'qs' strings. Plus, 'len' is updated

    if (*len == 0) {  // If the sequences only contained {N,n,*,$}, the sequences might have become of null length. Therefore, they shall be discarded
        return 4;
    }

    if (truth_bed != NULL) {
        if (pos < truth_bed[0] || pos > truth_bed[1])
            return 4;
    }

    // To decide whether it's obvious invariant, obvious variant or not obvious. "count" counts the number of differing reads
    first_base = seq[0];
    for (int i = 0; i < *len; i++) {
        if (seq[i] != '.' && seq[i] != ',') {  // In case there is a differing read, the counter is increased
            count++;
            if ((seq[i] - first_base) != 0 && abs(seq[i] - first_base) != 32) {  // Check independently of upper or lower case
                flag = 1;
            }
        }
    }
    if (count == 0) {                         // There is no disagreeing read
        return 3;                             // obv_invar - OBVIOUS INVARIANT
    } else if (count == *len && flag == 0) {  // All the reads are disagreeing, and equal with one another
        return 2;                             // obv_var - OBVIOUS VARIANT
    } else {                                  // Some reads (not all) are disagreeing
        return 1;                             // nobv - NOT OBVIOUS
    }
}

void clean_seq(char *seq, char *qs, int *len, int *fwd_ref, int *fwd_alt, int *rev_ref, int *rev_alt) {
    /*
    Retrieves '*', '$', 'N' and 'n' characters from the nucleotide sequence, with their corresponding quality score.
    The length is accordingly updated, by substracting the number of {*,N,n} detected.

    Other string segments appearing in 'seq' shall also be retrieved. That is for example "^K" or indels such as "+7AGTAAGT".

    Note that since 'seq' and 'qs' are pointers to a string (being a string an array of numbers ending in '\0'), the easiest way to do this task is by
    sliding all the working characters to the front, and then just assigning '\0' after the last working character. In this way, the string will have
    been effectively shortened, thus avoiding the creation of new pointers and strings.

    It returns the number of forward reads mapping REF (fwd_ref), the number of forward reads mapping an alternate (fwd_alt), the number of reverse
    reads mapping REF (rev_ref) and the number of reverse reads mapping an alternate (rev_alt). It will be useful later to compute the strand bias (SB)
    */
    long i = 0, j = 0, k = 0;
    int indel_length = 0, num_length = 0;
    *fwd_ref = 0;
    *fwd_alt = 0;
    *rev_ref = 0;
    *rev_alt = 0;
    while (i < strlen(seq)) {
        if (seq[i] == '$') {  // Retrieve '$', which wasn't already counted in 'len'
            i++;
            continue;
        } else if (seq[i] == '*' || seq[i] == 'N' || seq[i] == 'n') {
            i++;
            j++;
            continue;
        } else if (seq[i] == '^') {  // Retrieve the "^K" elements
            i = i + 2;
            continue;
        } else if (seq[i] == '+' || seq[i] == '-') {  // Retrieve indels
            num_length = 0;
            indel_length = 0;
            while ((int)(seq[i + 1 + num_length]) >= 48 && (int)(seq[i + 1 + num_length]) <= 57) {  // While it is a number
                num_length++;
            }
            for (int x = 0; x < num_length; x++) {
                indel_length = indel_length + (seq[i + 1 + x] - '0') * pow(10, num_length - x - 1);  // The length of the indel is calculated by turning into an integer the characters that
                                                                                                     // formed the number after '+' and '-'
            }
            i = i + 1 + num_length + indel_length;
            continue;
        } else if (seq[i] == '.') {  // mapping REF in forward
            *fwd_ref = *fwd_ref + 1;
        } else if (seq[i] == ',') {  // mapping REF in reverse
            *rev_ref = *rev_ref + 1;
        } else if (seq[i] == 'A' || seq[i] == 'C' || seq[i] == 'G' || seq[i] == 'T') {  // mapping an alternate in forward
            *fwd_alt = *fwd_alt + 1;
        } else if (seq[i] == 'a' || seq[i] == 'c' || seq[i] == 'g' || seq[i] == 't') {  // mapping an alternate in reverse
            *rev_alt = *rev_alt + 1;
        }
        seq[k] = seq[i];
        qs[k] = qs[j];
        k++;
        i++;
        j++;
    }
    seq[k] = '\0';  // shortening the strings
    qs[k] = '\0';   // shortening the strings
    *len = k;       // updating 'len'
}

void add_to_nobv(int pos, int len, int fwd_ref, int fwd_alt, int rev_ref, int rev_alt, char ref, char *seq, char *qs, FILE *nobv_file, int *s, int *nobv_indxs, int *nobv_cnt, char *nobv_file_name) {
    /*
    Since the positions (and their corresponding sequences and quality scores) have been classified as NOT OBVIOUS, the aim of this function is to:

        - Make the the nucleobase and quality score sequences strings easily understandable, by turning '.' and ',' characters into their corresponding
        upper case reference nucleobase letter {A,T,G,C} and turning the quality scores previously represented in ASCII into the Phred+33 standard.

        - Effectively add the position to the EM algorithm.

        - Add the all information of the base into a 9-column matrix, 'nobv_indxs'. It shall be later used to analize the decisions taken with the
        algorithm with the reference information.

        - Update the NOT OBVIOUS counter, 'nobv_cnt'.

        - In case of interest, the line shall also be written down in the corresponding 'nobv_file'.
    */
    int l = 0, m = 0, q = 0, ref_int = 4;
    ;
    for (int i = 0; i < len; i++) {
        if (seq[i] == '.' || seq[i] == ',') {  // If they are the same as the reference, assign them the reference
            seq[i] = ref;
        }
        if (seq[i] >= 97 && seq[i] <= 122) {  // If they are in lower case, turn into upper case
            seq[i] = seq[i] - 32;
        }

        q = ((int)qs[i]) - 33;  // Turn into the Phred+33 scale

        decide_m_l(&m, &l, pos, seq[i], q);  // Considering the nucleobase value and the corresponding quality score, the learner (m) and base (l) can be assigned
        s[offset_3d(*nobv_cnt, m, l)]++;     // Position of this read effectively added to the algorithm once appering in the s tensor
    }

    // Add to indxs vectors. For practical use and compatibility with other data types, ref shall be saved as an int {0,1,2,3}
    ref_int = nuc_to_num(ref);

    nobv_indxs[offset_2d(*nobv_cnt, 0, 9)] = pos;      // first column. Position.
    nobv_indxs[offset_2d(*nobv_cnt, 1, 9)] = ref_int;  // second column. Reference nucleobase in integer format.
    nobv_indxs[offset_2d(*nobv_cnt, 2, 9)] = len;      // third column. Length of the sequence.
    nobv_indxs[offset_2d(*nobv_cnt, 3, 9)] = fwd_ref;  // fourth column. Number of reads mapping REF in forward.
    nobv_indxs[offset_2d(*nobv_cnt, 4, 9)] = fwd_alt;  // fifth column. Number of reads mapping ALT in forward.
    nobv_indxs[offset_2d(*nobv_cnt, 5, 9)] = rev_ref;  // sixth column. Number of reads mapping REF in reverse.
    nobv_indxs[offset_2d(*nobv_cnt, 6, 9)] = rev_alt;  // seventh column. Number of reads mapping ALT in reverse.
                                                       // eigth column will remain not filled until the decisions have been obtained.
                                                       // ninth column will remain not filled until the QS of the probabilities are obtained.
    *nobv_cnt = *nobv_cnt + 1;                         // The number of NOT OBVIOUS positions is updated. Note that it is used as index for loading the nobv_indxs.

    // Write in the nobv_file
    if (strcmp(nobv_file_name, "-") != 0) {
        fprintf(nobv_file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%s\t", pos, ref, len, fwd_ref, fwd_alt, rev_ref, rev_alt, seq);
        for (int i = 0; i < len - 1; i++) {
            fprintf(nobv_file, "%d,", ((int)qs[i] - 33));
        }
        fprintf(nobv_file, "%d\n", ((int)qs[len - 1]) - 33);
    }
}

void add_to_obv_var(int pos, int len, int fwd_ref, int fwd_alt, int rev_ref, int rev_alt, char ref, char *seq, FILE *obv_var_indx_file, int *obv_var_indxs, int *obv_var_cnt, char *obv_var_file_name) {
    /*
    Since the positions (and their corresponding sequences and quality scores) have been classified as OBVIOUS VARIANT, the aim of this function is to:

        - Add all the information referring to this position into a 9-column matrix, 'obv_var_indxs'. It shall be used to store the obvious variant results.

        - Update the OBVIOUS VARIANT counter, 'obv_var_cnt'.

        - In case of interest, the line shall also be written down in the corresponding 'obv_var_indx_file'.
    */
    int nuc_int = 4;
    int ref_int = 4;

    // Add to indxs vectors. For practical use and compatibility with other data types, 'ref' shall be saved as an int {0,1,2,3}
    nuc_int = nuc_to_num(seq[0]);
    ref_int = nuc_to_num(ref);
    obv_var_indxs[offset_2d(*obv_var_cnt, 0, 9)] = pos;      // first column. Position.
    obv_var_indxs[offset_2d(*obv_var_cnt, 1, 9)] = ref_int;  // second column. Reference nucleobase in integer format.
    obv_var_indxs[offset_2d(*obv_var_cnt, 2, 9)] = len;      // third column. Length of the sequence.
    obv_var_indxs[offset_2d(*obv_var_cnt, 3, 9)] = fwd_ref;  // fourth column. Number of reads mapping REF in forward.
    obv_var_indxs[offset_2d(*obv_var_cnt, 4, 9)] = fwd_alt;  // fifth column. Number of reads mapping ALT in forward.
    obv_var_indxs[offset_2d(*obv_var_cnt, 5, 9)] = rev_ref;  // sixth column. Number of reads mapping REF in reverse.
    obv_var_indxs[offset_2d(*obv_var_cnt, 6, 9)] = rev_alt;  // seventh column. Number of reads mapping ALT in reverse.
    obv_var_indxs[offset_2d(*obv_var_cnt, 7, 9)] = nuc_int;  // eighth column. The class decision. No conflict.
    obv_var_indxs[offset_2d(*obv_var_cnt, 8, 9)] = 1000;     // ninth column. The decision QS. It has been assigned '1000' arbitrarily, since there has been no probability computation.

    *obv_var_cnt = *obv_var_cnt + 1;  // The number of OBVIOUS VARIANT positions is updated. Note that it is used as index for loading the obv_var_indxs

    // Write in the 'obv_var_indx_file'
    if (strcmp(obv_var_file_name, "-") != 0) {
        fprintf(obv_var_indx_file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%s\n", pos, ref, len, fwd_ref, fwd_alt, rev_ref, rev_alt, seq);
    }
}

void add_to_obv_invar(int pos, char ref, FILE *obv_invar_indx_file, int *obv_invar_indxs, int *obv_invar_cnt, char *obv_invar_file_name) {
    /*
    Since the positions (and their corresponding sequences and quality scores) have been classified as OBVIOUS INVARIANT, the aim of this function is to:

        - Add the position to a vector, 'obv_invar_indxs'. It shall be later used to analize the decisions taken with the
        algorithm with the reference information.

        - Update the OBVIOUS VARIANT counter, 'obv_var_cnt'.

        - In case of interest, the line shall also be written down in the corresponding listing 'obv_invar_indx_file'.
    */

    // Add to indxs vectors
    obv_invar_indxs[*obv_invar_cnt] = pos;

    *obv_invar_cnt = *obv_invar_cnt + 1;  // The number of OBVIOUS INVARIANT positions is updated. Note that it is used as index for loading the obv_invar_indxs

    // Write in the obv_invar_indx_file
    if (strcmp(obv_invar_file_name, "-") != 0) {
        fprintf(obv_invar_indx_file, "%d\t%c\n", pos, ref);
    }
}

void decide_m_l(int *m, int *l, int pos, char nuc, int qs) {
    /*
    Considering that 's' can be expressed as a tensor, this function is used to determine the corresponding m- and l-axis position
        - The read's nucleobase contribution to this position will be translated from {A,C,T,G} into l={0,1,2,3}
        - The read's quality score contribution to this position will be translated into the corresponding learner: from {0...40...}
        into m={0,1,2,3,4,5,6}
    */

    *l = nuc_to_num(nuc);

    if (qs < 0) {
        printf("QS smaller than 0 (qs = %c) at %d\n", qs, pos);
    } else if (qs >= 0 && qs < 10) {
        *m = 0;
    } else if (qs >= 10 && qs < 20) {
        *m = 1;
    } else if (qs >= 20 && qs < 25) {
        *m = 2;
    } else if (qs >= 25 && qs < 30) {
        *m = 3;
    } else if (qs >= 30 && qs < 35) {
        *m = 4;
    } else if (qs >= 35 && qs < 40) {
        *m = 5;
    } else {  // qs>=40
        *m = 6;
    }
}

int *update_matrix(int *matrix, int cnt, int size) {
    /*
    Since the length of the incoming mpileup is unknown, the memory has to be managed dynamically, by updating each memory chunk when it gets full.
    This is the case of the 's' tensor, the 'nobv_indxs' and 'obv_var_indxs' matrixes and the 'obv_invar_indxs' vector.

    The procedure to fill this memory units dynamically is to start by saving a P*x*sizeof(int) memory space, where x=M*L if it is the 's' tensor,
    x=9 if it is the 'nobv_indxs' and 'obv_var_indxs' matrixes and x=1 for 'obv_invar_indxs' vector. In every mpileup line read, the function
    'update_matrix' checks whether the memory unit is full: if it's not, it continues until it fills. If it is, then the contents of this memory chunk are
    reallocated to a new and bigger space, with P*x*size_of(int) extra free space. This way, every time that it gets full, the memory is expanded by
    P*x*size_of(int). Therefore, considering that the size of all created memory units will be multiple from P*x*sizeof(int), it is also easier to
    check whether the memory unit is full: one just has to check whether the counter (in this case, cnt+1, since array positions in C start at 0)
    is divisible by P.

    Note that in every expansion, the newly extra allocated space has no assigned data, not even 0. That will not be a problem for 'nobv_indxs',
    'obv_var_indxs' and 'obv_invar_indxs', since their only goal is to save information, and no space will be left blank. However, every (n,m,l)
    position of the 's' will be called when de EM algorithm is running. Since these (n,m,l) positions from de 's' matrix will only be filled if, for the
    n'th position, the m'th learner targets the l'th base at least once, many of the remaining (n',m',l') positions will be left blank, and it is
    important that this "blank" corresponds to a 0 (meaning that for the n''th position, the m''th learner did not target the l''th base). Therefore,
    when we are expanding the 's' tensor, the newly allocated memory is filled with 0s.
    */
    if ((cnt + 1) % P == 0) {
        matrix = (int *)realloc(matrix, (cnt + 1 + P) * size * sizeof(int));
        if (matrix == NULL) {
            printf("\nError! memory not extended.\n");
            exit(EXIT_FAILURE);
        }
        if (size == M * L) {                                             // If we are treating an 's' tensor
            for (int i = (cnt)*M * L; i < (cnt + 1 + P) * M * L; i++) {  // For the newly allocated space
                matrix[i] = 0;                                           // Fill with 0s
            }
        }
    }
    return matrix;
}

int *trim_matrix(int *matrix, int cnt, int size) {
    /*
    Using the procedure described before, the memory unit has always less than P*x*sizeof(int) free space. After the end of the mpileup output is
    reached, the 's', 'nobv_indxs', 'obv_var_indxs' and 'obv_invar_indxs' will have inevitably expanded to acommodate the incoming data, having
    in the end less than P*x*sizeof(int) of free space. Therefore, this memory units have to be trimmed in order to be completely full, thus making
    their corresponding counter the size indicator of the matrix.

    To do so, this memory units are one last time reallocated to perfectly fit their size.
    */
    matrix = (int *)realloc(matrix, (cnt + 1) * size * sizeof(int));
    if (matrix == NULL) {
        printf("\nError! memory not extended.\n");
        exit(EXIT_FAILURE);
    }
    return matrix;
}

void initialize_pi(int T, double log_pi[K][T + 1]) {
    /*
    Pi initialization: it has been seen to output acceptable outcomes:
          AA     CC      TT      GG      AC      AT      AG      CT      CG      TG
    pi = [0.205  0.205   0.205   0.205   0.03    0.03    0.03    0.03    0.03    0.03]

    Note that the values from different iterations are also saved, for later study if necessary.
    */
    for (int t = 0; t < T; t++) {
        for (int k = 0; k < K; k++) {
            if (k <= 3) {
                log_pi[k][t] = LOG(0.205);
            } else {
                log_pi[k][t] = LOG(0.03);
            }
        }
    }
}

void initialize_gamma(int T, double log_gamma[L][K][M][T + 1]) {
    /*
    Gamma initialization: it has been seen to output acceptable outcomes:

     _  AA    CC    TT    GG    AC    AT    AG    CT    CG    TG   _
    |   0.4   0.2   0.2   0.2   0.3   0.3   0.3   0.2   0.2   0.2   | A
    |   0.2   0.4   0.2   0.2   0.3   0.2   0.2   0.3   0.3   0.2   | C
    |   0.2   0.2   0.4   0.2   0.2   0.3   0.2   0.3   0.2   0.3   | T
    |_  0.2   0.2   0.2   0.4   0.2   0.2   0.3   0.2   0.3   0.3  _| G

    Note that the values from different iterations are also saved, for later study if necessary.
    */
    for (int t = 0; t <= T; t++) {
        for (int m = 0; m < M; m++) {
            log_gamma[0][0][m][t] = LOG(0.4);
            log_gamma[0][1][m][t] = LOG(0.2);
            log_gamma[0][2][m][t] = LOG(0.2);
            log_gamma[0][3][m][t] = LOG(0.2);
            log_gamma[0][4][m][t] = LOG(0.3);
            log_gamma[0][5][m][t] = LOG(0.3);
            log_gamma[0][6][m][t] = LOG(0.3);
            log_gamma[0][7][m][t] = LOG(0.2);
            log_gamma[0][8][m][t] = LOG(0.2);
            log_gamma[0][9][m][t] = LOG(0.2);

            log_gamma[1][0][m][t] = LOG(0.2);
            log_gamma[1][1][m][t] = LOG(0.4);
            log_gamma[1][2][m][t] = LOG(0.2);
            log_gamma[1][3][m][t] = LOG(0.2);
            log_gamma[1][4][m][t] = LOG(0.3);
            log_gamma[1][5][m][t] = LOG(0.2);
            log_gamma[1][6][m][t] = LOG(0.2);
            log_gamma[1][7][m][t] = LOG(0.3);
            log_gamma[1][8][m][t] = LOG(0.3);
            log_gamma[1][9][m][t] = LOG(0.2);

            log_gamma[2][0][m][t] = LOG(0.2);
            log_gamma[2][1][m][t] = LOG(0.2);
            log_gamma[2][2][m][t] = LOG(0.4);
            log_gamma[2][3][m][t] = LOG(0.2);
            log_gamma[2][4][m][t] = LOG(0.2);
            log_gamma[2][5][m][t] = LOG(0.3);
            log_gamma[2][6][m][t] = LOG(0.2);
            log_gamma[2][7][m][t] = LOG(0.3);
            log_gamma[2][8][m][t] = LOG(0.2);
            log_gamma[2][9][m][t] = LOG(0.3);

            log_gamma[3][0][m][t] = LOG(0.2);
            log_gamma[3][1][m][t] = LOG(0.2);
            log_gamma[3][2][m][t] = LOG(0.2);
            log_gamma[3][3][m][t] = LOG(0.4);
            log_gamma[3][4][m][t] = LOG(0.2);
            log_gamma[3][5][m][t] = LOG(0.2);
            log_gamma[3][6][m][t] = LOG(0.3);
            log_gamma[3][7][m][t] = LOG(0.2);
            log_gamma[3][8][m][t] = LOG(0.3);
            log_gamma[3][9][m][t] = LOG(0.3);
        }
    }
}

void compute_alpha(int t, int n, double *alpha, double *alpha_sum, int T, double log_pi[K][T + 1], double log_gamma[L][K][M][T + 1], int *s);
void compute_gamma(int m, int t, int k, int T, int N, double log_gamma[L][K][M][T + 1], double *alpha, int *s);

void run_EMVC(int N, int T, double *alpha, double log_pi[K][T + 1], double log_gamma[L][K][M][T + 1], int *s, int *nobv_indxs) {
    /*
    Core of the EM algorithm. It follows the exact same formulation and calculations as in the corresponding paper.

    The algorithm has T iterations. For each iteration, the calculations can be divided in 2 parts:
        - E-step:   For each position (n), for each class (k), the corresponding alpha (a posteriori probability) is computed.
        - M-step:   For each learner (m), the corresponding gamma (confusion matrix) is computed.
                    For each class (k), the corresponding pi (a priori probability) is computed.

    Finally, after convergence of the results (after T iterations), a MAP decision is applied to decide, for each position (n), the class (k) that
    maximizes the a posteriori probability (alpha).

    In the end, the decided class and its corresponding alpha (expressed as -10*log10(1-alpha)) are saved in the 8th and 9th column of 'nobv_indxs'.

    For more information on this algorithm, and the calculations taking place in it, check the corresponding paper.

    Variable explanation:

        -   int k_max   In the final part of the algorithm, the assigned class is decided by running a MAP decision. Therefore, the class having the
                        highest a posteriori probability (that is, the probability described with "alpha") is the one assigned to that particular
                        position. The integer 'k_max' saves the class with the highest probability for a given position.

        -   double numerator1    For each iteration (t), for each position (n) and for each class (k), this double corresponds to the numerator of the
                                a posteriori probability (alpha).

        -   double denominator1  For each iteration (t), for each position (n) and for each class (k), this double corresponds to the denominator of the
                                a posteriori probability (alpha).

        -   double numerator2    For each iteration (t), for each learner (m), for each class (k) and for each nucleobase (l), this double corresponds to
                                the numerator of gamma. This gamma corresponds to (l,k)th position of the m'th learner confusion matrix a the t'th
                                iteration.

        -   double denominator2  For each iteration (t), for each learner (m), for each class (k) and for each nucleobase (l), this double corresponds to
                                the denominator of gamma. This gamma corresponds to (l,k)th position of the m'th learner confusion matrix a the t'th
                                iteration.

        -   double sum_check1    For each iteration (t), for each position (n), this double corresponds to the summation of every k'th alpha (a posteriori
                                probability). Since the a posteriori probabilities represent mutually exclusive events, 'sum_check1' must be equal to 1.
                                It is used just to ensure the calculations proceed properly. Note that, considering that these a posteriori probabilities
                                are doubles, a +-0.001 error margin is allowed.

        -   double sum   For each iteration (t) and for each class (k), this double corresponds to the numerator of the a priorit probability (pi)

        -   double f1[K][M]  For each position (n), for each class (k) and for each learner (m), this double is used to calculate the corresponding alpha
                            (a posteriori probability), representing the L-element product. It is used to calculate f2[K]

        -   double f2[K]     For each position (n) and for each class (k), this double is used to calculate the corresponding alpha (a posteriori probability),
                            representing the M-element product. It appears in both the numerator and the denominator of the corresponding alpha.
    */

    int k_max;
    double alpha_sum[K];
    // Iterations
    for (int t = 0; t < T; t++) {
        // printf("Iteration number: %d\n", t + 1);
        // EXPECTATION STEP, E-Step
        // Compute alpha
        memset(alpha_sum, 0, K * sizeof(double));
        for (int n = 0; n < N; n++) {
            compute_alpha(t, n, alpha, alpha_sum, T, log_pi, log_gamma, s);
        }
        // MAXIMIZATION STEP, M-Step
        // Compute gamma
        for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) {
                compute_gamma(m, t, k, T, N, log_gamma, alpha, s);
            }
        }
        // Compute pi
        for (int k = 0; k < K; k++) {
            log_pi[k][t + 1] = LOG(alpha_sum[k]) - LOG(N);
        }
        // Printing next a priori probability
        // print_pi(T, log_pi, t + 1);
    }

    // MAP decision. Calculate y = argmax(alpha)
    k_max = 0;
    // #pragma omp parallel for
    for (int n = 0; n < N; n++) {
        // Compute argmax(alpha)
        for (int k = 0; k < K; k++) {
            if (alpha[offset_2d(n, k, K)] > alpha[offset_2d(n, k_max, K)]) {
                k_max = k;
            }
        }
        // Assign class decided through EM
        nobv_indxs[offset_2d(n, 7, 9)] = k_max;  // Finally, the third column of 'nobv_indxs' is filled with the algorithm decisions
        if (alpha[offset_2d(n, k_max, K)] != 1) {
            nobv_indxs[offset_2d(n, 8, 9)] = (int)-10 * log10(1.0 - alpha[offset_2d(n, k_max, K)]);
        } else {  // When working with very high 'alpha's (almost 0), we may find that we exceed the precision of the compiler, and it returns that the result of
                  //-10*log10(1-alpha) = infinite. In that case, an arbitrary value of 1000 is assigned to the QS, just as it is done with the OBVIOUS VARIANT cases.
            nobv_indxs[offset_2d(n, 8, 9)] = 1000;
        }
    }
}

inline void compute_alpha(int t, int n, double *alpha, double *alpha_sum, int T, double log_pi[K][T + 1], double log_gamma[L][K][M][T + 1], int *s) {
    double log_f[K], numerator[K], denominator;

    double max_log_f = -DBL_MAX;
    // Calculating log(f) = sums of s(n,m,l)*log_gamma(m,l,k,t) for each k, it is the logarithm of equation (10) of https://upcommons.upc.edu/bitstream/handle/2117/367991/1-s2.0-S0031320322002023-main.pdf
    for (int k = 0; k < K; k++) {
        log_f[k] = 0;
        for (int m = 0; m < M; m++) {
            for (int l = 0; l < L; l++) {
                if (s[offset_3d(n, m, l)] > 0)
                    log_f[k] += s[offset_3d(n, m, l)] * log_gamma[l][k][m][t];
            }
        }
        if (log_pi[k][t] + log_f[k] > max_log_f) 
            max_log_f = log_pi[k][t] + log_f[k];
    }

    // Calculating alpha(n, k, t) = pi(k, t) * f(k) / sum(pi(k, t) * f(k)) for each k
    // which is equal to 2^(log_pi(k, t) + log(f(k))/ sum(2^(log_pi(k, t) + log(f(k))
    denominator = 0;
    for (int k = 0; k < K; k++)
        denominator += numerator[k] = EXP(log_pi[k][t] + log_f[k] - max_log_f);

    assert(denominator > 0);
    for (int k = 0; k < K; k++)
        alpha_sum[k] += alpha[offset_2d(n, k, K)] = numerator[k] / denominator;
}

inline void compute_gamma(int m, int t, int k, int T, int N, double log_gamma[L][K][M][T + 1], double *alpha, int *s) {
    double aux, numerator[L], denominator;

    denominator = 0.001;

    for (int l = 0; l < L; l++)
        numerator[l] = 0.001;
    // #pragma omp parallel for
    for (int n = 0; n < N; n++) {
        for (int l = 0; l < L; l++) {
            if (s[offset_3d(n, m, l)] != 0) {
                aux = alpha[offset_2d(n, k, K)] * s[offset_3d(n, m, l)];
                numerator[l] += aux;
                denominator += aux;
            }
        }
    }
    assert(denominator > 0);
    for (int l = 0; l < L; l++)
        log_gamma[l][k][m][t + 1] = LOG(numerator[l]) - LOG(denominator);

}

void compare_to_ref(int N, int *nobv_var_indxs, int *nobv_invar_indxs, int *nobv_var_cnt, int *nobv_invar_cnt, int *nobv_indxs) {
    /*
    Until now, we had only listed those positions whose class had not been yet determined and consequently they had entered the EM algorithm. Now the
    algorithm has run, and those NOT OBVIOUS positions have been assigned a class. It is important for further analysis to compare these decisions
    with the reference file, in order to determine whether the algorithm considers that a variant has been detected. Both the reference class and the decided
    class are saved in the 2nd and 8rd column of the 'nobv_indxs' matrix.

    Therefore, after running the algorithm, we may find ourselves in 2 situations:
        -   The decided class corresponds to the one of the reference file (note that that will only happen if k=0,1,2,3: AA, CC, GG, TT). In that case,
            we may account it as NOT OBVIOUS INVARIANT (nobv_invar)
        -   The decided class does not correspond to the one of the reference file (for k=4,5,6,7,8,9). In that case, we shall count it as NOT OBVIOUS VARIANT (nobv_var)

    The same way as before, the index (aka position/locus) of the decision is stored accordingly:
        -   If no variant was decided (NOT OBVIOUS INVARIANT), the information shall be stored in the 2-column matrix nobv_invar_indxs
        -   If a variant was decided (NOT OBVIOUS VARIANT), the information shall be stored in the 9-column matrix nobv_var_indxs

    Also, note there is a counter for each matrix, indicating the next position to be filled (plus serving an indicator of the number of positions of each type)
    */
    for (int i = 0; i < N; i++) {
        if (nobv_indxs[offset_2d(i, 7, 9)] == nobv_indxs[offset_2d(i, 1, 9)]) {                   // nobv_invar
            nobv_invar_indxs[offset_2d(*nobv_invar_cnt, 0, 2)] = nobv_indxs[offset_2d(i, 0, 9)];  // Save the index (position) number in first column and the nobv_invar_cnt'th line of nobv_invar_indxs
            nobv_invar_indxs[offset_2d(*nobv_invar_cnt, 1, 2)] = nobv_indxs[offset_2d(i, 7, 9)];  // Save the invariant class in second column and the nobv_invar_cnt'th line of nobv_invar_indxs
            *nobv_invar_cnt = *nobv_invar_cnt + 1;                                                // Updating counter
        } else {                                                                                  // nobv_var
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 0, 10)] = nobv_indxs[offset_2d(i, 0, 9)];      // first column, saving position/locus
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 1, 10)] = nobv_indxs[offset_2d(i, 1, 9)];      // second column, saving reference nucleobase
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 2, 10)] = nobv_indxs[offset_2d(i, 2, 9)];      // third column, saving depth
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 3, 10)] = nobv_indxs[offset_2d(i, 3, 9)];      // fourth column, saving number of forward reads mapping REF
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 4, 10)] = nobv_indxs[offset_2d(i, 4, 9)];      // fifth column, saving number of forward reads not mapping REF
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 5, 10)] = nobv_indxs[offset_2d(i, 5, 9)];      // sixth column, saving number of reverse reads mapping REF
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 6, 10)] = nobv_indxs[offset_2d(i, 6, 9)];      // seventh column, saving number of reverse reads not mapping REF
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 7, 10)] = nobv_indxs[offset_2d(i, 7, 9)];      // eighth column, saving decided class
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 8, 10)] = nobv_indxs[offset_2d(i, 8, 9)];      // ninth column, saving QS of the decision
            nobv_var_indxs[offset_2d(*nobv_var_cnt, 9, 10)] = i;                                   // tenth column, saving the index of the position in the nobv_indxs matrix
            *nobv_var_cnt = *nobv_var_cnt + 1;  // Updating counter
        }
    }
}

void compare_to_ground_truth(int *tp, int *fp, int *fn, int obv_var_cnt, int nobv_var_cnt, int *obv_var_indxs, int *nobv_var_indxs, FILE *gt_file) {
    /*
    Right now, we have every supposed variant stored in 'obv_var_indxs' (currently saving 'obv_var_cnt' positions) and in 'nobv_var_indxs' (currently
    saving 'nobv_var_cnt' positions). Also, we have already created the variables in charge of counting the TRUE POSITIVES, FALSE POSITIVES and FALSE
    NEGATIVES: 'tp', 'fp' and 'fn', respectively.
    The job of this function is to compare all decided variant and invariant positions to the ones from the ground truth. Again, note that:
        -   TRUE POSITIVE: the same variant that has been decided appears in the ground truth. That is, a particular position appears in either
            'obv_var_indxs' or in 'nobv_var_indxs', and in the ground truth file (note that in the ground truth file is composed only by detected
            variants). Furthermore, the variant showing in in the ground truth must be the same that has been assigned through our algorithm.
        -   FALSE POSITIVE: a variant that has been decided does not appear in the ground truth. This may happen because either a variant was decided
            where there was none, or because there was after all a variant in that position, but the decided class is not the one established in the
            ground truth.
        -   TRUE NEGATIVE: no variant has been detected where, according to the ground truth, there's indeed none. Again, note that these ones will not
            be counted because of the little importance they have when computing precision and sensitivity.
        -   FALSE NEGATIVE: no variant has been detected at a position where there is actually one. That is, a variant/line appears in the ground truth
            but does not appear in the 'obv_var_indxs' nor in the 'nobv_var_indxs'.

    To tackle this problem, we may use the "int read_ground_truth()" function (see description below), whose main purpose is to determine the class
    of the ground-truth-variant, by returning {0,1,...9} accordingly, and write down (onto 'pos') the locus of the corresponding variant. Note that,
    as explained in the corresponding function, when a faulty line is read, it returns '-1', and when an INDEL is scanned, it shall return '10'.

    Since all lines in the ground truth file and the 'obv_var_indxs' and 'nobv_var_indxs' are ordered by position, it comes specially handy to have
    index integers that point to the next-to-be-considered position of de EMVC-decided-variants: 'i_obv_var' and 'i_nobv_var'.

    For clarity, further detailed comments are written inbetween the code.
    */

    int pos = 0, class = 0, i_obv_var = 0, i_nobv_var = 0;
    // To ensure the ground truth file has been properly pointed
    if (gt_file == NULL) {
        exit(EXIT_FAILURE);
    }

    class = read_ground_truth(&pos, gt_file);

    /*
    Using a while, we will read every line of the ground truth file, either until the end, or until the positions of the variants being read surpass
    those of the 'obv_var_indxs', or 'nobv_var_indxs', to avoid counting as FALSE NEGATIVE lines that were not initially appearing in the output of
    the mpileup function.
    */
    while (class != 11 && ((pos <= obv_var_indxs[offset_2d(obv_var_cnt - 1, 0, 9)]) || (pos <= nobv_var_indxs[offset_2d(nobv_var_cnt - 1, 0, 10)]))) {
        /*
        If the class returned by the "int read_ground_truth()" function is either '10' or '-1', the loop shall continue without altering any of the counters
        (either the line was faulty, or the line contained an INDEL, and we don't work with those).
        */
        if (class == -1 || class == 10) {
            class = read_ground_truth(&pos, gt_file);
            continue;
        }

        /*
        If the position 'pos' of the newly-read ground-truth-variant is higher than the next-in-line position of 'nobv_var_indxs' or 'obv_var_indxs', it
        means that variants have been detected NOT OBVIOUSly or OBVIOUSly (respectively) where no variants have been signaled by the ground truth. Therefore,
        we shall account this situations as FALSE POSITIVE positions. The corresponding 'i_obv_var' or 'i_nobv_var' will have to be increased accordingly,
        until the appointed position does coincide with 'pos' or surpasses 'pos'.
        */
        while ((pos > nobv_var_indxs[offset_2d(i_nobv_var, 0, 10)]) && (i_nobv_var < nobv_var_cnt)) {
            *fp = *fp + 1;
            i_nobv_var++;
        }
        while ((pos > obv_var_indxs[offset_2d(i_obv_var, 0, 9)]) && (i_obv_var < obv_var_cnt)) {
            *fp = *fp + 1;
            i_obv_var++;
        }

        // After that, we may find ourselves in 5 situations:
        /*
        1.- 'pos' equals the position pointed by 'i_nobv_var', as in NOT OBVIOUS VARIANT (the position of a variant has been correctly identified)
            and the class pointed by 'i_nobv_var' equals the 'class' designated by the ground truth (the class of a variant has been correctly
            identified). In this case, we are dealing with a TRUE POSITIVE.
        */
        if ((pos == nobv_var_indxs[offset_2d(i_nobv_var, 0, 10)]) && (class == nobv_var_indxs[offset_2d(i_nobv_var, 7, 10)])) {
            *tp = *tp + 1;
            i_nobv_var++;
        }
        /*
        2.- 'pos' equals the position pointed by 'i_obv_var', as in OBVIOUS variant (the position of a variant has been correctly identified) and
            the class pointed by 'i_obv_var' equals the 'class' designated by the ground truth (the class of a variant has been correctly identified).
            In this case, we are also dealing with a TRUE POSITIVE.
        */
        else if ((pos == obv_var_indxs[offset_2d(i_obv_var, 0, 9)]) && (class == obv_var_indxs[offset_2d(i_obv_var, 7, 9)])) {
            *tp = *tp + 1;
            i_obv_var++;
        }
        /*
        3.- 'pos' equals the position pointed by 'i_nobv_var', as in NOT OBVIOUS VARIANT (the position of a variant has been correctly identified),
            but the class pointed by 'i_nobv_var' does not equal the 'class' designated by the ground truth (the class of a variant has been incorrectly
            identified). In this case, we are dealing with a FALSE POSITIVE.
        */
        else if ((pos == nobv_var_indxs[offset_2d(i_nobv_var, 0, 10)]) && (class != nobv_var_indxs[offset_2d(i_nobv_var, 7, 10)])) {
            *fp = *fp + 1;
            i_nobv_var++;
        }
        /*
        4.- 'pos' equals the position pointed by 'i_obv_var', as in OBVIOUS VARIANT (the position of a variant has been correctly identified), but
            the class pointed by 'i_obv_var' does not equal the 'class' designated by the ground truth (the class of a variant has been incorrectly
            identified). In this case, we are dealing with a FALSE POSITIVE.
        */
        else if ((pos == obv_var_indxs[offset_2d(i_obv_var, 0, 9)]) && (class != obv_var_indxs[offset_2d(i_obv_var, 7, 9)])) {
            *fp = *fp + 1;
            i_obv_var++;
        }
        /*
        5.- In any other case, we will be dealing with a FALSE NEGATIVE. This cases will mainly happen positions appointed by 'i_nobv_var' or 'i_obv_var'
            are now both higher than 'pos' (we have already ensured that they may not be lower than 'pos'). Therefore, it means that the ground truth
            counts a variant position that does not appear as OBVIOUS VARIANT nor NOT OBVIOUS VARIANT.
        */
        else {
            *fn = *fn + 1;
        }
        // We read the class for the next loop
        class = read_ground_truth(&pos, gt_file);
    }
}

int read_ground_truth(int *position, FILE *gt_file) {
    /*
    This function's purpose is to read a line of the ground_truth.vcf file, and extract the position in which the variant is found (note that all
    lines in the file represent variants), and the actual class for this corresponding position.

    For our purpose, the ground truth offers different relevant information, presented in tab-separated columns:
        -   CHROM: 1st column, the chromosome name.
        -   POS: 2nd column, position/locus with the variant.
        -   REF: 4th column, reference nucleobase. For proper use, this one should be the same as the one used as reference in the pileup.
        -   ALT: 5th column, alternate alleles for this variant.
        -   FILTER: 7th column, if a filter is applied and the line/variant has passed it. We will only work with "PASS"ed lines.
        -   GENOTYPE: 10th column, specially the part with 0/1, 1/1, 1/2...
                ·   If we have a 0/1 or 0|1 in INFO, it means that the variant presents one allele like the one in reference, and one alternate allele.
                ·   If we have a 1/1 or 1|1 in INFO, it means that both alleles of the variant present the alternate form.
                ·   If we have a 1/2 or 1|2 in INFO, it means that both alleles are different from the reference one, and between them, taking two different
                    alternate forms (comma-separated).
    For more information on the vcf format, check the corresponding data sheets.

    Remember that this algorithm does not manage INDELS yet, so this kind
    of lines have to be treated specially:
        "20 67500   .   T   TTGGTATCTAG 50  PASS    platforms=2;platformnames=Illumina,CG;datasets=2;datasetnames=HiSeqPE300x,CGnormal;callsets=3;callsetnames=HiSeqPE300xfreebayes,HiSeqPE300xGATK,CGnormal;datasetsmissingcall=IonExome,SolidPE50x50bp,SolidSE75bp;lowcov=CS_IonExomeTVC_lowcov,CS_SolidPE50x50GATKHC_lowcov,CS_SolidSE75GATKHC_lowcov    GT:PS:DP:GQ 0/1:.:627:922"
        "20 106703  .   C   CT,CTTT 50  PASS    platforms=2;platformnames=Illumina,CG;datasets=2;datasetnames=HiSeqPE300x,CGnormal;callsets=3;callsetnames=HiSeqPE300xfreebayes,HiSeqPE300xGATK,CGnormal;datasetsmissingcall=IonExome,SolidPE50x50bp,SolidSE75bp;lowcov=CS_IonExomeTVC_lowcov,CS_SolidPE50x50GATKHC_lowcov,CS_SolidSE75GATKHC_lowcov;filt=CS_CGnormal_filt  GT:PS:DP:GQ 1|2:.:426:355"
    All insertions and deletions have the same first nucleobase in both the REF column and ALT column. Therefore, we can distinguish and discard this
    INDEL lines by taking the first letter of the appointed alternatives/reference and check if the resulting class equals the one of reference.
        - In this example, for the first line, we would check it's a 0/1, and take the "T" from REF and the first "T" from ALT, thus meaning a "TT" class
        (k=2), which is equal to the reference class REF.
        - In this example, for the second line, we would check it's a 1|2, and take the "C" from the first ALT and the first "C" from the second ALT,
        thus resulting in a "CC" class (k=1), which is equal to the one as reference REF.
    In all these cases, instead of returning the corresponding NOT VARIANT class, we will return k=10, to notify the function user of this situation.

    Therefore, we can always work by taking the first letters of the corresponding REF or ALT (indicated by INFO), obtaining its class using the
    "int class_to_num()" function, comparing it with the reference class and deciding whether it's an actual variant for us or no. If it's not,
    we shall return '10', and if it is, we shall return the calculated class.

    */
    char *line = NULL, *token = NULL, *pos=NULL, *ref=NULL, *alt=NULL, *alt1, *alt2, *filter=NULL, *gt=NULL, *ptr, a1='A', a2='A';
    int class, ref0;
    size_t len = 0, ret = 0;
    int token_cnt;

    // READING THE FILE
    token_cnt = 0;
    // Obtaining a line of the VCF file. Saving its length in 'len' and string in 'line'.
    ret = getline(&line, &len, gt_file);

    if (ret == -1) {
        return 11;
    }

    // VCF files usually present a header providing information of the file in general, made up of several lines starting by "##". We shall skip this header
    if ((line[0] == '#')) {  // If the first character of the line is '#'
        return -1;
    }

    // Extracts tab-separated elements of the line (for future references, "tokens").
    token = strtok(line, "\t");

    while (token != NULL) {  // As long as there are tokens to extract
        /*
        The 'while' will keep working as long as there are tab-separated element in the line. Nontheless, VCF files allow to present different samples
        for the same reference. In this case, we would have more than 10 columns, and an error should show up. Also, it would not make sense to use
        a ground truth with multiple samples.
        */
        if (token_cnt > 9) {  // A token counter is used to check the column from which the token has been extracted.
            printf("\nError! This software does not support more than one input samples. Only the first sample of the ground_truth.vcf file will be used\n");
        }
        switch (token_cnt) {  // Depending on the particular token we are extracting, we shall save (or not) the information in the corresponding variable
            case 0:
                break;  // chromosome, CHROM
            case 1:
                pos = token;
                break;  // position, POS
            case 2:
                break;  // identification, ID
            case 3:
                ref = token;
                break;  // reference, REF
            case 4:
                alt = token;
                break;  // alternate alleles, ALT
            case 5:
                break;  // quality value, QUAL
            case 6:
                filter = token;
                break;  // filtered, FILTER
            case 7:
                break;  // information, INFO
            case 8:
                break;  // format, FORMAT
            case 9:
                gt = token;
                break;  // genotype + other info
            default:
                printf("Error when reading. Should not be more than 9 tokens\n");
        }
        token = strtok(NULL, "\t");  // Extract the next token
        token_cnt++;                 // Update the token counter
    }

    // In case the corresponding ground truth line hasn't passed the filter, we will discard it, and return '-1'
    if ((filter[0] != 'P') || (filter[1] != 'A') || (filter[2] != 'S') || (filter[3] != 'S')) {
        return -1;
    }

    // To turn the 'pos' read, which is a string, into an integer, and already save it in the corresponding function parameter
    *position = strtol(pos, &ptr, 10);

    /*
    Extract the two possible alternate alleles. There can't be more, since we are not working with more than one sample. They are comma-separated,
    so they can also be extracted as tokens.
    */
    alt1 = strtok(alt, ",");
    alt2 = strtok(NULL, ",");

    // gt[0] indicates whether the first allele is the reference or an alternate. It is the first number from 0/1, 1/1, 1/2...
    switch (gt[0]) {
        // a1 saves the nucleobase letter of the first allele
        case '0':
            a1 = ref[0];
            break;
        case '1':
            a1 = alt1[0];
            break;
        case '2':
            a1 = alt2[0];
            break;
        default:
            printf("Error when reading. Should not have any other value\n");
    }
    // gt[2] indicates whether the second allele is the reference or an alternate. It is the second number from 0/1, 1/1, 1/2...
    switch (gt[2]) {
        // a2 saves the nucleobase letter of the second allele
        case '0':
            a2 = ref[0];
            break;
        case '1':
            a2 = alt1[0];
            break;
        case '2':
            a2 = alt2[0];
            break;
        default:
            printf("Error when reading. Should not have any other value\n");
    }

    class = class_to_num(a1, a2);

    // To turn the REF nucleobase into an integer: {A,C,T,G}={0,1,2,3}
    ref0 = nuc_to_num(ref[0]);
    /*switch(ref[0]){
        case 'A': ref0=0; break;
        case 'C': ref0=1; break;
        case 'T': ref0=2; break;
        case 'G': ref0=3; break;
        default: printf("ref[0] has an imporper value");
    }*/

    // In case it is an INDEL line
    if (class == ref0) {
        return 10;
    }

    return class;
}

int class_to_num(char a1, char a2) {
    /*
    It has been seen that a function like this can be very practical (and clean), when determining the corresponding class of a position given the
    nucleobases of each chromosome.

    AA      class 0
    CC      class 1
    TT      class 2
    GG      class 3
    AC, CA  class 4
    AT, TA  class 5
    AG, GA  class 6
    CT, TC  class 7
    CG, GC  class 8
    TG, TG  class 9
    */

    if (a1 == 'A' && a2 == 'A') {
        return 0;
    } else if (a1 == 'C' && a2 == 'C') {
        return 1;
    } else if (a1 == 'T' && a2 == 'T') {
        return 2;
    } else if (a1 == 'G' && a2 == 'G') {
        return 3;
    } else if ((a1 == 'A' && a2 == 'C') || (a1 == 'C' && a2 == 'A')) {
        return 4;
    } else if ((a1 == 'A' && a2 == 'T') || (a1 == 'T' && a2 == 'A')) {
        return 5;
    } else if ((a1 == 'A' && a2 == 'G') || (a1 == 'G' && a2 == 'A')) {
        return 6;
    } else if ((a1 == 'C' && a2 == 'T') || (a1 == 'T' && a2 == 'C')) {
        return 7;
    } else if ((a1 == 'C' && a2 == 'G') || (a1 == 'G' && a2 == 'C')) {
        return 8;
    } else if ((a1 == 'T' && a2 == 'G') || (a1 == 'G' && a2 == 'T')) {
        return 9;
    } else {
        printf("Error! Could not calculate class\n");
        return -1;
    }
}

int nuc_to_num(char nuc) {
    /*
    Translates a upper or lower case nucleobase character into the corresponding number.
    */
    if (nuc == 'a' || nuc == 'A') {
        return 0;
    } else if (nuc == 'c' || nuc == 'C') {
        return 1;
    } else if (nuc == 't' || nuc == 'T') {
        return 2;
    } else if (nuc == 'g' || nuc == 'G') {
        return 3;
    } else {
        return -1;
    }
    // else {
    //     printf("Error! None of A,C,T,G\n nuc = %c\n", nuc);
    // }
}

double compute_precision(int tp, int fp) {
    /*
    For clarity and neatness, this function has been defined to compute the precision given the number of true positives 'tp' and false positives 'fp'
    It is given as a percentage with 4 decimal figures.
    */
    if ((tp + fp) == 0) {  // To make sure the denominator is not 0. This shall only happen with very few reads (or positions). In that case, precision will
                           // be 0/0.
        printf("\nError! Denominator can't be 0\n");
        return 0;
    }
    double prec = (double)100 * tp / (tp + fp);
    return prec;
}

double compute_sensitivity(int tp, int fn) {
    /*
    For clarity and neatness, this function has been defined to compute the sensitivity given the number of true positives 'tp' and false negatives 'fn'
    It is given as a percentage with 4 decimal figures.
    */
    if ((tp + fn) == 0) {  // To make sure the denominator is not 0. This shall only happen with very few reads (or positions). In that case, sensitivity will
                           // be 0/0.
        printf("Error! Denominator can't be 0\n");
        return 0;
    }
    double sens = (double)100 * tp / (tp + fn);
    return sens;
}

#include <locale.h>

void generate_vcf(char *chr, int obv_var_cnt, int nobv_var_cnt, int *obv_var_indxs, int *nobv_var_indxs, double *alpha, char *file_name) {
    /*
    Firstly, the corresponding VCF file is created and the proper header is printed on top.

    This function works with two different variant matrixes ('obv_var_indxs' and 'nobv_var_indxs'), each one properly ordered by position. The goal is to
    print the information they contain in the VCF format, mixing NOT OBVIOUS and OBVIOUS variants (in the end, it is of little interest for the user to
    know whether the detected variant has been decided through the algorithm or not) and ordering them by position.

    To print all variants by order of appearance (by position), 'i_obv' and 'i_nobv' play a specially important role: they will point the next-to-be-printed
    variant of each matrix. Therefore, they will range from 0 to 'obv_var_cnt' and 'nobv_var_cnt' respectively. The rest of variables are integers, chars or
    strings used to store information that will be printed in each line (for more information, check VCF data sheets):
        *chr    -   Chromosome name. Parameter input to the function.
        pos     -   Position/locus. First column of the 'obv_var_indxs' and 'nobv_var_indxs' files.
        id      -   Identifiers. Will not apply in this technique. Will always be '.'.
        ref     -   Reference base. Second column of the 'obv_var_indxs' and 'nobv_var_indxs' files.
        alt[4]  -   Alternate alleles. Comma-separated list of alternate non-reference alleles. Calculated through the 'void calculate_paramVCF()' function. It's
                    a 4-position-long char array, allowing to store a 3-characters-long string: in our case, it will only have to store a maximum of 2 alleles + comma.
        qual    -   Quality score for the assertation made in ALT. Will not apply in this technique. Will always be '0'.
        filter  -   No filter has been applied. Therefore, it is always set to '.'.
        info    -   Additional information. In our case, no extra information is displayed.
        format  -   Informing about the displayed data in the corresponding column. In our case, only the genotype (0/1, 1/1 and 1/2) will be displayed.
        gt[4]   -   Genotype of the variant, which is all the information displayed about it. Phased and unphased pairs are not specified (no distinction between
                    '/' or '|'), so all genotypes will be: '0/1', '1/1' or '1/2'. In every case its a 3-positions-long char array, thus a 4-positions-long string.

    The function works as follows. A while loop is used to browse through both the 'obv_var_indxs' and 'nobv_var_indxs', as long as i_obv+i_nobv does not exceed
    the total number of variants detected (which is obv_var_cnt+nobv_var_cnt). The next variant to be printed will be from the OBVIOUS matrix instead of the NOT
    OBVIOUS matrix in the following cases:
        - The next-in-line variant from the 'obv_var_indxs' matrix presents a lower position/locus than the next-in-line variant from the 'nobv_var_indxs'
          matrix. Note that 'i_obv' must have not yet exceeded its maximum value.
        - All variant from the NOT OBVIOUS matrix have already been printed. In this case, all the remaining variants from 'obv_var_indxs' shall be printed
          all together. Note that 'i_obv' must have not yet exceeded its maximum value.

    In the rest of cases, the next variant to be printed will come from the NOT OBVIOUS matrix ('nobv_var_indxs').
    */
    FILE *f;
    int i_obv = 0, i_nobv = 0, qual, pos, len, fwd_ref, fwd_alt, rev_ref, rev_alt;
    char ref, alt[4], gt[4], *id = ".", *filter = ".", *info = ".", *format = "GT:DP:FR:FA:RR:RA:PROB_AA:PROB_CC:PROB_TT:PROB_GG:PROB_AC:PROB_AT:PROB_AG:PROB_CT:PROB_CG:PROB_TG";
    f = fopen(file_name, "w");

    fprintf(f, "##fileformat=VCFv4.2\n");
    fprintf(f, "##source=EMVC\n");
    fprintf(f, "##phasing=unphased\n");
    fprintf(f, "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n");
    fprintf(f, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description='Read Depth'>\n");
    fprintf(f, "##FORMAT=<ID=FR,Number=2,Type=Integer,Description='Number of REF mapped forward'>\n");
    fprintf(f, "##FORMAT=<ID=FA,Number=2,Type=Integer,Description='Number of ALT mapped forward'>\n");
    fprintf(f, "##FORMAT=<ID=RR,Number=2,Type=Integer,Description='Number of REF mapped reverse'>\n");
    fprintf(f, "##FORMAT=<ID=RA,Number=2,Type=Integer,Description='Number of ALT mapped reverse'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_AA,Number=1,Type=Float,Description='Probability of the AA class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_CC,Number=1,Type=Float,Description='Probability of the CC class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_TT,Number=1,Type=Float,Description='Probability of the TT class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_GG,Number=1,Type=Float,Description='Probability of the GG class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_AC,Number=1,Type=Float,Description='Probability of the AC class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_AT,Number=1,Type=Float,Description='Probability of the AT class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_AG,Number=1,Type=Float,Description='Probability of the AG class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_CT,Number=1,Type=Float,Description='Probability of the CT class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_CG,Number=1,Type=Float,Description='Probability of the CG class'>\n");
    fprintf(f, "##FORMAT=<ID=PROB_TG,Number=1,Type=Float,Description='Probability of the TG class'>\n");

    fprintf(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN12878\n");

    // setlocale(LC_NUMERIC, "de_DE"); // ".OCP" if you want to use system settings
    while ((i_obv + i_nobv) < (obv_var_cnt + nobv_var_cnt)) {
        if (((obv_var_indxs[offset_2d(i_obv, 0, 9)] < nobv_var_indxs[offset_2d(i_nobv, 0, 10)]) || i_nobv == nobv_var_cnt) && (i_obv < obv_var_cnt)) {
            calculate_paramVCF(f, obv_var_indxs, i_obv, &pos, &ref, &len, &fwd_ref, &fwd_alt, &rev_ref, &rev_alt, &qual, alt, gt, 9);  // pos, ref, alt and gt are read and computed
            
            fprintf(f, "%s\t%d\t%s\t%c\t%s\t%d\t%s\t%s\t%s\t%s:%d:%d:%d:%d:%d:", 
                    chr, pos, id, ref, alt, qual, filter, info, format, gt, len, fwd_ref, 
                    fwd_alt, rev_ref, rev_alt); 
            
            if (obv_var_indxs[offset_2d(i_obv, 7, 9)] == 0) 
                fprintf(f, "1");
            else 
                fprintf(f, "0");
            for (int k = 1; k < K; k++) {
                if (obv_var_indxs[offset_2d(i_obv, 7, 9)] == k) 
                    fprintf(f, ":1");
                else 
                    fprintf(f, ":0");
            }
            fprintf(f, "\n"); 
            i_obv++;
                                                                                                                           // Index updated
        } else if (((obv_var_indxs[offset_2d(i_obv, 0, 9)] > nobv_var_indxs[offset_2d(i_nobv, 0, 10)]) || i_obv == obv_var_cnt) && (i_nobv < nobv_var_cnt)) {
            calculate_paramVCF(f, nobv_var_indxs, i_nobv, &pos, &ref, &len, &fwd_ref, &fwd_alt, &rev_ref, &rev_alt, &qual, alt, gt, 10);  // pos, ref, alt and gt are read and computed  
            
            fprintf(f, "%s\t%d\t%s\t%c\t%s\t%d\t%s\t%s\t%s\t%s:%d:%d:%d:%d:%d:", 
                    chr, pos, id, ref, alt, qual, filter, info, format, gt, len, fwd_ref, 
                    fwd_alt, rev_ref, rev_alt); 

            int alpha_idx = nobv_var_indxs[offset_2d(i_nobv, 9, 10)];
            fprintf(f, "%.10f", alpha[offset_2d(alpha_idx, 0, K)]);
            for (int k = 1; k < K; k++) {
                fprintf(f, ":%.10f", alpha[offset_2d(alpha_idx, k, K)]);
            }

            fprintf(f, "\n");   
            i_nobv++;                                                                                                            // Index updated
        } else {
            printf("Error!\n");
        }
    }
    fclose(f);
    // printf("A file has been generated by the name: '%s'\n", file_name);
}

void calculate_paramVCF(FILE *f, int *var_indxs, int i, int *pos, char *ref, int *len, int *fwd_ref, int *fwd_alt, int *rev_ref, int *rev_alt, int *qual, char *alt, char *gt, int ms) {
    /*
    At the i'th line of the 'var_indxs' memory unit (remember that it shall only be 'nobv_var_indxs' or 'obv_var_indxs'), it reads the line position in the
    chromosome and assigns it to 'pos'. The same does for 'ref', 'len', 'fwd_ref', 'fwd_alt', 'rev_ref', 'rev_alt' and 'qual'.
    At the i'th line of the 'var_indxs', it reads the corresponding class and translates it to the usual VCF nomenclature, using a genotype (0/1, 1/1, 1/2,
    etc.) and alternates (comma-separated):
        -Firstly, the decided class is considered ('var_indxs[offset_2d(i,2,3)]').
            · If it is of the kind {AA,CC,TT,GG} (k={0,1,2,3}), it will only be an actual variant if the decided class is different from the reference (for
              example, having ref='A', k='CC'). In this case, the ALT string will only be 1 position long (in the previous example, ALT="C"), and the genotype
              will always be of the type gt="1/1".

            · If it is of the kind {AC,AT,AG,CT,CG,TG} (k={4,5,6,7,8,9}), it will always be a variant to us, independently from the reference value. Considering
              that, we may find ourselves in two different situations:
                - The reference value corresponds to one of the EMVC-decided-variant alleles: for example class='AC', and ref='A' or ref='C'. In this case,
                  ALT will be a 1-position-long string saving the class allele not apearing in the reference. If class='AC', ref='A', alt="C". Correspondingly,
                  the genotype is always of the type gt="0/1".

                - The reference value does not correspond to one of the EMVC-decided-variant alleles: for example, class='AC', and ref='T' or ref='G'. In this
                  case, ALT will be a 3-position-long string saving each allele character, comma-separated. If class='AC', and ref='T', alt="A,C". Correspondingly,
                  the genotype is always of the type gt="1/2"

    *** Note that this function does not distinguish between phased and unphased genotypes, so no distinction is made between '|' and '/': All of them are
        with '/'.
    */
    *pos = var_indxs[offset_2d(i, 0, ms)];

    // Turn nucleobases from {0,1,2,3}-->{A,C,T,G}
    switch (var_indxs[offset_2d(i, 1, ms)]) {
        case 0:
            *ref = 'A';
            break;
        case 1:
            *ref = 'C';
            break;
        case 2:
            *ref = 'T';
            break;
        case 3:
            *ref = 'G';
            break;
        default:
            printf("Error! Reference value impossible\n");
    }

    *len = var_indxs[offset_2d(i, 2, ms)];

    *fwd_ref = var_indxs[offset_2d(i, 3, ms)];
    *fwd_alt = var_indxs[offset_2d(i, 4, ms)];
    *rev_ref = var_indxs[offset_2d(i, 5, ms)];
    *rev_alt = var_indxs[offset_2d(i, 6, ms)];

    *qual = var_indxs[offset_2d(i, 8, ms)];

    switch (var_indxs[offset_2d(i, 7, ms)]) {  // Considering the decided EMVC-decided-variant class...
        case 0:                               // AA
            if (*ref != 'A') {
                alt[0] = 'A';
                alt[1] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
                
            } 
            // else {
            //     printf("Error! Should not be a variant\n");
            // }
            break;
        case 1:  // CC
            if (*ref != 'C') {
                alt[0] = 'C';
                alt[1] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            }
            //  else {
            //     printf("Error! Should not be a variant\n");
            // }
            break;
        case 2:  // TT
            if (*ref != 'T') {
                alt[0] = 'T';
                alt[1] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            }
            //  else {
            //     printf("Error! Should not be a variant\n");
            // }
            break;
        case 3:  // GG
            if (*ref != 'G') {
                alt[0] = 'G';
                alt[1] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            }
            //  else {
            //     printf("Error! Should not be a variant\n");
            // }
            break;
        case 4:  // AC
            if (*ref == 'A') {
                alt[0] = 'C';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else if (*ref == 'C') {
                alt[0] = 'A';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else {
                alt[0] = 'A';
                alt[1] = ',';
                alt[2] = 'C';
                alt[3] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '2';
                gt[3] = '\0';
            }
            break;
        case 5:  // AT
            if (*ref == 'A') {
                alt[0] = 'T';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else if (*ref == 'T') {
                alt[0] = 'A';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else {
                alt[0] = 'A';
                alt[1] = ',';
                alt[2] = 'T';
                alt[3] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '2';
                gt[3] = '\0';
            }
            break;
        case 6:  // AG
            if (*ref == 'A') {
                alt[0] = 'G';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else if (*ref == 'G') {
                alt[0] = 'A';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else {
                alt[0] = 'A';
                alt[1] = ',';
                alt[2] = 'G';
                alt[3] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '2';
                gt[3] = '\0';
            }
            break;
        case 7:  // CT
            if (*ref == 'C') {
                alt[0] = 'T';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else if (*ref == 'T') {
                alt[0] = 'C';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else {
                alt[0] = 'C';
                alt[1] = ',';
                alt[2] = 'T';
                alt[3] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '2';
                gt[3] = '\0';
            }
            break;
        case 8:  // CG
            if (*ref == 'C') {
                alt[0] = 'G';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else if (*ref == 'G') {
                alt[0] = 'C';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else {
                alt[0] = 'C';
                alt[1] = ',';
                alt[2] = 'G';
                alt[3] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '2';
                gt[3] = '\0';
            }
            break;
        case 9:  // TG
            if (*ref == 'T') {
                alt[0] = 'G';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else if (*ref == 'G') {
                alt[0] = 'T';
                alt[1] = '\0';
                gt[0] = '0';
                gt[1] = '/';
                gt[2] = '1';
                gt[3] = '\0';
            } else {
                alt[0] = 'T';
                alt[1] = ',';
                alt[2] = 'G';
                alt[3] = '\0';
                gt[0] = '1';
                gt[1] = '/';
                gt[2] = '2';
                gt[3] = '\0';
            }
            break;
        default:
            printf("Error! Reference value impossible. pos = %d\t ref = %d\t class = %d\n", var_indxs[offset_2d(i, 0, ms)], var_indxs[offset_2d(i, 1, ms)], var_indxs[offset_2d(i, 4, ms)]);
    }
}

void print_init() {
    /*
    Function to print current specifications of the program executed.
    */
    printf("\n-----------------------------------------------------------------------------\n");
    printf("Version 15\n");
    printf("14/09/2022\n");
    printf("-----------------------------------------------------------------------------\n");
}

void print_s(int *s, int cnt) {
    /*
    Practical way to easily print the s tensor, if necessary. It will print all M learners and L nucleobases from the first to the cnt'th position
    It prints the L-dimension (different nucleobases) in closer columns (space-bar-separated) and the M-dimension (different learners) in further
    columns

    Example:

    m =     0               1               2               3               4               5               6
    l =     0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3
            A C T G         A C T G         A C T G         A C T G         A C T G         A C T G         A C T G

    n = 0   2 0 0 0         0 0 0 0         0 0 0 0         0 0 1 0         0 0 4 0         0 0 0 0         0 0 0 0
    n = 1   0 0 0 0         0 0 1 0         0 0 0 0         2 0 0 0         7 0 0 0         0 0 0 0         0 0 0 0
    n = 2   1 1 0 0         0 0 0 0         0 0 0 0         0 3 0 0         0 6 0 0         0 1 0 0         0 0 0 0
    n = 3   0 0 2 1         0 0 1 0         0 0 0 0         0 0 1 0         0 0 7 0         0 0 2 0         0 0 0 0
    n = 4   0 0 3 1         0 0 0 0         0 0 0 0         0 0 3 0         0 0 6 0         0 0 1 0         0 0 0 0
    n = 5   1 3 0 0         0 1 0 0         0 0 0 0         0 1 0 0         0 7 0 0         0 1 0 0         0 0 0 0
    n = 6   4 0 0 1         0 0 0 0         0 0 0 0         1 0 0 0         8 0 0 0         0 0 0 0         0 0 0 0
    n = 7   0 1 2 0         0 0 0 0         0 0 0 0         0 0 1 0         0 0 9 0         0 0 0 0         0 0 0 0
    n = 8   0 1 0 1         0 0 0 0         0 0 0 0         0 0 0 1         0 0 0 8         0 0 0 0         0 0 0 0
    ...

    (Note: the m, l (and their corresponding character equivalents {A,C,T,G}), and n values are not included in the printout)
    */
    for (int i = 0; i < cnt; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < L; k++) {
                printf("%d ", s[offset_3d(i, j, k)]);
            }
            printf("\t");
        }
        printf("\n");
    }
}

void print_s_file(int *s, int cnt, char *s_file_name) {
    /*
    Practical way to easily print the s tensor, if necessary. It will print all M learners and L nucleobases from the first to the cnt'th position
    It prints the L-dimension (different nucleobases) in closer columns (space-bar-separated) and the M-dimension (different learners) in further
    columns

    Example:

    m =     0               1               2               3               4               5               6
    l =     0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3         0 1 2 3
            A C T G         A C T G         A C T G         A C T G         A C T G         A C T G         A C T G

    n = 0   2 0 0 0         0 0 0 0         0 0 0 0         0 0 1 0         0 0 4 0         0 0 0 0         0 0 0 0
    n = 1   0 0 0 0         0 0 1 0         0 0 0 0         2 0 0 0         7 0 0 0         0 0 0 0         0 0 0 0
    n = 2   1 1 0 0         0 0 0 0         0 0 0 0         0 3 0 0         0 6 0 0         0 1 0 0         0 0 0 0
    n = 3   0 0 2 1         0 0 1 0         0 0 0 0         0 0 1 0         0 0 7 0         0 0 2 0         0 0 0 0
    n = 4   0 0 3 1         0 0 0 0         0 0 0 0         0 0 3 0         0 0 6 0         0 0 1 0         0 0 0 0
    n = 5   1 3 0 0         0 1 0 0         0 0 0 0         0 1 0 0         0 7 0 0         0 1 0 0         0 0 0 0
    n = 6   4 0 0 1         0 0 0 0         0 0 0 0         1 0 0 0         8 0 0 0         0 0 0 0         0 0 0 0
    n = 7   0 1 2 0         0 0 0 0         0 0 0 0         0 0 1 0         0 0 9 0         0 0 0 0         0 0 0 0
    n = 8   0 1 0 1         0 0 0 0         0 0 0 0         0 0 0 1         0 0 0 8         0 0 0 0         0 0 0 0
    ...

    (Note: the m, l (and their corresponding character equivalents {A,C,T,G}), and n values are not included in the printout)
    */
    FILE *file;
    if (strcmp(s_file_name, "-") != 0) {
        file = fopen(s_file_name, "w");
        for (int i = 0; i < cnt; i++) {
            fprintf(file, "%d\t\t\t", i);
            for (int j = 0; j < M; j++) {
                for (int k = 0; k < L; k++) {
                    fprintf(file, "%d ", s[offset_3d(i, j, k)]);
                }
                fprintf(file, "\t");
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

void print_indxs(int *indxs, int cnt, int col) {
    /*
    Practical way to easily print columns of 'indxs' matrixes
    */
    for (int i = 0; i < cnt; i++) {
        printf("%d\t", i);
        if (col == 1) {
            printf("%d\n", indxs[i]);
        } else if (col == 2) {
            printf("%d\t%d\n", indxs[offset_2d(i, 0, 9)], indxs[offset_2d(i, 1, 9)]);
        } else {
            printf("ERROR! Wrong function!");
        }
    }
}

void print_gamma(int T, double log_gamma[L][K][M][T + 1], int t) {
    /*
    Practical way to easily print the gamma matrix (confusion matrix) of each learner (m) at a given t'th iteration.

    Example:


    Confusion matrix of learner 1. gamma1=
    0.441   0.176   0.134   0.163   0.472   0.557   0.426   0.008   0.028   0.024
    0.230   0.521   0.197   0.137   0.496   0.027   0.031   0.479   0.418   0.003
    0.136   0.169   0.446   0.181   0.019   0.403   0.021   0.478   0.013   0.354
    0.193   0.134   0.222   0.519   0.013   0.013   0.523   0.035   0.541   0.618

    Confusion matrix of learner 2. gamma2=
    0.726   0.160   0.097   0.109   0.434   0.429   0.336   0.019   0.000   0.000
    0.107   0.660   0.078   0.068   0.544   0.021   0.009   0.500   0.521   0.000
    0.091   0.115   0.721   0.152   0.021   0.550   0.005   0.481   0.000   0.488
    0.077   0.064   0.104   0.671   0.001   0.000   0.650   0.000   0.479   0.512

    Confusion matrix of learner 3. gamma3=
    0.934   0.025   0.021   0.017   0.418   0.555   0.282   0.000   0.000   0.000
    0.023   0.946   0.020   0.011   0.582   0.000   0.002   0.724   0.540   0.000
    0.023   0.019   0.938   0.024   0.000   0.445   0.005   0.276   0.000   0.444
    0.020   0.010   0.022   0.947   0.000   0.000   0.712   0.000   0.460   0.556

    Confusion matrix of learner 4. gamma4=
    0.984   0.006   0.006   0.007   0.464   0.497   0.338   0.000   0.000   0.005
    0.004   0.983   0.005   0.004   0.536   0.000   0.001   0.684   0.532   0.006
    0.007   0.007   0.985   0.006   0.000   0.501   0.000   0.315   0.000   0.462
    0.006   0.004   0.004   0.983   0.000   0.002   0.661   0.001   0.468   0.528

    Confusion matrix of learner 5. gamma5=
    0.998   0.001   0.001   0.001   0.516   0.482   0.548   0.000   0.000   0.000
    0.000   0.997   0.001   0.001   0.484   0.000   0.000   0.451   0.541   0.001
    0.001   0.001   0.998   0.001   0.000   0.518   0.000   0.549   0.001   0.521
    0.001   0.001   0.000   0.997   0.000   0.000   0.452   0.000   0.458   0.478

    Confusion matrix of learner 6. gamma6=
    0.999   0.001   0.000   0.002   0.581   0.525   0.652   0.000   0.000   0.000
    0.000   0.998   0.000   0.000   0.419   0.000   0.000   0.332   0.451   0.000
    0.000   0.001   0.999   0.001   0.000   0.475   0.000   0.668   0.000   0.712
    0.000   0.000   0.000   0.997   0.000   0.000   0.348   0.000   0.549   0.288

    Confusion matrix of learner 7. gamma7=
    0.996   0.001   0.002   0.001   0.602   0.448   0.487   0.000   0.000   0.000
    0.000   0.995   0.002   0.000   0.398   0.000   0.000   0.476   0.365   0.000
    0.002   0.002   0.997   0.001   0.000   0.552   0.000   0.524   0.000   0.457
    0.002   0.002   0.000   0.997   0.000   0.000   0.513   0.000   0.635   0.543
    */
    for (int m = 0; m < M; m++) {
        printf("Confusion matrix of learner %d. gamma%d=\n", m + 1, m + 1);
        for (int l = 0; l < L; l++) {
            for (int k = 0; k < K; k++) {
                printf("%.3lf\t", EXP(log_gamma[l][k][m][t]));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void print_pi(int T, double log_pi[K][T + 1], int t) {
    /*
    Practical way to easily print the pi vector (a priori probability for every class (k)) at a given t'th iteration.

    Example:

    A priori probabilities. pi=
    k = 0       1       2       3       4       5       6       7       8       9
        AA      CC      TT      GG      AC      AT      AG      CT      CG      TG

        0.299   0.189   0.298   0.189   0.002   0.002   0.008   0.009   0.002   0.002

    (Note: the k values and their corresponding character equivalents {AA,CC,TT...} are not included in the printout)
    */
    printf("A priori probabilities. pi=\n");
    for (int k = 0; k < K; k++) {
        printf("%.3lf (log: %.3lf)\t", EXP(log_pi[k][t]), log_pi[k][t]);
    }
    printf("\n\n");
}

void print_info(int tp, int fp, int fn, int obv_var_cnt, int obv_invar_cnt, int nobv_cnt, int nobv_var_cnt, int nobv_invar_cnt, int disc_cnt, double prec, double sens, char *gt_file_name, char *results_file_name, FILE *results_file) {
    /*
    To print all the final information. There are 4 possible outcomes:
        1.- Only the number and different types of positions (# of OBVIOUS VARIANT, # of OBVIOUS INVARIANT, etc.) are displayed only on-screen: if no "ground_truth_file" has
            been introduced and the "results_file" flag was not set in the initial command.
        2.- Only the number and different types of positions (# of OBVIOUS VARIANT, # of OBVIOUS INVARIANT, etc.) are displayed on-screen and written down in a "results_file"
            if: no "ground_truth_file" has been introduced and the "results_file" flag was set in the initial command.
        3.- Both the number and different types of positions (# of OBVIOUS VARIANT, # of OBVIOUS INVARIANT, etc.) and the precision & sensitivity computations are displayed
            only on-screen if: a "ground_truth_file" has been introduced and the "results_file" flag was not set in the initial command.
        4.- Both the number and different types of positions (# of OBVIOUS VARIANT, # of OBVIOUS INVARIANT, etc.) and the precision & sensitivity computations are displayed
            both on-screen and wirtten down in a "results_file" if: a "ground_truth_file" has been introduced and the "results_file" flag was set in the initial command.
    */
    printf("# of OBVIOUS VARIANT = %d\n", obv_var_cnt);
    printf("# of OBVIOUS INVARIANT = %d\n", obv_invar_cnt);
    printf("# of NOT OBVIOUS = %d\n", nobv_cnt);
    printf("  # of NOT OBVIOUS VARIANT = %d\n", nobv_var_cnt);
    printf("  # of NOT OBVIOUS INVARIANT = %d\n", nobv_invar_cnt);
    printf("# of discarded positions = %d\n", disc_cnt);
    printf("Sum = %d\n", (obv_var_cnt + obv_invar_cnt + nobv_var_cnt + nobv_invar_cnt + disc_cnt));  // To check the total makes sense with the incomming positions

    if (strcmp(gt_file_name, "-") != 0) {
        printf("# of TRUE POSITIVES = %d\n", tp);
        printf("# of FALSE POSITIVES = %d\n", fp);
        printf("# of TRUE NEGATIVES = %d\n", (obv_var_cnt + obv_invar_cnt + nobv_var_cnt + nobv_invar_cnt - tp - fp - fn));
        printf("# of FALSE NEGATIVES = %d\n", fn);
        printf("Precision = %lf%%\n", prec);
        printf("Sensitivity = %lf%%\n", sens);
    }

    if (strcmp(results_file_name, "-") != 0) {
        fprintf(results_file, "# of OBVIOUS VARIANT = %d\n", obv_var_cnt);
        fprintf(results_file, "# of OBVIOUS INVARIANT = %d\n", obv_invar_cnt);
        fprintf(results_file, "# of NOT OBVIOUS = %d\n", nobv_cnt);
        fprintf(results_file, "  # of NOT OBVIOUS VARIANT = %d\n", nobv_var_cnt);
        fprintf(results_file, "  # of NOT OBVIOUS INVARIANT = %d\n", nobv_invar_cnt);
        fprintf(results_file, "# of discarded positions = %d\n", disc_cnt);
        fprintf(results_file, "Sum = %d\n", (obv_var_cnt + obv_invar_cnt + nobv_var_cnt + nobv_invar_cnt + disc_cnt));  // To check the total makes sense with the incomming positions
    }

    if ((strcmp(gt_file_name, "-") != 0) && (strcmp(results_file_name, "-") != 0)) {
        fprintf(results_file, "# of TRUE POSITIVES = %d\n", tp);
        fprintf(results_file, "# of FALSE POSITIVES = %d\n", fp);
        fprintf(results_file, "# of TRUE NEGATIVES = %d\n", (obv_var_cnt + obv_invar_cnt + nobv_var_cnt + nobv_invar_cnt - tp - fp - fn));
        fprintf(results_file, "# of FALSE NEGATIVES = %d\n", fn);
        fprintf(results_file, "Precision = %lf%%\n", prec);
        fprintf(results_file, "Sensitivity = %lf%%\n", sens);
    }
}

inline int offset_3d(int x, int y, int z) {
    /*
    Since all memory units are basically just an amorphous chunk of memory space, a way to visualize how data is organized within the 's' memory unit
    is by picturing a tensor. Memory that will be dynamically modified has been defined with calloc and realloc, thus allocating a chunk of space
    accessible like a vector/array: s[k]. This function allows to access a dynamically modified memory unit as if it had been defined staticly:
        s[n][m][l] = s[offset_3d(n,m,l)]
    With 'x', 'y' and 'z' we specify the x-axis position, the y-axis position and the z-axis position.

    Size X = N
    Size Y = M
    Size Z = L
    */
    return (x * M * L) + (y * L) + z;
}

inline int offset_2d(int x, int y, int col) {
    /*
    Since all memory units are basically just an amorphous chunk of memory space, a way to visualize how data is organized within the 'nobv_indxs'
    and 'obv_var_indxs' memory unit is by picturing a 3-column matrix. Memory that will be dynamically modified has been defined with calloc and realloc,
    thus allocating a chunk of space accessible like a vector/array: nobv_indxs[k]. This function allows to access a dynamically modified memory
    unit as if it had been defined staticly:
        nobv_indxs[n][i] = nobv_indxs[offset_2d(n,i,3)]
        obv_var_indxs[n][i] = obv_var_indxs[offset_2d(n,i,3)]
    With 'col', we are specifying the number of columns de working matrix has (both in 'obv_var_indxs' and 'nobv_indxs' col=3), and with 'x' and 'y' we
    select the line and column of interest.

    Size X
    Size Y = col
    */
    return (x * col + y);
}