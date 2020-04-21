# Negative_data
An adaptive method of defining negative mutation status for multi-sample comparison using next-generation sequencing

CONTACT: Lei.Wei@roswellpark.org

Copyright (c) 2020 by Roswell Park Cancer Institue. #########################################################

#Background
Cancer-related somatic mutations can be measured using next-generation sequencing (NGS). For a specific mutation in a given sample, a NGS test may yield three possible statuses: positive, negative, or unknown due to low covearge. We propose an adaptive mutation-specific negative (MSN) method to improve the classification of negative and unknown mutation statuses in the comparison of multiple "related" tumor samples of the same patient. This document provides information on how to run MSN to differentiate "unknown" from "negative" statuses in a number of "related" tumor samples from the same patient. This method also works on measurement of other samples containing tumor cells or tumor DNA, such as circulating tumor cells (CTCs) and circulating tumor DNA (ctDNA).

#Input data
Before running this current program, we assume the users have finished mutation calling using any somatic mutation caller. Then the mutations from all "related" samples of the same patient need to be combined into a list of unique mutations, and each mutation has been measured in all "related" samples to extract the numbers of mutant and wildtype reads using any appropriate program such as the SAMtools pileup function. The user also need to identify positive samples using any standard that they selected (e.x. by a minimum number of mutant reads, or a minimum VAF etc). The final input data to the MSN program is a data matrix containing the counts of mutant and wildtype reads of every mutation in all "related" samples. 

An example is provided in the provided file "input_example.txt". 

Example input:

    Mutation                sample1     sample2     sample3      sample4    sample5

    chr1.43285942.A.T       10,10,+     0,20,       4,5,+           0,8,         0,5,

    chr17.5378413.T.C       40,32,+     0,15,       17,31,+         0,20,       3,20,

    chr17.5378249.C.A       5,100,+     0,20,       0,30,           0,100,      0,150,

The first column is a list of mutant unique id (for example, "chr.position.reference allele.mutant allele"), then each of the rest columns represents a sample. For any mutation in a given sample, the following information is included: 1) the number of mutant reads; 2) the number of reference/wildtype reads; 3) whether the sample has a positive status for that mutation, as indicated by a plus sign "+". No sign is need for any non-positive samples. For example, "5,100,+" means there are 5 mutant reads,  100 wildtype reads, and the mutation is considered as positive in that sample.  Please note, all included mutations should be positive in at least one of the samples.

#Running the code
    ./MSN.pl -i Example_data/input_example.txt -o output.txt -p_cutoff 0.05
options:
-i input file path (required)
-o output file path (optional, default is input.'_ReEval_NegStatus'
-p the p value cutoff for rejecting the null hypothesis that a non-positive sample may contain the same frequency of mutant reads as one of the positive samples (optional, default is 0.05)

#Output

The program will output a file containing the updated mutation status:
* positive ("+"): the program will not modify any existing positive or generating any new positive status;

* negative("-"): re-assigned by the program;

* unknown("unknown"): re-assigned by the program. 

Example output using the above input:

    Mutation                sample1     sample2         sample3         sample4         sample5

    chr1.43285942.A.T       10,10,+     0,20,-          4,5,+           0,8,Unknown     0,5,Unknown

    chr17.5378413.T.C       40,32,+     0,15,-          17,31,+         0,20,-          3,20,Unknown

    chr17.5378249.C.A       5,100,+     0,20,Unknown    0,30,Unknown    0,100,Unknown   0,150,-

In the above example, the first mutation "chr1.43285942.A.T" was determined to be negative in sample2 (0,20,-), but unknown in sample4(0,8,Unknown) and sample5(0,5,Unknown). For samples with "unknown" status, additional coverage is recommended to determine the mutation's actual status in that sample. 


