
# MEpurity
# What is MEpurity?
Next-generation sequencing has revolutionized the study of cancer genomes. However, the tumor sample we obtained always consist of a mixture of normal and tumor tissue and the reads we get can be of multiple clonal types. MEpurity is a C++ program for calculating Tumor Purity using DNA methylation 450k data without matched normal control sample. 

For more details about MEpurity, please visit the page [Introduction of MEpurity](https://github.com/lbw1995/MEpurity).

# Install
## Requirement
To compile MEpurity you need two things: G++ (which usually are already installed on Linux) and Boost C++ Libraries. The last is not installed on Linux by default, but it can be downloaded from:
[Boostlib](https://www.boost.org/users/history/version_1_69_0.html)
## Install from source
Download the compressed source file MEpurity.tar.gz and do as follows:

    $ tar -xzvf MEpurity.tar.gz
    $ cd ./MEpurity/src
    $ g++ bmm.cpp file.cpp help.cpp kmeans.cpp -o MEpurity -I [path-to-boostlib]/include/
# Usage
    Version 0.1
    ./MEpurity [options]
## Required parameters:
<<<<<<< HEAD
<<<<<<< HEAD
    -f:     The Illumina Infinium Human Methylation 450K (450k) input data.
    -m:     The normal mapfile under path of this software.
    -o:     The output file path that you would like to contain the results.
## Optional parameters:
    -h      Show this help message and exit.
    -s      The number of CpG sites that you want to use in the map file. <int> (Default:80000)
=======
    -i:     The Illumina Infinium Human Methylation 450K (450k) input data.
    -p:     The file under the path of the software (parameters.txt) that stores the parameters about the distribution of beta values at different CpG sites in normal samples.
    -o:     The output file path that you would like to contain the results.
## Optional parameters:
    -h      Show this help message and exit.
    -s      The number of CpG sites that you want to use in the map file. <int> (Default:70000)
>>>>>>> 7d1da275950b6a6294c04af9a18aa76a4c296c7c
=======
    -i:     The Illumina Infinium Human Methylation 450K (450k) input data.
    -p:     The file under the path of the software (parameters.txt) that stores the parameters about the distribution of beta values at different CpG sites in normal samples.
    -o:     The output file path that you would like to contain the results.
## Optional parameters:
    -h      Show this help message and exit.
    -s      The number of CpG sites that you want to use in the map file. <int> (Default:70000)
>>>>>>> MEpurity/master
    -t      The maximum iteration time of bmm algorithm. <int> (Default:10000)
    -c      The least percemtage of sites belonging to a cluster that would not be filted. <float> (Default:0.01)
    -n      The original number of clusters. <int> (Default:10)
    -v      Output progress in terms of mixing coefficient (expected) values if 1. <bool> (Default:False)
## Example
<<<<<<< HEAD
<<<<<<< HEAD
    ./MEpurity -f ../test/test.txt -m ../map.txt -o ./output.txt
=======
    tar -zxvf ../test/test.tar.gz -C ../test
    ./MEpurity -i ../test/test.txt -p ../parameters.txt -o ./output.txt
>>>>>>> 7d1da275950b6a6294c04af9a18aa76a4c296c7c
=======
    tar -zxvf ../test/test.tar.gz -C ../test
    ./MEpurity -i ../test/test.txt -p ../parameters.txt -o ./output.txt
>>>>>>> MEpurity/master
# Output
The output file contains 2 columns:

1.Sample id

2.Tumor purity value
