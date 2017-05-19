#!/bin/bash

# First - check we have required tools
if command -v s3cmd &>/dev/null; then
    S3CMD="s3cmd"
elif command -v aws &>/dev/null; then
    S3CMD="aws s3"
else
    echo "Error: cannot find AWS commands 's3cmd' or 'aws s3'" >&2
    echo "See http://docs.aws.amazon.com/cli/latest/userguide/installing.html for installation instructions." >&2
    exit 1
fi

function print_usage {
    echo -e "USAGE: aws-igenomes.sh\n" \
        "\t\t[-g <genome name>]\n" \
        "\t\t[-s <source name>]\n" \
        "\t\t[-b <build name>]\n" \
        "\t\t[-t <reference type>]\n" \
        "\t\t[-o <output directory>]\n" \
        "\t\t[-q (quiet mode)]\n" \
        "\t\t[-h (usage help)]\n\n" \
        "All command line flags are optional. If not specified,\n" \
        "the script will prompt for input.\n\n" \
        "Please note that this script requires the AWS command\n" \
        "line tools to be installed and configured for authenticated\n" \
        "access.\n" >&2
}

# Build arrays of available genomes, sources and builds
BUILD_OPTS=(
    "Arabidopsis_thaliana/Ensembl/TAIR10"
    "Arabidopsis_thaliana/Ensembl/TAIR9"
    "Arabidopsis_thaliana/NCBI/build9.1"
    "Arabidopsis_thaliana/NCBI/TAIR10"
    "Bacillus_cereus_ATCC_10987/NCBI/2004-02-13"
    "Bacillus_subtilis_168/Ensembl/EB2"
    "Bos_taurus/Ensembl/Btau_4.0"
    "Bos_taurus/Ensembl/UMD3.1"
    "Bos_taurus/NCBI/Btau_4.2"
    "Bos_taurus/NCBI/Btau_4.6.1"
    "Bos_taurus/NCBI/UMD_3.1"
    "Bos_taurus/NCBI/UMD_3.1.1"
    "Bos_taurus/UCSC/bosTau4"
    "Bos_taurus/UCSC/bosTau6"
    "Bos_taurus/UCSC/bosTau7"
    "Bos_taurus/UCSC/bosTau8"
    "Caenorhabditis_elegans/Ensembl/WBcel215"
    "Caenorhabditis_elegans/Ensembl/WBcel235"
    "Caenorhabditis_elegans/Ensembl/WS210"
    "Caenorhabditis_elegans/Ensembl/WS220"
    "Caenorhabditis_elegans/NCBI/WS190"
    "Caenorhabditis_elegans/NCBI/WS195"
    "Caenorhabditis_elegans/UCSC/ce10"
    "Caenorhabditis_elegans/UCSC/ce6"
    "Canis_familiaris/Ensembl/BROADD2"
    "Canis_familiaris/Ensembl/CanFam3.1"
    "Canis_familiaris/NCBI/build2.1"
    "Canis_familiaris/NCBI/build3.1"
    "Canis_familiaris/UCSC/canFam2"
    "Canis_familiaris/UCSC/canFam3"
    "Danio_rerio/Ensembl/GRCz10"
    "Danio_rerio/Ensembl/Zv9"
    "Danio_rerio/NCBI/GRCz10"
    "Danio_rerio/NCBI/Zv9"
    "Danio_rerio/UCSC/danRer10"
    "Danio_rerio/UCSC/danRer7"
    "Drosophila_melanogaster/Ensembl/BDGP5"
    "Drosophila_melanogaster/Ensembl/BDGP5.25"
    "Drosophila_melanogaster/Ensembl/BDGP6"
    "Drosophila_melanogaster/NCBI/build4.1"
    "Drosophila_melanogaster/NCBI/build5"
    "Drosophila_melanogaster/NCBI/build5.3"
    "Drosophila_melanogaster/NCBI/build5.41"
    "Drosophila_melanogaster/UCSC/dm3"
    "Drosophila_melanogaster/UCSC/dm6"
    "Enterobacteriophage_lambda/NCBI/1993-04-28"
    "Equus_caballus/Ensembl/EquCab2"
    "Equus_caballus/NCBI/EquCab2.0"
    "Equus_caballus/UCSC/equCab2"
    "Escherichia_coli_K_12_DH10B/Ensembl/EB1"
    "Escherichia_coli_K_12_DH10B/NCBI/2008-03-17"
    "Escherichia_coli_K_12_MG1655/NCBI/2001-10-15"
    "Gallus_gallus/Ensembl/Galgal4"
    "Gallus_gallus/Ensembl/WASHUC2"
    "Gallus_gallus/NCBI/build2.1"
    "Gallus_gallus/NCBI/build3.1"
    "Gallus_gallus/UCSC/galGal3"
    "Gallus_gallus/UCSC/galGal4"
    "Glycine_max/Ensembl/Gm01"
    "Homo_sapiens/Ensembl/GRCh37"
    "Homo_sapiens/NCBI/build36.3"
    "Homo_sapiens/NCBI/build37.1"
    "Homo_sapiens/NCBI/build37.2"
    "Homo_sapiens/NCBI/GRCh38"
    "Homo_sapiens/NCBI/GRCh38Decoy"
    "Homo_sapiens/UCSC/hg18"
    "Homo_sapiens/UCSC/hg19"
    "Homo_sapiens/UCSC/hg38"
    "Macaca_mulatta/Ensembl/Mmul_1"
    "Mus_musculus/Ensembl/GRCm38"
    "Mus_musculus/Ensembl/NCBIM37"
    "Mus_musculus/NCBI/build37.1"
    "Mus_musculus/NCBI/build37.2"
    "Mus_musculus/NCBI/GRCm38"
    "Mus_musculus/UCSC/mm10"
    "Mus_musculus/UCSC/mm9"
    "Mycobacterium_tuberculosis_H37RV/Ensembl/H37Rv.EB1"
    "Mycobacterium_tuberculosis_H37RV/NCBI/2001-09-07"
    "Oryza_sativa_japonica/Ensembl/IRGSP-1.0"
    "Oryza_sativa_japonica/Ensembl/MSU6"
    "Pan_troglodytes/Ensembl/CHIMP2.1"
    "Pan_troglodytes/Ensembl/CHIMP2.1.4"
    "Pan_troglodytes/NCBI/build2.1"
    "Pan_troglodytes/NCBI/build3.1"
    "Pan_troglodytes/UCSC/panTro2"
    "Pan_troglodytes/UCSC/panTro3"
    "Pan_troglodytes/UCSC/panTro4"
    "PhiX/Illumina/RTA"
    "PhiX/NCBI/1993-04-28"
    "Pseudomonas_aeruginosa_PAO1/NCBI/2000-09-13"
    "Rattus_norvegicus/Ensembl/RGSC3.4"
    "Rattus_norvegicus/Ensembl/Rnor_5.0"
    "Rattus_norvegicus/Ensembl/Rnor_6.0"
    "Rattus_norvegicus/NCBI/RGSC_v3.4"
    "Rattus_norvegicus/NCBI/Rnor_5.0"
    "Rattus_norvegicus/NCBI/Rnor_6.0"
    "Rattus_norvegicus/UCSC/rn4"
    "Rattus_norvegicus/UCSC/rn5"
    "Rattus_norvegicus/UCSC/rn6"
    "Rhodobacter_sphaeroides_2.4.1/NCBI/2005-10-07"
    "Saccharomyces_cerevisiae/Ensembl/EF2"
    "Saccharomyces_cerevisiae/Ensembl/EF3"
    "Saccharomyces_cerevisiae/Ensembl/EF4"
    "Saccharomyces_cerevisiae/Ensembl/R64-1-1"
    "Saccharomyces_cerevisiae/NCBI/build2.1"
    "Saccharomyces_cerevisiae/NCBI/build3.1"
    "Saccharomyces_cerevisiae/UCSC/sacCer2"
    "Saccharomyces_cerevisiae/UCSC/sacCer3"
    "Schizosaccharomyces_pombe/Ensembl/EF1"
    "Schizosaccharomyces_pombe/Ensembl/EF2"
    "Sorangium_cellulosum_So_ce_56/NCBI/2007-11-27"
    "Sorghum_bicolor/Ensembl/Sbi1"
    "Staphylococcus_aureus_NCTC_8325/NCBI/2006-02-13"
    "Sus_scrofa/Ensembl/Sscrofa10.2"
    "Sus_scrofa/Ensembl/Sscrofa9"
    "Sus_scrofa/NCBI/Sscrofa10"
    "Sus_scrofa/NCBI/Sscrofa10.2"
    "Sus_scrofa/NCBI/Sscrofa9.2"
    "Sus_scrofa/UCSC/susScr2"
    "Sus_scrofa/UCSC/susScr3"
    "Zea_mays/Ensembl/AGPv2"
    "Zea_mays/Ensembl/AGPv3"
)
GENOME_OPTS=()
SOURCE_OPTS=()
for i in ${BUILD_OPTS[@]}; do
    j=$(dirname $i)
    k=$(dirname $j)
    if [[ ! " ${GENOME_OPTS[@]} " =~ " $k " ]]; then
        GENOME_OPTS+=($k)
    fi
    if [[ ! " ${SOURCE_OPTS[@]} " =~ " $j " ]]; then
        SOURCE_OPTS+=($j)
    fi
done

# Array of optional reference types
TYPE_SUFFIXES=(
    "gtf::Annotation/Genes/" # appended to command: genes.gtf
    "bed12::Annotation/Genes/" # appended to command: genes.bed
    "bismark::Sequence/BismarkIndex/"
    "bowtie::Sequence/BowtieIndex/"
    "bowtie2::Sequence/Bowtie2Index/"
    "bwa::Sequence/BWAIndex/"
    "fasta::Sequence/WholeGenomeFasta/"
    "star::Sequence/STARIndex/"
    "chromosomes::Sequence/Chromosomes/"
    "smrna::Annotation/SmallRNA/"
    "variation::Annotation/Variation/"
)
TYPE_OPTS=()
for i in ${TYPE_SUFFIXES[@]}; do
    TYPE_OPTS+=(${i%%::*})
done

# Command line flags
while getopts ":g:s:b:t:o:qh" opt; do
    case $opt in
        g)
            GENOME=$OPTARG
            ;;
        s)
            SOURCE=$OPTARG
            ;;
        b)
            BUILD=$OPTARG
            ;;
        t)
            TYPE=$OPTARG
            ;;
        o)
            OUTPUT_DIR=$OPTARG
            ;;
        q)
            QUIET=1
            ;;
        h)
            print_usage
            exit
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            print_usage
            exit 1;
            ;;
    esac
done

if [ ! $QUIET ]; then echo "AWS-iGenomes s3 sync script ($(date))" >&2; fi

# Check that we can access the bucket (requires proper configuration)
TEST_CMD="$S3CMD --region eu-west-1 ls s3://ngi-igenomes/igenomes/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.gtf > /dev/null 2>&1"
if [ ! $QUIET ]; then echo "Testing AWS S3 connection..." 1>&2; fi
eval $TEST_CMD
if [ $? != 0 ]; then
    echo -e "\nError: could not contact S3 bucket. Is AWS authentication set up?" 1>&2
    echo "See http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html for instructions." 1>&2
    exit;
fi

# Get reference genome
if [ ! $GENOME ]; then
    echo "Please enter a reference genome: (leave blank to see options)" >&2
    read GENOME
fi
while [[ ! " ${GENOME_OPTS[@]} " =~ " ${GENOME} " ]]; do
    if [ $GENOME ]; then
        echo -e "\nError - genome '$GENOME' not recognised!\n" >&2
    fi
    echo "Available options:" >&2
    for i in ${GENOME_OPTS[@]}; do
        echo -e "\t - ${i}" >&2 ;
    done
    echo "Please enter a reference genome:" >&2
    read GENOME
done

# Get reference source
if [ ! $SOURCE ]; then
    NUM_OPTS=0
    LAST_MATCH=
    for i in ${SOURCE_OPTS[@]}; do
        if [[ $i == *"${GENOME}"* ]]; then
            NUM_OPTS=$((NUM_OPTS+1))
            LAST_MATCH=${i#${GENOME}/}
        fi
    done
    if [ $NUM_OPTS == 1 ]; then
        SOURCE=$LAST_MATCH
    else
        echo "Please enter a reference source: (leave blank to see options)" >&2
        read SOURCE
    fi
fi
while [[ ! " ${SOURCE_OPTS[@]} " =~ " ${GENOME}/${SOURCE} " ]]; do
    if [ $SOURCE ]; then
        echo -e "\nError - source '$SOURCE' not recognised!\n" >&2
    fi
    echo "Available options:" >&2
    for i in ${SOURCE_OPTS[@]}; do
        if [[ $i == *"${GENOME}"* ]]; then
            echo -e "\t - ${i#${GENOME}/}" >&2 ;
        fi
    done
    echo "Please enter a reference source:" >&2
    read SOURCE
done

# Get reference build
if [ ! $BUILD ]; then
    NUM_OPTS=0
    LAST_MATCH=
    for i in ${BUILD_OPTS[@]}; do
        if [[ $i == *"${GENOME}/${SOURCE}/"* ]]; then
            NUM_OPTS=$((NUM_OPTS+1))
            LAST_MATCH=${i#${GENOME}/${SOURCE}/}
        fi
    done
    if [ $NUM_OPTS == 1 ]; then
        BUILD=$LAST_MATCH
    else
        echo "Please enter a reference build: (leave blank to see options)" >&2
        read BUILD
    fi
fi
while [[ ! " ${BUILD_OPTS[@]} " =~ " ${GENOME}/${SOURCE}/${BUILD} " ]]; do
    if [ $BUILD ]; then
        echo -e "\nError - build '$BUILD' not recognised!\n" >&2
    fi
    echo "Available options:" >&2
    for i in ${BUILD_OPTS[@]}; do
        if [[ $i == *"${GENOME}/${SOURCE}/"* ]]; then
            echo -e "\t - ${i#${GENOME}/${SOURCE}/}" >&2 ;
        fi
    done
    echo "Please enter a reference build:" >&2
    read BUILD
done

# Get reference type
if [ ! $TYPE ]; then
    echo "Please enter a reference type: (leave blank to see options)" >&2
    read TYPE
fi
while [[ ! " ${TYPE_OPTS[@]} " =~ " ${TYPE} " ]]; do
    if [ $TYPE ]; then
        echo -e "\nError - reference type '$TYPE' not recognised!\n" >&2
    fi
    echo "Available options:" >&2
    for i in ${TYPE_OPTS[@]}; do
        echo -e "\t - ${i}" >&2 ;
    done
    echo "Please enter a reference type:" >&2
    read TYPE
done

REF_SUFFIX=
for i in "${TYPE_SUFFIXES[@]}" ; do
    k="${i%%::*}"
    v="${i##*::}"
    if [ "$k" == "$TYPE" ]; then
        REF_SUFFIX=$v
    fi
done
if [ ! $REF_SUFFIX ]; then
    echo "Error, could not find reference suffix!" >&2
    exit 1;
fi
S3PATH="s3://ngi-igenomes/igenomes/${GENOME}/${SOURCE}/${BUILD}/$REF_SUFFIX"

# Default save directory
if [ ! $OUTPUT_DIR ]; then
    OUTPUT_DIR="$(pwd)/references/$GENOME/$SOURCE/$BUILD/$SUFFIX"
fi

CMD="$S3CMD --region eu-west-1 sync ${S3PATH} ${OUTPUT_DIR}"

# Get single files for GTF and BED
if [ "$TYPE" == "gtf" ]; then CMD+=' --exclude "*" --include "genes.gtf"'; fi
if [ "$TYPE" == "bed" ]; then CMD+=' --exclude "*" --include "genes.bed"'; fi

# Hide command output if in quiet mode
if [ $QUIET ]; then CMD+=' > /dev/null 2>&1'; fi

if [ ! $QUIET ]; then
    echo -e "Fetching references for:\n"\
        "\t        Genome: ${GENOME}\n"\
        "\t        Source: ${SOURCE}\n"\
        "\t         Build: ${BUILD}\n"\
        "\t          Type: ${TYPE}\n"\
        "\t       S3 Path: ${S3PATH}\n"\
        "\tSync Directory: ${OUTPUT_DIR}\n"\
        "\t  Full Command: ${CMD}\n" >&2
fi

# Try to download some genomes!
eval $CMD

if [ $? != 0 ]; then
    echo "Error: Command exited with a non-zero exit code - $?" 1>&2
fi

if [ ! $QUIET ]; then echo "AWS iGenomes script completed ($(date))" >&2; fi
