#!/bin/bash

function print_usage {
    echo -e "USAGE: aws-igenomes.sh\n" \
        "\t\t[-g <genome name>]\n" \
        "\t\t[-s <source name>]\n" \
        "\t\t[-b <build name>]\n" \
        "\t\t[-t <reference type>]\n" \
        "\t\t[-o <output directory>]\n" \
        "\t\t[-d (dry run, no downloads)]\n" \
        "\t\t[-q (quiet mode, non-interactive)]\n" \
        "\t\t[-h (usage help)]\n\n" \
        "All command line flags are optional. If not specified\n" \
        "and not running in quiet-mode, the script will prompt\n"\
        "for input and show available options.\n\n" \
        "Please note that this script requires the AWS command\n" \
        "line tools to be installed and configured for authenticated\n" \
        "access.\n" >&2
}

# Build arrays of available genomes, sources and builds
BUILD_OPTS=(
    "Arabidopsis_thaliana/Ensembl/TAIR10::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Arabidopsis_thaliana/Ensembl/TAIR9::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Arabidopsis_thaliana/NCBI/build9.1::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Arabidopsis_thaliana/NCBI/TAIR10::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bacillus_cereus_ATCC_10987/NCBI/2004-02-13::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Bacillus_subtilis_168/Ensembl/EB2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Bos_taurus/Ensembl/Btau_4.0::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/Ensembl/UMD3.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Bos_taurus/NCBI/Btau_4.2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/NCBI/Btau_4.6.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/NCBI/UMD_3.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/NCBI/UMD_3.1.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/UCSC/bosTau4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/UCSC/bosTau6::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/UCSC/bosTau7::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Bos_taurus/UCSC/bosTau8::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/Ensembl/WBcel215::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/Ensembl/WBcel235::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/Ensembl/WS210::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/Ensembl/WS220::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/NCBI/WS190::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/NCBI/WS195::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/UCSC/ce10::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Caenorhabditis_elegans/UCSC/ce6::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Canis_familiaris/Ensembl/BROADD2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Canis_familiaris/Ensembl/CanFam3.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Canis_familiaris/NCBI/build2.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Canis_familiaris/NCBI/build3.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Canis_familiaris/UCSC/canFam2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Canis_familiaris/UCSC/canFam3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Danio_rerio/Ensembl/GRCz10::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Danio_rerio/Ensembl/Zv9::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Danio_rerio/NCBI/GRCz10::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Danio_rerio/NCBI/Zv9::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Danio_rerio/UCSC/danRer10::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Danio_rerio/UCSC/danRer7::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/Ensembl/BDGP5::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/Ensembl/BDGP5.25::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/Ensembl/BDGP6::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/NCBI/build4.1::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/NCBI/build5::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/NCBI/build5.3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/NCBI/build5.41::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/UCSC/dm3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Drosophila_melanogaster/UCSC/dm6::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Enterobacteriophage_lambda/NCBI/1993-04-28::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Equus_caballus/Ensembl/EquCab2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Equus_caballus/NCBI/EquCab2.0::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Equus_caballus/UCSC/equCab2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Escherichia_coli_K_12_DH10B/Ensembl/EB1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Escherichia_coli_K_12_DH10B/NCBI/2008-03-17::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Escherichia_coli_K_12_MG1655/NCBI/2001-10-15::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Gallus_gallus/Ensembl/Galgal4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Gallus_gallus/Ensembl/WASHUC2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Gallus_gallus/NCBI/build2.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Gallus_gallus/NCBI/build3.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Gallus_gallus/UCSC/galGal3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Gallus_gallus/UCSC/galGal4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Glycine_max/Ensembl/Gm01::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Homo_sapiens/Ensembl/GRCh37::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Homo_sapiens/GATK/b37::gatk"
    "Homo_sapiens/GATK/hg19::gatk"
    "Homo_sapiens/GATK/hg38::gatk"
    "Homo_sapiens/GATK/GRCh37::gatk"
    "Homo_sapiens/GATK/GRCh38::gatk"
    "Homo_sapiens/NCBI/build36.3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Homo_sapiens/NCBI/build37.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Homo_sapiens/NCBI/build37.2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Homo_sapiens/NCBI/GRCh38::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Homo_sapiens/NCBI/GRCh38Decoy::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Homo_sapiens/UCSC/hg18::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Homo_sapiens/UCSC/hg19::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Homo_sapiens/UCSC/hg38::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Macaca_mulatta/Ensembl/Mmul_1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Mus_musculus/Ensembl/GRCm38::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Mus_musculus/Ensembl/NCBIM37::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Mus_musculus/NCBI/build37.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Mus_musculus/NCBI/build37.2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Mus_musculus/NCBI/GRCm38::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Mus_musculus/UCSC/mm10::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Mus_musculus/UCSC/mm9::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Mycobacterium_tuberculosis_H37RV/Ensembl/H37Rv.EB1::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Mycobacterium_tuberculosis_H37RV/NCBI/2001-09-07::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Oryza_sativa_japonica/Ensembl/IRGSP-1.0::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Oryza_sativa_japonica/Ensembl/MSU6::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Pan_troglodytes/Ensembl/CHIMP2.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Pan_troglodytes/Ensembl/CHIMP2.1.4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Pan_troglodytes/NCBI/build2.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Pan_troglodytes/NCBI/build3.1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Pan_troglodytes/UCSC/panTro2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Pan_troglodytes/UCSC/panTro3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Pan_troglodytes/UCSC/panTro4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "PhiX/Illumina/RTA::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "PhiX/NCBI/1993-04-28::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Pseudomonas_aeruginosa_PAO1/NCBI/2000-09-13::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Rattus_norvegicus/Ensembl/RGSC3.4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Rattus_norvegicus/Ensembl/Rnor_5.0::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Rattus_norvegicus/Ensembl/Rnor_6.0::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Rattus_norvegicus/NCBI/RGSC_v3.4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Rattus_norvegicus/NCBI/Rnor_5.0::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Rattus_norvegicus/NCBI/Rnor_6.0::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Rattus_norvegicus/UCSC/rn4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Rattus_norvegicus/UCSC/rn5::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Rattus_norvegicus/UCSC/rn6::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Rhodobacter_sphaeroides_2.4.1/NCBI/2005-10-07::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Saccharomyces_cerevisiae/Ensembl/EF2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Saccharomyces_cerevisiae/Ensembl/EF3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Saccharomyces_cerevisiae/Ensembl/EF4::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Saccharomyces_cerevisiae/Ensembl/R64-1-1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Saccharomyces_cerevisiae/NCBI/build2.1::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Saccharomyces_cerevisiae/NCBI/build3.1::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Saccharomyces_cerevisiae/UCSC/sacCer2::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Saccharomyces_cerevisiae/UCSC/sacCer3::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Schizosaccharomyces_pombe/Ensembl/EF1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Schizosaccharomyces_pombe/Ensembl/EF2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Sorangium_cellulosum_So_ce_56/NCBI/2007-11-27::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Sorghum_bicolor/Ensembl/Sbi1::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Staphylococcus_aureus_NCTC_8325/NCBI/2006-02-13::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq"
    "Sus_scrofa/Ensembl/Sscrofa10.2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Sus_scrofa/Ensembl/Sscrofa9::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna::variation"
    "Sus_scrofa/NCBI/Sscrofa10::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Sus_scrofa/NCBI/Sscrofa10.2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Sus_scrofa/NCBI/Sscrofa9.2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Sus_scrofa/UCSC/susScr2::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Sus_scrofa/UCSC/susScr3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Zea_mays/Ensembl/AGPv2::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
    "Zea_mays/Ensembl/AGPv3::gtf::bed12::bismark::bowtie::bowtie2::bwa::star::fasta::chromosomes::abundantseq::smrna"
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

# Array of reference type suffixes
TYPE_SUFFIXES=(
    "gtf::Annotation/Genes/" # appended to command: genes.gtf
    "bed12::Annotation/Genes/" # appended to command: genes.bed
    "bismark::Sequence/BismarkIndex/"
    "bowtie::Sequence/BowtieIndex/"
    "bowtie2::Sequence/Bowtie2Index/"
    "bwa::Sequence/BWAIndex/"
    "star::Sequence/STARIndex/"
    "fasta::Sequence/WholeGenomeFasta/"
    "chromosomes::Sequence/Chromosomes/"
    "abundantseq::Sequence/AbundantSequences/"
    "smrna::Annotation/SmallRNA/"
    "variation::Annotation/Variation/"
    "gatk::"
)

# Command line flags
while getopts ":g:s:b:t:o:dqh" opt; do
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
        d)
            DRYRUN=1
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

# First - check we have required tools
if command -v s3cmd &>/dev/null; then
    S3CMD="s3cmd"
elif command -v aws &>/dev/null; then
    S3CMD="aws s3"
else
    echo "Error: cannot find AWS commands 's3cmd' or 'aws s3'" >&2
    echo "See http://docs.aws.amazon.com/cli/latest/userguide/installing.html for installation instructions." >&2
    if [ $DRYRUN ]; then
        echo "Continuing as a dry run..">&2
        S3CMD="aws s3"
    else
        exit 1
    fi
fi

if [ ! $QUIET ]; then echo "AWS-iGenomes s3 sync script ($(date))" >&2; fi

# Check that we can access the bucket (requires proper configuration)
TEST_CMD="$S3CMD --no-sign-request --region eu-west-1 ls s3://ngi-igenomes/igenomes/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.gtf > /dev/null 2>&1"
if [ ! $QUIET ]; then echo "Testing AWS S3 connection..." 1>&2; fi
eval $TEST_CMD
if [ $? != 0 ]; then
    echo -e "\nError: could not contact S3 bucket. Is AWS authentication set up?" 1>&2
    echo "See http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html for instructions." 1>&2
    if [ $DRYRUN ]; then
        echo "Continuing as a dry run..">&2
    else
        exit 1
    fi
fi

# Get reference genome
if [ ! $GENOME ]; then
    if [ $QUIET ]; then
        echo "No reference genome specified. Disable quiet mode (-q) to see options. Exiting." >&2
        exit 1
    fi
    echo "Please enter a reference genome: (leave blank to see options)" >&2
    read GENOME
fi
while [[ ! " ${GENOME_OPTS[@]} " =~ " ${GENOME} " ]]; do
    if [ $GENOME ]; then
        if [ $QUIET ]; then
            echo "Error - genome '$GENOME' not recognised! Disable quiet mode (-q) to see options. Exiting." >&2
            exit 1
        fi
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
        echo -e "\nOnly one option for source type, setting to $SOURCE \n" >&2
    else
        if [ $QUIET ]; then
            echo "No reference source specified. Disable quiet mode (-q) to see options. Exiting." >&2
            exit 1
        fi
        echo "Please enter a reference source: (leave blank to see options)" >&2
        read SOURCE
    fi
fi
while [[ ! " ${SOURCE_OPTS[@]} " =~ " ${GENOME}/${SOURCE} " ]]; do
    if [ $SOURCE ]; then
        if [ $QUIET ]; then
            echo "Error - source '$SOURCE' not recognised! Disable quiet mode (-q) to see options. Exiting." >&2
            exit 1
        fi
        echo -e "\nError - source '$SOURCE' not recognised!\n" >&2
    fi
    echo "Available options:" >&2
    for i in ${SOURCE_OPTS[@]}; do
        if [[ $i == "${GENOME}"* ]]; then
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
        if [[ $i == "${GENOME}/${SOURCE}/"* ]]; then
            NUM_OPTS=$((NUM_OPTS+1))
            LAST_MATCH=${i#${GENOME}/${SOURCE}/}
            LAST_MATCH="${LAST_MATCH%%::*}"
        fi
    done
    if [ $NUM_OPTS == 1 ]; then
        BUILD=$LAST_MATCH
        echo -e "\nOnly one option for build type, setting to $BUILD \n" >&2
    else
        if [ $QUIET ]; then
            echo "No reference build specified. Disable quiet mode (-q) to see options. Exiting." >&2
            exit 1
        fi
        echo "Please enter a reference build: (leave blank to see options)" >&2
        read BUILD
    fi
fi
while [[ ! " ${BUILD_OPTS[@]} " =~ "${GENOME}/${SOURCE}/${BUILD}::"* ]]; do
    if [ $BUILD ]; then
        if [ $QUIET ]; then
            echo "Error - build '$BUILD' not recognised! Disable quiet mode (-q) to see options. Exiting." >&2
            exit 1
        fi
        echo -e "\nError - build '$BUILD' not recognised!\n" >&2
    fi
    echo "Available options:" >&2
    for i in ${BUILD_OPTS[@]}; do
        if [[ $i == "${GENOME}/${SOURCE}/"* ]]; then
            BOPTS=${i#${GENOME}/${SOURCE}/}
            BOPTS="${BOPTS%%::*}"
            echo -e "\t - ${BOPTS}" >&2 ;
        fi
    done
    echo "Please enter a reference build:" >&2
    read BUILD
done

# Now that we know the build, we can find possible reference types
TYPE_OPTS=()
for i in ${BUILD_OPTS[@]}; do
    if [[ $i == "${GENOME}/${SOURCE}/${BUILD}"* ]]; then
        TYPE_OPTS=(${i//::/ })
        TYPE_OPTS=("${TYPE_OPTS[@]:1}") # Remove first element (path)
    fi
done

# Get reference type
if [ ! $TYPE ]; then
    if [ ${#TYPE_OPTS[@]} == 1 ]; then
        TYPE=${TYPE_OPTS[0]}
        echo -e "\nOnly one option for reference type, setting to $TYPE \n" >&2
    else
        if [ $QUIET ]; then
            echo "No reference type specified. Disable quiet mode (-q) to see options. Exiting." >&2
            exit 1
        fi
        echo "Please enter a reference type: (leave blank to see options)" >&2
        read TYPE
    fi
fi
while [[ ! " ${TYPE_OPTS[@]} " =~ " ${TYPE} " ]]; do
    if [ $TYPE ]; then
        if [ $QUIET ]; then
            echo "Error - type '$TYPE' not recognised! Disable quiet mode (-q) to see options. Exiting." >&2
            exit 1
        fi
        echo -e "\nError - reference type '$TYPE' not recognised!\n" >&2
    fi
    echo "Available options:" >&2
    for i in ${TYPE_OPTS[@]}; do
        echo -e "\t - ${i}" >&2 ;
    done
    echo "Please enter a reference type:" >&2
    read TYPE
done

REF_SUFFIX='unset'
for i in "${TYPE_SUFFIXES[@]}" ; do
    k="${i%%::*}"
    v="${i##*::}"
    if [ "$k" == "$TYPE" ]; then
        REF_SUFFIX=$v
    fi
done
if [ "$REF_SUFFIX" == 'unset' ]; then
    echo "Error, could not find reference suffix!" >&2
    exit 1;
fi
S3PATH="s3://ngi-igenomes/igenomes/${GENOME}/${SOURCE}/${BUILD}/${REF_SUFFIX}"

# Default save directory
if [ ! $OUTPUT_DIR ]; then
    OUTPUT_DIR="$(pwd)/references/${GENOME}/${SOURCE}/${BUILD}/${REF_SUFFIX}"
fi

CMD="$S3CMD --no-sign-request --region eu-west-1 sync ${S3PATH} ${OUTPUT_DIR}"

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
if [ $DRYRUN ]; then
    echo "Skipping download as running in dry-run mode." >&2
else
    eval $CMD
fi

if [ $? != 0 ]; then
    echo "Error: Command exited with a non-zero exit code - $?" 1>&2
fi

if [ ! $QUIET ]; then echo "AWS iGenomes script completed ($(date))" >&2; fi
