#!/bin/bash

# First - check we have what we need
if command -v s3cmd &>/dev/null; then
    S3CMD="s3cmd"
elif command -v aws &>/dev/null; then
    S3CMD="aws s3"
else
    echo "Error: cannot find AWS commands 's3cmd' or 'aws s3'. These are required."
    exit 1
fi

# Get requirements from command line or prompt
if $1; then
    $GENOME=$1
else
    read $GENOME
fi

if $2; then
    $SOURCE=$2
else
    read $SOURCE
fi

if $3; then
    $BUILD=$3
else
    read $BUILD
fi

if $4; then
    $TYPE=$4
else
    read $TYPE
fi

# Get path suffix from $TYPE
case $TYPE in
    gtf)      $SUFFIX='Annotation/Genes/genes.gtf';;
    bed12)    $SUFFIX='Annotation/Genes/genes.gtf';;
    bismark)  $SUFFIX='Sequence/BismarkIndex/';;
    bowtie2)  $SUFFIX='Sequence/Bowtie2Index/';;
    bowtie1)  $SUFFIX='Sequence/BowtieIndex/';;
    bwa)      $SUFFIX='Sequence/BWAIndex/';;
    star)     $SUFFIX='Sequence/STARIndex/';;
    fasta)    $SUFFIX='Sequence/WholeGenomeFasta/';;
    *) echo "Error - download type not recognised! [gtf|bed12|bismark|bowtie2|bowtie1|bwa|star|fasta]"; exit 1;;
esac

# Try to download some genomes!
`$S3CMD --region eu-west-1 sync s3://ngi-igenomes/igenomes/$GENOME/$SOURCE/$BUILD/$SUFFIX references/$GENOME/$SOURCE/$BUILD/$SUFFIX`
