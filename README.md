# AWS-iGenomes

### Common reference genomes hosted on AWS S3
#### Download script & command builder: [https://ewels.github.io/AWS-iGenomes/](https://ewels.github.io/AWS-iGenomes/)

![Amazon Web Services](images/AWS_logo.png)

## Introduction
In NGS bioinformatics, a typical analysis run involves aligning raw DNA sequencing reads against a known reference genome. A different reference is needed for every species, and many species have several references to choose from. Each tool then builds its own indices against these references. As such, one analysis run typically requires a number of different files. For example: raw underlying DNA sequence, annotation (GTF files) and index file for use the chosen alignment tool.

These files are quite large and take time to generate. Downloading and building them for each AWS run often takes a significant of the total run time and resources, which is very wasteful. To help with this, we have created an AWS S3 bucket containing the [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) references, with a few additional indices for a extra tools on top of this base dataset. The iGeomes initiative aims to collect and standardise a number of common species, references and tool indices.

This data is hosted in an S3 bucket (~5TB) and crucially is uncompressed (unlike the `.tar.gz` files held on the illumina iGenomes FTP servers). AWS runs can by pull just the required files to their local file storage before running.  This has the advantage of being faster, cheaper and more reproducible.

## Download Script
To make usage easier, this repository contains a script ([`aws-igenomes.sh`](aws-igenomes.sh)) which can sync the AWS-iGenomes for you. It requires the AWS command line tools to be installed and configured with authentication. Required references can be supplied on the command line or given through prompts when running the script.

This repository is hosted using GitHub pages, so the script can be run in a single command as follows:

```bash
curl -fsSL https://ewels.github.io/AWS-iGenomes/aws-igenomes.sh | bash
```

For more details, see [https://ewels.github.io/AWS-iGenomes/](https://ewels.github.io/AWS-iGenomes/)

### Command Builder
If you'd prefer to just get a sync command for the files you need, you can use the web-based command builder that's available at [https://ewels.github.io/AWS-iGenomes/](https://ewels.github.io/AWS-iGenomes/)

## Instructions
### Bucket details
The details of the S3 bucket are as follows:

* Bucket Name: `ngi-igenomes`
* Bucket ARN: `arn:aws:s3:::ngi-igenomes`
* Region: _EU (Ireland)_

### Description of Files
A full list of available files can be seen in [`ngi-igenomes_file_manifest.txt`](ngi-igenomes_file_manifest.txt)

The following species have reference builds available:

* _Arabidopsis thaliana_
* _Bacillus cereus_ ATCC 10987
* _Bacillus subtilis_ 168
* _Bos taurus_
* _Caenorhabditis elegans_
* _Canis familiaris_
* _Danio rerio_
* _**Drosophila melanogaster**_
* _Enterobacteriophage lambda_
* _Equus caballus_
* _Escherichia coli_ K 12 DH10B
* _Escherichia coli_ K 12 MG1655
* _Gallus gallus_
* _Glycine max_
* _**Homo sapiens**_
* _Macaca mulatta_
* _**Mus musculus**_
* _Mycobacterium tuberculosis_ H37RV
* _Oryza sativa japonica_
* _Pan troglodytes_
* _PhiX_
* _Pseudomonas aeruginosa_ PAO1
* _Rattus norvegicus_
* _Rhodobacter sphaeroides_ 2.4.1
* _**Saccharomyces cerevisiae**_
* _Schizosaccharomyces pombe_
* _Sorangium cellulosum_ So ce 56
* _Sorghum bicolor_
* _Staphylococcus aureus_ NCTC 8325
* _Sus scrofa_
* _Zea mays_

Most of these species then have references from multiple sources and builds. For example, _Mus musculus_
has the following:

* Ensembl
    * `GRCm38`, `NCBIM37`
* NCBI
    * `build37.1`, `build37.2`, `GRCm38`
* UCSC
    * `mm10`, `mm9`

Within each reference build, the following resources are typically available (with a few exceptions):

* Gene annotation in `GTF` and `BED` format
* Sequence `FASTA` files:
    * Whole genome files
    * Separate chromosomes
    * Abundant sequences
* Alignment indices for the following tools:
    * [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
    * [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
    * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    * [BWA](http://bio-bwa.sourceforge.net/)
    * [STAR](https://github.com/alexdobin/STAR)
* For some genomes:
    * smRNA ([miRBase](http://www.mirbase.org/))
    * Variation

An additional special-case is the GATK bundles, available for _Homo sapiens_ (`b37`, `hg19`, `hg38`, `GRCh37` and `GRCh38`).

See [Data origin](#data-origin) below for more details of how these files were generated.

### Costs, billing and authentication
The S3 bucket is currently set to be completely open access (there were problems with the previous _Requester Pays_
policy). This will remain the case until the credits awarded to fund this project from Amazon run out or expire
(hopefully stable for some time yet).

Note that if if possible, it's best for us if you run in the same region as this S3 bucket (`eu-west`, Ireland).
Then there should be no data transfer fees and the resource should stay around for longer.
From the [EC2 FAQ](https://aws.amazon.com/ec2/faqs/):

> There is no Data Transfer charge between two Amazon Web Services within the same region (i.e. between Amazon EC2 US West and another AWS service in the US West). Data transferred between AWS services in different regions will be charged as Internet Data Transfer on both sides of the transfer.

### Basic Usage
How you use this resource largely depends on how you're using AWS. Very generally however, you can retrieve your required data by using the [AWS Command Line Interface](https://aws.amazon.com/cli/).

For example, using the `aws sync` command:

```bash
aws s3 sync s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/ ./my_refs/
```

If the `aws` tool isn't installed, probably the easiest way to get it is using `pip`:

```bash
pip install --upgrade --user awscli
```

Remember that you must configure the tool with some kind of AWS authentication to access the contents of the s3 bucket.

For more information and help, see the [AWS CLI user guide](http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-set-up.html).

### Usage with Nextflow
[Nextflow](https://www.nextflow.io/) is a powerful workflow manager allowing the creation of bioinformatics analysis pipelines. It was created to help the transition from traditional academic HPC systems to cloud computing. As such, it has extensive built-in support for a number of AWS features. One such feature is native integration with s3. This means that you can specify paths to required reference files in your pipeline which are stored in s3 and Nextflow will automatically retrieve them.

The repository contains an example Nextflow config file containing common paths and a suggested usage example: [`nextflow.config`](nextflow.config)

For an example of this in action, see our [NGI-RNAseq pipeline](https://github.com/SciLifeLab/NGI-RNAseq/). The `aws` profile [config](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/conf/aws.config#L61-L193) contains s3 paths and our [regular HPC config](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/conf/uppmax.config#L113-L245) contains comparable regular file paths. This allows us to run the pipeline on either our HPC system or AWS with the same command and no extra setup.

## Data origin
This resource is based on the [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) references. These were downloaded and unpacked in April 2016.

After unpacking, references were added for [STAR](https://github.com/alexdobin/STAR), [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) and [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). A new reference directory was contained for each reference and the index built (see commands below).

**A full list of available files can be seen in this repository:  [`ngi-igenomes_file_manifest.txt`](ngi-igenomes_file_manifest.txt)**

#### STAR
```bash
module load star/2.5.1b

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ./ --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 100
```
(if no GTF file available, `--sjdbGTFfile genes.gtf --sjdbOverhang 100` was not specified).

#### Bismark
```bash
module load bowtie/1.1.2
module load bowtie2/2.2.6
module load bismark/0.14.5

bismark_genome_preparation ./
bismark_genome_preparation --bowtie2 ./
```

#### BED12
BED12 files were generated using the [`gtf2bed`](https://github.com/ExpressionAnalysis/ea-utils/blob/master/clipper/gtf2bed) tool from [ea-utils](https://expressionanalysis.github.io/ea-utils/).

#### GATK Bundles
The [GATK Resource Bundles](https://software.broadinstitute.org/gatk/download/bundle) for builds
`b37`, `hg19` and `hg38` were downloaded from the Broad FTP server on 2017-05-25. For more information
about their contents, please see
[this article](https://software.broadinstitute.org/gatk/documentation/article.php?id=1213).

Please note that `b37/CEUTrio.HiSeq.WGS.b37.NA12878.bam` and associated files are not included.
This file is ~355GB and with the FTP download limiting from Broad it was going to take nearly
a year to transfer.

#### More details
STAR, Bismark and BED12 additions were kindly done by the UPPMAX team. Full details and exactly scripts used for this can be found at [github.com/UPPMAX/bio-data](https://github.com/UPPMAX/bio-data/tree/master/igenomes/index).


## The Future
We are currently in discussion with the Open Data team at Amazon about making this into a Public Data resource. Until this happens, Amazon have been kind enough to provide us with a grant to cover the expenses of hosting this data on S3 for one year (April 2017 until April 2018).

Any updates to this arrangement will be posted here. If you have any questions please get in touch with [Phil Ewels](http://phil.ewels.co.uk) (phil.ewels@scilifelab.se, [@ewels](https://github.com/ewels)) or create an issue on this repository.


## Changelog

#### Version v0.3 (dev)

* Made a web interface for generating aws s3 sync commands _(not everyone likes random command line scripts..)_
* Now that Amazon are taking the cost of the hosting, everything is fully public
    * Added `--no-sign-request` to the commands so that they work without authentication
* Added new GRCh37 and GRCh38 builds for GATK
    * Different to the existing hg18 and hg19 builds only in that the file organisation is cleaner and consistent with the rest of iGenomes (old builds left for backwards-compatibility)
    * Contain new indexes for BWA. More to be added in the future.

#### [Version v0.2](https://github.com/ewels/AWS-iGenomes/releases/tag/v0.2) - 2016-05-25

* Added GATK bundles `b37`, `hg19` and `hg38` from the Broad FTP download
* Minor download script updates

#### [Version v0.1](https://github.com/ewels/AWS-iGenomes/releases/tag/v0.1) - 2016-05-23
Initial released. Repository created with file-list of the iGenomes resource, with added
BED12, STAR and Bismark indices. Download bash script written and basic website created
at https://ewels.github.io/AWS-iGenomes/

## Credits
The [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) resource was created by illumina. All credit for the collection and standardisation of this data should go to them!

This S3 resource was set up and documented by Phil Ewels ([@ewels](https://github.com/ewels)). The additional references not found in the base iGenomes resource were created with the help of Wesley Schaal ([@wschaal](https://github.com/wschaal)) - a system administrator at [UPPMAX](https://www.uppmax.uu.se/) (Uppsala Multidisciplinary Center for Advanced Computational Science).

The resource was initially developed for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
