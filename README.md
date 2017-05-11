# AWS-iGenomes

**Common reference genomes hosted on AWS S3**

## Introduction
In NGS bioinformatics, a typical analysis run involves aligning raw DNA sequencing reads against a known reference genome. A different reference is needed for every species, and many species have several references to choose from. Each tool then builds its own indices against these references. As such, one analysis run typically requires a number of different files. For example: raw underlying DNA sequence, annotation (GTF files) and index file for use the chosen alignment tool.

These files are quite large and take time to generate. Downloading and building them for each AWS run often takes a significant of the total run time and resources, which is very wasteful. To help with this, we have created an AWS S3 bucket containing the [illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.htm) references, with a few additional indices for a extra tools on top of this base dataset. The iGeomes initiative aims to collect and standardise a number of common species, references and tool indices.

This data is hosted in an S3 bucket (~5TB) and crucially is uncompressed (unlike the `.tar.gz` files held on the illumina iGenomes FTP servers). AWS runs can by pull just the required files to their local file storage before running.  This has the advantage of being faster, cheaper and more reproducible.

## Costs, billing and authentication
The S3 bucket is set to use the _Requester Pays_ policy. This means that our account won't be charged if lots of other people use the resource. Unfortunately, this means that access has to be limited to authenticated requests only. Usually this shouldn't be a problem, and full read-access is granted to any authenticated AWS user.

Note that if you are running in the same region as this S3 bucket (`eu-west`, Ireland) then there should be no data transfer fees and the resource should be free to use. From the [EC2 FAQ](https://aws.amazon.com/ec2/faqs/):

> There is no Data Transfer charge between two Amazon Web Services within the same region (i.e. between Amazon EC2 US West and another AWS service in the US West). Data transferred between AWS services in different regions will be charged as Internet Data Transfer on both sides of the transfer.

## Instructions
### Bucket details
The details of the S3 bucket are as follows:

* Bucket Name: `ngi-igenomes`
* Bucket ARN: `arn:aws:s3:::ngi-igenomes`
* Region: _EU (Ireland)_

### Description of Files
A full list of available files can be seen in [`ngi-igenomes_file_manifest.txt`](ngi-igenomes_file_manifest.txt).

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

For an example of this in action, see our [NGI-RNAseq pipeline](https://github.com/SciLifeLab/NGI-RNAseq/). The `aws` profile [config](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/conf/aws.config#L61-L193) contains s3 paths and our [regular HPC config](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/conf/uppmax.config#L113-L245) contains comparable regular file paths. This allows us to run the pipeline on either our HPC system or AWS with the same command and no extra setup.

## The Future
We are currently in discussion with the Open Data team at Amazon about making this into a Public Data resource. Until this happens, Amazon have been kind enough to provide us with a grant to cover the expenses of hosting this data on S3 for one year (April 2017 until April 2018).

Any updates to this arrangement will be posted here. If you have any questions please get in touch with [Phil Ewels](http://phil.ewels.co.uk) (phil.ewels@scilifelab.se, [@ewels](https://github.com/ewels)) or create an issue on this repository.


## Credits
The [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) resource was created by illumina. All credit for the collection and standardisation of this data should go to them!

This S3 resource was set up and documented by Phil Ewels ([@ewels](https://github.com/ewels)). The additional references not found in the base iGenomes resource were created with the help of Wesley Schaal ([@wschaal](https://github.com/wschaal)) - a system administrator at [UPPMAX](https://www.uppmax.uu.se/) (Uppsala Multidisciplinary Center for Advanced Computational Science).

The resource was initially developed for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

---

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
