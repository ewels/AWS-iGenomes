# CLAUDE.md

Guidance for Claude Code (and other AI assistants) working in this repository.

## What this repo is

**AWS-iGenomes** distributes the [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)
reference genomes (plus extra tool indices) as an uncompressed, openly accessible S3
bucket, so cloud bioinformatics pipelines can pull just the files they need instead of
rebuilding references each run. It is an [AWS Open Data](https://registry.opendata.aws/aws-igenomes/)
resource.

- **Bucket:** `s3://ngi-igenomes/` (region `eu-west-1`, ~9 TB), open/public read.
- The repo itself is mostly **tooling and documentation**, not the data.

## Key files

| File | Purpose |
|---|---|
| `aws-igenomes.sh` | Bash helper that builds/runs `aws s3 sync` commands to fetch references. |
| `index.html` | Web-based command/path builder (published via GitHub Pages). |
| `ngi-igenomes_file_manifest.txt` | Full file listing (regenerate from the live bucket with `aws s3 ls s3://ngi-igenomes/igenomes/ --recursive`). |
| `nextflow.config` | iGenomes path config for use in Nextflow / nf-core pipelines. |
| `README.md` | User-facing documentation. |

## Path layout in the bucket

References follow the iGenomes convention:
`igenomes/<Species>/<Source>/<Build>/Sequence|Annotation/<IndexType>/...`
e.g. `igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/`. The manifest file is the
source of truth for what exists; prefer reading it over guessing paths.

## Private analysis (not committed)

Egress / cost / download analysis of the bucket lives entirely under **`.claude/analysis/`**,
which is **gitignored and for private review only** — it is never committed or published.
It contains the AWS billing report, its methodology (account IDs, profiles, CLI commands),
and the CloudWatch metrics-filter tooling. If you are asked to work on any of that, look in
`.claude/analysis/` (start with `egress-report-methodology.md`). Do not surface its contents
in committed files.

## Conventions

- **Never commit anything under `.claude/`** (analysis, methodology, AWS specifics) or any
  generated screenshots / scratch data. Public repo files must contain no account IDs,
  credentials, profiles, or internal infrastructure detail.
- When editing the file manifest or path logic, keep `aws-igenomes.sh`, `index.html` and
  `nextflow.config` consistent with each other.
- This is a public, widely-used community resource — be conservative with changes to the
  sync script and documented paths.
