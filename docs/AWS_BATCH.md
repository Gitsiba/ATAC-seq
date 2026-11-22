# Run the ATAC-seq workflow on AWS Batch

This guide describes the cloud execution path for the existing paired-end Nextflow workflow. It is an infrastructure starter guide, not a claim that an AWS deployment has already been provisioned or tested.

## Architecture

- Amazon S3 stores input FASTQ files, the reference FASTA, Nextflow work files, and published results.
- Amazon ECR stores the pinned ATAC-seq tool container.
- AWS Batch runs each Nextflow process as a containerised job.
- IAM roles grant the Nextflow launcher and Batch compute environment only the permissions they need.
- CloudWatch provides job and compute-environment logs.

## Prerequisites

1. An AWS account with billing alerts enabled.
2. AWS CLI v2 configured locally using an IAM Identity Center profile or another approved credential mechanism.
3. An S3 bucket in the selected region.
4. An ECR repository for the workflow image.
5. An AWS Batch managed compute environment, job queue, and job definition/permissions configured for the selected VPC.
6. A Nextflow installation and Java runtime on the machine that launches the workflow.
7. IAM permissions for the launcher to submit and monitor Batch jobs and access the S3 work directory. Batch compute roles must be able to pull the ECR image and read/write the required S3 objects.

Do not commit credentials, account IDs tied to private infrastructure, patient data, or unpublished sequencing data to Git.

## 1. Build and publish the container

From the repository root:

```bash
docker build -t atacseq-tools:0.1.0 -f containers/atacseq/Dockerfile .
```

Create an ECR repository in your chosen region, authenticate Docker using the AWS CLI, tag the image with your ECR repository URI, and push it. Use your own account ID and repository URI; do not copy example credentials into source files.

The image contains pinned versions of FastQC, fastp, Bowtie 2, SAMtools, MACS2, and MultiQC. Rebuild and retest the image when changing tool versions.

## 2. Configure the AWS Batch profile

The configuration in `conf/aws-batch.config` expects these runtime parameters:

- `--aws_region`: AWS region (default: `ap-south-1`).
- `--aws_batch_queue`: existing AWS Batch job queue name or ARN.
- `--aws_workdir`: S3 URI used for Nextflow work files.
- `--container_image`: ECR image URI and tag.

Example invocation (replace all placeholders):

```bash
nextflow run main.nf \
  -c conf/aws-batch.config \
  --samplesheet samplesheet.csv \
  --reference_fasta reference.fa \
  --genome_size hs \
  --outdir results \
  --aws_region ap-south-1 \
  --aws_batch_queue YOUR_BATCH_QUEUE \
  --aws_workdir s3://YOUR_BUCKET/atacseq/work \
  --container_image YOUR_ACCOUNT.dkr.ecr.ap-south-1.amazonaws.com/atacseq-tools:0.1.0
```

For remote execution, make sure the input files are accessible to the AWS Batch jobs. Stage FASTQ files and the reference FASTA in S3 and use S3 paths in the samplesheet and reference parameter when running in the cloud.

## 3. Security and cost controls

- Prefer IAM roles and short-lived credentials; do not use long-lived access keys in scripts.
- Restrict S3 access to the project bucket and required prefixes.
- Keep S3 buckets private and enable encryption.
- Set AWS Budgets and billing alerts before launching jobs.
- Start with a tiny synthetic dataset and a small compute environment.
- Set sensible vCPU, memory, retry, and job-concurrency limits.
- Review and remove unused Batch compute environments, job queues, ECR images, logs, and test objects when no longer needed.
- Configure S3 lifecycle policies for temporary work data only after confirming retention requirements.
- Never use real patient or confidential research data for a public demonstration.

## Important implementation notes

- The current workflow was originally authored with per-process Conda directives. The AWS Batch profile supplies a shared container image for cloud execution; validate the container runtime and process environment on a small test before production use.
- The workflow's input staging and output publishing must be checked with S3 paths in the target AWS account.
- AWS Batch networking, service roles, execution roles, job queue, and compute environment are account-level resources and are not created by this repository configuration.
- No AWS resources are created by CI. The GitHub Actions workflow remains a lightweight software validation workflow.
- This starter profile is not a validated clinical pipeline or a production deployment. Validate scientific outputs and cloud permissions before using research data.
