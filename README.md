# Helper tools for popscle

Collection of tools to make [popscle](https://github.com/statgen/popscle) easier to use.



## Filter BAM file for usage with popscle dsc-pileup

Filter BAM file for usage with popscle dsc-pileup by keeping reads:
  - which overlap with SNPs in the VCF file
  - and which have a cell barcode (default: "CB" tag) contained in the cell barcode list
Keeping only relevant reads for popscle dsc-pileup can speedup it up quite significantly
(depending on the reduction of the number of reads in the filtered BAM file vs original).

```
$ ./filter_bam_file_for_popscle_dsc_pileup.sh
Usage:   filter_bam_file_for_popscle_dsc_pileup input_bam_filename barcodes_tsv_filename vcf_filename output_bam_filename [barcode_tag]

Purpose: Filter BAM file for usage with popscle dsc-pileup by keeping reads:
           - which overlap with SNPs in the VCF file
           - and which have a cell barcode (default: "CB" tag) contained in the cell barcode list
         Keeping only relevant reads for popscle dsc-pileup can speedup it up quite significantly
         (depending on the reduction of the number of reads in the filtered BAM file vs original).

```

###  Example

```bash
# Create filtered BAM with only the reads dsc-pileup needs.
./filter_bam_file_for_popscle_dsc_pileup.sh \
    ./samples_to_demultiplex/outs/possorted_genome_bam.bam \
    ./samples_to_demultiplex/outs/filtered_feature_bc_matrix/barcodes.tsv \
    samples.vcf \
    /tmp/samples_to_demultiplex.filter_bam_file_for_popscle_dsc_pileup.bam

# Use filtered BAM file for dsc-pileup.
popscle dsc-pileup \
    --sam /tmp/samples_to_demultiplex.filter_bam_file_for_popscle_dsc_pileup.bam \
    --vcf samples.vcf \
    --group-list ./samples_to_demultiplex/outs/filtered_feature_bc_matrix/barcodes.tsv \
    --out samples_to_demultiplex.pileup
```



## Sort VCF file in the same order as the BAM file

Sort VCF file in the same order as the BAM file so the following `popscle dsc-pileup`
error can be solved easily:

```
[E:%s] Your VCF/BCF files and SAM/BAM/CRAM files have different ordering of chromosomes. SAM/BAM/CRAM file has %s before %s, but VCF/BCF file has %s after %s"
```

```
$ ./sort_vcf_same_as_bam.sh
Usage: sort_vcf_same_as_bam BAM_file VCF_file [VCF_type]

Arguments:
  - BAM_file: BAM file from which to get the contig order to sort the VCF file.
  - VCF_file: VCF file to sort by contig order as defined in the BAM file.
  - VCF_type: VCF ouput file type (default: same as input VCF file type):
              v: uncompressed VCF, z: compressed VCF,
              u: uncompressed BCF, b: compressed BCF

Purpose:
  Sort VCF file in the same order as the BAM file, so it can be used with popscle.
```

### Examples

```bash
# Sort VCF file in the same order as the BAM file, so it can be used with popscle.
./sort_vcf_same_as_bam.sh \
    ./samples_to_demultiplex/outs/possorted_genome_bam.bam \
    samples.vcf \
  > /tmp/samples.sorted_as_in_bam.vcf

# Sort gzipped VCF file in the same order as the BAM file and write compressed VCF file.
./sort_vcf_same_as_bam.sh \
    ./samples_to_demultiplex/outs/possorted_genome_bam.bam \
    samples.vcf.gz \
  > /tmp/samples.sorted_as_in_bam.vcf.gz

# Sort gzipped VCF file in the same order as the BAM file and write uncompressed VCF file.
./sort_vcf_same_as_bam.sh \
    ./samples_to_demultiplex/outs/possorted_genome_bam.bam \
    samples.vcfi.gz \
    v \
  > /tmp/samples.sorted_as_in_bam.vcf
```



## Create filtered VCF files.

[BCFtools](https://www.htslib.org) can be used for filtering VCF files.

Looking at the [BCFtools filtering expressions manual](https://www.htslib.org/doc/bcftools.html#expressions)
gives an idea how to create your own filters for mutations.



### Import functions

Import functions in current shell.

```bash
# Import functions.
source filter_vcf_file_for_popscle.sh

# Check if all needed programs are installed.
check_if_programs_exists
```



### Get number of samples in BCF/VCF file.

```bash
get_number_of_samples_in_vcf [VCF_file]
```

Example:

```
$ get_number_of_samples_in_vcf DGRP2.source_NCSU.dm6.final.bcf
205
```



### Get samples names in BCF/VCF file.

Get all the sample names available in the BCF/VCF file (after the `FORMAT` column).

```
get_samples_names_in_vcf [VCF_file]
```

Example:

```
$ get_samples_names_in_vcf DGRP2.source_NCSU.dm6.final.bcf | head
DGRP-021
DGRP-026
DGRP-028
DGRP-031
DGRP-032
DGRP-038
DGRP-040
DGRP-041
DGRP-042
DGRP-045
```



### Subset samples from BCF/VCF file.

Extract only certain samples from a VCF file with multiple samples.

```bash
subset_samples_from_vcf comma_separated_samples_names [VCF_file]
```

Example:

```
subset_samples_from_vcf DGRP-032,DGRP-026,DGRP-042 DGRP2.source_NCSU.dm6.final.bcf | get_samples_names_in_vcf
DGRP-032
DGRP-026
DGRP-042
```



### Only keep SNPs from BCF/VCF file.

Only keep SNPs from VCF file (filter out INDELs and other mutations).

```bash
only_keep_snps [VCF_file]
```



### Filter out mutations missing genotype info for one or more samples.

Filter out mutations missing genotype (`./.`) info for one or more samples.

For those mutations no info is available if the sample has the reference and/or mutations,
so it might be better to skip this mutation in `popscle dsc-pileup`.

```bash
filter_out_mutations_missing_genotype_for_one_or_more_samples [VCF_file]
```



### Filter out mutations heterozygous for one or more samples.

Filter out mutations that are heterozygous for one or more samples.

This can be useful to reduce the number of mutations for `popscle dsc-pileup` when working with inbred lines
(all mutations are supposed to be homozygous).
In combination with `filter_out_mutations_not_unique_for_one_sample`, the number of mutations can be reduced
even further.

```bash
filter_out_mutations_heterozygous_for_one_or_more_samples [VCF_file]
```



### Filter out mutations homozygous reference in all samples.

Filter out mutations that have homozygous reference calls in all samples.

If the mutation position contains the reference for both alleles in all samples,
the mutation is not informative and can be skipped for `popscle dsc-pileup`.

```bash
filter_out_mutations_homozygous_reference_in_all_samples [VCF_file]
```



### Filter out mutations heterozygous in all samples.

Filter out mutations that are heterozygous in all samples.

If all samples are inbred lines, you might want to remove all non-homozygous SNPs.

```bash
filter_out_mutations_heterozygous_in_all_samples [VCF_file]
```



### Filter out mutations homozygous in all samples.

Filter out mutations that are homozygous in all samples.

If the mutation position contains the mutation for both alleles in all samples,
the mutation is not informative and can be skipped for `popscle dsc-pileup`.

```bash
filter_out_mutations_homozygous_in_all_samples [VCF_file]
```



### Only keep mutations heterozygous or homozygous in one sample.

Only keep mutations (heterozygous/homozygous) which are found only in
one sample, but not at all (heterozygous/homozygous) in other samples.

```bash
only_keep_mutations_heterozygous_or_homozygous_in_one_sample [VCF_file]
```



### Calculate allele frequency, allele count and total number of alleles.

Calculate allele frequency (`AF`), allele count (`AC`) and total number of alleles (`AN`) from genotype info of each sample.

This will add `AF`, `AC` and `AN` info fields or update those fields based on the genotype info of each sample in case they
were set incorrectly.

It is recommended to run this function before running:
  - `filter_out_mutations_homozygous_reference_in_all_samples`: needs correct value for `AC`.
  - `filter_out_mutations_homozygous_in_all_samples`: needs correct value for `AC` and `AN`.
  - `only_keep_mutations_homozygous_in_one_sample`: needs correct value for `AC`.
  - `only_keep_mutations_heterozygous_in_one_sample`: needs correct value for `AC`.
  - `only_keep_mutations_heterozygous_or_homozygous_in_one_sample`: needs correct value for `AC`.

Running it after `subset_samples_from_vcf` is also recommended as that function only updates 'AF' but not `AC` and `AN`.

```bash
calculate_AF_AC_AN_values_based_on_genotype_info [VCF_file]
```



### Examples

Create (minimal) VCF file for `popscle dsc-pileup` for 3 inbread lines (homozygous genotype SNPs are very common):
  - Only keep mutations for 3 selected samples
  - Only keep SNPs
  - (Re)calculate allele frequency (`AF`), allele count (`AC`), total number of alleles (`AN`).
  - Remove all SNPs which are missing genotype information for at least one sample (not useful to call those positions in `popscle dsc-pileup`).
  - Remove all SNPs which are homozygous reference in all samples (not useful to call those positions in `popscle dsc-pileup`).
  - Remove all SNPs which are homozygous in all samples (not useful to call those positions in `popscle dsc-pileup`).
  - Remove all SNPs which are heterozygous in at least one sample (those mutations shouldn't exist in inbred lines).
  - Only keep SNPs heterozygous or homozygous in one sample (but as heterozygous mutations are already filtered out, only keep homozygous ones).

```bash
subset_samples_from_vcf DGRP-032,DGRP-026,DGRP-042 DGRP2.source_BCM-HGSC.dm6.final.bcf \
  | only_keep_snps \
  | calculate_AF_AC_AN_values_based_on_genotype_info \
  | filter_out_mutations_missing_genotype_for_one_or_more_samples \
  | filter_out_mutations_homozygous_reference_in_all_samples \
  | filter_out_mutations_homozygous_in_all_samples \
  | filter_out_mutations_heterozygous_for_one_or_more_samples \
  | only_keep_mutations_heterozygous_or_homozygous_in_one_sample \
  > output.vcf
```

Create (minimal) VCF file for `popscle dsc-pileup` (heterozygous genotype SNPs are very common):
  - Only keep SNPs.
  - (Re)calculate allele frequency (`AF`), allele count (`AC`), total number of alleles (`AN`).
  - Remove all SNPs which are missing genotype information for at least one sample (not useful to call those positions in `popscle dsc-pileup`).
  - Remove all SNPs which are homozygous reference in all samples (not useful to call those positions in `popscle dsc-pileup`).
  - Remove all SNPs which are homozygous in all samples (not useful to call those positions in `popscle dsc-pileup`).

```bash
only_keep_snps input.vcf \
  | calculate_AF_AC_AN_values_based_on_genotype_info \
  | filter_out_mutations_missing_genotype_for_one_or_more_samples \
  | filter_out_mutations_homozygous_reference_in_all_samples \
  | filter_out_mutations_homozygous_in_all_samples \
  > output.vcf
```

