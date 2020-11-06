#!/bin/bash
#
# Copyright (C): 2020 - Gert Hulselmans
#
# Purpose: Functions for filtering VCF files for usage with popscle by removing mutations which are not informative.
#
#
# BCFtools filtering expressions manual:
#   https://www.htslib.org/doc/bcftools.html#expressions



# Function to check if any of the programs in a pipe failed.
check_exit_codes () {
    local GET_PIPESTATUS="${PIPESTATUS[@]}";
    local exit_code;

    for exit_code in ${GET_PIPESTATUS} ; do
        if [ ${exit_code} -ne 0 ] ; then
             return ${exit_code};
        fi
    done

    return 0;
}



# Check if necessary programs are installed.
check_if_programs_exists () {
    local exit_code=0;

    # Check if awk is installed.
    if ! type awk > /dev/null 2>&1 ; then
        printf 'Error: "awk" could not be found in PATH.\n';
        exit_code=2;
    fi

    # Check if bcftools is installed.
    if ! type bcftools > /dev/null 2>&1 ; then
        printf 'Error: "bcftools" could not be found in PATH.\n';
        exit_code=2;
    fi

    return ${exit_code};
}



get_number_of_samples_in_vcf () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Only look at the VCF header and count the number of samples in the "#CHROM" line.
    bcftools view -h "${vcf_input_file}" \
        | awk \
            -F '\t' \
            '
            {
                if ( $1 == "#CHROM" ) {
                    if ( NF > 9 ) {
                        nbr_samples = NF - 9;
                        print nbr_samples;
                    } else {
                        print "0";
                    }
                }
            }'

    check_exit_codes;

    return $?;
}



get_samples_names_in_vcf () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Only look at the VCF header and print the sample names listed in the "#CHROM" line.
    bcftools view -h "${vcf_input_file}" \
        | awk \
            -F '\t' \
            -v "vcf_input_file=${vcf_input_file}" \
            '
            {
                if ( $1 == "#CHROM" ) {
                    if  ( NF > 9 ) {
                        # Print all sample names.
                        for (sample_column_idx=10 ; sample_column_idx <= NF; sample_column_idx++) {
                            print $sample_column_idx;
                        }

                        exit(0);
                    } else {
                        printf("Error: No sample names found in VCF file \"%s\".\n", vcf_input_file);

                        exit(1);
                    }
                }
            }'

    check_exit_codes;

    return $?;
}



subset_samples_from_vcf () {
    # Comma separated list of samples to extract from VCF file.
    local samples="${1}";

    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${2:-/dev/stdin}";
    
    if [ ${#@} -lt 1 ] ; then
        printf 'Usage: subset_samples_from_vcf comma_separated_samples_names [VCF_file]\n';
        return 1;
    fi

    # Extract specific samples from VCF file.
    bcftools view --samples "${samples}" "${vcf_input_file}";

    return $?;
}



only_keep_snps () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out all non SNPs mutations.
    bcftools view --types 'snps' "${vcf_input_file}";

    return $?;
}



filter_out_mutations_missing_genotype_for_one_or_more_samples () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out mutations which have missing genotypes for one or more samples
    # as those mutations are not very informative.
    bcftools view --genotype '^miss' "${vcf_input_file}";

    return $?;
}



filter_out_mutations_heterozygous_for_one_or_more_samples () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out mutations which are heterozygous for one or more samples.
    bcftools view --genotype '^het' "${vcf_input_file}";

    return $?;
}



filter_out_mutations_homozygous_reference_in_all_samples () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out mutation which are homozygous reference in all samples.
    bcftools view --exclude 'AC=0' "${vcf_input_file}";

    return $?;
}



filter_out_mutations_heterozygous_in_all_samples () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out mutation which are heterozygous in all samples.
    bcftools view --exclude 'COUNT(GT="het")=N_SAMPLES' "${vcf_input_file}";

    return $?;
}



filter_out_mutations_homozygous_in_all_samples () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out mutations which are homozygous in all samples.
    bcftools view \
        --exclude 'COUNT(GT="AA") = N_SAMPLES' \
        "${vcf_input_file}";

    return $?;
}



only_keep_mutations_homozygous_in_one_sample () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Only keep mutations (homozygous) which are found only in one sample,
    # but not at all (heterozygous/homozygous) in other samples.
    #bcftools view --include 'AC=2 && ( GT = "1|1" | GT = "1/1")' "${vcf_input_file}";
    bcftools view \
        --include 'COUNT(GT="AA") = 1 && COUNT(GT="RR") = (N_SAMPLES - 1)' \
        "${vcf_input_file}";

    return $?;
}



only_keep_mutations_heterozygous_or_homozygous_in_one_sample () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Only keep mutations (heterozygous/homozygous) which are found only in
    # one sample, but not at all (heterozygous/homozygous) in other samples.
    bcftools view \
        --include '( COUNT(GT="AA") = 1 || COUNT(GT="AR") = 1 ) && COUNT(GT="RR") = (N_SAMPLES - 1)' \
        --include '( AC=2 && ( GT = "1|1" | GT = "1/1") || AC=1 && ( GT = "0|1" | GT = "1|0" | GT = "0/1" | GT = "0/1") )' \
        "${vcf_input_file}";

    return $?;
}



calculate_AF_AC_AN_values_based_on_genotype_info () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # (Re)calculate AF, AC, AN values bases on the genotype info provided for each sample.
    bcftools plugin fill-tags "${vcf_input_file}" -- --tags 'AF,AC,AN';

    return $?;
}

