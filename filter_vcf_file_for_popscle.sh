#!/bin/bash
#
# Copyright (C): 2020 - Gert Hulselmans
#
# Purpose: Functions for filtering VCF files for usage with popscle by removing mutations which are not informative.



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



get_contig_order_from_bam () {
    local bam_input_file="${1}";
    local output_type="${2}";

    if [ ${#@} -ne 2 ] ; then
        printf 'Usage: get_contig_order_from_bam BAM_file output_type\n\n';
        printf 'Arguments:\n';
        printf '  - BAM_file: BAM file from which to get the contig order and contig lengths.\n';
        printf '  - output_type:\n';
        printf '      - "names":        Return contig names.\n';
        printf '      - "chrom_sizes":  Return contig names and contig lengths.\n';
        printf '      - "vcf":          Return VCF header section for contigs.\n\n';
        return 1;
    fi

    case "${output_type}" in
        'names')
            ;;
        'chrom_sizes')
            ;;
        'vcf')
            ;;
        *)
            printf 'Error: output_type "%s" is not supported.\n' "${output_type}";
            return 1;
            ;;
    esac

    # Get the order of the contigs from the BAM header.
    samtools view -H "${bam_input_file}" \
      | awk \
            -F '\t' \
            -v output_type="${output_type}" \
            '
            {
                # Only look at sequence header fields.
                if ($1 == "@SQ") {
                    contig_idx += 1;
                    contig_name = "";
                    contig_length = "";

                    # Extract contig (chromosome) name and contig (chromosome) length.
                    for (i = 2; i <= NF; i++) {
                        if ($i ~ /^SN:/) {
                            contig_name = substr($i, 4);
                        }

                        if ($i ~ /^LN:/) {
                            contig_length = substr($i, 4);
                        }

                        # Create contig order to name and contig order to length and vcf contig appings.
                        contig_idx_to_name[contig_idx] = contig_name;
                        contig_idx_to_length[contig_idx] = contig_length;
                        contig_idx_to_vcf_contig[contig_idx] = sprintf("##contig=<ID=%s,length=%s>", contig_name, contig_length);
                    }
                }
            } END {
                if (output_type == "names") {
                    contig_names = "";

                    for (contig_idx = 1; contig_idx <= length(contig_idx_to_name); contig_idx++) {
                        contig_names = contig_names " " contig_idx_to_name[contig_idx];
                    }

                    # Print all contig names (without leading space).
                    print substr(contig_names, 2);
                } else if (output_type == "chrom_sizes") {
                    # Print all contig names with their length in a TAB separated fashion.
                    for (contig_idx = 1; contig_idx <= length(contig_idx_to_name); contig_idx++) {
                        print contig_idx_to_name[contig_idx] "\t" contig_idx_to_length[contig_idx];
                    }
                } else if (output_type == "vcf") {
                    # Print VCF header section for contigs.
                    for (contig_idx = 1; contig_idx <= length(contig_idx_to_vcf_contig); contig_idx++) {
                        print contig_idx_to_vcf_contig[contig_idx];
                    }
                }
            }'

      check_exit_codes;

      return $?;
}



sort_vcf_same_as_bam () {
    local bam_input_file="${1}";
    local vcf_input_file="${2}";

    if [ ${#@} -ne 2 ] ; then
        printf 'Usage: sort_vcf_same_as_bam BAM_file VCF_file\n\n';
        printf 'Arguments:\n';
        printf '  - BAM_file: BAM file from which to get the contig order to sort the VCF file.\n';
        printf '  - VCF_file: VCF file to sort by contig order as defined in the BAM file.\n\n';
        printf 'Purpose:\n';
        printf '  Sort VCF file in the same order as the BAM file, so it can be used with popscle.\n\n.';
        return 1;
    fi

    # Sort VCF file by same chromosome order as BAM file.
    bcftools reheader \
        -h <(
              # Create new VCF header:
              #   - Get VCF header of VCF input file.
              #   - Remove all contig header lines and "#CHROM" line from the VCF header.
              #   - Append contig headers in the order they appear in the input BAM file.
              #   - Add "#CHROM" line as last line of the new VCF header.
              bcftools view -h "${vcf_input_file}" \
                | awk \
                    '
                    {
                        if ($1 !~ /^##contig=/ && $1 !~ /^#CHROM/) {
                            # Remove all contig header lines and "#CHROM" line.
                            print $0;
                        }
                    }' \
                      | cat \
                            - \
                            <(get_contig_order_from_bam "${bam_input_file}" 'vcf') \
                            <(bcftools view -h "${vcf_input_file}" | tail -n 1)
            ) \
        "${vcf_input_file}" \
      | bcftools sort;

    check_exit_codes;

    return $?;
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



filter_out_mutations_homozygous_in_all_samples () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out mutations which are homozygous in all samples.
    bcftools view --exclude 'AC=AN' "${vcf_input_file}";

    return $?;
}



filter_out_mutations_not_unique_for_one_sample () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # Filter out mutations which are found homozygous in more than one sample.
    bcftools view --include 'AC=2' "${vcf_input_file}";

    return $?;
}



calculate_AF_AC_AN_values_based_on_genotype_info () {
    # VCF input file to use or stdin when no VCF input file is given.
    local vcf_input_file="${1:-/dev/stdin}";

    # (Re)calculate AF, AC, AN values bases on the genotype info provided for each sample.
    bcftools plugin fill-tags "${vcf_input_file}" -- --tags 'AF,AC,AN';

    return $?;
}

