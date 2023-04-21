#!/usr/bin/env python3

from __future__ import annotations

import argparse
import glob
import io
import os.path
import shutil
import sys

import polars as pl

try:
    from isal import igzip as gzip_mod  # type: ignore[import]
except ImportError:
    import gzip as gzip_mod


def write_popscle_pileup_cel_full_filename(
    popscle_dsc_output_prefix: str,
    popscle_dsc_output_full_prefix: str,
) -> pl.DataFrame | None:
    popscle_pileup_cel_dfs = []

    popscle_pileup_cel_full_filename = f"{popscle_dsc_output_full_prefix}.pileup.cel.gz"

    if os.path.exists(popscle_pileup_cel_full_filename):
        print(
            f'Error: popscle pileup CEL full file "{popscle_pileup_cel_full_filename}" already exists.',
            file=sys.stderr,
        )
        return None

    for i, popscle_pileup_cel_part_filename in enumerate(
        sorted(glob.glob(f"{popscle_dsc_output_prefix}.*.pileup.cel.gz"))
    ):
        print(
            f'Reading partial popscle pileup CEL file "{popscle_pileup_cel_part_filename}" ...',
            file=sys.stderr,
        )

        popscle_pileup_cel_dfs.append(
            # Read partial popscle pileup CEL file.
            pl.read_csv(
                popscle_pileup_cel_part_filename,
                separator="\t",
                has_header=True,
                dtypes={
                    "#DROPLET_ID": pl.Int64,
                    "BARCODE": pl.Utf8,
                    "NUM.READ": pl.Int64,
                    "NUM.UMI": pl.Int64,
                    "NUM.UMIwSNP": pl.Int64,
                    "NUM.SNP": pl.Int64,
                },
            )
            .rename({"#DROPLET_ID": "DROPLET_ID_PARTITIONED"})
            .with_columns(
                # Add current partition as a column.
                pl.lit(i).alias("PARTITION")
            )
        )

    # Combine partial popscle pileup CEL files and add real "DROPLET_ID" for the full
    # dataset.
    popscle_pileup_cel_df = pl.concat(popscle_pileup_cel_dfs).with_row_count(
        name="DROPLET_ID",
        offset=0,
    )

    with gzip_mod.open(popscle_pileup_cel_full_filename, "w") as fh_full:
        print(
            f'Writing popscle pileup PLP full file "{popscle_pileup_cel_full_filename}" ...',
            file=sys.stderr,
        )

        # Create BytesIO object to temporarily write the corrected popscle pileup CEL
        # file to.
        bytes_io_tsv = io.BytesIO()

        # Remove "DROPLET_ID_PARTITIONED" and "PARTITION" columns before writing
        # corrected popscle pileup CEL file.
        popscle_pileup_cel_df.select(
            [
                pl.col("DROPLET_ID").alias("#DROPLET_ID"),
                pl.col("BARCODE"),
                pl.col("NUM.READ"),
                pl.col("NUM.UMI"),
                pl.col("NUM.UMIwSNP"),
                pl.col("NUM.SNP"),
            ],
        ).write_csv(
            bytes_io_tsv,
            has_header=True,
            separator="\t",
        )

        # Write BytesIO object with corrected popscle pileup CEL output
        # to full popscle pileup CEL file.
        fh_full.write(bytes_io_tsv.getbuffer())

    # Return corrected popscle pileup CEL output with "DROPLET_ID_PARTITIONED" and
    # "PARTITION" columns as Polars DataFrame.
    return popscle_pileup_cel_df


def write_popscle_pileup_plp_full_filename(
    popscle_pileup_cel_df: pl.DataFrame,
    popscle_dsc_output_prefix: str,
    popscle_dsc_output_full_prefix: str,
) -> bool:
    popscle_pileup_plp_full_filename = f"{popscle_dsc_output_full_prefix}.pileup.plp.gz"

    if os.path.exists(popscle_pileup_plp_full_filename):
        print(
            f'Error: popscle pileup PLP full file "{popscle_pileup_plp_full_filename}" already exists.',
            file=sys.stderr,
        )
        return False

    with gzip_mod.open(popscle_pileup_plp_full_filename, "w") as fh_full:
        print(
            f'Writing popscle pileup PLP full file "{popscle_pileup_plp_full_filename}" ...',
            file=sys.stderr,
        )

        for i, popscle_pileup_plp_part_filename in enumerate(
            sorted(glob.glob(f"{popscle_dsc_output_prefix}.*.pileup.plp.gz"))
        ):
            print(
                f'Reading partial popscle pileup PLP file "{popscle_pileup_plp_part_filename}" ...',
                file=sys.stderr,
            )

            # Create BytesIO object to temporarily write the corrected popscle pileup
            # PLP file to.
            bytes_io_tsv = io.BytesIO()

            (
                # Read partial popscle pileup PLP file.
                pl.read_csv(
                    popscle_pileup_plp_part_filename,
                    separator="\t",
                    has_header=True,
                    dtypes={
                        "#DROPLET_ID": pl.Int64,
                        "SNP_ID": pl.Int64,
                        "ALLELES": pl.Utf8,
                        "BASEQS": pl.Utf8,
                    },
                )
                .lazy()
                .rename({"#DROPLET_ID": "DROPLET_ID_PARTITIONED"})
                .with_columns(
                    # Add current partition as a column.
                    pl.lit(i).alias("PARTITION")
                )
                # Correct "DROPLET_ID" column from partial popscle pileup
                # PLP file with real "DROPLET_ID" for the full dataset.
                .join(
                    popscle_pileup_cel_df.lazy(),
                    on=["PARTITION", "DROPLET_ID_PARTITIONED"],
                    how="inner",
                )
                .select(
                    pl.col("DROPLET_ID").alias("#DROPLET_ID"),
                    pl.col("SNP_ID"),
                    pl.col("ALLELES"),
                    pl.col("BASEQS"),
                )
                .collect(streaming=True)
                .write_csv(
                    bytes_io_tsv,
                    # Write header only for first partial popscle pileup PLP file.
                    has_header=i == 0,
                    separator="\t",
                )
            )

            # Write/append BytesIO object with corrected popscle pileup PLP output
            # to full popscle pileup PLP file.
            fh_full.write(bytes_io_tsv.getbuffer())

    return True


def write_popscle_pileup_umi_full_filename(
    popscle_dsc_output_prefix: str,
    popscle_dsc_output_full_prefix: str,
) -> bool:
    line_count = 0
    popscle_pileup_umi_full_filename = f"{popscle_dsc_output_full_prefix}.pileup.umi.gz"

    if os.path.exists(popscle_pileup_umi_full_filename):
        print(
            f'Error: popscle pileup UMI full file "{popscle_pileup_umi_full_filename}" already exists.',
            file=sys.stderr,
        )
        return False

    with gzip_mod.open(popscle_pileup_umi_full_filename, "wt") as fh_full:
        print(
            f'Writing popscle pileup UMI full file "{popscle_pileup_umi_full_filename}" ...',
            file=sys.stderr,
        )

        for popscle_pileup_umi_part_filename in sorted(
            glob.glob(f"{popscle_dsc_output_prefix}.*.pileup.umi.gz")
        ):
            print(
                f'Reading partial popscle pileup UMI file "{popscle_pileup_umi_part_filename}" ...',
                file=sys.stderr,
            )
            with gzip_mod.open(popscle_pileup_umi_part_filename, "rt") as fh:
                for line in fh:
                    print(
                        str(line_count),
                        line.split("\t", 1)[1],
                        sep="\t",
                        end="",
                        file=fh_full,
                    )
                    line_count += 1

    return True


def write_popscle_pileup_var_full_filename(
    popscle_dsc_output_prefix: str,
    popscle_dsc_output_full_prefix: str,
) -> bool:
    popscle_pileup_var_full_filename = f"{popscle_dsc_output_full_prefix}.pileup.var.gz"
    popscle_pileup_var_most_complete_df = None
    popscle_pileup_var_most_complete_filename = None

    if os.path.exists(popscle_pileup_var_full_filename):
        print(
            f'Error: popscle pileup VAR full file "{popscle_pileup_var_full_filename}" already exists.',
            file=sys.stderr,
        )
        return False

    # Read each partial popscle dsc pileup var file and take the one that contains the
    # most mutations as it seems that popscle will always include all mutations with
    # the same SNP_ID.
    for popscle_pileup_var_part_filename in sorted(
        glob.glob(f"{popscle_dsc_output_prefix}.*.pileup.var.gz")
    ):
        print(
            f'Reading partial popscle pileup VAR file "{popscle_pileup_var_part_filename}" ...',
            file=sys.stderr,
        )

        popscle_pileup_var_df = pl.read_csv(
            popscle_pileup_var_part_filename,
            separator="\t",
            has_header=True,
            dtypes={
                "#SNP_ID": pl.Int64,
                "CHROM": pl.Utf8,
                "POS": pl.Int64,
                "REF": pl.Utf8,
                "ALT": pl.Utf8,
                "AF": pl.Utf8,
            },
        ).rename({"#SNP_ID": "SNP_ID"})

        if popscle_pileup_var_most_complete_df is None:
            popscle_pileup_var_most_complete_df = popscle_pileup_var_df
            popscle_pileup_var_most_complete_filename = popscle_pileup_var_part_filename
        elif popscle_pileup_var_df.height > popscle_pileup_var_most_complete_df.height:
            # Check if the last element of the smallest file is found exactly in the
            # bigger. If this assert fails, our assumption is wrong.
            assert (
                popscle_pileup_var_most_complete_df[
                    popscle_pileup_var_most_complete_df.height - 1
                ].to_struct("LAST_SNP")
                == popscle_pileup_var_df[
                    popscle_pileup_var_most_complete_df.height - 1
                ].to_struct("LAST_SNP")
            ).item(), (
                "SNP_IDs do not match between different popscle dsc pileup var files."
            )

            popscle_pileup_var_most_complete_df = popscle_pileup_var_df
            popscle_pileup_var_most_complete_filename = popscle_pileup_var_part_filename

    if popscle_pileup_var_most_complete_filename:
        print(
            f'Writing popscle pileup VAR full file "{popscle_pileup_var_full_filename}" ...',
            file=sys.stderr,
        )
        shutil.copy2(
            popscle_pileup_var_most_complete_filename, popscle_pileup_var_full_filename
        )

        return True

    return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge popscle dsc pileup outputs from multiple popscle dsc "
        "pileup runs with the same VCF file as input, but with a different list "
        "of cell barcodes."
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="popscle_dsc_output_prefix",
        action="store",
        type=str,
        required=True,
        help="popscle pileup dsc output prefix for partial popscle dsc pileup output.",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="popscle_dsc_output_full_prefix",
        action="store",
        type=str,
        required=True,
        help="popscle pileup dsc output prefix for full popscle dsc pileup output.",
    )

    args = parser.parse_args()

    popscle_pileup_cel_full_filename = (
        f"{args.popscle_dsc_output_full_prefix}.pileup.cel.gz"
    )
    popscle_pileup_plp_full_filename = (
        f"{args.popscle_dsc_output_full_prefix}.pileup.plp.gz"
    )
    popscle_pileup_umi_full_filename = (
        f"{args.popscle_dsc_output_full_prefix}.pileup.umi.gz"
    )
    popscle_pileup_var_full_filename = (
        f"{args.popscle_dsc_output_full_prefix}.pileup.var.gz"
    )

    output_exists = False

    if os.path.exists(popscle_pileup_cel_full_filename):
        print(
            f'Error: popscle pileup CEL full file "{popscle_pileup_cel_full_filename}" already exists.',
            file=sys.stderr,
        )
        output_exists = True

    if os.path.exists(popscle_pileup_plp_full_filename):
        print(
            f'Error: popscle pileup PLP full file "{popscle_pileup_plp_full_filename}" already exists.',
            file=sys.stderr,
        )
        output_exists = True

    if os.path.exists(popscle_pileup_umi_full_filename):
        print(
            f'Error: popscle pileup UMI full file "{popscle_pileup_umi_full_filename}" already exists.',
            file=sys.stderr,
        )
        output_exists = True

    if os.path.exists(popscle_pileup_var_full_filename):
        print(
            f'Error: popscle pileup VAR full file "{popscle_pileup_var_full_filename}" already exists.',
            file=sys.stderr,
        )
        output_exists = True

    if output_exists:
        sys.exit(1)

    popscle_pileup_cel_df = write_popscle_pileup_cel_full_filename(
        popscle_dsc_output_prefix=args.popscle_dsc_output_prefix,
        popscle_dsc_output_full_prefix=args.popscle_dsc_output_full_prefix,
    )

    if popscle_pileup_cel_df is None:
        sys.exit(1)

    write_popscle_pileup_plp_full_filename(
        popscle_pileup_cel_df=popscle_pileup_cel_df,
        popscle_dsc_output_prefix=args.popscle_dsc_output_prefix,
        popscle_dsc_output_full_prefix=args.popscle_dsc_output_full_prefix,
    )

    write_popscle_pileup_umi_full_filename(
        popscle_dsc_output_prefix=args.popscle_dsc_output_prefix,
        popscle_dsc_output_full_prefix=args.popscle_dsc_output_full_prefix,
    )
    write_popscle_pileup_var_full_filename(
        popscle_dsc_output_prefix=args.popscle_dsc_output_prefix,
        popscle_dsc_output_full_prefix=args.popscle_dsc_output_full_prefix,
    )


if __name__ == "__main__":
    main()
