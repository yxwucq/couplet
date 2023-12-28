#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
"""
The script resolves two SingleShot reads
Usage: python run_couplet.py  --fq1=TEST_R1.fq.gz  --fq2=TEST_R2.fq.gz --mismatch_threshold=0.05 --phred=prob

"""

import subprocess
import copy
import multiprocessing
import argparse
import logging
import re
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import sys
import os
import functools

from couplet.utils import (
    get_base_name,
    get_outfile_names,
    check_python_version,
    load_quality_table,
)
from couplet.export import (
    log_core_stats,
    log_additional_stats,
    log_and_plot_additional_stats_single,
)
from couplet.resolve import (
    resolve_phred_min,
    resolve_phred_prob,
    resolve_phred_with_qtable,
)
from couplet.rules import (
    five_bp_rule,
    five_bp_error_codes,
    five_bp_modifications,
    six_bp_rule,
    six_bp_error_codes,
    six_bp_modifications,
    ResolutionRule,
)
from couplet.core import resolve_read_pair, update_stats
from couplet.split_fastq_file import split_fastq_file, merge_fastq_files, shell_split_fastq_file, shell_merge_fastq_files
from couplet.postprocess_stats import postprocess_stats

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fq1",
        help="forward read names",
    )

    parser.add_argument(
        "--fq2",
        help="reverse read names",
    )
    parser.add_argument(
        "--out_dir",
        help="output directory",
        default=".",
    )
    parser.add_argument(
        "--num_cores",
        help="Number of cores to use",
        default=4,
        type=int,
    )

    parser.add_argument(
        "--mismatch-threshold",
        help="Enter a value between 0 and 1 for allowed fraction of mismatches",
        default=0.05,
        type=float,
    )
    parser.add_argument(
        "--phred",
        help="Enter a mode: min or prob",
        default="prob",
        choices=["min", "prob", "qtable"],
        type=str,
    )
    parser.add_argument(
        "--quality-table",
        help="Enter the path to a quality table",
        default=None,
        type=str,
    )

    parser.add_argument(
        "--rule",
        help="Enter a rule name: 5bp or 6bp",
        default="5bp",
        choices=["5bp", "6bp"],
        type=str,
    )

    parser.add_argument(
        "--gap-penalty",
        default=-2,
        help="Gap penalty parameter for NW alignment default is -2",
        type=float,
    )
    parser.add_argument(
        "--end-gap-penalty",
        default=-5,
        help="gap penalty for begining and end gaps  for NW alignment default is -5",
        type=float,
    )
    parser.add_argument(
        "--match-award",
        default=1,
        help="Match award parameter for NW alignment default is 1",
        type=float,
    )
    parser.add_argument(
        "--mismatch-penalty",
        default=-1,
        help="Mismatch  penalty parameter for NW alignment default is -1",
        type=float,
    )
    parser.add_argument(
        "--no-mismatch-aware-trimming",
        dest="mismatch_aware_trimming",
        action="store_false",
        help="Flag to turn off mismatch aware trimming.",
    )
    parser.add_argument(
        "--min-read-length",
        default=50,
        help="Minimum read length required to not discard a read.",
        type=float,
    )
    parser.add_argument(
        "--additional-stats",
        dest="additional_stats",
        action="store_true",
        help="Flag to turn on logging of additional statistics. This the default behaviour.",
    )
    parser.add_argument(
        "--no-additional-stats",
        dest="additional_stats",
        action="store_false",
        help="Flag to turn off logging of additional statistics.",
    )
    parser.add_argument(
        "--generate-plots",
        dest="generate_plots",
        action="store_true",
        help="Flag to indicate that plotting should take place here, rather than in post-processing.",
    )
    parser.add_argument(
        "--episeq-in-qname",
        dest="episeq_in_qname",
        action="store_true",
        help="Flag to export the epigenetic sequence as part of the fastq read name (or QNAME in .sam specification).",
    )
    parser.add_argument(
        "--orig-quals-in-sam-tag",
        dest="orig_quals_in_sam_tag",
        action="store_true",
        help="Flag to turn on storing of the original quality scores (Q1, Q2) in a SAM tag (XQ and YQ respectively).",
    )

    parser.add_argument(
        "--xe-tag",
        dest="xe_tag",
        action="store_true",
        help="Flag to turn on storing of modification information in a CEGX tag (XE).",
    )

    parser.set_defaults(
        episeq_in_qname=False,
        generate_plots=False,
        additional_stats=True,
        mismatch_aware_trimming=True,
        xe_tag=False,
        orig_quals_in_sam_tag=False,
    )

    args = parser.parse_args()
    if len(sys.argv) == 1:
        # display help message when no args are passed.
        parser.print_help()
        sys.exit()
    
    if not args.fq1.endswith('_R1.fq.gz'):
        raise ValueError("fq1 file name should end with _R1.fq.gz")
    if not args.fq2.endswith('_R2.fq.gz'):
        raise ValueError("fq2 file name should end with _R2.fq.gz")
    
    return args

def main():
    global_args = parse_args()
    # Check Python version is above 3.9.0
    check_python_version([3, 8, 0])
    
    logging.info("Splitting fastq files to " + str(global_args.num_cores*2) + " parts")
    
    split_fastq_R1 = shell_split_fastq_file(global_args.fq1)
    split_fastq_R2 = shell_split_fastq_file(global_args.fq2)
    print(split_fastq_R1)
    print(split_fastq_R2)
    assert len(split_fastq_R1) == len(split_fastq_R2)
    split_file_length = len(split_fastq_R1)
    
    # Run couplet on each split file
    pool = multiprocessing.Pool(processes=global_args.num_cores)
    arg_list = []
    for i in range(1, split_file_length+1):
        local_args = copy.deepcopy(global_args) # should not contain list or dict
        local_args.fq1 = global_args.fq1.replace('_R1.fq.gz', f"_SPLIT_{i}_R1.fq.gz")
        local_args.fq2 = global_args.fq2.replace('_R2.fq.gz', f"_SPLIT_{i}_R2.fq.gz")
        print("Pushing " + local_args.fq1 + " and " + local_args.fq2)
        arg_list.append(local_args)
    
    print(len(arg_list))
    pool.map(run_couplet, arg_list)
    
    pool.close()
    pool.join()
    
    resolve_file_list = [global_args.fq1.replace('_R1.fq.gz', f"_SPLIT_{i}_resolved.fq.gz") for i in range(1, split_file_length+1)]
    discard_r1_list = [global_args.fq1.replace('_R1.fq.gz', f"_SPLIT_{i}_discarded_R1.fq.gz") for i in range(1, split_file_length+1)]
    discard_r2_list = [global_args.fq1.replace('_R1.fq.gz', f"_SPLIT_{i}_discarded_R2.fq.gz") for i in range(1, split_file_length+1)]
    
    for file_list in [resolve_file_list, discard_r1_list, discard_r2_list]:
        for file in file_list:
            if not os.path.exists(file):
                file_list.remove(file)
        if len(file_list) == 0:
            raise ValueError("No files in " + str(file_list))
    
    for file_list in [resolve_file_list, discard_r1_list, discard_r2_list]:
        shell_merge_fastq_files(file_list, os.path.join(global_args.out_dir, file_list[0].split('/')[-1].replace('_SPLIT_1_', '_')))

    for temp_file in resolve_file_list + discard_r1_list + discard_r2_list:
        subprocess.run(["rm", temp_file])

    stats_file_list = [global_args.fq1.replace('_R1.fq.gz', f"_SPLIT_{i}_couplet.yaml") for i in range(1, split_file_length+1)]
    
    postprocess_stats(stats_file_list, None, global_args.out_dir, global_args.fq1.split('/')[-1].replace('_R1.fq.gz', ''))
    
    # remove all split files
    subprocess.run(f"rm {global_args.fq1.replace('_R1.fq.gz', '_SPLIT_*')}", shell=True)
    
def run_couplet(args):
    print("Running couplet on " + args.fq1 + " and " + args.fq2)
    # Setup logging
    log_output_file = get_base_name(args) + "_couplet.log"
    logging.basicConfig(
        filename=log_output_file,
        format="%(asctime)-15s %(filename)s:%(lineno)d %(message)s",
        force=True,
    )
    LOGGER = logging.getLogger("root")
    LOGGER.setLevel(logging.INFO)
    LOGGER.addHandler(logging.StreamHandler())

    # Log all input arguments
    LOGGER.info("Command line inputs:")
    for k, v in vars(args).items():
        LOGGER.info(k + ": " + str(v))

    # Choosing mode phred calc
    if args.phred == "min":
        resolve_phred_fn = resolve_phred_min
    elif args.phred == "qtable":
        if args.quality_table is None:
            LOGGER.error(
                "Requested 'qtable' mode for Phred score resolution, but no quality table was provided. Exiting."
            )
            sys.exit(1)
        elif os.path.exists(args.quality_table):
            LOGGER.info("Using quality table: " + args.quality_table)
            quality_table = load_quality_table(args.quality_table)
            resolve_phred_fn = functools.partial(
                resolve_phred_with_qtable, quality_table=quality_table
            )
        else:
            LOGGER.error(
                f"Requested 'qtable' mode for Phred score resolution, but quality table '{args.quality_table}' does not exist. Exiting."
            )
            sys.exit(1)
    else:
        resolve_phred_fn = resolve_phred_prob

    if args.rule == "5bp":
        rule = ResolutionRule(
            five_bp_rule,
            five_bp_error_codes,
            five_bp_modifications,
            match_award=args.match_award,
            gap_penalty=args.gap_penalty,
            mismatch_penalty=args.mismatch_penalty,
            end_gap_penalty=args.end_gap_penalty,
        )
    elif args.rule == "6bp":
        print("using 6bp rule")
        rule = ResolutionRule(
            six_bp_rule,
            six_bp_error_codes,
            six_bp_modifications,
            match_award=args.match_award,
            gap_penalty=args.gap_penalty,
            mismatch_penalty=args.mismatch_penalty,
            end_gap_penalty=args.end_gap_penalty,
        )
    else:
        print("Not defined !")
        sys.exit()

    (
        out_discard_file1,
        out_discard_file2,
        resolved_path,
        stats_file,
        additional_stats_file,
    ) = get_outfile_names(args)

    LOGGER.info(f"Writing resolved read to: {resolved_path}")
    stats = {}
    with gzip.open(out_discard_file1, mode="wt") as out_discard_file1F, gzip.open(
        out_discard_file2, mode="wt"
    ) as out_discard_file2F, gzip.open(resolved_path, mode="wt") as resolved_handle:
        read1_iter = SeqIO.parse(gzip.open(args.fq1, "rt"), "fastq")
        read2_iter = SeqIO.parse(gzip.open(args.fq2, "rt"), "fastq")
        for r1, r2 in zip(read1_iter, read2_iter):
            result, read_stats = resolve_read_pair(
                r1,
                r2,
                rule,
                resolve_phred_fn=resolve_phred_fn,
                mismatch_threshold=args.mismatch_threshold,
                use_mismatch_aware_trimming=args.mismatch_aware_trimming,
                min_read_length=args.min_read_length,
                episeq_in_qname=args.episeq_in_qname,
                log_additional_stats=args.additional_stats,
                orig_quals_in_sam_tag=args.orig_quals_in_sam_tag,
                xe_tag=args.xe_tag,
            )
            read_stats["num_reads"] = 1
            stats = update_stats(stats, read_stats)
            if result is not None:
                # print(f"writing to {resolved_handle}")
                SeqIO.write(result, resolved_handle, "fastq")
            else:
                # print(f"writing to {out_discard_file1F}")
                SeqIO.write(r1, out_discard_file1F, "fastq")
                SeqIO.write(r2, out_discard_file2F, "fastq")

    # Log and plot statistics
    log_core_stats(stats, stats_file)
    if args.additional_stats:
        if args.generate_plots:
            log_and_plot_additional_stats_single(stats, rule, additional_stats_file)
        else:
            log_additional_stats(stats, rule, additional_stats_file)


if __name__ == "__main__":
    main()
