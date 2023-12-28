#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#

from typing import List, Optional

import argparse
import yaml
import logging
from couplet.export import (
    log_core_stats_merged,
    log_and_plot_additional_stats_merged,
)


def postprocess_stats(input_stats_files: List[str], 
                      input_additional_stats_files: Optional[List[str]], 
                      output_folder: str, 
                      output_prefix: str ) -> None:

    # Setup LOGGER
    log_output_file = (
        output_folder + "/" + output_prefix + "_post_process.log"
    )
    logging.basicConfig(
        filename=log_output_file,
        format="%(asctime)-15s %(filename)s:%(lineno)d %(message)s",
    )
    LOGGER = logging.getLogger("root")
    LOGGER.setLevel(logging.INFO)

    # Standard stats
    log_core_stats_merged(
        input_stats_files, output_folder, output_prefix
    )

    # Mismatch stats
    if input_additional_stats_files is not None:
        log_and_plot_additional_stats_merged(
            input_additional_stats_files, output_folder, output_prefix
        )


