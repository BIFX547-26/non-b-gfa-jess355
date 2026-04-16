from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from non_b_gfa import cli


class CliTests(unittest.TestCase):
    def test_main_builds_expected_run_config(self) -> None:
        argv = [
            "gfa",
            "-seq",
            "input.fa",
            "-out",
            "result",
            "-chrom",
            "chr1",
            "-skipWGET",
            "-skipGQ",
            "-maxDRrep",
            "42",
        ]

        with patch.object(sys, "argv", argv), patch("non_b_gfa.cli.run_analysis") as run_analysis:
            cli.main()

        config = run_analysis.call_args.args[0]
        self.assertEqual(Path("input.fa"), config.sequence_file)
        self.assertEqual(Path("result"), config.output_prefix)
        self.assertEqual("chr1", config.chrom)
        self.assertFalse(config.do_wget)
        self.assertFalse(config.findGQ)
        self.assertEqual(42, config.maxDRrep)

    def test_main_enables_wget_only_when_requested(self) -> None:
        argv = [
            "gfa",
            "-seq",
            "input.fa",
            "-out",
            "result",
            "-doWGET",
        ]

        with patch.object(sys, "argv", argv), patch("non_b_gfa.cli.run_analysis") as run_analysis:
            cli.main()

        config = run_analysis.call_args.args[0]
        self.assertTrue(config.do_wget)


if __name__ == "__main__":
    unittest.main()
