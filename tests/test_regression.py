from __future__ import annotations

import tarfile
import tempfile
import unittest
from pathlib import Path

from non_b_gfa.defaults import RunConfig
from non_b_gfa.runner import run_analysis


class RegressionTests(unittest.TestCase):
    def test_sample_outputs_match_reference_files(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        fixture_tar = repo_root / "test_files.tar"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            reference_dir = tmp_path / "reference"
            output_dir = tmp_path / "output"
            reference_dir.mkdir()
            output_dir.mkdir()

            with tarfile.open(fixture_tar) as archive:
                for member in archive.getmembers():
                    if not member.isfile():
                        continue
                    destination = reference_dir / member.name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    extracted = archive.extractfile(member)
                    self.assertIsNotNone(extracted, member.name)
                    destination.write_bytes(extracted.read())

            config = RunConfig(
                sequence_file=reference_dir / "gfa_test.fasta",
                output_prefix=output_dir / "gfa_test",
                do_wget=False,
            )
            run_analysis(config)

            for reference_file in sorted(reference_dir.glob("gfa_test_*.gff")) + sorted(reference_dir.glob("gfa_test_*.tsv")):
                produced_file = output_dir / reference_file.name
                self.assertTrue(produced_file.exists(), reference_file.name)
                self.assertEqual(reference_file.read_text(encoding="ascii"), produced_file.read_text(encoding="ascii"))


if __name__ == "__main__":
    unittest.main()
