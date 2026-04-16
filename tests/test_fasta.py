from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from non_b_gfa.fasta import read_fasta_records


class ReadFastaRecordsTests(unittest.TestCase):
    def test_reads_multiple_records_and_normalizes_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "sample.fa"
            fasta_path.write_text(
                ">seq one\n"
                "ACgt-NN12\n"
                ">seq two extra words\n"
                "ttaa\n",
                encoding="ascii",
            )

            records = read_fasta_records(fasta_path)

        self.assertEqual(2, len(records))
        self.assertEqual("seq one", records[0].title)
        self.assertEqual("acgtnn", records[0].sequence)
        self.assertEqual("seq two extra words", records[1].title)
        self.assertEqual("ttaa", records[1].sequence)

    def test_rejects_non_fasta_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "bad.fa"
            fasta_path.write_text("not fasta\nacgt\n", encoding="ascii")

            with self.assertRaisesRegex(ValueError, "FASTA format"):
                read_fasta_records(fasta_path)


if __name__ == "__main__":
    unittest.main()
