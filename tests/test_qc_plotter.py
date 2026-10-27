import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "07_qc_metrics_and_plots.py"
spec = importlib.util.spec_from_file_location("atac_qc_plotter", SCRIPT)
plotter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plotter)


class InputParsingTests(unittest.TestCase):
    def test_read_numeric_column_with_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sizes.tsv"
            path.write_text("fragment_length\n120\n180\n", encoding="utf-8")
            np.testing.assert_array_equal(plotter.read_numeric_column(path), [120.0, 180.0])

    def test_read_peak_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "peaks.narrowPeak"
            path.write_text("chr1\t10\t40\tpeak1\t100\t.\t2\t3\t4\t15\n", encoding="utf-8")
            np.testing.assert_array_equal(plotter.read_peak_lengths(path), [30])

    def test_invalid_peak_interval_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.bed"
            path.write_text("chr1\t40\t10\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                plotter.read_peak_lengths(path)

    def test_empty_input_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.tsv"
            path.write_text("# no observations\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                plotter.read_numeric_column(path)


if __name__ == "__main__":
    unittest.main()
