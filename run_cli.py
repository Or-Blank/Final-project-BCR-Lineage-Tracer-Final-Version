"""
run_cli.py
==========
Command-line entry point. Run from your project folder:

    python run_cli.py --input data.xlsx --output-dir results/
    python run_cli.py --help

Same as:  python -m bcr_lineage_tracer --input data.xlsx ...
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from BCR_lineage_tracer.main import main

if __name__ == "__main__":
    main()
