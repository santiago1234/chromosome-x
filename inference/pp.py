import sys
from pathlib import Path

sys.path.append('.')

from tracts.driver import run_tracts

script_path = Path(sys.argv[0]).resolve()
script_directory = script_path.parent

script_path = '/Users/santiago/Reasearch/projects/chromosome-x/inference/pp.py'

driver_filename = "one_pulse.yaml"

run_tracts(driver_filename = driver_filename, script_dir = script_path)
