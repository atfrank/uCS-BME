# uCS-BME: Framework for Using Unassigned NMR Chemical Shifts to Model RNA Secondary Structure

# Requires:

* BME code:
	** slightly modified version of BME is included with this repo
	** original BME code can be found [here](https://github.com/sbottaro/BME.git)

* Python modules
```
conda create -n ucsbme
conda activate ucsbme
pip install sklearn seaborn matplotlib tqdm
```

# Get chemical shifts (CS)
```
usage: sh/get_chemical_shifts.sh <id>  <#models> <path-to-ct> <ss2cs-path>
```
* Example
```
export CSBME=path/to/this/repo
bash sh/get_chemical_shifts.sh 2LU0 2 data/ ${CSBME}/SS2CS
```

# Reweight conformational library using unassigned data

```
usage: ucsbme.py [-h] -e EXPERIMENTAL -s SIMULATED -e1 ERROR_ONE -e2 ERROR_TWO [-o OUTPUT] [-t TMPDIR] [-d]

optional arguments:
  -h, --help            show this help message and exit
  -e EXPERIMENTAL, --experimental EXPERIMENTAL
                        Experimental peaks (peaks those peaks to have names '(F1) [ppm]' (for C or N nuclei) and '(F2) [ppm]' (for H
                        nuclei))
  -s SIMULATED, --simulated SIMULATED
                        Simulated chemical shift table
  -e1 ERROR_ONE, --error_one ERROR_ONE
                        Simulated chemical shift error for dimenison 1
  -e2 ERROR_TWO, --error_two ERROR_TWO
                        Simulated chemical shift error for dimenison 2
  -o OUTPUT, --output OUTPUT
                        Output prefix for generated files
  -t TMPDIR, --tmpdir TMPDIR
                        Location used to store auxillary files
  -d, --degub           run in debug mode
```
* Example:
```
python ucsbme.py  -e SARS-CoV-2/5_UTR/iminos_experimental.csv -s SARS-CoV-2/5_UTR/iminos_simulated_test.csv -e1 1.89 -e2 0.39 -o data/tmp/test -t data/tmp/
```
