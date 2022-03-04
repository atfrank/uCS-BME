# uCS-BME: Framework for Using Unassigned NMR Chemical Shifts to Model RNA Secondary Structure

# Requires:

* BME code:
	* slightly modified version of BME is included with this repo
	* original BME code can be found [here](https://github.com/KULL-Centre/BME/)

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
export uCSBME=path/to/this/repo
bash sh/get_chemical_shifts.sh 2LU0 2 data/ ${uCSBME}/SS2CS
```

# Reweight conformational library using unassigned data

```
usage: ucsbme.py [-h] -e EXPERIMENTAL -s SIMULATED -e1 ERROR_ONE -e2 ERROR_TWO -n1 NAME_ONE -n2 NAME_TWO [-o OUTPUT] [-t TMPDIR]

optional arguments:
  -h, --help            show this help message and exit
  -e EXPERIMENTAL, --experimental EXPERIMENTAL
                        Experimental peaks (peaks those peaks to have names
  -s SIMULATED, --simulated SIMULATED
                        Simulated chemical shift table
  -e1 ERROR_ONE, --error_one ERROR_ONE
                        Simulated chemical shift error for dimenison 1
  -e2 ERROR_TWO, --error_two ERROR_TWO
                        Simulated chemical shift error for dimenison 2
  -n1 NAME_ONE, --name_one NAME_ONE
                        Column name for dimenison 1
  -n2 NAME_TWO, --name_two NAME_TWO
                        Column name for dimenison 2
  -o OUTPUT, --output OUTPUT
                        Output prefix for generated files
  -t TMPDIR, --tmpdir TMPDIR
                        Location used to store auxillary files
```
* Example:
```
python ucsbme.py  -e SARS-CoV-2/5_UTR/iminos_experimental.csv -s SARS-CoV-2/5_UTR/iminos_simulated_test.csv -e1 1.89 -e2 0.39 -n1 "(F1) [ppm]" -n2 "(F2) [ppm]" -o data/tmp/test -t data/
```

* Input Format:
	* Unassigned experimental data
		```
		Peak,Region,Type,Index (F2),Index (F1),(F2) [ppm],(F1) [ppm],(F2) [Hz],(F1) [Hz],Intensity [abs],Annotation,
		1,1,AUTOMATIC,195.0,352.0,13.6096,148.1867,10889.4361,12014.4555,53341.56,,
		2,1,AUTOMATIC,204.0,68.0,13.5171,161.8697,10815.4541,13123.8305,55225.48,,
		3,1,AUTOMATIC,145.0,51.0,14.1233,162.6888,11300.4472,13190.2368,55865.04,,
		4,1,AUTOMATIC,214.0,352.0,13.4144,148.1867,10733.2518,12014.4555,57653.94,,
		5,1,AUTOMATIC,213.0,350.0,13.4247,148.2830,10741.4721,12022.2680,58536.78,,
		6,1,AUTOMATIC,201.0,63.0,13.5479,162.1106,10840.1147,13143.3618,58555.99,,
		7,1,AUTOMATIC,132.0,65.0,14.2568,162.0143,11407.3101,13135.5493,59100.34,,
		8,1,AUTOMATIC,227.0,341.0,13.2808,148.7167,10626.3889,12057.4243,62112.55,,
		9,1,AUTOMATIC,308.0,371.0,12.4487,147.2713,9960.5509,11940.2368,64530.53,,
		...
		...
		```
	* Simulated chemical shift data
		```
		model,resid,resname,nucleus,simcs,id
		1,1,ADE,N1,223.7364120000003,.
		1,2,URA,N1,164.06234799999976,.
		1,3,URA,N1,163.26905999999977,.
		1,4,ADE,N1,214.66633199999896,.
		1,5,ADE,N1,220.41685066666628,.
		1,6,ADE,N1,220.40126114285718,.
		1,7,GUA,N1,146.8135439999999,.
		1,8,GUA,N1,147.89474566666672,.
		1,9,URA,N1,166.7628119999997,.
		1,10,URA,N1,167.26346799999976,.
		...
		...
		```
