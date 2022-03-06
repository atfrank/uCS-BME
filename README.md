# uCS-BME: Framework for Using Unassigned 2D NMR Chemical Shifts to Model RNA Secondary Structure

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

# Example: 5-UTR of the SAR-CoV-2 RNA

## Get chemical shifts (CS)
```
usage: ss2cs_batch.py [-h] -i INPUT -n NUMBER_OF_STRUCTURES -o OUTPUT -s SS2CS_PATH

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input prefix for secondary structures stored as CT files
  -n NUMBER_OF_STRUCTURES, --number_of_structures NUMBER_OF_STRUCTURES
                        number of secondary structures
  -o OUTPUT, --output OUTPUT
                        path to save output CSV
  -s SS2CS_PATH, --ss2cs_path SS2CS_PATH
                        path to SS2CS
```
* Use SS2CS to predicted chemical shifts from a collection of secondary structures
	* In this example, the script expects to find:
		* ```user_all_1.ct```, ```user_all_2.ct```, ```user_all_3.ct```, ..., ```user_all_12.ct``` in  ```data/SARS-CoV-2/5_UTR/``` 
		* the SS2CS model and required data in ```${uCSBME}/SS2CS/```
	* Predicted chemical shifts will be stored in ```test/simulated_cs.csv```
	* NOTE: To reproduce modeling in paper, change ```-n 12``` to ```-n 500```
	```
	export uCSBME=path/to/this/repo
	rm -rfv test && mkdir test
	python SS2CS/ss2cs_batch.py -i data/SARS-CoV-2/5_UTR/user_all -n 12 -o test/simulated_cs.csv -s ${uCSBME}/SS2CS/ &> /dev/null	
	# uncomment to predict chemical shifts for all 500 5'-UTR structures
	# python SS2CS/ss2cs_batch.py -i data/SARS-CoV-2/5_UTR/user_all -n 500 -o test/simulated_cs.csv -s ${uCSBME}/SS2CS/ &> /dev/null
	```
## Reweighting conformational library using 2D unassigned CS data

```
usage: ucsbme.py [-h] -ex EXPERIMENTAL -si SIMULATED -n1 NAME_ONE -n2 NAME_TWO -e1 ERROR_ONE -e2 ERROR_TWO [-ou OUTPUT] [-im] [-tm TMPDIR] [-se SEPARATION]

optional arguments:
  -h, --help            show this help message and exit
  -ex EXPERIMENTAL, --experimental EXPERIMENTAL
                        CSV file with experimental peaks (peaks those peaks to have names
  -si SIMULATED, --simulated SIMULATED
                        CSV file with simulated chemical shift table
  -n1 NAME_ONE, --name_one NAME_ONE
                        Column name for dimenison 1 (C or N) in the experimental peak file
  -n2 NAME_TWO, --name_two NAME_TWO
                        Column name for dimenison 2 (H) in the experimental peak file
  -e1 ERROR_ONE, --error_one ERROR_ONE
                        Estimated chemical shift prediction error for dimenison 1
  -e2 ERROR_TWO, --error_two ERROR_TWO
                        Estimated chemical shift prediction error for dimenison 2
  -ou OUTPUT, --output OUTPUT
                        Output prefix for generated files
  -im, --imino_only     Use only imino simulated chemical shifts (GUA: N1/H1 and URA: N3/H3)
  -tm TMPDIR, --tmpdir TMPDIR
                        Location used to store auxillary files
  -se SEPARATION, --separation SEPARATION
                        Separation character in CSV
```
* Run uCS-BME on the 5-UTR of the SAR-CoV-2 RNA
	* Use imino chemical shifts (```-im``` flag)
	* Store result in ```test/test``` as specified (```-ou``` flag)
		```		
		python ucsbme.py -ex data/SARS-CoV-2/5_UTR/iminos_experimental.csv -si test/simulated_cs.csv -e1 1.89 -e2 0.39 -n1 "(F1) [ppm]" -n2 "(F2) [ppm]" -ou test/test -tm test/tmp/ -im
		```

* Input Format:
	* Unassigned experimental data:
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
	* Simulated chemical shift data:
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
# Citations:
* If you us SS2CS to compute chemical shifts from secondary structures of RNA we request that you cite:

		@article{zhang2021probabilistic,
		  title={Probabilistic Modeling of RNA Ensembles Using NMR Chemical Shifts},
		  author={Zhang, Kexin and Frank, Aaron T},
		  journal={The Journal of Physical Chemistry B},
		  volume={125},
		  number={35},
		  pages={9970--9978},
		  year={2021},
		  publisher={ACS Publications}
		}
* If you use uCS-BME, we request that you cite:
		@article{moudgal2021probabilistic,
		  title={Using Unassigned NMR Chemical Shifts to Model RNA Secondary Structure},
		  author={Moudgal, Neel and Arhin, Grace and Frank, Aaron T},
		  journal={The Journal of Physical Chemistry B} [In Review],
		  volume={},
		  number={},
		  pages={},
		  year={2021},
		  publisher={ACS Publications}
		}
		as well as:	
		@incollection{bottaro2020integrating,
		title={Integrating molecular simulation and experimental data: a Bayesian/maximum entropy reweighting approach},
		author={Bottaro, Sandro and Bengtsen, Tone and Lindorff-Larsen, Kresten},
		booktitle={Structural Bioinformatics},
		pages={219--240},
		year={2020},
		publisher={Springer}
		}

