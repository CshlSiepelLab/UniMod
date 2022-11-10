---
TITLE: "Unimod"
AUTHORS: "Yixin Zhao; Lingjie Liu; Adam Siepel"
---

# Model-based characterization of the equilibrium dynamics of transcription initiation and promoter-proximal pausing in human cells

## Overview
In our manuscript, we described our initial model and two types of extensions.
One allows pause sites to vary across cells, and the other allows for both varied pause sites and steric hindrance of initiation at steady state. Here, we provide scripts of the implementation, and illustrate how we can use them to estimate initiation rates, pause release rates
and landing pad occupancy for both synthetic and experimental data.

<p align="center">
  <img src="figures/figure1AB.png" alt="unimod" width="500"/>
</p>

<p align = "center">
	Fig. 1 The initial probabilistic model for transcription initiation, promoter-proximal pausing, and elongation
</p>

## Dependencies

The unified model is implemented in the statistical programming language [R](https://www.r-project.org/), and depends on a couple of packages. One of the easiest ways to install them is via [conda](https://docs.conda.io/en/latest/).

```
conda create -n unimod --file environment.yml
```

Once installed, you can activate the environment then run the examples within it,

```
conda activate unimod
```

Test data could be downloaded from [here](http://compgen.cshl.edu/yizhao/unimod/data/), and assumed to be placed within the data directory.

## Examples

### Estimate rates based on simulated data

```
Usage: ./estimate_rates_simulation.R [options]
Estimate transcription rates based on simulated data

Options:
	-h, --help
		Show this help message and exit

	-r CHARACTER, --rds=CHARACTER
		Input file produced by SimPol [default NULL]

	-s LOGICAL, --steric=LOGICAL
		Infer landing-pad occupancy or not [default FALSE]

	-d CHARACTER, --outputDir=CHARACTER
		Directory for saving results [default .]
```

The input data is produced by [SimPol](https://github.com/CshlSiepelLab/SimPol), a simulator we developed for simulating the dynamics of RNA Polymerase (RNAP) on DNA template. One of the outputs from SimPol, pos.RDS, records the last 100 steps of the simulation, containing the information of RNAP positions in every cell. Therefore, we can utilize this information to sample cells, then sample read counts conditional on local RNAP frequency. We can later use this synthetic read counts to infer the transcription rates. The whole process is finished by doing

```
./estimate_rates_simulation.R -r ../data/k50ksd25kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t40n20000s33h17_pos.RDS -d ../outputs/simulation
```

The prefix "k50ksd25kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t40n20000s33h17" of the input file
indicates the parameters we used in this test data set, which is further explained [here](https://github.com/CshlSiepelLab/SimPol#usage). Under the given parameters, we simulated 20,000 cells in total for the equivalent of 40 min. (400,000 time slices). We then randomly sampled 5,000 of the 20,000 cells 50 times for each run. The output csv file contains the following columns:

1. trail, refers to the number of run, from 1 to 50
2. chi, the $\chi$ estimates
3. beta_org, the $\beta$ estimates from the initial model
4. beta_adp, the $\beta$ estimates from the adapted model which allows pause sites to vary across cells

Details of the model and the simulation could be found in the method section [here](https://www.biorxiv.org/content/10.1101/2022.10.19.512929v1.full).

We can also use the same R script to infer landing-pad occupancy,

```
./estimate_rates_simulation.R -r ../data/k50ksd25kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t40n20000s33h17_pos.RDS -s T -d ../outputs/simulation
```

And in this case, a fifth column phi will show up in the result. These are the $\phi$ estimates referring to the occupancy.  

## Citation
Zhao, Y., Liu, L. & Siepel, A. Model-based characterization of the equilibrium dynamics of transcription initiation and promoter-proximal pausing in human cells. 2022.10.19.512929 Preprint at [bioRxiv](https://doi.org/10.1101/2022.10.19.512929) (2022).
