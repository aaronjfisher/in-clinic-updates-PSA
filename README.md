Prostate Surveillance
===========

Prediction modeling for active surveillance of prostate cancer

Analysts can reproduce our results using a sun grid engine computing cluster, by running

`bash run-MCMC-IS.sh`

This will ultimately result in plots being saved to the `plots` folder.

### Contents

#### Folders
* `simulation-data` - contains code for simulation data based on the Johns Hopkins Active Surveillance cohort (JHAS), and results of these simulations. 
* `plots` - folder to save final plots of the analysis

#### Files
* `run-MCMC-IS.sh` - run the pipeline for MCMC and for IS
* `call-jags-cluster.R` - run MCMC using jags. This function calls the script
    - call-jags-functions-setup.R, which in turn calls either
        + `model-for-jags-IOP_BX-IOP_SURG.txt`,
        + `model-for-jags-NIOP_BX-IOP_SURG.txt`,
        + `model-for-jags-IOP_BX-NIOP_SURG.txt`, or
        + `model-for-jags-NIOP_BX-NIOP_SURG.txt`, depending on variables that tell whether informative observation processes for biopsies and surgeries should be included in the model.
* `combine-jags-results.R` - concatenates the results of parallel runs of `call-jags-cluster.R`
* `IS-fitting.R` - run IS on individual patients. This function calls the script
    - `IS-function-setup.R`
* `IS-combine-results.r` - aggregate results from `IS-fitting.R` and generate plots







