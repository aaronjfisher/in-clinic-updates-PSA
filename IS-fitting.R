


######################### 
##     Workflow:

# Define Scenario
# Load data
# Load functions from `IS-functions-setup.R`
# Generate particles or candidate draws for the posterior for each subject. (`gen_particles`)
# These particles are weighted or accepted/rejected to get an estimated posterior for each subject (`posteriors_all_subjects`).  
# All new subjects share the same candidate draws, but have different weights (for importance sampling (IS)) or have different values accepted (for acceptance/rejection sampling (RS)).
# Compare results to assess performance.
#########################



taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))

########################
# Packages and basic internal functions
library(dplyr)
library(MASS)
library(ggplot2)
library(tensor)

invLogit <- function(x)
	return(exp(x)/(1+exp(x)))
########################


#IOP (Informative observation process)
IOP_BX <- TRUE #Informative observation process for biopsy
IOP_SURG <- TRUE #Informative observation process for surgery
leave_one_out <- FALSE
crop <- FALSE

IOPs<-paste0(
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG')
batch_path<-paste0(
	'batches/',
	IOPs,
	'/leave_one_out_',leave_one_out,
	'/crop_',crop,
	'/batch-1/')
posterior_path<-paste0(batch_path,
	'concatenated_posterior.rds')

###############
# Load Data


psa_data_full<-read.csv("simulation-data/psa-data-sim.csv")
pt_data_full<-read.csv("simulation-data/pt-data-sim.csv") #true eta is now contained here, in patient data (pt data)
bx_data_full <- read.csv("simulation-data/bx-data-sim.csv")

#Splines for surg and bx
couldve_had_biopsy_all <- !is.na(bx_data_full$bx_here)
did_have_biopsy_all <- couldve_had_biopsy_all & bx_data_full$bx_here

if(IOP_SURG|IOP_BX){
	sum(couldve_had_biopsy_all)
	sum(did_have_biopsy_all)
}






oo <- readRDS(posterior_path)
	# We currently use the output from the full model as an
	# approximation for the output of a leave-one-out fit
	# (oo = leave-one-out output)
of <- readRDS(posterior_path)
	# output from full model

n_reps <- 10 # How many draws of subj sepecific variables to get for each proposal of population variables.

n_post<-length(oo$p_eta) #posterior size for population variables.
(P<-n_post*n_reps) #Number of particles in our set.



missing_etas <- which(is.na(pt_data_full$obs_eta))
	# obs_eta=1 if observed and aggressive, 0 if observed
	# and not aggressive, or NA if not observed
eta_true_or_jags<-pt_data_full$obs_eta
eta_true_or_jags[missing_etas]<-colMeans(of$eta_hat_means)
	#indexed by subject
###############



##############################
####### Run Functions
##############################

source('IS-functions-setup.R') # contains `posteriors_all_subjects` function

##### Generate particles:

seed<-101
set.seed(seed)


ps1 <- gen_particles(oo,n_reps=n_reps,verbose=TRUE) #particle set 1
runifs <- runif(P) #used for rejection sampling

N <- max(psa_data_full$subj) #number of subjects
e_ss_threshold <- 1000 #minimum expected sample size
n_draws_init <- n_post/10
	#number of proposals to use in first 
	#iteration of "dynamic" approach

if(taskID==1){
	dynamic <- posteriors_all_subjects(
		missing_etas=missing_etas,
		ps=ps1,
		runifs=runifs, 
		e_ss_threshold=e_ss_threshold,
		n_draws_init=n_draws_init,
		get_ZW_approach=FALSE,
		psa_data_full=psa_data_full,
		bx_data_full=bx_data_full,
		verbose=TRUE
		)

	small <- posteriors_all_subjects(
		missing_etas=missing_etas,
		ps=ps1,
		runifs=runifs,
		e_ss_threshold=0,
		n_draws_init=50000,
		get_ZW_approach=FALSE,
		psa_data_full=psa_data_full,
		bx_data_full=bx_data_full,
		verbose=TRUE
		)


	save('seed', 'n_reps','posterior_path',
		'small',
		'dynamic',
		file=paste0(batch_path,Sys.Date(),'_seed_',seed,'_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_P-',P,'_online_fit_results.RData'))
}

if(taskID>1){

	#For "big" approach, run IS in parallel for 20 partitions of the subjects.
	cc <- cut(1:length(missing_etas),breaks=20)
	taskInds <- which(as.numeric(cc)==taskID-1)

	big <- posteriors_all_subjects(
		missing_etas=missing_etas[taskInds],
		ps=ps1,
		runifs=runifs, 
		e_ss_threshold=0,
		n_draws_init=P,
		get_ZW_approach=FALSE,
		psa_data_full=psa_data_full,
		bx_data_full=bx_data_full,
		verbose=TRUE
		)

	saveRDS(big,
		file=paste0(batch_path,Sys.Date(),'_seed_',seed,'_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_P-',P,'_online_fit_big_taskID-',taskID,'.rds'))
}
