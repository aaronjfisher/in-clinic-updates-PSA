# Jobs for fitting MCMC on training set,
# running IS on specific patients,
# and plotting agreement.

# qsub commands contain parallel jobs
# each qsub call is run sequentially



# There are four logical arguments passed to R scripts from the command line
# these correspond respectively to:
	# * IOP_BX - whether we model an informative observation process for biopsies
	# * IOP_SURG - whether we model an informative observation process for surgeries
	# * leave_one_out - whether we leave out a subject from the training dataset
	# * crop - if leave_one_out=TRUE, crop determines whether we should remove all data from the left out subject (crop=FALSE), or only the most recent data from that subject (crop=TRUE)





###########################
# Run MCMC on entire sample

qsub -N jagsIOPall -t 1-400 -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args TRUE TRUE FALSE FALSE" call-jags-cluster.R jagsIOPall.Rout'
	#Example alternate call with leave one out
		# qsub -N IOPleave1out -t 204-1000 -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args TRUE TRUE TRUE TRUE" call-jags-cluster-master.R IOPleave1out.Rout'


# qstat | wc -l

# rm IOPall.o*
# rm IOPall.e*
# tail IOPall.Rout -n 30

# Combine parallel results
qsub -N combineJAGS -hold_jid jagsIOPall -V -l mf=4G,h_vmem=4G -cwd -b y 'R CMD BATCH --no-save "--args TRUE TRUE FALSE FALSE" combine-jags-results.R combine-jags-results.Rout'



###########################
# Importance sampling


#Fit IS
# task ID ranges from 1-21 (1 for small & dynamic approach, 2-21 for big approach)
qsub -N fitIS -hold_jid combineJAGS -t 1-21 -V -l mf=10G,h_vmem=10G -cwd -b y 'R CMD BATCH --no-save IS-fitting.R fitIS.Rout'

# qstat | wc -l
# rm fitIS.e*
# rm fitIS.o*
# tail fitIS.Rout -n 10


qsub -N combIS -hold_jid fitIS -V -l mf=3G,h_vmem=3G -cwd -b y 'R CMD BATCH --no-save IS-combine-results.R IS-combine-results.Rout'
# This will save plots of the results in the 'plots' folder


