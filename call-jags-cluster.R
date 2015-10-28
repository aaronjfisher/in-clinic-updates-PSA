



#All values passed to `args` are logical
args<-as.logical(commandArgs(TRUE))
	#The base case is:
	# args = c(TRUE, TRUE, FALSE, FALSE)

IOP_BX <- args[1]
IOP_SURG <- args[2]
leave_one_out <- args[3]
crop <- args[4]

#import environment variable, used for running multiple chains in parallel
taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if(is.na(taskID)) taskID <- 0 #for testing with qrsh & interactive jobs

#star is the index of a subject to not include in the training data
#this subject can then be estimated with IS.
#if star=0, no subjects are excluded from the training set.
#Setting star=0 can serve as an approximation to a leave-one-out simulation.
if( leave_one_out) {SEED <- 0		; star <- taskID}
if(!leave_one_out) {SEED <- taskID 	; star <- 0; 		crop <- FALSE} #here, parallel is just used to extend the number of MCMC draws.


(save_path <- paste0(
	'batches/',
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG',
	'/leave_one_out_',leave_one_out,
	'/crop_',crop,'/batch-1/'))

if(!dir.exists(save_path)) dir.create(save_path,recursive=TRUE)

# MCMC settings
# ni <- 100; nb <- 20; nt <- 5; nc <- 1
# ni <- 1000; nb <- 20; nt <- 5; nc <- 1
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1
# ni <-  250000; nb <- 25000; nt <- 10; nc <- 1 ##Do not use this last setting. It results in vectors too large for memory and crashe.
	

(len.sim <- round((ni - nb)/nt))

#save file to show that job has started
cat('',file=paste0(save_path,Sys.Date(),'_started_SEED-',SEED,'_star-',star,'_crop-',crop,'_len.sim-',len.sim,'.txt'))

source('call-jags-functions-setup.R') #includes function `call_jags_psa`
outJAGS <- call_jags_psa(seed=SEED) #final output of this script


#To reduce file size of results, make new object outJAGS2

#Ignore subject specific parameters for most or all subjects
params2ignore<-c('b_vec','rho_int_slope','eta_hat','p_rc')
if(IOP_BX) params2ignore <- unique(c(params2ignore, 'p_bx'))
if(IOP_SURG) params2ignore <- unique(c(params2ignore, 'p_surg'))


#first store diagnostic and general information
outJAGS2 <- outJAGS[!names(outJAGS) %in% c('sims.matrix','sims.list','sims.array')] 

#then define the sims.list (_sl) for outJAGS2
outJAGS2_sl <- outJAGS$sims.list[!names(outJAGS$sims.list) %in% params2ignore]
outJAGS2_sl$eta_hat_means <- apply(outJAGS$sims.list$eta_hat,2,mean)
missing_etas <- which(is.na(pt_data$obs_eta))
n_observed_eta <- sum(!is.na(pt_data$obs_eta)) 
	#these by design are front-loaded into obs_eta
if(crop){ #if `crop`, save subj-specific posterior for subj star. otherwise do not.
	outJAGS2_sl$b_vec_star <- outJAGS$sims.list$b_vec[,star,]
	outJAGS2_sl$eta_hat_star <- NULL
	if(star %in% missing_etas)
		outJAGS2_sl$eta_hat_star <- outJAGS$sims.list$eta_hat[,star-n_observed_eta]
}

outJAGS2$sims.list <- outJAGS2_sl

#Function analogous to `head()`, but for arrays
corner <- function(x, n=6 ){
	d <- dim(x)
	ld<-length(d)

	if(ld==0) return(x[1:min(n,length(x))])

	#If x is an array:
	text_n <- paste0('1:',pmin(n,dim(x)))
	text_ind <- paste0('x[', paste(text_n,collapse=','), ']')
	eval(parse(text=text_ind))
}
sizeMb<-function(x) {object.size(x)/1000000}

lapply(outJAGS2$sims.list,corner)
lapply(outJAGS2$sims.list,sizeMb)
sizeMb(outJAGS2)

# str(outJAGS2$sims.list)

#########
#Save results
saveRDS(outJAGS2,file=paste0(save_path,Sys.Date(),'_posterior_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_SEED-',SEED,'_star-',star,'_crop-',crop,'_nsim-',len.sim,'.rds'))
saveRDS(lapply(outJAGS$sims.list,corner),file=paste0(save_path,Sys.Date(),'_full_posterior_corner_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_SEED-',SEED,'_star-',star,'_crop-',crop,'_nsim-',len.sim,'.rds'))
	#note, this is not outJAGS2, it's a preview of outJAGS

