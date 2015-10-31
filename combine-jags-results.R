# This file contains code for combining output from parallel JAGS scripts
# results in saved file 'concatenated_posterior.rds'

#All values passed to `args` are logical
args<-as.logical(commandArgs(TRUE))
	#The base case is:
	# args = c(TRUE, TRUE, FALSE, FALSE)


#Navigate to results folder
IOP_BX <- args[1]
IOP_SURG <- args[2]
leave_one_out <- args[3]
crop <- args[4]

(save_path <- paste0(
	'batches/',
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG',
	'/leave_one_out_',leave_one_out,
	'/crop_',crop,'/batch-1/'))

setwd(save_path)

#Highlight posterior files of interest
files<-dir()[grep('posterior',dir())]
files<- files[!grepl('concatenated_posterior',files)]
files<- files[!grepl('full_posterior_corner',files)]
files<- files[!grepl('SEED-0',files)]
nruns<-length(files)


####################
####################
# Calculate dimension needed for concatenated candidate set.
# Initialize the concatenated posterior list (ol).

o<-readRDS(files[1])
#str(o$sims.list)
o<-o[names(o)!='sims.matrix']
d1<-dim(o$sims.list$p_eta)[1] #len.sim (size of MCMC posterior) for each JAGS run
ol<-o$sims.list
	#output-long, to store concatenated results
	#initializing (with NAs) in this way preserves naming structure
for(i in 1:length(o$sims.list)){
	#For vector element (eta_hat_means)
	if(class(o$sims.list[[i]])=='numeric'){
		ol[[i]]<-matrix(NA,nruns,length(o$sims.list[[i]]))
		next
	}

	#For other elements
	dOli<-dim(o$sims.list[[i]])
	dOli[1]<-dOli[1]*nruns 
	ol[[i]] <- array(NA,dim=dOli)
}
print(object.size(ol),units='Mb')




####################
####################
# Fill in ol
	
pb <- txtProgressBar(min = 1, max = nruns, char = "=", 
        style = 3)
for(r in 1:nruns){
	ind_r <- ((r-1)*d1+1):(d1*r)
	o_r <-  readRDS(files[r])$sims.list

	for(i in 1:length(ol)){
		if(class(o_r[[i]])=='numeric'){
			ol[[i]][r,] <- o_r[[i]]
		}else{
			text_olr <- paste0(
				'ol[[i]][ind_r', 
				paste(rep(',',length(dim(o_r[[i]]))-1),collapse=''), 
				']',
				' <- o_r[[i]]')
			eval(parse(text=text_olr))
		}
	}

	setTxtProgressBar(pb, r)
}
print(object.size(ol),units='Mb')

saveRDS(ol,file='concatenated_posterior.rds')


