



#' Generate candidate particles to weight, or accept/reject. The same candidates can be used for multiple new subjects.
#'
#' @param oo output from leave-one-out JAGS call
#' @param n_reps number of subject-specific parameters to draw for each draw from the population-level parameters.
#' @param verbose whether progress bars are shown.
gen_particles<-function(oo,n_reps,verbose=TRUE){

	K<-2  # two latent classes
	n_post<-length(oo$p_eta) # number of particles from JAGS output
	P<-n_post*n_reps # total number of particles we'll be weighting over

	# Only expand the subj-specific params
	# To avoid repeating redundant information, store subj-specific variables as arrays with:
		#first index denoting 1:n_post
		#second index denoting 1:n_reps
	# Also reshape the `mu` to matrix form

	mu<-array(NA,dim=c(n_post,2,K))
	mu[,1,1]<-oo$mu_int[,1]
	mu[,1,2]<-oo$mu_int[,2]
	mu[,2,1]<-oo$mu_slope[,1]
	mu[,2,2]<-oo$mu_slope[,2]
	
	out<-oo[names(oo)!='eta_hat_means'] #Remove extra summarizing element from posterior.
	out$mu_int <- out$mu_slope <- NULL
	out$mu <- mu

	eta <- array( 
		rbinom(n_post*n_reps,1,prob=rep(oo$p_eta,times=n_reps)),dim=c(n_post,n_reps)
		) #array is filled in first by columns.
	
	
	(expand_time<-system.time({
	#random effects proposals
	b_vec_star <- array(NA,dim=c(n_post,n_reps,2))
	if(verbose) pb_sim <- txtProgressBar(min = 0, max = P, char = "=", style=3)
	for(r in 1:n_reps){ # n_reps = number of times we expand the particles from oo.
	for(oo_ind in 1:n_post){ #index for oo
		#

		#define covariance matrix
		cov_for_this_bvec <- diag(c(oo$sigma_int[oo_ind]^2,oo$sigma_slope[oo_ind]^2))
		cov_for_this_bvec[1,2]<-
		cov_for_this_bvec[2,1]<-oo$cov_int_slope[oo_ind]

		#generate random effects
		b_vec_star[oo_ind,r,]<-mvrnorm(1,mu=mu[oo_ind,,eta[oo_ind,r]+1], Sigma=cov_for_this_bvec)

		if(verbose) setTxtProgressBar(pb_sim,(r-1)*n_post+oo_ind)
	}}}))

	out$b_vec_star <- b_vec_star
	out$eta <- eta
	return(out)
}




#' Get likelihood of each particle,
#' given subject-star's data,
#' and a proposed set of subject and population-level parameters
#'
#' @param ps A particle set outputed from get_particles
#' @param psa_data_star psa data for subject of interest
#' @param bx_data_star biopsy data for subject of interest
get_likelihood<-function(ps, psa_data_star, bx_data_star){

	n_post<-dim(ps$eta)[1]
	n_reps<-dim(ps$eta)[2]

	#Setup subject data
	Y_star <- psa_data_star$log_psa
	X_star <- psa_data_star$vol_std 
	Z_star <- cbind(1, psa_data_star$age_std)

	RC_star <- dplyr::filter(bx_data_star, bx_here==1)$rc
	BX_star <- dplyr::filter(bx_data_star, !is.na(bx_here))$bx_here
	SURG_star <- bx_data_star$surg
	prev_pos_biopsy <- dplyr::filter(bx_data_star, !is.na(bx_here))$prev_G7
		##bx_here=NA means that the patient has reclassified. 
		#prev_G7=1 only in intervals with bx_here=1 & rc=1
		# (that is, the interval where reclassification occurred)
		# or later intervals where bx_here=NA


	#Only take the covariates for which we have actual biopsies,
		#(this can lead to a silent error otherwise)
	couldve_had_biopsy_star <- !(is.na(bx_data_star$bx_here))
	did_have_biopsy_star <- couldve_had_biopsy_star & bx_data_star$bx_here

	V_RC_star <- dplyr::mutate(bx_data_star, intercept=1) %>%
		dplyr::select(intercept, 
			contains("rc_time_ns"), 
			contains("rc_date_ns"), 
			rc_age_std ) %>%
		dplyr::filter(did_have_biopsy_star)

	d_V_RC<-dim(V_RC_star)[2] #Note, these don't include eta yet
	d_Z<-dim(Z_star)[2]

	if(IOP_BX){
		U_BX_star <- bx_data_star %>% 
			dplyr::filter(couldve_had_biopsy_star) %>%#covariates just for BX
			dplyr::mutate(intercept=1) %>%
			dplyr::select(intercept, 
					contains("bx_time_ns"),
					contains("bx_date_ns"),
					contains("bx_age_ns"),
					contains("bx_num_prev_bx_ns") )
	}

	if(IOP_SURG){
		W_SURG_star<- 	bx_data_star %>%
		dplyr::mutate(intercept=1) %>%
		dplyr::select(intercept,
					contains("surg_time_ns"),
					contains("surg_date_ns"),
					contains("surg_age_ns"),
					surg_num_prev_bx_std,
					prev_G7)
	}




	#########
	#likelihood of PSA data
	#To vectorize likelihood fits, we use expanded vectors
		# Expanded vectors are have suffix `_exp`
	# 1) organize parameters in an array (see further details on dimension below)
	# 2) get likelihood of each visit
	# 3) aggregate over visits to get joint log likelihood
			

	if(length(Y_star)==0){
		LL_Y<-0
	}else{
		# arrays in this section generally are of the form:
			# 1:n_visits x 1:n_post x 1:n_reps
			# arrays in `ps` have first two dimensions: 1:n_post x 1:n_reps x [dimension of parameter]

		#effect of eta has already been added in gen_particles function. Eta determines the b_vec that are selected.

		n_visits <- length(Y_star)

		mu_obs_psa <- aperm(
			tensor(#Z_star_X_bvec - result is 1:n_post x 1:n_reps x 1:n_visits
				ps$b_vec_star, Z_star,
				alongA=3, alongB=2
			) + 
			tensor( #beta_X_star
				array(rep(ps$beta, times=n_reps),dim=c(n_post,n_reps,1)),# beta_tensor 
				as.matrix(X_star), 
				alongA=3, alongB=2
			),c(3,1,2)) #need to get n_visits on first dim, so Y_star (length=n_visits) recycles correctly.

		sigma_obs_psa <- aperm( #recycle over n_visits
			array(ps$sigma_res, dim=c(n_post,n_reps,n_visits)),
			c(3,1,2)
			)

		LL_Y_j <- log(dnorm(Y_star,mean=mu_obs_psa, sd=sigma_obs_psa))
		rm(mu_obs_psa, sigma_obs_psa)

		LL_Y<-apply(LL_Y_j,MARGIN=c(2,3),FUN=sum)
		rm(LL_Y_j)
	}

	

	#'  The remaining likelihoods are all binomial, with a
	#'  common form, so we can calculate them more easily with a
	#'  common function
	#'
	#' For either SURG, BX, or RC, get the likelihood for all visits from a subject.
	#' This returns a vector of length P, with the joint log-likelihood of all visits, conditional on each particle
	#' It works under the model that for any particle p,
	#' we have a bernoulli outcome with:
	#'  logit(mean)= W %*% coeffs[p,1:dim(W)] + eta * coeffs[p,dim(W)+1] + interact_SURG * eta * W[dim(W)] * coeffs[p,dim(W)+2]
	#' Where `interact_SURG` is an indicator that adds an interaction term for the SURG regression
	get_joint_LL_measurements<-function(W, outcomes, coeffs, eta, interact_SURG=FALSE){

		if(length(outcomes)==0) return(0)

		n_visits <- length(outcomes)
		if(n_visits != dim(W)[1]) error('Outcome length does not match covariate length')

		d_W <- dim(W)[2]
		n_post <- dim(eta)[1]
		n_reps <- dim(eta)[2]

		#Calculate the sum of 3 terms:
			# coeffs_W; coef_eta; and (optionally) coef_eta_interact_SURG
			# Each will have dimension 1:n_post x 1:n_reps x 1:n_visits

		coeffs_W <-  tensor(
			array(coeffs[,1:d_W],dim=c(n_post,1,d_W)),
			as.matrix(W),
			alongA=3,
			alongB=2
			)
			#results in 1:n_post x 1 x 1:n_visits
				# same for all 1:n_reps because no person-specific effects here.
				# thus, the 2nd dimension is 1.

		#
		eta_coeffs_dWp1 <- coeffs[,d_W+1] * array(eta,dim=c(dim(eta),1))
			#dWp1 indicates d_W+1 
			#results in 1:n_post x 1:n_reps x 1
				#same for all 1:n_visits because eta is constant over time.
				#thus, the third dimesion is 1
			
		eta_coeffs_dWp2_G7 <- 0
		if(interact_SURG){
			eta_coeffs_dWp2 <- coeffs[,d_W+2] * array(eta,dim=c(dim(eta),1))
				#n_post x n_reps x 1 (same for all visits)
			eta_coeffs_dWp2_G7 <- tensor(
				eta_coeffs_dWp2,
				as.matrix(W[,d_W]),
				alongA=3,
				alongB=2
				) #n_post x n_reps x n_visits
		}

		logit_p_exp <- aperm(  
			#Need to put n_visits on the first dimension for use in dbinom.
			aperm(tensor(coeffs_W,matrix(1,nrow=n_reps),2,2), c(1,3,2)) + 
			tensor(eta_coeffs_dWp1,matrix(1,nrow=n_visits),3,2) +
			eta_coeffs_dWp2_G7,
			c(3,1,2))

		p_exp<- invLogit(logit_p_exp)


		LL_j <- log(dbinom(x=outcomes,size=1,prob=invLogit(logit_p_exp)))

		LL <- apply(LL_j,MARGIN=c(2,3),sum)

		return(LL)
	}


	LL_RC <- get_joint_LL_measurements(
		W=V_RC_star,
		outcomes=RC_star,
		coeffs=ps$gamma_RC,
		eta=ps$eta,
		interact_SURG=FALSE)


	LL_BX <-LL_SURG <- 0 #zero if not included in model
	if(IOP_BX){
		LL_BX <- get_joint_LL_measurements(
		W=U_BX_star,
		outcomes=BX_star,
		coeffs=ps$nu_BX,
		eta=ps$eta,
		interact_SURG=FALSE)
		
	}
	if(IOP_SURG){
		LL_SURG <- get_joint_LL_measurements(
			W=W_SURG_star,
			outcomes=SURG_star,
			coeffs=ps$omega_SURG,
			eta=ps$eta,
			interact_SURG=TRUE)
	}

	W <- exp(LL_Y + LL_BX + LL_RC + LL_SURG)

	return(W)
}





##############################
# Estimate posterior means for each subject.

# For each subject, calculate and store:
	# for Importance sampling
		# * the weighted eta posterior mean
		# * the effective sample size for the number of posterior draws
	# for Rejection sampling
		# * average eta over accepted draws
		# * number of accepted draws from the posterior



#' @param data_star a list of dataframes for subject star
#' @param ps particle set from get_particles
#' @param runifs for use in rejection sampling
#' @param rej_const for use in rejection sampling
#' @param e_ss_threshold for dynamic method, when should we increase the number of particles used.
#' @param n_draws_init how many particle draws to use on first attempt at estimating risk
#' @param get_ZW_approach currently in development
posterior_star<-function( data_star, ps, runifs, rej_const=NULL,	e_ss_threshold=800,	n_draws_init = min(length(ps$eta), e_ss_threshold*2),get_ZW_approach=TRUE){

	P <- length(ps$eta)
	n_post <- dim(ps$eta)[1]
	n_reps <- dim(ps$eta)[2]
	
	if(e_ss_threshold > P) stop("Effective sample size exceeds the number of particles")
	if(n_draws_init > P) stop('n_draws_init must be less than or equal to the number of particles')

	#Want to approximately double the size of the particle set for each time we do not
	# exceed the minimum effective sample size.
	#Get seq starting at 0, ending at P, and roughly increasing by a factors of 2 at each step, from n_draws_init
	

	l2P <- log(P,base=2)
	l2n1 <- log(n_draws_init,base=2)
	l2n1_big <- log(ceiling(n_draws_init/n_post)*n_post,base=2)
	l2n_post <- log(n_post,base=2)

	seq1<-c()
	seq2<-seq(from=l2n1_big,to=l2P,by=1) #if(n_draws_init > n_post)
	if(l2n1 <= l2n_post){
		seq1 <- c(seq(from=l2n1,to=l2n_post,by=1),l2n_post)
		seq2 <- seq(from=l2n_post,to=l2P,by=1)
	}
	#seq1 goes up to the first n_post entries.
	#after that, seq2 goes up by multiples of n_post (by doubling the size),
	#so that we always get clean slices of the array. 

	log_breaks <- unique( c(-Inf, seq1, seq2, l2P) )
	breaks <- round(2^log_breaks)

	# Double checks -- these errors shouldn't ever be triggered, but are here as a redundancy.
	if(any(breaks > n_post & breaks %% n_post!=0)) stop('error occured')

	#Function to select a subset (inds) of the particle set
	#1:n_post index is stored in the first dimension of ps
	select_nth_dim <- function(x, n=1, inds, verbose=FALSE, dim_intact=TRUE){

	    d <- dim(x)
	    ld <- length(d)
	    if(is.null(d)) ld <- 1 #if we have a vector

	    if(n>ld) stop('n must be less than dimension length of x')

	    if (ld == 1) return(x[inds])

	    prefix <-
	    suffix <- ''

	    if(dim_intact){
	    	dOut <- d
	    	dOut[n] <- length(inds)

	    	prefix <- 'array('
	    	suffix <- ',dim=dOut)'
	    }

	    text_output <- 
	    	paste0(prefix,"x[",
	    		paste0(rep(",",n-1),collapse=''),
	    		'inds',
	    		paste0(rep(",",ld-n),collapse=''), 
	    		"]", suffix
	    	)
	    if(verbose) message(paste('returning',text_output))
	    eval(parse(text = text_output))
	}




	subject_specific <- c('eta','b_vec_star')

	#Calculate likelihood of particles.
	#If effective sample size not met, double size of particle set and recalculate.
	#Initialize likelihood as full array based on using all 
	#elements in ps, and fill this array as needed based
	#on effective sample size
	time_IS<-system.time({
	likelihood <- matrix(NA,n_post,n_reps)
	
	for(i in 2:(length(breaks))){

		#inds_1 and inds_2 select the first and 2nd dimensions of the arrays in ps

		these_inds_vector<-(breaks[i-1]+1):breaks[i]
		inds_1 <- unique(these_inds_vector %% n_post)
		inds_1[which(inds_1==0)] <- n_post
			# replace "zero" entries with n_post
			# to adjust for how %%  sets these to zero.
		inds_2 <- unique(ceiling(these_inds_vector/n_post))


		if(max(inds_2)>1 & length(inds_1)!=n_post) stop('dimension mismatch has occured') #this should never happen, but putting in a redundancy just in case.


		psi <- list()
		for(j in 1:length(ps)){
			#subset by n_post
			part_j <- select_nth_dim(ps[[j]],n=1,inds=inds_1,dim_intact=TRUE)
			psi[[j]]<-part_j
			if(names(ps)[j] %in% subject_specific){
				#Further subset by n_reps
				psi[[j]] <- select_nth_dim(part_j,n=2,inds=inds_2,dim_intact=TRUE)
			}
			names(psi)[j]<-names(ps)[j]
		}

		likelihood_i <- get_likelihood(
			ps=psi, 
			psa_data_star=data_star$PSA,
			bx_data_star=data_star$BX
			)
		
		likelihood[inds_1,inds_2] <- likelihood_i

		# Below, save copies in W without altering `likelihood`, so 
		# likelihood can be appended if needed in next stage.
		W <- c(likelihood) / sum(likelihood,na.rm=TRUE)
		W <- W[!is.na(W)]
		effective_ss <- 1/crossprod(W)
		
		last_ind_used <- breaks[i]
		if(effective_ss >= e_ss_threshold) break

	}

	if(effective_ss < e_ss_threshold) warning('Even with full particle set, the maximum effective sample size required was not met.')

	inds_cumulative <- which(!is.na(likelihood),arr.ind=TRUE)
	ps_cumulative_eta <- select_nth_dim(
				ps[[j]],
				n=1,
				inds=unique(inds_cumulative[,1]),
				dim_intact=TRUE) %>%
			select_nth_dim(.,n=2,
				inds=unique(inds_cumulative[,2]),
				dim_intact=TRUE)
	
	## Importance Weighting ##
	etas_IS_star <- crossprod(W,c(ps_cumulative_eta))

	})['elapsed'] #end of system.time()
	names(time_IS) <- c()

	## Rejection Sampling ##
	if(is.null(rej_const)) rej_const <- max(likelihood,na.rm=TRUE)
	accept_ind <- which(
		(c(likelihood[!is.na(likelihood)])/rej_const) >= runifs[1:last_ind_used]
		)
	num_accepted_star <- length(accept_ind)
	etas_RS_star <- mean(ps_cumulative_eta[accept_ind])

	## Zhenke Wu's 2015 Approach ##
	# (still in development)
	n_missing_ZW <- 
	etas_ZW_star <- NA
	if(get_ZW_approach){
		psi_eta_1 <-
		psi_eta_0 <- psi
		psi_eta_1$eta[] <- 1
		psi_eta_0$eta[] <- 0

		likelihood_eta_1 <- get_likelihood( #denote as l1
				ps=psi_eta_1,
				psa_data_star=data_star$PSA,
				bx_data_star=data_star$BX
				) 
		likelihood_eta_0 <- get_likelihood( #denote as l0
				ps=psi_eta_0,
				psa_data_star=data_star$PSA,
				bx_data_star=data_star$BX
				)

		#Want to calculate 
			# l1 * p_eta / (l0 *(1-p_eta) + l1 *p_eta)
			# 1 / [ l0 *(1-p_eta)/(l1 * p_eta) + 1 ]
		# and then take the mean
		# Can also implement on log scale
		W_ZW <- likelihood_eta_1 * psi$p_eta / 
			(likelihood_eta_1 * psi$p_eta + likelihood_eta_0 * (1-psi$p_eta))

		n_missing_ZW <- sum(is.na(W_ZW)) #for some particles, we end up with 0 / 0, because likelihood_eta_1 & likelihood_eta_2 both end up rounding down to zero.
		etas_ZW_star <- mean(W_ZW,na.rm=TRUE)
	}

	return(list(
		#Importance Sampling:
		W=W,
		etas_IS_star = etas_IS_star,
		time_IS= time_IS,
		effective_ss_star = effective_ss,
		particle_draws = last_ind_used,
		#Rejection Sampling:
		etas_RS_star = etas_RS_star,
		num_accepted_star = num_accepted_star,
		#ZW Sampling
		etas_ZW_star= etas_ZW_star,
		ns_missing_ZW = n_missing_ZW
	))
}








###### Run loop over subjects
posteriors_all_subjects <- function(missing_etas,
	ps, runifs, e_ss_threshold, n_draws_init, get_ZW_approach=FALSE,
	psa_data_full=psa_data_full,
	bx_data_full=bx_data_full,
	verbose=TRUE
	){


	N <- length(missing_etas)

	#names for output of this function
	namesOut <-c(
		'subj',
		'time_fit_total', #time spent on subject star
		'time_IS', #time spent total doing IS (possibly dynamically)
		'effective_ss', #effective sample for importance weighting
		'particle_draws', #number of proposals used
		'num_accepted', #number of proposals accepted by rejection sampling
		'etas_RS', #mean among proposals accepted by RS
		'etas_IS', #posterior mean estimates from IS
		'etas_ZW', #currently still in development
		'ns_missing_ZW') #diagnostic/debugging for ZW method (still in development)
	out<-data.frame(matrix(NA,N,length(namesOut)))
	names(out)<-namesOut

	time_fit_total<-rep(NA,N)


	if(verbose) pb <- txtProgressBar(min = 1, max = length(missing_etas), char = "=", style=3)
	for(i in 1:length(missing_etas)){

		#To test for star = a random number, use
		#star <- sample(missing_etas,1)

		data_star <- list(
			PSA=filter(psa_data_full, subj == missing_etas[i]),
			BX=filter(bx_data_full, subj == missing_etas[i])
		)

		rej_const=NULL
		time_fit_total<-system.time({
			post_star <- posterior_star(
				data_star,ps,
				runifs=runifs,
				rej_const=rej_const,
				e_ss_threshold = e_ss_threshold,
				n_draws_init = n_draws_init,
				get_ZW_approach=get_ZW_approach)
		})['elapsed']

		#remove '_star' from names
		names(post_star)<-sapply(names(post_star),function(x){
			strsplit(x,'_star')[[1]][1]
		})
		
		vec_star<-c('subj'=missing_etas[i],
			unlist(post_star[names(post_star) %in% namesOut])
			)
		out[i,]<- vec_star[namesOut]
		out[i,'time_fit_total']<- time_fit_total
		

		if(verbose) setTxtProgressBar(pb, i)
	}


	return(out)

}



