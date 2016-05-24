
#########################################################
# MCMC sampling for DLpar for UN estimates
#########################################################

tfr.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {

	if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
    nr_simu <- mcmc$iter
    nr_countries <- mcmc$meta$nr_countries_estimation
    nr_countries_all <- mcmc$meta$nr_countries
    has_extra_countries <- nr_countries != nr_countries_all
    Triangle_c4.low <- mcmc$meta$Triangle_c4.low
    Triangle_c4.up <- mcmc$meta$Triangle_c4.up
    lower_d <- mcmc$meta$d.low
    upper_d <- mcmc$meta$d.up

    id_DL <- mcmc$meta$id_DL[mcmc$meta$id_DL <= nr_countries]
    id_DL_all <- mcmc$meta$id_DL
    id_early <- mcmc$meta$id_early[mcmc$meta$id_early <= nr_countries]
    id_early_all <- mcmc$meta$id_early
    id_notearly <- mcmc$meta$id_notearly[mcmc$meta$id_notearly <= nr_countries]
    id_notearly_all <- mcmc$meta$id_notearly
    nr_notearly <- length(id_notearly)
    nr_DL <- length(id_DL)
    sqrt.nr.dl.plus1 <- sqrt(nr_DL + 1)
    CoverCplus1 <- nr_DL/(nr_DL+1)

    chi0 <- mcmc$meta$chi0
    psi0 <- mcmc$meta$psi0
    nu_psi0 <- mcmc$meta$nu.psi0
    nu_psiD <- nu_psi0 + nr_DL
    prod_psi0 <- nu_psi0*psi0^2

    mean_eps_tau_0 <- mcmc$meta$mean.eps.tau0
    sd_eps_tau_0 <- mcmc$meta$sd.eps.tau0
    nu_tau0 <- mcmc$meta$nu.tau0
    nu_tauD <- nu_tau0 + nr_notearly
    prod_tau0 <- nu_tau0*sd_eps_tau_0^2

    delta4_0 <- mcmc$meta$delta4.0
    Triangle4_0  <- mcmc$meta$Triangle4.0
    nu4_D <- mcmc$meta$nu4 +  nr_DL
    prod_delta4_0 <- mcmc$meta$nu4*delta4_0^2
    
    nu_deltaD <- mcmc$meta$nu.delta0 +  nr_DL
    a_delta <- nu_deltaD/2
    prod_delta0 <- mcmc$meta$nu.delta0*mcmc$meta$delta0^2 
    
    suppl.T <- if(!is.null(mcmc$meta$suppl.data$regions)) mcmc$meta$suppl.data$T_end else 0
    
    mcenv <- as.environment(mcmc) # Create an environment for the mcmc stuff in order to avoid 
					   			  # copying of the mcmc list 
    
    mcenv$thin <- thin
    ones <- matrix(1, ncol=nr_DL, nrow=3)

	if(is.null(mcenv$eps_Tc)) mcenv$eps_Tc <- get_eps_T_all(mcmc)
	
	if(has_extra_countries) {
		if (verbose)
			cat('\nCountries', mcmc$meta$regions$country_code[(nr_countries+1):nr_countries_all], 
							'not included in the estimation.\n')
	}

  	# sd_Tc with sigma0 and sd_tau_eps
    # the non-constant variance is sum of sigma0 and add_to_sd_Tc
    # matrix with each column one country
    mcenv$add_to_sd_Tc <- matrix(NA, mcenv$meta$T_end-1+suppl.T, nr_countries_all)
    mcenv$data.list <- list()
    for (country in 1:nr_countries_all){
    	# could exclude 1:(tau_c-1) here
    	this.data <- get.observed.tfr(country, mcenv$meta)
    	this.data <- this.data[1:(mcenv$meta$T_end_c[country]-1)]
    	mcenv$data.list[[country]] <- this.data
        mcenv$add_to_sd_Tc[1:(mcenv$meta$T_end_c[country]-1),country] <- (
        				this.data - mcenv$S_sd)*ifelse(this.data > mcenv$S_sd, -mcenv$a_sd, mcenv$b_sd)
	}
	mcenv$const_sd_dummie_Tc <- matrix(0, mcenv$meta$T_end-1+suppl.T, nr_countries_all)
	mid.years <- as.integer(c(if(suppl.T > 0) rownames(mcenv$meta$suppl.data$tfr_matrix) else c(), rownames(mcenv$meta$tfr_matrix)))
	mcenv$const_sd_dummie_Tc[mid.years[-length(mid.years)] < 1975,] <- 1

  	mcenv$mean_eps_Tc <- matrix(0, mcenv$meta$T_end -1+suppl.T, nr_countries_all)
    idx.tau_c.id.notearly.all <- matrix(c(mcenv$meta$tau_c[id_notearly_all], id_notearly_all), ncol=2)
	mcenv$mean_eps_Tc[idx.tau_c.id.notearly.all] <- mcenv$mean_eps_tau
	
	idx.tau_c.id.notearly <- matrix(c(mcenv$meta$tau_c[id_notearly], id_notearly), ncol=2)
    ################################################################### 
    # Start MCMC
	############
    for (simu in start.iter:nr_simu) {
    	if(verbose.iter > 0 && (simu %% verbose.iter == 0))
        	cat('\nIteration:', simu, '--', date())
        unblock.gtk('bDem.TFRmcmc')
        #################################################################
        ## a_sd, b_sd, f_sd and sigma0
        #################################################################
        # updates sd_Tc
        # start with this to get the right sd_Tc in the next steps!!
        mcmc.update.abSsigma0const(mcenv, idx.tau_c.id.notearly)

       #################################################################### 
        # 2. mean_eps_tau sd_eps_tau: gibbs step
        ##################################################################

        eps_taus <- mcenv$eps_Tc[idx.tau_c.id.notearly] # only for estimation countries
        mcenv$mean_eps_tau <- rnorm(1, mean = (sum(eps_taus)/(mcenv$sd_eps_tau^2) +
                        mean_eps_tau_0/(sd_eps_tau_0^2))/
                        (1/sd_eps_tau_0^2 + nr_notearly/(mcenv$sd_eps_tau^2)),
                             sd = 1/sqrt(1/sd_eps_tau_0^2 + nr_notearly/(mcenv$sd_eps_tau^2)) )
        sum_dt2 = sum(  (eps_taus - mcenv$mean_eps_tau)^2)
        mcenv$sd_eps_tau <- sqrt(1/rgamma(1,nu_tauD/2, rate = 1/2*(prod_tau0 + sum_dt2)))

 		#update all not-early countries
		mcenv$sd_Tc[idx.tau_c.id.notearly.all] <- mcenv$sd_eps_tau
		mcenv$mean_eps_Tc[idx.tau_c.id.notearly.all] <- mcenv$mean_eps_tau
        #################################################################### 
        # 2. chi's and psi's: gibbs step
        ##################################################################
        # this is phi_c, the logit transformed d_c's
        d_transformed_declinecountries <- 
                        log((mcenv$d_c[id_DL] - lower_d)/(
                                upper_d-mcenv$d_c[id_DL])) 
        mcenv$chi <- rnorm(1, mean = (sum(d_transformed_declinecountries)/(mcenv$psi^2) +chi0/(psi0^2))/(1/psi0^2 + nr_DL/(mcenv$psi^2)),
                             sd = 1/sqrt(1/psi0^2 + nr_DL/(mcenv$psi^2)))
        sum_dt2 = sum(  (d_transformed_declinecountries - mcenv$chi)^2)
        mcenv$psi <- sqrt(1/rgamma(1,nu_psiD/2, rate = 1/2*(prod_psi0 + sum_dt2)))

        #################################################################### 
        # 2. Triangle4 and delta4_star: gibbs step
        ##################################################################
        # logit transformed Triangle_c4's
        delta4.squared <- mcenv$delta4^2
        delta4.0.squared <- delta4_0^2
        Triangle_c4_transformed <- 
                        log((mcenv$Triangle_c4[id_DL] - Triangle_c4.low)/
                                        (Triangle_c4.up-mcenv$Triangle_c4[id_DL])) 
        mcenv$Triangle4 <- rnorm(1, mean = (sum(Triangle_c4_transformed)/(delta4.squared) + 
        										Triangle4_0/delta4.0.squared)/
                                                  (1/delta4.0.squared + nr_DL/(delta4.squared)),
                             sd = 1/sqrt(1/delta4.0.squared + nr_DL/delta4.squared))

        sum_dt2 = sum(  (Triangle_c4_transformed - mcenv$Triangle4)^2)
        mcenv$delta4 <- sqrt(1/rgamma(1,nu4_D/2, rate = 1/2*(prod_delta4_0 + sum_dt2)))

        #################################################################### 
        # country-specific parameters: d_c, gamma's, U_c and Triangle_c4
        ##################################################################
        
        for (country in id_DL_all){
                mcmc.update.d(country, mcenv)
                mcmc.update.gamma(country, mcenv)
                mcmc.update.Triangle_c4(country, mcenv)
                #print(c(country, mcmc$Triangle_c4[country]))
        } 
 
         # U_c updated only for countries with early decline
         for (country in id_early_all){
                mcmc.update.U(country, mcenv)
         } 

         ##################################################################
         #mcenv# alpha_i's and delta_i's, with gibbs step
         ##################################################################

         mcenv$alpha <- rnorm(3, mean = (apply(mcenv$gamma_ci[id_DL,], 2, sum)/(mcenv$delta^2) + 
         										mcenv$meta$alpha0.p/(mcenv$meta$delta0^2))/
                                                        (1/(mcenv$meta$delta0^2) + nr_DL/(mcenv$delta^2)),
                                sd = 1/sqrt(1/(mcenv$meta$delta0^2) + nr_DL/(mcenv$delta^2)))

         sum_gammas <- apply( (mcenv$gamma_ci[id_DL,] - t(ones*mcenv$alpha) )^2, 2, sum)
         mcenv$delta <- sqrt(1/rgamma(3,a_delta , rate = 1/2*(prod_delta0 + sum_gammas) ) )

         ################################################################### 
         # write samples simu/thin to disk
         ##################################################################
         mcenv$finished.iter <- simu
         mcenv$rng.state <- .Random.seed
         if (simu %% thin == 0){
         	mcenv$length <- mcenv$length + 1
         	flush.buffer <- FALSE
            if (simu + thin > nr_simu) flush.buffer <- TRUE
            store.mcmc(mcenv, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
         }
	}       # end simu loop MCMC
	.cleanup.mcmc(mcenv)
	resmc <- as.list(mcenv)
	class(resmc) <- class(mcmc)
    return(resmc)
}

.cleanup.mcmc <- function(mcmc) {
	if(is.environment(mcmc)) {
		rm(list=mcmc$dontsave[mcmc$dontsave != 'meta'], envir=mcmc)
		return(NULL)
	}
	for(rmitem in mcmc$dontsave[mcmc$dontsave != 'meta']) mcmc[[rmitem]] <- NULL
	return(mcmc)
}

unblock.gtk <- function(option, options.list=NULL) {
	if(!getOption(option, default=FALSE)) return()
	if(!is.null(options.list)) options(options.list)
	# This is to unblock the GUI, if the run is invoked from bayesDem
	# In such a case the gtk libraries are already loaded
	while(do.call('gtkEventsPending', list()))
		do.call('gtkMainIteration', list())

}

tfr.mcmc.sampling.extra <- function(mcmc, mcmc.list, countries, posterior.sample,
											 iter=NULL, burnin=2000, verbose=FALSE, verbose.iter=100) {
	#run mcmc sampling for countries given by the index 'countries'
	nr_simu <- iter
	if (is.null(iter))
    	nr_simu <- mcmc$length
    nr_countries_all <- mcmc$meta$nr_countries
    nr_countries <- length(countries)
    Triangle_c4.low <- mcmc$meta$Triangle_c4.low
    Triangle_c4.up <- mcmc$meta$Triangle_c4.up
    lower_d <- mcmc$meta$d.low
    upper_d <- mcmc$meta$d.up

    id_DL <- mcmc$meta$id_DL[is.element(mcmc$meta$id_DL, countries)]
    id_early <- mcmc$meta$id_early[is.element(mcmc$meta$id_early, countries)]
    id_notearly <- mcmc$meta$id_notearly[is.element(mcmc$meta$id_notearly, countries)]
	    
	suppl.T <- if(!is.null(mcmc$meta$suppl.data$regions)) mcmc$meta$suppl.data$T_end else 0
	const_sd_dummie_Tc_extra <- matrix(0, mcmc$meta$T_end-1+suppl.T, nr_countries)
	mid.years <- as.integer(c(if(suppl.T > 0) rownames(mcmc$meta$suppl.data$tfr_matrix) else c(), rownames(mcmc$meta$tfr_matrix)))
	const_sd_dummie_Tc_extra[mid.years[-length(mid.years)] < 1975,] <- 1
	
	# get values of the hyperparameters (sample from the posterior)
    hyperparameter.names <- tfr.parameter.names(trans=FALSE)
    hyperparameters <- list()
    sampled.index <- sample(posterior.sample, nr_simu, replace=TRUE)
    th.burnin <- get.thinned.burnin(mcmc, burnin)
    for (par in hyperparameter.names) {
    	hyperparameters[[par]] <- c()
    	for(mc in mcmc.list) {
    		if (no.traces.loaded(mc)  || th.burnin < mc$traces.burnin) {
    			traces <- bdem.parameter.traces(mc, par, burnin=th.burnin)
        	} else {
          		traces <- get.burned.tfr.traces(mc, par, th.burnin)
       		}
       		hyperparameters[[par]] <- rbind(hyperparameters[[par]], traces)
       	}
    	hyperparameters[[par]] <- hyperparameters[[par]][sampled.index,]
    }
	mcmc.orig <- mcmc
	mcenv <- as.environment(mcmc) # Create an environment for the mcmc stuff in order to avoid 
					              # copying of the mcmc list 

	updated.var.names <- c('gamma_ci', 'd_c', 'Triangle_c4', 'U_c')
	idx.tau_c.id.notearly <- matrix(c(mcmc$meta$tau_c[id_notearly], id_notearly), ncol=2)
    ################################################################### 
    # Start MCMC
	############
    for (simu in 1:nr_simu) {
        if(verbose.iter > 0 && (simu %% verbose.iter == 0))
			cat('\nIteration:', simu, '--', date())
			unblock.gtk('bDem.TFRmcmcExtra')
        # set hyperparameters for this iteration
        for (par in hyperparameter.names) {
        	if(is.null(dim(hyperparameters[[par]]))) {
        		mcenv[[par]] <- hyperparameters[[par]][simu]
        	} else {
        		mcenv[[par]] <- hyperparameters[[par]][simu,]
        	}
        }
        # compute eps_T, mean_eps_t and sd_Tc
        if(is.null(mcenv$eps_Tc)) mcenv$eps_Tc <- get_eps_T_all(mcenv)
         	
        add_to_sd_Tc_extra <- matrix(NA, mcenv$meta$T_end-1 + suppl.T, nr_countries)
        mcenv$data.list <- list()
    	for (icountry in 1:length(countries)){
    		country <- countries[icountry]
    		this.data <- get.observed.tfr(country, mcenv$meta)
    		this.data <- this.data[1:(mcenv$meta$T_end_c[country]-1)]
    		mcenv$data.list[[country]] <- this.data
			add_to_sd_Tc_extra[1:(mcenv$meta$T_end_c[country]-1),icountry] <- (
						this.data - mcenv$S_sd)*ifelse(this.data > mcenv$S_sd, -mcenv$a_sd, mcenv$b_sd)
		}
		mcenv$mean_eps_Tc <- matrix(0, mcenv$meta$T_end -1 + suppl.T, nr_countries_all)
        mcenv$sd_Tc <- matrix(NA, mcenv$meta$T_end -1 + suppl.T, nr_countries_all)
       	mcenv$sd_Tc[,countries] <- ifelse(const_sd_dummie_Tc_extra==1, 
         						mcenv$const_sd, 1)*
            			   ifelse((mcenv$sigma0 + add_to_sd_Tc_extra)>0, mcenv$sigma0 + add_to_sd_Tc_extra, 
            						mcenv$meta$sigma0.min)
         	
 		mcenv$sd_Tc[idx.tau_c.id.notearly] <- mcenv$sd_eps_tau
 		mcenv$mean_eps_Tc[idx.tau_c.id.notearly] <- mcenv$mean_eps_tau
        #################################################################### 
        # country-specific parameters: d_c, gamma's, U_c and Triangle_c4
        ##################################################################
        for (country in id_DL) {
        	mcmc.update.d(country, mcenv)
            mcmc.update.gamma(country, mcenv)
            mcmc.update.Triangle_c4(country, mcenv)
        } 
 
         # U_c updated only for countries with early decline
         for (country in id_early){
                mcmc.update.U(country, mcenv)
         } 

         ################################################################### 
         # write samples to disk
         ##################################################################
         #update the original mcmc with the new values
         for(var in updated.var.names) {
         	mcmc.orig[[var]] <- mcenv[[var]]
         }
         flush.buffer <- FALSE
         append <- TRUE
		 if (simu == 1) {
			append <- FALSE
			flush.buffer <- TRUE
		 } else {
			if (simu + 1 > nr_simu) flush.buffer <- TRUE
		 }
         store.mcmc(mcmc.orig, append=append, flush.buffer=flush.buffer, countries=countries, 
         				verbose=verbose)
         
	}       # end simu loop MCMC
	mcmc.orig <- .cleanup.mcmc(mcmc.orig)
    return(mcmc.orig)
}
