
DLcurve <- function(DLpar, tfr, p1, p2){ 
 # gives the DL-decrement
 # DLpar is c(Delta1, Delta2, Delta3, Delta4, d_c)
 # tfr is a vector for which the decrements for this curve need to  	be calculated
 	dlvalue <- rep(0.0, length(tfr))
	res <- .C("doDLcurve", as.numeric(DLpar), as.numeric(tfr), p1, p2, length(tfr), dl_values=dlvalue)
	return(res$dl_values)
#    t_mid1 <- DLpar[4] + DLpar[3] + DLpar[2] + 0.5 * DLpar[1]
#    t_mid3 <- DLpar[4] + 0.5 * DLpar[3]
#    DLcurve <- 5 * DLpar[5] * (-1/(1 + exp(-log(p1^2)/DLpar[1] * 
#        (tfr - t_mid1))) + 1/(1 + exp(-log(p2^2)/DLpar[3] * (tfr - 
#        t_mid3))))
#    return(ifelse((DLcurve < 0)|(tfr <= 1), 0, DLcurve))

}

##################################################################
# function to get the distortion for country for a given set of DLparameters
# note: this function gives only the eps's in (tau, lambda-1)
# (because rest is NA!!)
get_eps_T = function (DLpar, country, tfr_matrix, start, lambda, p1, p2) 
{
#    eps_T <- NULL
    tfr <- tfr_matrix[start:lambda, country]
    ldl <- length(tfr)-1
    dl <- DLcurve(DLpar, tfr[1:ldl], p1, p2)
#    for (t in start:(lambda - 1)) {
#        eps_T <- c(eps_T, tfr_matrix[t + 1, country] - tfr_matrix[t, 
#            country] + DLcurve(DLpar, tfr_matrix[t, country], 
#            p1, p2))
#    }
    eps_T <- tfr[2:(ldl+1)] - tfr[1:ldl] + dl
    return(eps_T)
}

#get.dl.index <- function(country, meta) {
#	index <- meta$start_c[country] : meta$lambda_c[country]
#	if(is.null(meta$suppl.data) || meta$has.suppl.data[country]) return(index)
#	return(index + mcmc$meta$suppl.data$T_end)
#}

get.eps.T <- function (DLpar, country, meta) 
{
    tfr <- get.observed.tfr(country, meta)[meta$start_c[country]:meta$lambda_c[country]]
    ldl <- length(tfr)-1
    dl <- DLcurve(DLpar, tfr[1:ldl], meta$dl.p1, meta$dl.p2)
    return(tfr[2:(ldl+1)] - tfr[1:ldl] + dl)
}

get_eps_T_all <- function (mcmc) {
	suppl.T <- if(!is.null(mcmc$meta$suppl.data$regions)) mcmc$meta$suppl.data$T_end else 0
	eps_Tc <- matrix(NA, mcmc$meta$T_end-1 + suppl.T, mcmc$meta$nr_countries)
    for (country in mcmc$meta$id_DL){
    	theta <- c((mcmc$U_c[country]-mcmc$Triangle_c4[country])*exp(mcmc$gamma_ci[country,])/                                     
                                sum(exp(mcmc$gamma_ci[country,])), mcmc$Triangle_c4[country], mcmc$d_c[country])
        eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country] <- get.eps.T(theta, country, mcmc$meta)
        }
    return(eps_Tc)
}

find.lambda.for.one.country <- function(tfr, T_end) {
	lambda <- T_end
	if ( sum(tfr<2, na.rm=TRUE)>2 ){
		period <- .get.T.start.end(tfr)[1]+2
		while (period<=T_end){
			if (( (tfr[period] - tfr[period-1]) >0 )& 
                    ( (tfr[period-1] - tfr[period-2]) >0) & 
                    (prod(tfr[(period-2):period]<2)==1) 
             	) {
				lambda = period-1
				period = T_end+1
			} else { 
				period = period +1
            }
         }
	}
	return(lambda)
}



.get.T.start.end <- function(tfr) {
	# Return first index after NAs at the beginning of the time series
	isna.tfr <- is.na(tfr)
	start <- if(sum(isna.tfr) > 0 && isna.tfr[1]) which(diff(cumsum(isna.tfr))==0)[1]+1 else 1
	return(c(start, sum(!isna.tfr) + start - 1))
} 

get.observed.tfr <- function(country.index, meta, matrix.name='tfr_matrix', matrix.name.suppl=matrix.name)
	return(get.observed.with.supplemental(country.index, meta[[matrix.name]], meta$suppl.data, matrix.name.suppl))

get.observed.with.supplemental <- function(country.index, matrix, suppl.data, matrix.name='tfr_matrix') {
	data <- matrix[,country.index]
	if(!is.null(suppl.data[[matrix.name]])) {
    	supp.c.idx <- suppl.data$index.from.all.countries[country.index]
    	if(is.na(supp.c.idx)) {sdata <- rep(NA, nrow(suppl.data[[matrix.name]])); names(sdata) <- rownames(suppl.data[[matrix.name]])}
    	else sdata <- suppl.data[[matrix.name]][,supp.c.idx]
    	data <- c(sdata, data)
    }
	return(data)
}

find.tau.lambda.and.DLcountries <- function(tfr_matrix, min.TFRlevel.for.start.after.1950 = 5, #5.5, 
												max.diff.local.and.global.max.for.start.at.loc = 0.53, #0.5,
												delta.for.local.max = 0.001,
												suppl.data=NULL) 
# gets tau_c and puts NAs before tau_c
# gets ids of DL (where decline has been observed)
# and divides those into early and not early
# find lambda_c's based on definition tfr increased twice, below 2
{
	post.v2 <- TRUE
	if(!is.null(getOption("TFRphase2.pre.v2", NULL)) && getOption("TFRphase2.pre.v2")==TRUE) {
		# This is for backward-compatibility to make some publications reproducible.
		min.TFRlevel.for.start.after.1950 = 5.5
		max.diff.local.and.global.max.for.start.at.loc = 0.5
		delta.for.local.max = 0
		post.v2 <- FALSE
		warning("Phase II is searched for using pre-v2.0 methodology. To switch to the current method set 'options(TFRphase2.pre.v2=FALSE)'.")
	}
    T_end <- dim(tfr_matrix)[1]
    nr_countries <- dim(tfr_matrix)[2]
    T_end_c <- lambda_c <-rep(T_end, nr_countries)
    #has.suppl <- rep(FALSE, nr_countries)
    T.suppl <- if(is.null(suppl.data$regions)) 0 else dim(suppl.data$tfr_matrix)[1]
    tau_c <- start_c <- rep(NA, nr_countries)
    for (country in 1:nr_countries) {
    	data <- get.observed.with.supplemental(country, tfr_matrix, suppl.data)
    	has.suppl <- length(data) > T_end && !is.na(suppl.data$index.from.all.countries[country])
    	# ignoring NAs at the beginning
    	T.start.end <- .get.T.start.end(data)
    	T.start <- T.start.end[1]
    	T_end_c[country] = T.start.end[2]
    	lT <- T_end_c[country] - T.start + 1
    	local_max_indices <- rep(NA, lT)
    	d <- diff(data[T.start:T_end_c[country]])
        does_tfr_decrease <- ifelse(d < delta.for.local.max, 1, 0)
        local_max_indices[1] <- does_tfr_decrease[1]
   		# in middle only a local max if increase is followed by decrease
        local_max_indices[-c(1, lT)] = diff(does_tfr_decrease)
   		# at end local max if preceded by an increase 
        local_max_indices[lT] = 1 - does_tfr_decrease[lT - 1]
 		value_global_max = max(data, na.rm = TRUE)
 		max_index <- max(seq(T.start, T_end_c[country]) * (local_max_indices > 0) * 
 						ifelse(data[T.start:T_end_c[country]] >
            				value_global_max - max.diff.local.and.global.max.for.start.at.loc, 1, 0))
        # move the point to the right if there are more recent points with the same values
        if(post.v2) {
        	is.same <- c(0, ifelse(abs(d) < delta.for.local.max, 1, 0), 0)
        	cs.same <- cumsum(is.same[(max_index-T.start+1):(lT+1)])
        	max_index <- max_index + cs.same[min(which(diff(cs.same)==0))]
        }
        tau_c[country] <- max_index
        start_c[country] <- tau_c[country]

		if(data[tau_c[country]] < min.TFRlevel.for.start.after.1950) {
        	if ((post.v2 && as.integer(names(data)[T.start]) > 1855) || !post.v2) {
        		tau_c[country] <- -1
        		start_c[country] <- which(!is.na(data))[1] # first data point that is not NA
        	}
        }
        if (tau_c[country] > 1 && T.suppl > 0 && has.suppl) 
        	suppl.data$tfr_matrix[1:min(tau_c[country] - 1, T.suppl),suppl.data$index.from.all.countries[country]] <- NA
        if(tau_c[country] > T.suppl + 1)
            tfr_matrix[1:(tau_c[country] - 1 - T.suppl), country] <- NA

        lambda_c[country] <- find.lambda.for.one.country(data, T_end_c[country])
        if (lambda_c[country] < T_end_c[country]) { # set NA all values between lambda_c and T_c_end
         	if(lambda_c[country] < T.suppl) {
         		suppl.data$tfr_matrix[(lambda_c[country] + 1):min(T.suppl, T_end_c[country]),
         										suppl.data$index.from.all.countries[country]] <- NA
         		if(T_end_c[country] > T.suppl) tfr_matrix[1:(T_end_c[country]-T.suppl),country] <- NA
         	} else tfr_matrix[(lambda_c[country] - T.suppl + 1):(T_end_c[country]-T.suppl),country] <- NA
        }
    }

    id_Tistau <- seq(1, nr_countries)[tau_c == T_end_c]
        # not needed in fit
    id_DL <- seq(1, nr_countries)[(tau_c != T_end_c)
    #& (apply(tfr_matrix, 2, min, na.rm = TRUE) >1) # excludes Macao and Hong Kong
                ]
        # for par in BHM
   id_early <- seq(1, nr_countries)[(tau_c == -1)
    #& (apply(tfr_matrix, 2, min, na.rm = TRUE) >1)
    			]
    # needed for which U_c's to update, this is updated in mcmc
   id_notearly = seq(1, nr_countries)[(tau_c != -1)  & (tau_c != T_end_c)
    #& (apply(tfr_matrix, 2, min, na.rm = TRUE) >1)
    			]
    return(list(tau_c = tau_c,  id_Tistau = id_Tistau, id_DL = id_DL, id_early = id_early,
        		id_notearly = id_notearly, tfr_matrix = tfr_matrix, T_end_c=T_end_c, 
        		lambda_c=lambda_c, start_c=start_c, suppl.matrix=suppl.data$tfr_matrix
         ))
}

mcmc.meta.ini <- function(...,
						U.c.low,
						start.year=1950, present.year=2015, 
						wpp.year=2015, my.tfr.file = NULL, my.locations.file = NULL,
						proposal_cov_gammas = NULL, # should be a list with elements 'values' and 'country_codes'
						verbose=FALSE
					 ) {
	# Initialize meta parameters - those that are common to all chains.
	args <- list(...)
	mcmc.input <- list()
	for (arg in names(args)) mcmc.input[[arg]] <- args[[arg]]
	mcmc.input$U.c.low.base <- U.c.low
	mcmc.input$start.year <- start.year
	mcmc.input$present.year <- present.year
	mcmc.input$wpp.year <- wpp.year
	if(present.year-3 > wpp.year) warning("present.year is much larger then wpp.year. Make sure WPP data for present.year are available.")
	tfr.with.regions <- set_wpp_regions(start.year=start.year, present.year=present.year, wpp.year=wpp.year, 
										my.tfr.file = my.tfr.file, my.locations.file=my.locations.file, verbose=verbose)

	meta <- do.meta.ini(mcmc.input, tfr.with.regions,  
						proposal_cov_gammas=proposal_cov_gammas, verbose=verbose)
	return(structure(c(mcmc.input, meta), class='bayesTFR.mcmc.meta'))
}
	
	
do.meta.ini <- function(meta, tfr.with.regions, proposal_cov_gammas = NULL, 
						use.average.gammas.cov=FALSE, burnin=200, verbose=FALSE) {
	results_tau <- find.tau.lambda.and.DLcountries(tfr.with.regions$tfr_matrix, suppl.data=tfr.with.regions$suppl.data)
	tfr_matrix_all <- tfr.with.regions$tfr_matrix_all
	tfr_matrix_observed <- tfr.with.regions$tfr_matrix
	updated.tfr.matrix <- results_tau$tfr_matrix
	suppl.data <- tfr.with.regions$suppl.data
	if(!is.null(suppl.data$regions)) {
		suppl.data$tfr_matrix_all <- tfr.with.regions$suppl.data$tfr_matrix
		suppl.data$tfr_matrix <- results_tau$suppl.matrix
		suppl.data$T_end <- dim(suppl.data$tfr_matrix)[1]
		suppl.data$nr_countries <- dim(suppl.data$tfr_matrix)[2]
	}
    lambda_c = results_tau$lambda_c
    nr_countries = length(lambda_c)
    nr_countries_estimation <- tfr.with.regions$nr_countries_estimation
                                                   
    # uniform prior for U_c, make lower bound country specific
    if(any(apply(updated.tfr.matrix, 2, function(x) all(is.na(x))))) { # some countries start Phase III before 1950
    	tfr_min_c <- c()
 		# loop over countries to find minimum
 		for (country in 1:nr_countries){
    		data <- get.observed.with.supplemental(country, updated.tfr.matrix, suppl.data)
    		tfr_min_c <- c(tfr_min_c, min(data, na.rm=TRUE))
    	}
    } else
    	tfr_min_c <- apply(updated.tfr.matrix, 2, min, na.rm = TRUE)
    lower_U_c <- ifelse(tfr_min_c > meta$U.c.low.base, tfr_min_c, meta$U.c.low.base)
	prop_cov_gammas <- array(NA, c(nr_countries,3,3))
	if(use.average.gammas.cov) {
		cov.to.average <- get.cov.gammas(sim.dir=meta$output.dir, burnin=burnin)$values
		if (all(is.na(cov.to.average))) {
			warning('Covariance of gamma is NA for all countries. Average from default covariance will be used.', 
						immediate.=TRUE)
			e <- new.env()
			data('proposal_cov_gammas_cii', envir=e)
			cov.to.average <- e$proposal_cov_gammas_cii$values
		}
	} else {
		# get default proposal_cov_gammas_cii and match with country codes of this run
		e <- new.env()
    	data('proposal_cov_gammas_cii', envir=e)
    	current.country.codes <- tfr.with.regions$regions$country_code
    	matched.index <- match(e$proposal_cov_gammas_cii$country_codes, current.country.codes)
    	is.notNA <- !is.na(matched.index)
    	prop_cov_gammas[matched.index[is.notNA],,] <- e$proposal_cov_gammas_cii$values[is.notNA,,]
    
		if (!is.null(proposal_cov_gammas)) { #user-specified, overwrites defaults for given countries
			matched.index <- match(proposal_cov_gammas$country_codes, current.country.codes)
			is.notNA <- !is.na(matched.index)
			prop_cov_gammas[matched.index[is.notNA],,] <- proposal_cov_gammas$values[is.notNA,,]
		}		
		cov.to.average <- prop_cov_gammas
	}
	# where NAs, put averages
	isNA <- apply(is.na(prop_cov_gammas), 1, any)
	if (any(isNA)) {
		avg <- matrix(NA, 3, 3)
		for(i in 1:3)
			avg[,i] <- apply(cov.to.average[,,i], 2, mean, na.rm=TRUE)
		for(is.na.country in 1:sum(isNA)) 
			prop_cov_gammas[(1:nr_countries)[isNA][is.na.country],,] <- avg
	}

	return(list(
			tfr_matrix=updated.tfr.matrix, 
			tfr_matrix_all=tfr.with.regions$tfr_matrix_all,
			tfr_matrix_observed=tfr_matrix_observed,
            tau_c = results_tau$tau_c, lambda_c = lambda_c,
            proposal_cov_gammas_cii = prop_cov_gammas,
            start_c = results_tau$start_c, 
            id_Tistau = results_tau$id_Tistau, 
            id_DL = results_tau$id_DL, 
            id_early = results_tau$id_early,
            id_notearly = results_tau$id_notearly,
            id_notearly_estimation = results_tau$id_notearly[results_tau$id_notearly <= nr_countries_estimation],
			U.c.low=lower_U_c,
            nr_countries=nr_countries,
            nr_countries_estimation=nr_countries_estimation,
            T_end=dim(tfr.with.regions$tfr_matrix)[1], T_end_c=results_tau$T_end_c, 
            regions=tfr.with.regions$regions,
            suppl.data=suppl.data
            ))

}

# ini MCMC for UN estimates

mcmc.ini <- function(chain.id, mcmc.meta, iter=100,
					 S.ini=5, 
					 a.ini=0, 
					 b.ini=a.ini, 
					 sigma0.ini=0.1, 
					 const.ini=1, 				 
					 gamma.ini=1, Triangle_c4.ini = 1.85,
					 d.ini=0.17,
					 save.all.parameters=FALSE,
					 verbose=FALSE
					 ) {
				 		 	
	nr_countries <- mcmc.meta$nr_countries

    ############################################
    # U_c, the starting levels of the decline
    U_c <- runif(nr_countries, mcmc.meta$U.c.low, mcmc.meta$U.up)
    for (country in mcmc.meta$id_notearly){
		U_c[country] = get.observed.tfr(country, mcmc.meta)[mcmc.meta$tau_c[country]]
	}
	
	##############################################
	# NON-CONST SD
    ##############################################
    # for non-constant sd: sd decreases linearly from TFR of f_sd (S in report) to both sides
    ### notation in report:
    ### a = a_sd, b = b_sd, S = f_sd, sigma0 = sigma0, c = const_sd
    S_sd <- S.ini # max SD
    a_sd <- a.ini
    b_sd <- b.ini
    sigma0 <- sigma0.ini
    const_sd <- const.ini

    sd_eps_tau <- mcmc.meta$sd.eps.tau0
    mean_eps_tau <- mcmc.meta$mean.eps.tau0


    ######################################################
    # for tranformed d: dt ~ N(chi, psi^2)
    ######################################################
    # initiate psi and chi
    chi <- mcmc.meta$chi0
    psi <- mcmc.meta$psi0

	# initiate d:
	d_c= rep(d.ini, nr_countries)

	##################################################################
	# alpha and deltas
    alpha <- mcmc.meta$alpha0.p
    delta <- rep(mcmc.meta$delta0, 3)

    ##################################################################
    # gammas
    gamma_ci <- matrix(gamma.ini, nrow=nr_countries, ncol=3)
                
	### for Triangle4's:
    Triangle4 <- mcmc.meta$Triangle4.0
    delta4 <- mcmc.meta$delta4.0
    T_end = mcmc.meta$T_end
   	Triangle_c4 <- rep(Triangle_c4.ini, nr_countries)
    for (country in mcmc.meta$id_DL) {
    	data <- get.observed.tfr(country, mcmc.meta)
    	minf <- min(data, na.rm = TRUE)
        if (minf < mcmc.meta$Triangle_c4.up) {
        	Triangle_c4[country] = max(mcmc.meta$Triangle_c4.low+0.0001, minf)
          }
    }
   	dontsave.pars <- c('add_to_sd_Tc', 'const_sd_dummie_Tc', 'meta')
    if (!save.all.parameters) dontsave.pars <- c(dontsave.pars, 'eps_Tc')
    if (!exists(".Random.seed")) runif(1)	    	
	mcmc <- structure(list(
						meta=mcmc.meta,
                        U_c=U_c, d_c=d_c, gamma_ci=gamma_ci, 
                        Triangle_c4 = Triangle_c4,
                        delta4 = delta4, Triangle4 = Triangle4,
                        alpha=alpha, delta=delta, psi=psi, chi=chi, 
                        a_sd=a_sd, b_sd=b_sd, const_sd=const_sd,
                        S_sd=S_sd, sigma0=sigma0,sd_eps_tau = sd_eps_tau, 
                        mean_eps_tau = mean_eps_tau,
                        d.ini=d.ini, gamma.ini=gamma.ini,Triangle_c4.ini=Triangle_c4.ini,
                        iter=iter, finished.iter=1, length = 1,
                        id=chain.id,
                        output.dir=paste('mc', chain.id, sep=''),
                        traces=0, traces.burnin=0, 
                        rng.state = .Random.seed,
                        compression.type=mcmc.meta$compression.type,
                        dontsave=dontsave.pars
                        ),
                   class='bayesTFR.mcmc')
                   
	##################################################################
	# distortions
	##################################################################
	# note: the eps will always be NA outside (tau_c, lambda-1)!!
	# ini the epsilons
	mcmc$eps_Tc <- get_eps_T_all(mcmc)
	return(mcmc)
}

mcmc.meta.ini.extra <- function(mcmc.set, countries=NULL, my.tfr.file = NULL, 
									my.locations.file=NULL, burnin = 200, verbose=FALSE) {
	update.regions <- function(reg, ereg, id.replace, is.new, is.old) {
		nreg <- list()
		for (name in c('code', 'area_code', 'country_code')) {
			reg[[name]][id.replace] <- ereg[[name]][is.old]
			nreg[[name]] <- c(reg[[name]], ereg[[name]][is.new])
		}
		for (name in c('name', 'area_name', 'country_name')) {
			reg[[name]][id.replace] <- as.character(ereg[[name]])[is.old]
			nreg[[name]] <- c(as.character(reg[[name]]), 
									  as.character(ereg[[name]])[is.new])
		}
		return(nreg)
	}
	meta <- mcmc.set$meta
	#create tfr matrix only for the extra countries
	tfr.with.regions <- set.wpp.extra(meta, countries=countries, 
									  my.tfr.file = my.tfr.file, my.locations.file=my.locations.file, verbose=verbose)
	if(is.null(tfr.with.regions)) return(list(meta=meta, index=c()))
	has.mock.suppl <- FALSE
	if(is.null(tfr.with.regions$suppl.data$regions) && !is.null(meta$suppl.data$regions)) {
		# create mock suppl.data in order to get the right data indices
		nrc <- length(tfr.with.regions$regions$code)
		mock.suppl <- list(regions=tfr.with.regions$regions, 
							tfr_matrix=matrix(NA, nrow=nrow(meta$suppl.data$tfr_matrix), ncol=nrc,
												dimnames=list(rownames(meta$suppl.data$tfr_matrix), NULL)),
							index.to.all.countries=1:nrc,
							index.from.all.countries=1:nrc)
		tfr.with.regions$suppl.data <- mock.suppl
		has.mock.suppl <- TRUE
	}
	Emeta <- do.meta.ini(meta, tfr.with.regions=tfr.with.regions, 
								use.average.gammas.cov=TRUE, burnin=burnin,
						 		verbose=verbose)
			 		
	# join the new meta with the existing one
	is.old <- tfr.with.regions$is_processed
	is.new <- !tfr.with.regions$is_processed
	nold <- sum(is.old)
	nr_countries.all <- meta$nr_countries + Emeta$nr_countries - nold
	if (nold > 0) {
		codes.replace <- tfr.with.regions$regions$country_code[is.old]
		id.replace <- unlist(sapply(codes.replace, get.country.object, meta=meta)['index',])
	} else {id.replace <- c()}
	proposal_cov_gammas_cii <- array(NA, c(nr_countries.all, 3, 3))
	for (i in 1:3) {
		meta$proposal_cov_gammas_cii[id.replace,i,] <- matrix(
								Emeta$proposal_cov_gammas_cii[is.old,i,], ncol=3)
		proposal_cov_gammas_cii[,i,] <- rbind(meta$proposal_cov_gammas_cii[,i,], 
							matrix(Emeta$proposal_cov_gammas_cii[is.new,i,], ncol=3))
	}
	new.meta <- list(proposal_cov_gammas_cii = proposal_cov_gammas_cii,
					 nr_countries=nr_countries.all
					)
					
	for (name in c('tfr_matrix', 'tfr_matrix_all', 'tfr_matrix_observed')) {
		meta[[name]][,id.replace] <- Emeta[[name]][,is.old]
		new.meta[[name]] <- cbind(meta[[name]], Emeta[[name]][,is.new])
	}
	for (name in c('tau_c', 'lambda_c', 'start_c', 'U.c.low', 'T_end_c')) {
		meta[[name]][id.replace] <- Emeta[[name]][is.old]
		new.meta[[name]] <- c(meta[[name]], Emeta[[name]][is.new])
	}
	idx.old <- (1:length(is.old))[is.old]
	idx.new.wo.old <- cumsum(is.new)
	for (name in c('id_Tistau', 'id_DL', 'id_early', 'id_notearly')) {
		is.inold <- is.element(Emeta[[name]], idx.old)
		if(any(is.inold)) {
			codes <- unlist(sapply(Emeta[[name]][is.inold], 
								get.country.object, meta=Emeta, index=TRUE)['code',])
			idx <- unlist(sapply(codes, get.country.object, meta=meta)['index',])
			not.incl <- (1:length(idx))[!is.element(meta$regions$country_code[idx], codes)]
			if (length(not.incl) > 0)
				meta[[name]] <- meta[[name]][-not.incl]
		}
		remove.from.old <- !is.element(idx.old, Emeta[[name]])
		if(any(remove.from.old)) {
			idx <- which(is.element(meta[[name]], id.replace[remove.from.old]))
			meta[[name]] <- meta[[name]][-idx]
		}
		new.meta[[name]] <- meta[[name]]
		idx2 <- !is.inold
		if(any(idx2))
			new.meta[[name]] <- c(new.meta[[name]], idx.new.wo.old[Emeta[[name]][idx2]] + meta$nr_countries)
	}
	new.meta[['regions']] <- update.regions(meta$regions, Emeta$regions, id.replace, is.new, is.old)
	if(!is.null(Emeta$suppl.data$regions) && !has.mock.suppl) {
		suppl.id.replace <- meta$suppl.data$index.from.all.countries[id.replace]
		suppl.id.replace <- suppl.id.replace[!is.na(suppl.id.replace)]
		suppl.is.old <- which(is.old)[which(is.element(meta$suppl.data$index.from.all.countries[id.replace], suppl.id.replace))]
		suppl.old <- Emeta$suppl.data$index.from.all.countries[suppl.is.old]
		suppl.is.new <- which(is.new & !is.na(Emeta$suppl.data$index.from.all.countries))
		suppl.new <- Emeta$suppl.data$index.from.all.countries[suppl.is.new]
		for (name in c('tfr_matrix', 'tfr_matrix_all')) {
			meta$suppl.data[[name]][,suppl.id.replace] <- Emeta$suppl.data[[name]][,suppl.old]
			new.meta$suppl.data[[name]] <- cbind(meta$suppl.data[[name]], Emeta$suppl.data[[name]][,suppl.new])
		}
		suppl.is.old.tmp <- rep(FALSE, Emeta$suppl.data$nr_countries)
		suppl.is.old.tmp[suppl.is.old] <- TRUE
		new.meta$suppl.data$regions <- update.regions(meta$suppl.data$regions, Emeta$suppl.data$regions, 
												suppl.id.replace, suppl.new, suppl.old)
		n.new <- ncol(new.meta$suppl.data$tfr_matrix) - ncol(meta$suppl.data$tfr_matrix)
		new.meta$suppl.data$nr_countries <- ncol(new.meta$suppl.data$tfr_matrix)
		new.meta$suppl.data$T_end <- nrow(new.meta$suppl.data$tfr_matrix)
		new.meta$suppl.data$index.from.all.countries <- meta$suppl.data$index.from.all.countries
		new.meta$suppl.data$index.to.all.countries <- meta$suppl.data$index.to.all.countries
		if (n.new > 0) {
			new.meta$suppl.data$index.from.all.countries <- c(new.meta$suppl.data$index.from.all.countries, rep(NA, sum(is.new)))
			new.meta$suppl.data$index.from.all.countries[meta$nr_countries + suppl.is.new] <- seq(meta$suppl.data$nr_countries + 1, 
												length=n.new)
			new.meta$suppl.data$index.to.all.countries <- c(new.meta$suppl.data$index.to.all.countries, 
											seq(meta$nr_countries+1, new.meta$nr_countries)[suppl.is.new])
		} 
	}
	index <- id.replace
	if (new.meta$nr_countries > meta$nr_countries) 
		index <- c(index, seq(meta$nr_countries+1, new.meta$nr_countries))
	for (item in names(new.meta)) {
		meta[[item]] <- new.meta[[item]]
	}

	return(list(meta=meta, index=index, index.replace=id.replace, 
				index_DL=index[is.element(index, new.meta$id_DL)]))
}

mcmc.ini.extra <- function(mcmc, countries, index.replace=NULL) {
	nr.countries.extra <- length(countries)
	nreplace <- length(index.replace)
	Uc <- runif(nr.countries.extra, mcmc$meta$U.c.low, mcmc$meta$U.up)
	if(nreplace > 0) {
		mcmc$U_c[index.replace] <- Uc[1:nreplace]
		mcmc$d_c[index.replace] <- mcmc$d.ini
		mcmc$Triangle_c4[index.replace] <- mcmc$Triangle_c4.ini
		mcmc$gamma_ci[index.replace,] <- matrix(mcmc$gamma.ini, nrow=nreplace, ncol=3)
	}
	U_c <- mcmc$U_c
	if(nr.countries.extra > nreplace)
		U_c <- c(U_c, Uc[(nreplace+1):nr.countries.extra])
    for (country in mcmc$meta$id_notearly[is.element(mcmc$meta$id_notearly, countries)]){
		U_c[country] = get.observed.tfr(country, mcmc$meta)[mcmc$meta$tau_c[country]]
	}
	mcmc.update <- list(U_c=U_c, d_c=c(mcmc$d_c, rep(mcmc$d.ini, nr.countries.extra-nreplace)), 
						gamma_ci=rbind(mcmc$gamma_ci, 
									   matrix(mcmc$gamma.ini, nrow=nr.countries.extra-nreplace, ncol=3)), 
                        Triangle_c4 = c(mcmc$Triangle_c4, rep(mcmc$Triangle_c4.ini, nr.countries.extra-nreplace))
                        )
	for (item in names(mcmc.update)) {
		mcmc[[item]] <- mcmc.update[[item]]
	}
	return(mcmc)
}