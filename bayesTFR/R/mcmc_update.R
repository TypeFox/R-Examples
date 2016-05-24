
mcmc.update.abS <- function(what, eps_Tc_temp, mcmc) {
        # 'what' is one of ('a', 'b', 'S')
        
        var.index <- (1:3)[what == c('a', 'b', 'S')]
        abS.values <- list(a=mcmc$a_sd, b=mcmc$b_sd, S=mcmc$S_sd)
        var.value <- abS.values[[var.index]]
        var.low <- mcmc$meta[[paste(what, 'low', sep='.')]]
        var.up <- mcmc$meta[[paste(what, 'up', sep='.')]]
        var.width <- mcmc$meta[[paste(what, 'width', sep='.')]]
        var.name <- paste(what, 'sd', sep='_')
        
        ##############################################################
        # slice sampling starts with getting the z = log_cond full likelihood - exp(1)
        # and then finding a new estimate such that its log_cond_full_likelihood is higher than z
        # the new estimate is found within interval around current estimate
        ##############################################################

        z <- log_cond_abf_sd(mcmc$add_to_sd_Tc, mcmc$const_sd, mcmc$sigma0, eps_Tc_temp, 
                                                mcmc$const_sd_dummie_Tc, mcmc$meta$sigma0.min) - rexp(1)
        v <- runif(1)   
        interval <- c(max(var.value - v*var.width, var.low),
                                  min(var.value + (1-v)*var.width,var.up))

        add_to_sd_Tc_prop <- matrix(NA, nrow(mcmc$add_to_sd_Tc), ncol(mcmc$add_to_sd_Tc))
        
        #while (TRUE){
        for(i in 1:50) {
                var_prop <- runif(1,interval[1], interval[2])
                abS.values[[what]] <- var_prop
                for (country in 1:mcmc$meta$nr_countries){
                        add_to_sd_Tc_prop[1:length(mcmc$data.list[[country]]), country] <- (mcmc$data.list[[country]] - abS.values$S)*
                                ifelse(mcmc$data.list[[country]] > abS.values$S, 
                                -abS.values$a, abS.values$b)
                }  
                like <- log_cond_abf_sd(add_to_sd_Tc_prop, mcmc$const_sd, mcmc$sigma0, eps_Tc_temp, 
                                         mcmc$const_sd_dummie_Tc, mcmc$meta$sigma0.min)     
                if (like >= z) {
                        mcmc[[var.name]] <- var_prop
                        mcmc$add_to_sd_Tc <- add_to_sd_Tc_prop
                        return()
                } else {
                        # shrink interval
                        if (var_prop < var.value){
                                interval[1] <- var_prop
                        } else {
                                interval[2] <- var_prop
                        }
                        if(abs(interval[1]-interval[2]) < 1e-10) break
                } # end else
        }
        mcmc[[var.name]] <- var_prop
        mcmc$add_to_sd_Tc <- add_to_sd_Tc_prop
        warning("Cannot find likelihood increase for ", what, ".\n Final interval: [", interval[1], ',', interval[2], ']\n',
        			"New likelihood: ", like, ", original likelihood: ", z, .immediate=TRUE)
}

mcmc.update.sigma0const <- function(what, log.like.func, eps_Tc_temp, mcmc) {
        # 'what' is one of ('sigma0', 'const')
        var.index <- (1:2)[what == c('sigma0', 'const')]
        var.values <- list(sigma0=mcmc$sigma0, const=mcmc$const_sd)
        var.name <- c('sigma0', 'const_sd')[var.index]
        var.value <- var.values[[var.index]]
        var.low <- mcmc$meta[[paste(what, 'low', sep='.')]]
        var.up <- mcmc$meta[[paste(what, 'up', sep='.')]]
        var.width <- mcmc$meta[[paste(what, 'width', sep='.')]]
        
        z <- eval(call(log.like.func, mcmc$add_to_sd_Tc, var.values[['const']], var.values[['sigma0']], 
                                                eps_Tc_temp, mcmc$const_sd_dummie_Tc, mcmc$meta$sigma0.min)) - rexp(1)
        v <- runif(1)
        interval <- c(max(var.value - v*var.width, var.low),
                                  min(var.value + (1-v)*var.width,var.up))

        #while (TRUE){
        for(i in 1:50) {
                var_prop <- runif(1,interval[1], interval[2])
                var.values[[what]] <- var_prop
                like <- eval(call(log.like.func, mcmc$add_to_sd_Tc, var.values[['const']], var.values[['sigma0']], 
                                 eps_Tc_temp, mcmc$const_sd_dummie_Tc, mcmc$meta$sigma0.min))
                if (like >= z) {
                        mcmc[[var.name]] <- var_prop
                        return()
                } else {
                        # shrink interval
                        if (var_prop < var.value){
                                interval[1] <- var_prop
                        } else {
                                interval[2] <- var_prop
                        }
                        if(abs(interval[1]-interval[2]) < 1e-10) break
                } # end else
        }
        mcmc[[var.name]] <- var_prop
        warning("Cannot find likelihood increase for ", what, ".\n Final interval: [", interval[1], ',', interval[2], ']\n',
        			"New likelihood: ", like, ", original likelihood: ", z, .immediate=TRUE)
}

mcmc.update.abSsigma0const <- function(mcmc, id.not.early.index) {
# updates a_sd, b_sd, f_sd, sigma0, and const_sd
#################################
        
# propose a new a, b, f sigma0 and c (each separately)
# then calculate the new distortions to get its log_cond_like
# then check if you want to accept the new parameter, or shrink the interval and repeat
# if accepted, update the parameters itself and add_to_sd_Tc matrix
# and at the end, update sd_Tc (as this is needed when updating other parameters)

# temp put NAs in eps_ctau's
	eps_Tc_temp = mcmc$eps_Tc
	# put NAs for countries that are not included in the estimation
	eps_Tc_temp[(mcmc$meta$nr_countries_estimation+1):mcmc$meta$nr_countries] <- NA
	
	# put NAs on tau_c spots for id_not_early countries
	eps_Tc_temp[id.not.early.index] <- NA
	
    mcmc.update.abS('a', eps_Tc_temp, mcmc)
    mcmc.update.abS('b', eps_Tc_temp, mcmc)
    mcmc.update.abS('S', eps_Tc_temp, mcmc)
    mcmc.update.sigma0const('sigma0', 'log_cond_sigma0', eps_Tc_temp, mcmc)
    mcmc.update.sigma0const('const', 'log_cond_const_sd', eps_Tc_temp, mcmc)
    mcmc$sd_Tc <- ifelse(mcmc$const_sd_dummie_Tc==1, mcmc$const_sd, 1)*
            ifelse((mcmc$sigma0 + mcmc$add_to_sd_Tc)>0, mcmc$sigma0 + mcmc$add_to_sd_Tc, 
            mcmc$meta$sigma0.min)
}

##############################################################################
mcmc.update.Triangle_c4 <- function(country, mcmc) {
# within country loop, updates Triangle_c4
# thus also the eps
# (thus also lambda_c and the NAs in the eps!)
        Triangle_c4_trans <- log(max(mcmc$Triangle_c4[country] - mcmc$meta$Triangle_c4.low, 1e-20)/
                        max(mcmc$meta$Triangle_c4.up - mcmc$Triangle_c4[country], 1e-20))
        epsT.idx <- mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1)
        lepsT.idx <- length(epsT.idx)
#        z <- (log_cond_Triangle_c4_trans(Triangle_c4_trans, 
#                mcmc$eps_Tc[epsT.idx,country],
#                mcmc$sd_Tc[epsT.idx,country],
#                mcmc$mean_eps_Tc[epsT.idx,country], mcmc$Triangle4, mcmc$delta4) - rexp(1))
        log_cond <- 0.0
        logcondt <- .C("log_cond_Triangle_c4_trans", Triangle_c4_trans, 
                mcmc$eps_Tc[epsT.idx,country],
                mcmc$sd_Tc[epsT.idx,country],
                mcmc$mean_eps_Tc[epsT.idx,country],
                lepsT.idx, mcmc$Triangle4, mcmc$delta4, log_cond=log_cond)
        z <-  logcondt$log_cond - rexp(1)
#		stop('')
        v <- runif(1)

        interval <- c(Triangle_c4_trans - v*mcmc$meta$Triangle_c4.trans.width, 
                                  Triangle_c4_trans + (1-v)*mcmc$meta$Triangle_c4.trans.width)
        
        #while (TRUE){
        for(i in 1:50) {
                Triangle_c4_trans_prop <- runif(1,interval[1], interval[2])
                Triangle_c4_prop <- (mcmc$meta$Triangle_c4.up*exp(Triangle_c4_trans_prop) +mcmc$meta$Triangle_c4.low)/
                                                (1+exp(Triangle_c4_trans_prop))
                        theta_prop <-  c((mcmc$U_c[country]-Triangle_c4_prop)*
                                                exp(mcmc$gamma_ci[country,])/
                                                sum(exp(mcmc$gamma_ci[country,])), 
                                                                Triangle_c4_prop, 
                                                mcmc$d_c[country])
                  eps_T_prop <- get.eps.T(theta_prop,country, mcmc$meta)
                  like <- .C("log_cond_Triangle_c4_trans", Triangle_c4_trans_prop, eps_T_prop,
                			mcmc$sd_Tc[epsT.idx,country],
                			mcmc$mean_eps_Tc[epsT.idx,country], lepsT.idx, mcmc$Triangle4, mcmc$delta4,
                			log_cond=log_cond)$log_cond
                  if (like >= z) {
                        mcmc$eps_Tc[epsT.idx, country] <- eps_T_prop
                        mcmc$Triangle_c4[country] <- Triangle_c4_prop
                        return()
                } else {
                        # shrink interval
                        if (Triangle_c4_prop < mcmc$Triangle_c4[country]){
                                interval[1] <- Triangle_c4_trans_prop
                        } else {
                                interval[2] <- Triangle_c4_trans_prop
                        }
                        if(abs(interval[1]-interval[2]) < 1e-10) break
                } # end else
        }
        mcmc$eps_Tc[epsT.idx, country] <- eps_T_prop
        mcmc$Triangle_c4[country] <- Triangle_c4_prop
        warning("Cannot find likelihood increase for Triangle_c4.\n Final interval: [", interval[1], ',', interval[2], ']\n',
        			"New likelihood: ", like, ", original likelihood: ", z, .immediate=TRUE)
}



############################################
mcmc.update.gamma <- function(country, mcmc) {
# within country loop, update gamma_c's
#################################
# MH-step: propose new vector of gamma's, 
# calculate the new distortions and log-full-cond-likelihood
# accept or reject, if accepted update gamma's and the distortions

#################################
# block update, propose new gamma_ci
#################################
	gamma_prop <- rmvnorm(1, mcmc$gamma_ci[country,], mcmc$meta$proposal_cov_gammas_cii[country,,])
    pci_prob <- exp(gamma_prop)/sum(exp(gamma_prop))
    theta_prop <- c(pci_prob*(mcmc$U_c[country] - mcmc$Triangle_c4[country]), 
                    mcmc$Triangle_c4[country], mcmc$d_c[country]) 
    eps_T_prop <- get.eps.T(theta_prop, country, mcmc$meta)
    idx <- mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1)
    prob_accept <- exp(log_like_gammas(gamma_prop,eps_T_prop,
                                       mcmc$sd_Tc[idx, country], mcmc$mean_eps_Tc[idx, country], 
                                       mcmc$alpha, mcmc$delta) - 
                       log_like_gammas(mcmc$gamma_ci[country,,drop=FALSE],
                                       mcmc$eps_Tc[idx, country], mcmc$sd_Tc[idx, country], 
                                       mcmc$mean_eps_Tc[idx, country], 
                                       mcmc$alpha, mcmc$delta) )
                                
	if (runif(1) < prob_accept){
		mcmc$gamma_ci[country,] <- gamma_prop
		mcmc$eps_Tc[idx, country] <- eps_T_prop
    }
    return()
}



mcmc.update.d <- function(country, mcmc) {
# if accepted, update d_c and the distortions
        d_trans <- log((mcmc$d_c[country] - mcmc$meta$d.low)/(mcmc$meta$d.up - mcmc$d_c[country]))
        z <- (log_cond_d_trans(d_trans, 
                        mcmc$eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country],
                  mcmc$sd_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country],
                        mcmc$mean_eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country],
                         mcmc$psi, mcmc$chi)
                - rexp(1))

        v <- runif(1)
        interval <- c(d_trans - v*mcmc$meta$d.trans.width, 
                                  d_trans + (1-v)*mcmc$meta$d.trans.width)
        
        theta_prop <-  c((mcmc$U_c[country]-mcmc$Triangle_c4[country])*
                                                exp(mcmc$gamma_ci[country,])/
                                                sum(exp(mcmc$gamma_ci[country,])), mcmc$Triangle_c4[country], 
                                                mcmc$d_c[country])
        while (TRUE){
                d_trans_prop <- runif(1,interval[1], interval[2])
                d_prop <- (mcmc$meta$d.up*exp(d_trans_prop) +mcmc$meta$d.low)/(1+exp(d_trans_prop))
                eps_T_prop <- get.eps.T(c(theta_prop[-5], d_prop),country, mcmc$meta)
                if (log_cond_d_trans(d_trans_prop,eps_T_prop, 
                          mcmc$sd_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country],
                          mcmc$mean_eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country],
                                mcmc$psi, mcmc$chi) >= z) {
                        mcmc$eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country] <- eps_T_prop
                        mcmc$d_c[country] <- d_prop
                        return()
                } else {
                        # shrink interval
                        if (d_prop < mcmc$d_c[country]){
                                interval[1] <- d_trans_prop
                        } else {
                                interval[2] <- d_trans_prop
                        }
                } # end else
        }
}


mcmc.update.U <- function(country, mcmc) {
        z  <- (log_cond_U(mcmc$eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country],
                mcmc$sd_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1),country],
                mcmc$mean_eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1),country])
                        - rexp(1))
          
        # find U_prop with g(theta_i*)>z
        # first define interval in which to search
        v <- runif(1)
        interval <- c(max(mcmc$U_c[country] - v*mcmc$meta$U.width, mcmc$meta$U.c.low[country]), 
                                  min(mcmc$U_c[country] + (1-v)*mcmc$meta$U.width,mcmc$meta$U.up))

        theta_current <-  c((mcmc$U_c[country]-mcmc$Triangle_c4[country])*exp(mcmc$gamma_ci[country,])/
                sum(exp(mcmc$gamma_ci[country,])), mcmc$Triangle_c4[country], mcmc$d_c[country])
        theta_prop <- theta_current
        while(TRUE) {
                U_prop <- runif(1,interval[1], interval[2])
                # keep proportions the same, just update the deltas
                theta_prop[1:3] <- theta_current[1:3]/(mcmc$U_c[country] - 
                                                        mcmc$Triangle_c4[country])*(U_prop - mcmc$Triangle_c4[country])
                eps_T_prop <- get.eps.T(theta_prop, country, mcmc$meta)
                if (log_cond_U(eps_T_prop, mcmc$sd_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1),country],
                                                mcmc$mean_eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1),country])
                                 >= z) {
                        mcmc$eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country] <- eps_T_prop
                        mcmc$U_c[country] <- U_prop
                        return()
                } else {
                        # shrink interval
                        if (U_prop < mcmc$U_c[country]){
                                interval[1] <- U_prop
                        } else {
                                interval[2] <- U_prop
                        }
                } # end else, interval was still ok
        }
}
