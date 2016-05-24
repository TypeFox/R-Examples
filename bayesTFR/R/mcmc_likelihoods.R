####################################################
# Leontine Alkema
####################################
# get the full conditionals for the parameters for which 
# slicing sampling is used, or MH

############################################
# 3. a_sd, b_sd and f_sd, and sigma and const_sd
############################################
log_cond_abf_sd <- function (add_to_sd_Tc, const_sd, sigma0, eps_Tc, const_sd_dummie_Tc, sigma0_min) {
    sd_Tc_prop = ifelse(const_sd_dummie_Tc == 1, const_sd, 1) * 
        ifelse((sigma0 + add_to_sd_Tc) > 0, sigma0 + add_to_sd_Tc, 
            sigma0_min)
    return(sum(dnorm(eps_Tc, mean = 0, sd =  sd_Tc_prop, log = TRUE), na.rm = TRUE))
}


log_cond_sigma0 <- function (add_to_sd_Tc, const_sd, sigma0, eps_Tc, const_sd_dummie_Tc, sigma0_min) {
    sd_Tc_prop = ifelse(const_sd_dummie_Tc == 1, const_sd, 1) * 
        ifelse((sigma0 + add_to_sd_Tc) > 0, sigma0 + add_to_sd_Tc, 
            sigma0_min)
    return(sum(dnorm(eps_Tc, mean = 0, sd =  sd_Tc_prop, log = TRUE), na.rm = TRUE))
}

log_cond_const_sd <- function (add_to_sd_Tc, const_sd, sigma0, eps_Tc, const_sd_dummie_Tc, sigma0_min) {
    sd_Tc_prop = ifelse(const_sd_dummie_Tc == 1, const_sd, 1) * 
        ifelse((sigma0 + add_to_sd_Tc) > 0, sigma0 + add_to_sd_Tc, 
            sigma0_min)
    return(sum(dnorm(eps_Tc, mean = 0, sd =  sd_Tc_prop, log = TRUE), na.rm = TRUE))
}

############################################
# 3. d_c_transformed
############################################
log_cond_d_trans <- function(d_trans, eps_T, sd_eps_T, mean_eps_T, psi, chi){ 
 log_cond_d_trans <- (-1/(2*psi^2) *(d_trans-chi)^2 +
                                sum(dnorm(eps_T, mean = mean_eps_T, sd = sd_eps_T, log = TRUE)))
 return(log_cond_d_trans)
}


############################################
# 4. U_c
############################################
log_cond_U <- function( eps_T, sd_eps_T, mean_eps_T){ 
# early decline countries, thus 
# note: gets eps_t and sd_t only for (1, lambda-1)
# when proposing new U_c,
# keep proportions fixed
# change the deltas such that prop are still the same and add up to 1
  # use normal prio
  # -1/(2*sigmaf^2) *(U- F)^2 +
# these are the early decline countries
return(sum(dnorm(eps_T, mean = mean_eps_T, sd = sd_eps_T, log = TRUE)))
}


############################################
# 5. gammas
############################################

# for block updates
log_like_gammas <- function(gamma, eps_T, sd_eps_T, mean_eps_T, alpha, delta) {
        #log_like_gamma <- (sum(dnorm(eps_T, mean = mean_eps_T, sd = sd_eps_T, log = TRUE)) + 
        #                dmvnorm(gamma, alpha, diag(delta^2), log = TRUE))
        log_like_gamma <- (sum(dnorm(eps_T, mean = mean_eps_T, sd = sd_eps_T, log = TRUE)) + 
                        fastdmvnorm(gamma, alpha, diag(delta^2)))
        return(log_like_gamma)
}

#### the standard dmvnorm is slow because it checks a lot of thing (like is Sigma symmetric)
#### so I replace it with a faster version without checks (it's dangerous though!!!)
# Source: Implementation in R of the Parallel Adaptive Wang-Landau algorithm (rpawl)

fastdmvnorm <- function(x, mu, Sigma){
	# assumed log=TRUE
    distval <- mahalanobis(x, center = mu, cov = Sigma)
    logdet <- sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    return(logretval)
}


###########################
log_cond_sigma0.tau = function(eps_Tc, sd_Tc, mean_eps_Tc){
 return(sum(dnorm(eps_Tc, mean = mean_eps_Tc, sd = sd_Tc ,log = TRUE), na.rm = TRUE))
}

##############################################################################
log_cond_Triangle_c4_trans <- function(Triangle_c4_trans, eps_T, sd_eps_T, mean_eps_T,Triangle4, delta4){ 
 # eps are non-NAs
# stop('')
 log_cond_Triangle_c4_trans <- (
                -1/(2*delta4^2) *(Triangle_c4_trans-Triangle4)^2 +
                sum(dnorm(eps_T, mean = mean_eps_T, sd = sd_eps_T, log = TRUE))         )
 return(log_cond_Triangle_c4_trans)
}

############################################
