# # Unit tests for iGWAS
#
#
# set.seed(2)
# J <- 500
# mu <- - 2*abs(rnorm(J))
# frac_sig <- 0.1
# X <- rbinom(J, 1, frac_sig)
# t1 <- rnorm(J, X*mu,1)
# t2 <- rnorm(J, X*mu,1)
# P_current <- pnorm(t2)
# P_prior <- pnorm(t1)
#
# N_current <- 1
# N_prior <- 1
#
# res_unw <- iGWAS(P_current, N_current, P_prior, N_prior, weighting_method="unweighted")
# #sum(res_unw$sig_ind)
# res_w <- iGWAS(P_current, N_current, P_prior, N_prior, sides=1)
#
#
#
