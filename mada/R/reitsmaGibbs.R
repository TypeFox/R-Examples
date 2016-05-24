# #################################################
# ## t_alpha-Reitsma model with Bayesian methods ##
# ## linear mixed model version                  ##
# #################################################
# ## Philipp Doebler, 2012                       ##
# #################################################
# ## philipp.doebler@googlemail.com              ##
# #################################################
# 
# # library(boot)
# # library(glmx)
# # library(mvtnorm)
# # library(MCMCpack)
# # #library(RcppGSL)
# 
# ### Notation from Psychological Methods paper
# 
# ## y_i: number of TP
# ## m_i: numer of positives
# ## z_i: number of FP
# ## n_i: number of negatives
# 
# ## p_i: "true" sensitivity
# ## q_i: "true" fpr
# 
# ## Sigma: covariance matrix
# ## beta: vector of parameters, first two are the mean, 
# ##       then covariate parameters follow
# 
# ## ap: alpha_p transformation parameter for sensitivities
# ## aq: alpha_q transformation parameter for fprs
# 
# ## N: number of primary studies
# 
# ## d2p, d2q: variances from delta method for p_i, q_i
# 
# ## non trivial priors 
# ## for Sigma: inverser Wishart with identity matrix and nu = 3
# 
# reitsmamcmc <- function(data, covariates = NULL, cc = 0.5, 
#                         b= 250, R = 250, thin = 1,
#                         alphasens = 1, alphafpr = 1, nu_wish = 3, V_wish = diag(2),
#                         update_alpha = FALSE, tune = 0.05){
# 
# data <- data + cc
# y <- data$TP
# m <- y + data$FN
# z <- data$FP
# n <- z + data$TN
#   
# N <- length(y)
# 
# ## values for p_i, q_i from standard estimators
# p <- y/m
# q <- z/n
# 
# ap <- alphasens ## initial values
# aq <- alphafpr
# 
# thetap <- talpha(ap)$linkfun(p)
# thetaq <- talpha(aq)$linkfun(q)
# 
# d2p <- (ap*(1-p) - (2-ap))^2/(m*p*(1-p))
# d2q <- (aq*(1-q) - (2-aq))^2/(n*q*(1-q))
# 
# 
# ## start values for parameters with linear model
# mui <- cbind(thetap, thetaq)
# 
# ## puzzle together design matrix for lm.fit
# X.lm <- cbind(rep(1,N)); K <- 1  
# if(!is.null(covariates)){X.lm <- cbind(X.lm,covariates);
#                          K <- 1 + ncol(covariates)}
# ## fit lm
# fit <- lm.fit(X.lm, mui)    
# 
# beta <- as.vector(t(fit$coefficients))
# Sigma <- cov(mui)
# 
# ## setup design matrices for bivariate approach
# X <- array(0, dim = c(N, 2, 2*K))
# X[,1,2*(1:K)-1] <- X.lm
# X[,2,2*(1:K)] <- X.lm
# 
# ## setup arrays for storing
# betaA <- array(NA, dim = c(R,2*K))
# colnames(betaA) <- c("tsens", "tfpr", colnames(covariates))
# SigmaA <- array(NA, dim = c(R,2,2))
# dimnames(SigmaA) <- list(NULL, c("tsens", "tfpr"), c("tsens", "tfpr")) 
# alphaA <- array(NA, dim = c(R,2))
# colnames(alphaA) <- c("alphasens", "alphafpr")
# accept <- 0
# deviance <- array(NA, dim = c(R,1))
# colnames(deviance) <- "deviance"
# 
# ## prepare likelihood quotient for alpha updates
# ## for this note that the likelihood only depends on
# ## the lower part of the hierarchy, i.e. the likelihood of mui
# ## cancels out in the likelihood quotient
# 
# loglik <- function(ap,aq, d2p, d2q){
#   sum(log(ap/p + (2-ap)/(1-p))) + sum(log(aq/q + (2-aq)/(1-q))) +
#     sum(dnorm(talpha(ap)$linkfun(p), mui[,1], sqrt(d2p), log = TRUE)) + 
#     sum(dnorm(talpha(aq)$linkfun(q), mui[,2], sqrt(d2q), log = TRUE))
# }
# 
# ## TODO: need a different log lik for deviance calculations
# 
# ## sample in turn from conditional distributions
# for(j in 1:(b+R)){
#   SigmaInv <- solve(Sigma)
#   
#   ## sample from conditional posterior for mu_i
#   for(i in 1:N){
#     DiInv <- diag(1/c(d2p[i], d2q[i]))
#     Omegai <- solve(DiInv + SigmaInv)
#     nui <- Omegai %*% (DiInv %*% c(thetap[i], thetaq[i]) + SigmaInv %*% as.vector(X[i,,] %*% beta))
#     mui[i,] <- rmvnorm(1, mean = nui, sigma = Omegai)
#   }
#   ## sample from conditional posterior for beta
#   OmegaInv <- matrix(0, ncol = 2*K, nrow = 2*K)  
#   for(i in 1:N){
#   OmegaInv <- OmegaInv + t(X[i,,]) %*% SigmaInv %*% X[i,,]
#   }
# #  Xsum <- apply(X, c(2,3), sum)
#   Omega <- solve(OmegaInv)
#   nu <- numeric(2*K)
#   for(i in 1:N){
#     nu <- nu + t(X[i,,]) %*% SigmaInv %*% mui[i,]
#   }
#   nu <- Omega %*% nu  
#   beta <- as.vector(rmvnorm(1, mean = nu, sigma = Omega))
#   
#   ## sample from conditional posterior for Sigma
#   a <- mui - t(apply(X, 1, function(x){x %*% beta}))
#   A <- t(a) %*% a 
#   Sigma <- riwish(N + nu_wish, A + V_wish)
#   
#   ## perform MH steps for alpha parameters
#   if(update_alpha){    
#     ## try naive uniform proposals, is a symmetric proposal
#     apstar <- runif(1,ap-tune,ap+tune)
#     aqstar <- runif(1,aq-tune,aq+tune)
#     ## calculate likelihood quotient
#     d2pstar <- (apstar*(1-p) - (2-apstar))^2/(m*p*(1-p))
#     d2qstar <- (aqstar*(1-q) - (2-aqstar))^2/(n*q*(1-q))
#     LR <- exp(loglik(apstar, aqstar, d2pstar, d2qstar) - 
#                 loglik(ap,aq,d2p,d2q))  
#     ## if u < LR, then update alphas and derived quantities.
#     if(runif(1) < LR){ap <- apstar
#                       aq <- aqstar
#                       
#                       thetap <- talpha(ap)$linkfun(p)
#                       thetaq <- talpha(aq)$linkfun(q)
#                       
#                       d2p <- (ap*(1-p) - (2-ap))^2/(m*p*(1-p))
#                       d2q <- (aq*(1-q) - (2-aq))^2/(m*q*(1-q))
#                       accept <- accept + 1
#     }  
#   } # end of if(update_alpha)
#   
# if(j %% 100 == 0){cat("current iteration: ", j, "\n")
#                   if(update_alpha){cat("current acceptance rate", 
#                                        accept/j, "\n")}
# }
# ## save current values in Arrays
#   if(j > b){
#     betaA[(j-b),] <- beta
#     SigmaA[(j-b),,] <- Sigma
#     alphaA[(j-b),] <- c(ap,aq)
#     
#   }
# }# end of loop over 1:R+b
# 
# SigmaA <- matrix(SigmaA,ncol =4)
# colnames(SigmaA) <- c("VARtsens","COVtsenstfpr", "COVtsenstfpr", "VARtfpr")
# 
# return(list(betaA = mcmc(betaA, start = b+1, end = R+b, thin = thin),
#             SigmaA = mcmc(SigmaA, 
#                           start = b+1, end = R+b, thin = thin),
#             alphaA = mcmc(alphaA, start = b+1, end = R+b, thin = thin),
#             accept = accept))
# }# end of function