rbayesBLP=function(Data, Prior, Mcmc){
#
# Keunwoo Kim 02/06/2014
#
# Purpose: 
#      draw theta_bar and Sigma via hybrid Gibbs sampler (Jiang, Manchanda, and Rossi, 2009)
#
# Arguments:
#    Data
#      X: J*T by H (if IV is used, the last column is endogeneous variable.)
#      share: vector of length J*T
#      J: number of alternatives (excluding outside option)
#      Z: instrumental variables (optional)
#
#    Prior
#      sigmasqR
#      theta_hat
#      A
#      deltabar
#      Ad
#      nu0
#      s0_sq
#      VOmega
#
#    Mcmc
#      R: number of MCMC draws
#      H: number of draws for Monte-Carlo integration
#
#      s: scaling parameter of MH increment
#      cand_cov: var-cov matrix of MH increment
#      (minaccep: lower bound of target range of acceptance rate)
#      (maxaccep: upper bound of target range of acceptance rate)
#
#      theta_bar_initial
#      r_initial
#      tau_sq_initial
#      Omega_initial
#      delta_initial
#
#      tol: convergence tolerance for the contraction mapping
#
# Output:
#      a List of tau_sq (or Omega and delta), 
#      theta_bar, r (equivalent to Sigma) draws, Sigma draws,
#      relative numerical efficiency of r draws, tunned parameters for MH, and
#      acceptance rate
#
  
pandterm=function(message) {stop(message,call.=FALSE)}  
#
# check for data
#
if(missing(Data)) {pandterm("Requires Data argument -- list of X and share")}
if(is.null(Data$X)) {pandterm("Requires Data element X")} else {X=Data$X}
if(is.null(Data$share)) {pandterm("Requires Data element share")} else {share=Data$share}
if(is.null(Data$J)) {pandterm("Requires Data element J")} else {J=Data$J}
if(is.null(Data$Z)) {IV=FALSE; Z=matrix(0); I=1} else {IV=TRUE; I=ncol(Z)}
K=ncol(X)
  
if (length(share) != nrow(X)) {pandterm("Mismatch in the number of observations in X and share")}
T=length(share)/J

#
# check for prior 
#
if(missing(Prior)) {  
  c=50
  sigmasqRoff=1
  sigmasqRdiag=log((1+sqrt(1-4*(2*(c(1:K)-1)*sigmasqRoff-c)))/2)/4
  sigmasqR=c(sigmasqRdiag, rep(1, K*(K-1)/2))
  A=BayesmConstant.A*diag(K)
  theta_hat=rep(0,K)
  nu0=K+1
  s0_sq=1
  deltabar=rep(0,I)
  Ad=BayesmConstant.A*diag(I)
  VOmega=BayesmConstant.BLPVOmega
}
else {
  if(is.null(Prior$sigmasqR)) {
    c=50
    sigmasqRoff=1
    sigmasqRdiag=log((1+sqrt(1-4*(2*(c(1:K)-1)*sigmasqRoff-c)))/2)/4
    sigmasqR=c(sigmasqRdiag, rep(1, K*(K-1)/2))
  } else {
    sigmasqR=Prior$sigmasqR
  }
  if(is.null(Prior$A)) {A=BayesmConstant.A*diag(K)} else {A=Prior$A}
  if(is.null(Prior$theta_hat)) {theta_hat=rep(0,K)} else {theta_hat=Prior$theta_hat}
  if(is.null(Prior$nu0)) {nu0=K+1} else {nu0=Prior$nu0}
  if(is.null(Prior$s0_sq)) {s0_sq=1} else {s0_sq=Prior$s0_sq}
  if(is.null(Prior$deltabar)) {deltabar=rep(0,I)} else {deltabar=Prior$deltabar}
  if(is.null(Prior$Ad)) {Ad=BayesmConstant.A*diag(I)} else {Ad=Prior$Ad}
  if(is.null(Prior$VOmega)) {VOmega=BayesmConstant.BLPVOmega} else {VOmega=Prior$VOmega}
}
  
if(length(sigmasqR) != K*(K+1)/2) pandterm("sigmasqR is of incorrect dimension")
if(sum(dim(A)==c(K,K)) != 2) pandterm("A is of incorrect dimension")
if(length(theta_hat) != K) pandterm("theta_hat is of incorrect dimension")
if((length(nu0) != 1) | (nu0 <=0)) pandterm("nu0 should be a positive number")
if((length(s0_sq) != 1) | (s0_sq <=0)) pandterm("s0_sq should be a positive number")
if(length(deltabar) != I) pandterm("deltabar is of incorrect dimension")
if(sum(dim(Ad)==c(I,I)) != 2) pandterm("Ad is of incorrect dimension")
if(sum(dim(VOmega)==c(2,2)) != 2) pandterm("VOmega is of incorrect dimension")

#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R and H")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$H)) {pandterm("Requires element H of Mcmc")} else {H=Mcmc$H}
if(is.null(Mcmc$initial_theta_bar)) {initial_theta_bar=rep(0,K)} else {initial_theta_bar=Mcmc$initial_theta_bar}
if(is.null(Mcmc$initial_r)) {initial_r=rep(0,K*(K+1)/2)} else {initial_r=Mcmc$initial_r}
if(is.null(Mcmc$initial_tau_sq)) {initial_tau_sq=0.1} else {initial_tau_sq=Mcmc$initial_tau_sq}
if(is.null(Mcmc$initial_Omega)) {initial_Omega=diag(2)} else {initial_Omega=Mcmc$initial_Omega}
if(is.null(Mcmc$initial_delta)) {initial_delta=rep(0,I)} else {initial_delta=Mcmc$initial_tau_sq}
if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}

if(is.null(Mcmc$s)+is.null(Mcmc$cand_cov)==0){
  s=Mcmc$s
  cand_cov=Mcmc$cand_cov
  tuning_auto=FALSE
}
if(is.null(Mcmc$s)+is.null(Mcmc$cand_cov)==1) pandterm("If you want to control tuning parameters, write both parameters.")
if(is.null(Mcmc$s)+is.null(Mcmc$cand_cov)==2){
  s=BayesmConstant.RRScaling/sqrt(K*(K+1)/2)
  cand_cov=diag(c(rep(0.1,K),rep(1,K*(K-1)/2)))
  tuning_auto=TRUE
}

if(is.null(Mcmc$tol)) {tol=BayesmConstant.BLPtol} else {tol=Mcmc$tol}
minaccep=0.3
maxaccep=0.5

if(length(initial_theta_bar)!=K) pandterm("initial_theta_bar is of incorrect dimension")
if(length(initial_r)!=(K*(K+1)/2)) pandterm("initial_r is of incorrect dimension")
if(initial_tau_sq<0) pandterm("initial_tau_sq should be positive")
if(sum(dim(initial_Omega)==c(2,2))!=2) pandterm("initial_Omega is of incorrect dimension")
if(length(initial_delta)!=I) pandterm("initial_delta is of incorrect dimension")
if(nprint<0) { pandterm('nprint must be >=0') }
  
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Data Dimensions:",fill=TRUE)
cat(" ",T," market(time); ",J+1," alternatives (including outside option); ",fill=TRUE)
cat(" ",fill=TRUE)
if (IV==TRUE){
  cat(" ",I," instrumental variable(s) ",fill=TRUE)
  cat(" ",fill=TRUE)
}
cat("Prior Parameters:",fill=TRUE)
cat("  thetahat",fill=TRUE)
print(theta_hat)
cat("  A",fill=TRUE)
print(A)
cat("  sigmasqR",fill=TRUE)
print(sigmasqR)
cat("  nu0",fill=TRUE)
print(nu0)
if (IV==TRUE){
  cat("  VOmega",fill=TRUE)
  print(VOmega)
  cat("  deltabar",fill=TRUE)
  print(deltabar)
  cat("  Ad",fill=TRUE)
  print(Ad)
}
if (IV==FALSE){
  cat("  s0_sq",fill=TRUE)
  print(s0_sq)
}
cat(" ",fill=TRUE)
cat("MCMC Parmameters: ",fill=TRUE)
cat(" ",R," reps; keeping every ",keep,"th draw; printing every ",nprint,"th draw",fill=TRUE)
cat(" ",H," draws for Monte-Carlo integration",fill=TRUE)
cat(" ",fill=TRUE)
cat("Contraction Mapping Tolerance: ",fill=TRUE)
cat(" until max(abs((mu1-mu0)/mu0)) <",tol,fill=TRUE)
cat(" ",fill=TRUE)

if (tuning_auto){
  cat("  automatically tuning parameters of RW M-H increment",fill=TRUE)
  cat(" ",fill=TRUE)
  cat("  target acceptance rate is between ",minaccep*100,"% and ",maxaccep*100,"%",fill=TRUE)
  cat(" ",fill=TRUE)
} else{
  cat("  scaling parameter of RW M-H increment is given as",fill=TRUE)
  print(s)
  cat(" ",fill=TRUE)
  cat("  var-cov matrix of RW M-H increment is given as",fill=TRUE)
  print(cand_cov)
  cat(" ",fill=TRUE)
}

# draw for MC integration
draw <- matrix(rnorm(K*H), K, H)

#
# tuning RW Metropolis-Hastings  
#

# if auto-tuning
complete1 <- 0

initial_theta_bar2 <- initial_theta_bar
initial_r2 <- initial_r
initial_tau_sq2 <- initial_tau_sq
initial_Omega2 <- initial_Omega
initial_delta2 <- initial_delta
rdraws <- NULL
if (tuning_auto){  
  cat("Tuning RW Metropolis-Hastings Increment...",fill=TRUE)
  cat(" ",fill=TRUE)   
  cat("-If acceptance rate < ",minaccep*100,"% => s/5",fill=TRUE)
  cat("-If acceptance rate > ",maxaccep*100,"% => s*3",fill=TRUE)
  cat("-If acceptance rate is ",minaccep*100,"~",maxaccep*100,"% => complete tuning",fill=TRUE)
  cat(" ",fill=TRUE)
  while (complete1==0){
    
    cat("  try s=",s,fill=TRUE)    
    out1 <- bayesBLP_rcpp_loop(IV, X, Z, share, J, T, draw, 500, sigmasqR, A, theta_hat, deltabar, Ad, nu0, s0_sq, VOmega,
                               s^2, cand_cov, initial_theta_bar2, initial_r2, initial_tau_sq2, initial_Omega2, initial_delta2, tol, 1, 0)
    initial_theta_bar2 <- as.vector(out1$thetabardraw[,500])
    initial_r2 <- as.vector(out1$rdraw[,500])
    initial_tau_sq2 <- out1$tausqdraw[500]
    if (IV==TRUE) {initial_Omega2 <- matrix(out1$Omegadraw[,500],2,2)}
    if (IV==TRUE) {initial_delta2 <- as.vector(out1$deltadraw[,500])}
    cat("    acceptance rate is ",out1$acceptrate,fill=TRUE)
    
    if (out1$acceptrate>0.20 & out1$acceptrate<0.80){
      rdraws <- cbind(rdraws, out1$rdraw)   
      cat("    (r draws stored)",fill=TRUE)
    }    
    if (out1$acceptrate<minaccep){      
      s <- s/5
    }else if (out1$acceptrate>maxaccep){      
      s <- s*3
    }else{      
      complete1 <- 1      
      cat(" ",fill=TRUE)
      cat("    (tuning completed.)",fill=TRUE)        
    }
  }  
  
  # scaling tunned var-cov matrix from r draws
  scale_opt <- s*sqrt(diag(cand_cov))
  
  Omega <- cov(t(rdraws))
  scale_Omega <- sqrt(diag(Omega))
  corr_opt <- Omega / (scale_Omega%*%t(scale_Omega))
  
  s <- 1
  cand_cov <- corr_opt * (scale_opt%*%t(scale_opt))  
  
  cat(" ",fill=TRUE)
  cat("Tuning Completed:",fill=TRUE)
  cat("  s=",s,fill=TRUE)    
  cat("  var-cov=",fill=TRUE)
  print(cand_cov)  
  cat(" ",fill=TRUE)
}

#
# main run  
#
cat("Starting Random Walk Metropolis-Hastings Sampler for BLP",fill=TRUE)
out <- bayesBLP_rcpp_loop(IV, X, Z, share, J, T, draw, R, sigmasqR, A, theta_hat, deltabar, Ad, nu0, s0_sq, VOmega,
                          s^2, cand_cov, initial_theta_bar, initial_r, initial_tau_sq, initial_Omega, initial_delta, tol, keep, nprint) 
out$s <- s
out$cand_cov <- cand_cov

attributes(out$tausqdraw)$class=c("bayesm.mat","mcmc")
attributes(out$tausqdraw)$mcpar=c(1,R,keep)
attributes(out$thetabardraw)$class=c("bayesm.mat","mcmc")
attributes(out$thetabardraw)$mcpar=c(1,R,keep)
attributes(out$rdraw)$class=c("bayesm.mat","mcmc")
attributes(out$rdraw)$mcpar=c(1,R,keep)
attributes(out$Sigmadraw)$class=c("bayesm.mat","mcmc")
attributes(out$Sigmadraw)$mcpar=c(1,R,keep)
attributes(out$Omegadraw)$class=c("bayesm.mat","mcmc")
attributes(out$Omegadraw)$mcpar=c(1,R,keep)
attributes(out$deltadraw)$class=c("bayesm.mat","mcmc")
attributes(out$deltadraw)$mcpar=c(1,R,keep)

return(out)
}
