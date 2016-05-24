# This is code from van der Laan, Mark J. and Gruber, Susan, "Targeted Minimum Loss Based Estimation of an Intervention Specific Mean Outcome" (August 2011). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 290.
#biostats.bepress.com/ucbbiostat/paper290
# It is here for testing only. It handles a subset of the cases handled by the package.

# Susan Gruber, sgruber@berkeley.edu
# October 27, 2011
# TMLE to estimate an intervention-specific mean outcome
#-----------------------------------------------------------------------------
# function: getEstimates (corrected IC calculation Jan 12, 2011)
# purpose: IPTW, Parametric MLE, TMLE estimates of tmt-specific mean outcome
# arguments:
#   d: dataset, wide format, following the time-ordering of the nodes
#   Inodes: tmt and censoring node  columns in d
#   Lnodes: time-dependent covariate and outcome columns in d
#   Ynodes: intermediate and final event node columns  in d (subset of Lnodes)
#   Qform:  regression formulas for Q_1 through Q_K+1
#   gform:  regression formulas for each treatment and censoring event
#   gbd:   lower bound on estimated probabilities for g-factors
#   family: regression family for initial estimates of Q nodes
#   move.to.weight: [not part of Gruber's original code, added to match ltmle - moves the 1/g from into the weight in the update step]
#-----------------------------------------------------------------------------
ltmle.sg <- function(d, Inodes, Lnodes, Ynodes, Qform, gform, gbd=0, family="quasibinomial", move.to.weight){
  #-----------------------------------------------------------------------------
  # function: estQ
  # purpose: parametric estimation of  Q_k, targeted when h is supplied
  #-----------------------------------------------------------------------------
  estQ <- function(Q.kplus1,d, Qform, uncensored,deterministic, h=NULL, family){ 
    Qform <- update.formula(Qform, Q.kplus1 ~ .)
    m <- glm(as.formula(Qform), data=data.frame(Q.kplus1, d)[uncensored&!deterministic,], 
             family=family)
    Q1W <- predict(m, newdata=d, type="response")
    if(!is.null(h)){
      off <- qlogis(Bound(Q1W, c(.0001, .9999)))
      if (move.to.weight) {
        m <- glm(Q.kplus1 ~ offset(off), data=data.frame(Q.kplus1, h, off),  
                 subset=uncensored & !deterministic, family="quasibinomial", weights=h)
        Q1W  <- plogis(off + coef(m))
      } else {
        m <- glm(Q.kplus1 ~ -1 + h + offset(off), data=data.frame(Q.kplus1, h, off),  
                 subset=uncensored & !deterministic, family="quasibinomial")
        Q1W  <- plogis(off + coef(m)*h)
      }
    }
    Q1W[deterministic] <- 1
    return(Q1W)
  }
  
  #-----------------------------------------------------------------------------
  # function: estg
  # purpose: parametric estimation of each g-factor
  #-----------------------------------------------------------------------------
  estg <- function(d, form, Inodes, Ynodes){
    n <- nrow(d)
    n.g <- length(form)
    gmat <- matrix(NA, nrow=n, ncol=n.g)
    uncensored <- rep(TRUE, n)
    deterministic <- rep(FALSE, n)
    for (i in 1:n.g) {
      if(any(Ynodes < Inodes[i])){
        Ynode.prev <- max(Ynodes[Ynodes < Inodes[i]])
        deterministic <- d[,Ynode.prev]==1
      }       
      m <- glm(as.formula(form[i]), data=d,subset= uncensored & !deterministic,
               family="binomial")
      gmat[,i] <- predict(m, newdata=d, type="response")
      gmat[deterministic,i] <- 1
      uncensored <- d[,Inodes[i]] == 1
    }       
    return(gmat)
  }
  
  #-----------------------------------------------------------------------------
  # function: bound
  # purpose: truncate values within supplied bounds
  #-----------------------------------------------------------------------------
  Bound <- function(x, bounds){
    x[x<min(bounds)] <- min(bounds)
    x[x>max(bounds)] <- max(bounds)
    return(x)
  }
  
  n <- nrow(d)
  n.Q <- length(Lnodes)
  n.g <- length(Inodes)
  
  g1W <- estg(d, gform, Inodes, Ynodes)
  cum.g1W <- Bound(t(apply(g1W, 1, cumprod)), c(gbd,1))
  empirical.meanwt <- mean(1/cum.g1W[,n.g], na.rm=TRUE)
  
  cum.g1W[is.na(cum.g1W)] <- Inf 
  
  
  iptw <- mean(d[,Lnodes[n.Q]] * d[,Inodes[n.g]] / cum.g1W[,n.g])
  wt.n <- 1 / cum.g1W[,n.g] / empirical.meanwt
  iptw.wt.n <- mean(d[,Lnodes[n.Q]] * d[,Inodes[n.g]] * wt.n)
  
  # Gcomp and TMLE
  Qstar <- Qinit <- d[, Lnodes[n.Q]]
  IC <- rep(0, n)
  for (i in n.Q:1){
    Inode.cur  <- which.max(Inodes[Inodes < Lnodes[i]])
    uncensored <- d[,Inodes[Inode.cur]] == 1
    if(any(Ynodes < Lnodes[i])){
      Ynode.prev <- max(Ynodes[Ynodes < Lnodes[i]])
      deterministic <- d[,Ynode.prev]==1
    } else {
      deterministic <- rep(FALSE, n)
    } 	
    Qinit <- estQ(Qinit, d, Qform[i], uncensored, deterministic, 
                  family = family)
    Qstar.kplus1 <- Qstar				
    Qstar <- estQ(Qstar.kplus1, d, Qform[i], uncensored, 
                  deterministic, h = 1/cum.g1W[,Inode.cur], family=family)
    IC[uncensored] <- (IC + (Qstar.kplus1 - Qstar)/
                         cum.g1W[,Inode.cur])[uncensored]
  }
  IC <- IC + Qstar - mean(Qstar)  
  return(c(iptw = iptw, iptw.wt.n=iptw.wt.n,  Gcomp = mean(Qinit), tmle = mean(Qstar), 
           var.tmle = var(IC)/n))
}


