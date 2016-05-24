penMSM <- function(type="fused", d, X, PSM1, PSM2, lambda1, lambda2, w, betastart,
                   nu = 0.5, tol = 1e-10, max.iter = 50, trace=TRUE, diagnostics=TRUE, 
                   family = "coxph", poissonresponse = NULL, poissonoffset = NULL, 
                   constant.approx = 1e-8){
  if(type=="lasso"){
    PSM <- PSM1
    lambda <- lambda1
  }
  if(type=="fused"){
    PSM <- rbind(PSM1, PSM2)
    lambda <- c(lambda1, lambda2)
  }
  betahatold <- B <- betastart
  diff <- Inf
  count <- 1
  lh <- l <- 0
  if(family == "coxph"){
    risksetlist <- buildrisksets(entry = d$entry, exit = d$exit, trans = d$trans, event = d$event)
    risksetlist <- risksetlist$Ri
  }
  if(diagnostics){
    Flist <- NULL
    slist <- NULL
    Alist <- NULL}
  if(trace){
    cat("start estimation:\n")
    cat("  .")
  }    
  while( (diff > tol) & (count < (max.iter+1)) ){
    if(family == "coxph"){
      F <- fisherinfo(beta = betahatold, X = X, risksetlist = risksetlist, event = d$event)
      s <- scorevector(beta = betahatold, X = X, risksetlist = risksetlist, event = d$event)
    }
    if(family == "poisson"){
      mu <- as.numeric(poissonoffset * exp(X %*% betahatold))
      F <- fishercpp(Xcpp = X, mucpp = mu)
      s <- scorevectorP(mu = mu, X = X, event = poissonresponse)
    }
    A <- penaltymatrix(lambda = lambda, PSM = PSM, beta = betahatold, w = w, constant = constant.approx)
    if(diagnostics){
      Flist[[count]] <- F
      slist[[count]] <- s
      Alist[[count]] <- A
    }
    M <- -F - A
    M <- svd(M)
    M <- M$v %*% diag(1/M$d) %*% t(M$u)
    betahatnew <- betahatold - nu * M %*% (s - A %*% betahatold)
    ## betahatnew <- betahatold - nu * solve(-F - A) %*% (s - A %*% betahatold)
    diff <- as.numeric(sum(abs(betahatnew - betahatold))/sum(abs(betahatnew)))
    B <- cbind(B, betahatnew)
    betahatold <- betahatnew
    if(trace){
      if(count%%20 != 0){
        cat(".")
      }else{
        cat(paste(" coef: ", paste(round(betahatnew, 2), collapse=", "), sep=""))
        cat("\n  .")
      }
    }
    count <- count+1
  }
  if(trace){
    if(count >= max.iter){
      cat(" estimation stopped because max. number of iterations was reached.")
      cat("\n")
      cat(paste(" relative change in last iteration: ", round(diff, 4), sep=""))
    }
    cat("\n")
  }
  beta <- B[, ncol(B)]
  M <- F + A
  M <- svd(M)
  M <- M$v %*% diag(1/M$d) %*% t(M$u)
  df <- sum(diag(F %*% M))
  if(family == "coxph"){
    nlpl <- -lpl(beta = beta, X = X, risksetlist = risksetlist, event = d$event)
  }
  if(family == "poisson"){
    nlpl <- -llP(beta = beta, X = X, event = poissonresponse, offset = poissonoffset)
  }
  aic <- 2*(nlpl + df)
  N <- nrow(d)
  gcv <- (1/N) * nlpl/(N*((1- df/N)^2))
  returnlist <- list(B = B, aic = aic, gcv = gcv, df = df)
  if(diagnostics){
    returnlist$F <- Flist
    returnlist$s <- slist
    returnlist$A <- Alist
  }
  return(returnlist)
}