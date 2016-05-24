penmodelEM <- function(parms, vbeta, data, design="pop", base.dist="Weibull", method="data", mode="dominant", q=0.02){

  agemin <- attr(data, "agemin")
  if(is.null(agemin)) stop("specify agemin for data")
  
  newdata <- carrierprobgeno(data, method=method, mode=mode, q=q)
  theta = theta0 = c(log(parms), vbeta)
  est0 <- est <- theta
  dd <- lval0 <- lval <- 1
  i <- 0
  lval <- loglikem(est, est0, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=FALSE)

  while(dd>0.00001){
    i <- i+1
    est0 <- est
    lval0 <- lval
    nlm.est <- nlm(loglikem, est0, theta0=est0, data=newdata, design=design, base.dist=base.dist, agemin=agemin, 
                   vec=FALSE, hessian=T)
    lval <- nlm.est$minimum
    est <- nlm.est$estimate
    dd <- abs(lval0-lval)
    #dd <- abs(sum(est-est0))    
    #print(c(i, dd, est))
  }
cat("Iterations = ", i, "\n")
  
  EST <- c(exp(nlm.est$estimate[1:2]), nlm.est$estimate[3:4])
  Var <- try(solve(nlm.est$hessian), TRUE)
  if(!is.null(attr(Var,"class"))) { stop("Model didn't converge.\nTry again with different initial values")
   }
  else{  
  se <- sqrt(diag(Var))
  se.exp <-exp(nlm.est$estimate)*se
  SE <- c(se.exp[1:2], se[3:4])
    
  grad <- numericGradient(loglikem, nlm.est$estimate, theta0=nlm.est$estimate, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=TRUE)
  Jscore <- t(grad)%*%grad
  H <- nlm.est$hessian
  #H <- numericNHessian(llik.retro.NF.noasc.vector, est, dat=cc.dat)
  RobustVar <- Var%*%(Jscore)%*%Var
  RobustSE <- sqrt(diag(RobustVar))
  RobustSE[1:2] <- RobustSE[1:2]*exp(nlm.est$estimate[1:2]) 
  }
    
  
  parms.cov <- Var
  parms.se <- SE
  parms.rcov <- RobustVar
  parms.rse <- RobustSE
  
  names(EST)<- names(parms.se)<-  names(parms.rse) <- c("lambda","rho" , "beta.sex","beta.gene")
  rownames(parms.cov) <- colnames(parms.cov) <-c("lambda","rho" , "beta.sex","beta.gene")
  rownames(parms.rcov) <- colnames(parms.rcov) <-c("lambda","rho" , "beta.sex","beta.gene")
  
  ageonset <- agemin:90
  
  p1 <- penf(nlm.est$estimate, ageonset, sex=1, mut=1, base.dist, agemin)  
  p2 <- penf(nlm.est$estimate, ageonset, sex=0, mut=1, base.dist, agemin)  
  p3 <- penf(nlm.est$estimate, ageonset, sex=1, mut=0, base.dist, agemin)  
  p4 <- penf(nlm.est$estimate, ageonset, sex=0, mut=0, base.dist, agemin)  
  
  pen.est <- pen.ci(nlm.est$estimate, RobustVar, age=70, base.dist, agemin)
  pen70.est <- pen.est[1,]
  pen70.se <- pen.est[2,]
  pen70.ci <- pen.est[3:4,] 
  rownames(pen70.ci) <- c("lowerlimit", "upperlimit")
  
  out <- list(  parms.est=EST, parms.cov=parms.cov, parms.se=parms.se, parms.rse=parms.rse,
                pen70.est=pen70.est, pen70.se=pen70.se, pen70.ci=pen70.ci,
                ageonset=ageonset,  
                pen.maleCarr=p1, pen.femaleCarr=p2, 
                pen.maleNoncarr=p3, pen.femaleNoncarr=p4)
  
  class(out) <- "penmodel"
  attr(out, "design") <- design
  attr(out, "base.dist") <- base.dist
  attr(out, "agemin") <- agemin
  attr(out, "data") <- data
  attr(out, "iterations") <- i
  
  
return(out)
  }
