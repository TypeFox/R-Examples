penmodel <- function(parms, vbeta, data, design="pop", base.dist="Weibull"){
  
  agemin <- attr(data, "agemin")
  if(is.null(agemin)) stop("specify agemin for data")

  if(any(is.na(data$mgene))) stop("data include missing genetic information, use penmodelEM function to fit")
  est1 <- nlm(loglik, c(log(parms), vbeta), data=data, design=design, base.dist=base.dist, agemin=agemin, hessian=TRUE)
  
  EST <- c(exp(est1$estimate[1:2]), est1$estimate[3:4])
  Var <- solve(est1$hessian)
  se <- sqrt(diag(Var))
  se.exp <-exp(est1$estimate)*se
  SE <- c(se.exp[1:2], se[3:4])
  
  grad <- numericGradient(loglik.vec, est1$estimate, data=data, design=design, base.dist=base.dist, agemin=agemin)
 Jscore <- t(grad)%*%grad
 H <- est1$hessian
 #H <- numericNHessian(llik.retro.NF.noasc.vector, est, dat=cc.dat)
 RobustVar <- solve(H)%*%(Jscore)%*%solve(H)
 RobustSE <- sqrt(diag(RobustVar))
 RobustSE[1:2] <- RobustSE[1:2]*exp(est1$estimate[1:2]) 
  

  parms.cov <- Var
  parms.se <- SE
  parms.rcov <- RobustVar
  parms.rse <- RobustSE
  
  names(EST)<- names(parms.se)<-  names(parms.rse) <- c("lambda","rho" , "beta.sex","beta.gene")
  rownames(parms.cov) <- colnames(parms.cov) <-c("lambda","rho" , "beta.sex","beta.gene")
  rownames(parms.rcov) <- colnames(parms.rcov) <-c("lambda","rho" , "beta.sex","beta.gene")

  ageonset <- agemin:90

  p1 <- penf(est1$estimate, ageonset, sex=1, mut=1, base.dist, agemin)  
  p2 <- penf(est1$estimate, ageonset, sex=0, mut=1, base.dist, agemin)  
  p3 <- penf(est1$estimate, ageonset, sex=1, mut=0, base.dist, agemin)  
  p4 <- penf(est1$estimate, ageonset, sex=0, mut=0, base.dist, agemin)  

  pen.est <- pen.ci(est1$estimate, RobustVar, age=70, base.dist, agemin)
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
return(out)
  
}