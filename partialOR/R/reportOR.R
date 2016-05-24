# 12/31/2012 11:48:53 AM

reportOR <- function(fit,dd,ci=0.95) {

# reports deviances for the F and H model, LR test, 
# OR, SE of log(OR) and CI based on parameters of both models
  m <- dim(dd)[2] - 2
  za <- qnorm(1/2 +ci/2 )
  fit0 <- fit$fit0
  fitH <- fit$fitH
  fitF <- fit$fitF

  LR <- round(  -2*fitH$value - fitF$dev, 3)  
  P.LR <- round(1-pchisq(LR,m),4)           #  LR-test P-value
  cat("\n"," Partial Odds Ratio estimation","\n")
  cat("\n"," Model      "," Deviance")
  cat("\n"," Null       ",format(round(summary(fit0)$dev,2),nsmall=2,width=9),
      "\n"," Full       ",format(round(fitF$dev,2),nsmall=2,width=9),
      "\n"," Homogen.   ",format(round(-2*fitH$value,2),nsmall=2,width=9),
      "\n"," LR-test of homogeneity: LR =",format(LR,nsmall=2),
             ",df =", m, ", P-value =", format(P.LR,nsmall=3),"\n")

# unadjusted OR
  lcomb <- c(-1,-1,1)  # linear combination coef's
  lor0 <- as.numeric(drop( t(coef(fit0) )%*% lcomb )) 
  t1 <- table(dd$x,dd$y)
  t1[t1==0] <- 0.5
  selor0 <- sqrt(sum(1/t1))
  or0.ci <- round(exp(c(lor0-za*selor0,lor0,lor0+za*selor0)),3)
  
# H-model
  lorh   <- drop(t(lcomb) %*% fitH$par[1:3])  # log OR H-model
  vparh  <- solve( -fitH$hessian )[1:3,1:3]
  selorh <- sqrt(drop (t(lcomb) %*% vparh %*% lcomb ))
  orh.ci <- round(exp(c(lorh-za*selorh,lorh,lorh+za*selorh)),3)

# MH-type estimator under the full model
  vparf <- solve( fitF$Hessian )
  eparf <- fitF$coefficients
  zz <- cbind(1,dd[,-(1:2)])
  h1 <- exp( cbind(0,t(eparf %*% t(zz))) )  # exp(b'zi)
  p  <- h1 / apply(h1,1,sum)                # pi's 00,01,10,11
  tt <- sum(p[,1]*p[,4])
  nn <- sum(p[,2]*p[,3])
  or.mh <- tt/nn
  d01 <- apply(    -2*p[,1]*p[,2]*p[,4]*zz, 2,sum)/tt +
         apply(-p[,2]*p[,3]*(1-2*p[,2])*zz, 2,sum)/nn 
  d10 <- apply(    -2*p[,1]*p[,3]*p[,4]*zz, 2,sum)/tt +
         apply(-p[,2]*p[,3]*(1-2*p[,3])*zz, 2,sum)/nn
  d11 <- apply( p[,1]*p[,4]*(1-2*p[,4])*zz, 2,sum)/tt +
         apply(     2*p[,2]*p[,3]*p[,4]*zz, 2,sum)/nn
  dlogorf <- c(d01,d10,d11)
  selorf <- sqrt(dlogorf %*% vparf %*% dlogorf)
  orf.ci <- round(c(or.mh/exp(za*selorf),or.mh,or.mh*exp(za*selorf)),3)

  cat("\n"," Unadjusted OR:", 
      "\n","   log(OR) =", round(lor0,4),", SE =", round(selor0,4), 
      "\n","   OR =", format(or0.ci[2],nsmall=3),
      ", ", paste(round(100*ci,0),"%-CI: ",sep=""), 
      format(or0.ci[1],nsmall=3),"to", format(or0.ci[3],nsmall=3))      

  cat("\n"," MH-type Full model estimate:", 
      "\n","   log(OR) =", round(log(or.mh),4),", SE =", round(selorf,4), 
      "\n","   OR =", format(orf.ci[2],nsmall=3),
      ", ", paste(round(100*ci,0),"%-CI: ",sep=""), 
      format(orf.ci[1],nsmall=3),"to", format(orf.ci[3],nsmall=3))      

  cat("\n"," Homogeneity-model estimate:",
      "\n","   log(OR) =", round(lorh,4),", SE =", round(selorh,4), 
      "\n","   OR =", format(orh.ci[2],nsmall=3),
      ", ", paste(round(100*ci,0),"%-CI: ",sep=""),
      format(orh.ci[1],nsmall=3),"to ",format(orh.ci[3],nsmall=3),"\n")
  return(invisible(NULL))
} # end reportOR

