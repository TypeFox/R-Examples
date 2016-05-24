pen.ci <- function(est, cov, age=70, base.dist, agemin){
   
  sest <- mvrnorm(n=10000, est, cov) 
  
  
  out <- cbind( 
    apply(sest, 1, penf, age=70, base.dist=base.dist, agemin=agemin, sex=1, mut=1), 
    apply(sest, 1, penf, age=70, base.dist=base.dist, agemin=agemin, sex=0, mut=1), 
    apply(sest, 1, penf, age=70, base.dist=base.dist, agemin=agemin, sex=1, mut=0), 
    apply(sest, 1, penf, age=70, base.dist=base.dist, agemin=agemin, sex=0, mut=0))
  
  est.pen <- 100*c(
    penf(est, age=70, base.dist=base.dist, agemin=agemin, sex=1, mut=1), 
    penf(est, age=70, base.dist=base.dist, agemin=agemin, sex=0, mut=1), 
    penf(est, age=70, base.dist=base.dist, agemin=agemin, sex=1, mut=0), 
    penf(est, age=70, base.dist=base.dist, agemin=agemin, sex=0, mut=0))

  
  se.pen <- sqrt(apply(100*out, 2, var))
  re <- round(rbind(
        Estimate=est.pen, SE=se.pen, apply(out*100, 2, quantile, prob=c(0.025, 0.975), na.rm=T)),2)
  colnames(re) <- c("Male Carrier", "Female Carrier", "Male Noncarrier", "Female Noncarrier")
  
  return(re)
}