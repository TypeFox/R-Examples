###
### This function calculates the CI of bounds using the B-method
### and Bonferroni given the bootstrap draws 
###

boundsCI <- function(lb.rep, ub.rep, lb.est, ub.est, alpha) {

  reps <- length(lb.rep)
  
  ## bonferroni
  bon.lower <- quantile(lb.rep, alpha/2)
  bon.upper <- quantile(ub.rep, 1-alpha/2)
  
  ## b-method
  bmin.dif <- lb.rep-lb.est
  bmax.dif <- ub.est-ub.rep
  bmin.max.dif <- c(bmin.dif, bmax.dif)
  b.sup <- rep(NA, 2*reps)
  
  for (i in 1:(2*reps)) {
    b.sup[i]<-max(sum(bmin.dif<=bmin.max.dif[i])/reps, #emp. dis func
                    sum(bmax.dif<=bmin.max.dif[i])/reps)
  }
  beta <- quantile(b.sup, 1-alpha)
  b.lower <- lb.est-quantile(bmin.dif, beta)
  b.upper <- ub.est+quantile(bmax.dif, beta)
  names(b.lower) <- names(bon.lower)
  names(b.upper) <- names(bon.upper)
  
  ##output
  return(list(bonferroni = c(bon.lower, bon.upper),
              bmethod = c(b.lower, b.upper)))
}
