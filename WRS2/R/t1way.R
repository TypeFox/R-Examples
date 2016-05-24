t1way <- function(formula, data, tr = 0.2) {

  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   
  
  if (tr==0.5) warning("Comparing medians should not be done with this function!")

  grp <- 1:length(x) 
  
  J <-length(grp)
  h <- vector("numeric",J)
  w <- vector("numeric",J)
  xbar <- vector("numeric",J)
  nv <- NA
  for(j in 1:J){
    xx <- !is.na(x[[j]])
    val <- x[[j]]
    x[[j]] <- val[xx]  # Remove missing values
    nv[j] <- length(x[[j]])
    h[j] <- length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
    w[j] <- h[j]*(h[j]-1)/((length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr))
    xbar[j] <- mean(x[[grp[j]]],tr)
 }
 u <- sum(w)
 xtil <- sum(w*xbar)/u
 A <- sum(w*(xbar-xtil)^2)/(J-1)
 B <- 2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
 TEST <- A/(B+1)
 nu1 <- J-1
 nu2 <- 1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
 sig <- 1-pf(TEST,nu1,nu2)
 result <- list(test = TEST, df1 =nu1, df2 = nu2, p.value = sig, call = cl)
 class(result) <- c("t1way")
 result
}
