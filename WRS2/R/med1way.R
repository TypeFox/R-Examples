med1way <-
function(formula, data, iter = 1000) {
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  alpha <- 0.05
  crit <- NA
  SEED <- FALSE
  x <- split(model.extract(mf, "response"), mf[,2])   
  
  grp <- 1:length(x)      
  J <- length(grp)  # The number of groups to be compared
  n <- vector("numeric",J)
  w <- vector("numeric",J)
  xbar <- vector("numeric",J)
  
  for(j in 1:J){
    xx <- !is.na(x[[j]])
    val <- x[[j]]
    x[[j]] <- val[xx]  # Remove missing values
    w[j] <- 1/msmedse(x[[grp[j]]], sewarn = FALSE)^2
    xbar[j]<-median(x[[grp[j]]])
    n[j]<-length(x[[grp[j]]])
 }
 pval <- NA
 u <- sum(w)
 xtil <- sum(w*xbar)/u
 TEST <- sum(w*(xbar-xtil)^2)/(J-1)
 
 if(is.na(crit)){
   temp <- med1way.crit(n,alpha,SEED=SEED,iter=iter,TEST=TEST)
   crit.val <- temp$crit.val
 }
 if(!is.na(crit)) crit.val <- crit
 result <- list(test = TEST, crit.val = crit.val, p.value = temp$p.value, call = cl)
 class(result) <- c("med1way")
 result
}
