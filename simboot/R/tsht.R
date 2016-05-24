tsht <- 
function(X, f, theta, cmat, conf.level, alternative, R, args)
{ 
  bargs <- args
  XOBS <- as.data.frame(X)
  
  estindsum <- function(X, f, cmat, theta)
    {
      estsum <- theta(X = X, f = f)
      SE <- sqrt(estsum$varest)
      estC <- (cmat %*% estsum$estimate)
      varC <- (cmat^2) %*% (estsum$varest)
      teststat <- estC/ sqrt(varC)
      return(
             list(
                  teststat = teststat,
                  estC = estC,
                  varC = varC,
                  cmat = cmat
                  )
             )
    }
  
  EST <- estindsum(X = XOBS, f = f, cmat = cmat, theta = theta)

  teststat.org <- EST$estC / sqrt(EST$varC)
  
  OBS <- EST$estC
  
  BTeststat <- function(X, i, f, cmat, obs)
    {
      XNEW <- as.data.frame(X[i, ])
      est <- estindsum(X = XNEW, f = f, cmat = cmat, theta = theta)
      Teststat <- (est$estC - obs)/ sqrt(est$varC)
      return(Teststat)
    }
  
  bargs$data <- as.data.frame(X)
  bargs$statistic <- BTeststat
  bargs$strata = f
  bargs$f <- f
  bargs$cmat <- cmat
  bargs$obs <- OBS
  bargs$R <- R
  if(is.null(bargs$R))
    {
      bargs$R <- 999
    }
  if(is.null(bargs$sim))
    {
      bargs$sim <- "ordinary"
    }
  if(is.null(bargs$stype))
    {
      bargs$stype <- "i"
    }
  
  bootout <- do.call("boot", bargs)

  matraw <- matrix( c( teststat.org, bootout$t ), byrow = TRUE, ncol = ncol( bootout$t ), dimnames = NULL)
  
  # teststat<-bootout$t
  alpha <- 1 - conf.level
  switch(alternative, 
         "two.sided" =
         {
           maxabsT <- apply(X = bootout$t, MARGIN = 1, FUN = function(x){
             max(abs(x), na.rm = TRUE)
           })
           count <- sapply( lapply( X = teststat.org, FUN = function( x ){
             maxabsT >= abs( x )
           }), FUN = sum )
           
           countraw <- apply( apply( X = matraw, MARGIN = 2, FUN = function( x ){
             abs( x[2:length( x )] ) >= abs( x[1] )
           }), MARGIN = 2, FUN = sum)
             
           pval <- count / R
           pvalraw <- countraw / R
           
           quant <- quantile(maxabsT, probs = 1 - alpha, na.rm = TRUE)
           LOWER <- EST$estC - quant * sqrt(EST$varC)
           UPPER <- EST$estC + quant * sqrt(EST$varC)
         },
         "less" =
         {
           maxT <- apply(X = bootout$t, MARGIN = 1, FUN = function(x){
             max(x, na.rm = TRUE)
           })

           count <- sapply( lapply( X = teststat.org, FUN = function( x ){
             maxT >= x
           }), FUN = sum )
           
           countraw <- apply( apply( X = matraw, MARGIN = 2, FUN = function( x ){
             x[2:length( x )] >=  x[1]
           }), MARGIN = 2, FUN = sum)
             
           pval <- count / R
           pvalraw <- countraw / R
                        
           quant <- quantile(maxT, probs = 1-alpha, na.rm = TRUE)
           LOWER <- NA
           UPPER <- EST$estC + quant * sqrt(EST$varC)
         },
         "greater" =
         {
           minT <- apply(X = bootout$t, MARGIN = 1, FUN = function(x){
             min(x, na.rm = TRUE)
           })

           count <- sapply( lapply( X = teststat.org, FUN = function( x ){
             minT <= x
           }), FUN = sum )
           
           countraw <- apply( apply( X = matraw, MARGIN = 2, FUN = function( x ){
             x[2:length( x )] <=  x[1]
           }), MARGIN = 2, FUN = sum)
           
           pval <- count / R
           pvalraw <- countraw / R
           
           quant <- quantile(minT, probs = alpha, na.rm = TRUE)
           LOWER <- EST$estC + quant * sqrt(EST$varC)
           UPPER <- NA
         })
  conf.int <- cbind(EST$estC, LOWER, UPPER)
  colnames(conf.int) <- cbind("estimate", "lower", "upper")
  
  p.value <- matrix(c( pval, pvalraw ), ncol = 2, dimnames = list(dimnames(cmat)[[1]], c("adj. p", "raw p")))
  
  return(list(conf.int = conf.int, p.value = p.value, conf.level = conf.level, alternative = alternative))
}

