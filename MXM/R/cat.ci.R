####################
#### G^2 (and X^2) test of (un)conditional independence
####################

cat.ci <- function(xi, yi, cs, dataset, type = NULL, rob = FALSE, R = 1) {
  ## the xi and yi are two numbers, 1 and 2 for example
  ## indicating the two variables whose conditional independence 
  ## will be tested
  ## xi, yi and cs must be different, non onverlapping numbers
  ## cs is one or more numbers indicating the conditioning variable(s)
  ## it is et to 0 by default. In this case an uncodntional test of 
  ## independence is  performed
  ## dataset is the whole dataset, and is expected to be a matrix
  ## type is set to NULL be default. This argument is not taken into consideration anywhere
  ## rob is FALSE by default, even if it is TRUE it is not taken into cosideration
  ## the type and rob arguments are put here so as to have the same signature as condi
  
  if ( sum(cs == 0) > 0 ) {  ## There are no conditioning variables
    a1 <- chisq.test(dataset[, xi], dataset[, yi], correct = FALSE)  ## chi-square test, but faster than PLL
    stat <- as.numeric( a1$statistic )
    dof <- as.numeric( a1$parameter )
    pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    res <- c( as.numeric(stat), pval, dof )  

  } else {   ## There are conditioning variables
    dat <- cbind( dataset[, c(xi, yi, cs)] ) 
    pa <- ncol(dat)
    colnames(dat) <- paste("V", 1:pa, sep = "")
    xnam <- colnames(dat)[3:pa]
    form <- as.formula( paste("~ V1 + V2 ", paste(xnam, collapse = "+"), sep = "+") )
    mod <- xtabs(form , dat)  ## creates all the contingency tables 
    forma <- as.formula(paste( paste("~", "V1*", paste(xnam, collapse= "*"),
    sep = ""), paste("V2*", paste(xnam, collapse = "*"), sep = ""), sep = "+" ) )
    b1 <- summary( MASS::loglm(forma, mod) )$tests[1, 1:2]  ## PLL model
    
    if ( nrow(dat) > 5 * b1[2] ) {  ## condition to perform the test
      res <- as.numeric( c( b1[1], pchisq(b1[1], b1[2], lower.tail = FALSE, log.p = TRUE), b1[2] ) )
    } else  res <- c( 0, 0, 1 )  
  }
  
  names(res) <- c("Chi-squared test", "logged p-value", "df")
  res
}