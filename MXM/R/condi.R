#####################
#####################
##### Conditional indpendence test for continuous variables
#####
#####################
#####################

condi <- function(ind1, ind2, cs, dat, type = "pearson", rob = FALSE, R = 1) {
  ## ind1 and ind2 are the two indices of the two variables whose correlation is of interest
  ## cs is a vector with the indices of of variable(s),  over which the condition takes place
  ## dat is the data, a matrix form
  ## type is either "pearson" or "spearman"
  ## For a robust estimation of the PEarson correlation set rob = TRUE or FALSE otherwise

  n <- nrow(dat) ## sample size
  d <- sum( cs>0 )  ## dimensionality of cs
  
  ## NOTE: if you set test = "spearman", you must have the ranks of the data in 
  ## the dat argument and not the data themselves. This is to speed up the computations
  if (type == "spearman") { 
    rob = FALSE   ## if spearman is chosen, robust is set to FALSE  
  }

  if (rob == FALSE) {
    if ( d == 0 ) {
      r <- cor(dat[, ind1], dat[, ind2])
    } else {
      tmpm <- dat[, c(ind1, ind2, cs) ]
      tmpm <- as.matrix(tmpm)
      corrMatrix <- cor(tmpm)
      xyIdx <- 1:2
      csIdx <- 3:( d + 2 ) # or csIdx = 3
      residCorrMatrix <- ( corrMatrix[xyIdx, xyIdx] ) - as.matrix( corrMatrix[xyIdx, csIdx] ) %*% 
      ( solve( as.matrix( corrMatrix[csIdx, csIdx] ) , rbind( corrMatrix[csIdx, xyIdx] ) ) )
      r <- abs( residCorrMatrix[1, 2] / sqrt( residCorrMatrix[1, 1] * residCorrMatrix[2, 2]) )
    }
    
  } else {  ## robust estimation using M estimation
    # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
    if ( d == 0 ) {
      # b1 <- coef( robust::lmRob( dat[, ind1] ~ dat[, ind2], control = cont) )[2]
      # b2 <- coef( robust::lmRob( dat[, ind2] ~ dat[, ind1], control = cont) )[2]
       b1 <- coef( MASS::rlm( dat[, ind1] ~ dat[, ind2], maxit = 2000 ) )[2]
       b2 <- coef( MASS::rlm( dat[, ind2] ~ dat[, ind1], maxit = 2000 ) )[2]
      
      r <- sqrt( abs(b1 * b2) )
    } else {
      # e1 <- resid( robust::lmRob( dat[, ind1] ~  dat[, ind2] + dat[, cs], control = cont ) )
      # e2 <- resid( robust::lmRob( dat[, ind2] ~  dat[, ind1] + dat[, cs], control = cont ) )
      e1 <- resid( MASS::rlm( dat[, ind1] ~.,  data = data.frame( dat[, c(ind2, cs) ] ), maxit = 2000 ) )
      e2 <- resid( MASS::rlm( dat[, ind2] ~., data = data.frame( dat[, c(ind1, cs) ] ), maxit = 2000 ) )
      r <- cor(e1, e2)
    }
  }
  
  if (type == "pearson") {
    stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(n - d - 3) )  ## absolute of the test statistic
  } else {
    stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(n - d - 3) ) /  1.029563  ## absolute of the test statistic
  }
  
  pvalue <- log(2) + pt(stat, n - d - 3, lower.tail = FALSE, log.p = TRUE)  ## logged p-value
  result <- c(stat, pvalue, n - d - 3)
  names(result) <- c('test', 'logged.p-value', 'df')
  result
}
