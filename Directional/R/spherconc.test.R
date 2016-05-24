################################
#### ANOVA for spherical data (Test for equality of concentration parameters)
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 226-227
################################

spherconc.test <- function(x, ina) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ni <- as.vector(table(ina))
  x <- as.matrix(x)
  x <- x/sqrt(rowSums(x^2))  ## makes sure x are unit vectors
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  if (p == 3) {
    S <- aggregate(x , by = list(ina), mean)
    Rbi <- sqrt( rowSums(S[, -1]^2) )   ## the mean resultant length of each group
    S <- colMeans(x)
    Rb <- sqrt(sum(S^2))  ## the mean resultant length of all the data
    if (Rb < 0.44) {
      ## case 1
      g1 <- wi <- numeric(g)
      wi <- (5 * (ni - 5))/3
      g1 <- asin(3 * Rbi/sqrt(5))
      U1 <- sum(wi * g1^2) - (sum(wi * g1))^2/sum(wi)
      stat <- U1
      pvalue <- 1 - pchisq(stat, g - 1)
      mess <- paste('The mean resultant length is less than 0.44. U1 was calculated')
    }
    if (Rb >= 0.44 & Rb <= 0.67) {
      ## case 2
      g2 <- wi <- numeric(g)
      wi <- (ni - 4)/0.394
      g2 <- asin((Rbi + 0.176)/1.029)
      U2 <- sum(wi * g2^2) - (sum(wi * g2))^2/sum(wi)
      stat <- U2
      pvalue <- 1 - pchisq(stat, g - 1)
      mess <- paste('The mean resultant length is between 0.44 and 0.67. 
     U2 was calculated')  
    }
    if (Rb > 0.67) {
      ## case 3
      Ri <- Rbi * ni
      vi <- 2 * (ni - 1)
      v <- 2 * (n - g)
      d <- 1/(3 * (g - 1)) * (sum(1/vi) - 1/v)
      U3 <- 1/(1 + d) * (v * log((n - sum(Ri))/v) - sum(vi * log((ni - Ri)/vi)))
      stat <- U3
      pvalue <- 1 - pchisq(U3, g - 1)
      mess <- paste('The mean resultant length is more than 0.67. U3 was calculated')
    }
  } else {
    stat <- NA
    pvalue <- NA
    mess <- paste("This test is valid only for spherical data")
  }
  res <- c(stat, pvalue)
  names(res) <- c('test', 'p-value')
  list(mess = mess, res = res)
}
