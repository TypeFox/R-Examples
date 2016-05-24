qfcharact <- function(loa, flagged, zsc, nfactors, floa, av_rel_coef=0.8) {
  nqsorts <- nrow(loa)
  loa_sq <- loa^2
  #number of loading q-sorts
  nload <- colSums(flagged)
  #Eigenvalues
  eigenvals <- colSums(loa_sq)
  #Total explained variance
  expl_var <- 100*(eigenvals/nqsorts)
  #Reliability
  reliability <- av_rel_coef*nload/(1+(nload-1)*av_rel_coef)
  #Standard Error of Factor Scores
  se_fscores <- sapply(zsc, sd)*sqrt(1-reliability)
  #FACTOR MATRIXES
  #correlation among factors
  f_cor <- cor(zsc)
  #SE of differences
  sed <- matrix(data = NA, nrow = nfactors, ncol = nfactors)
  colnames(sed) <- paste("f", 1:nfactors, sep="")
  row.names(sed) <- paste("f", 1:nfactors, sep="")
  f <- 1
  while (f <= ncol(floa)) {
    g <- 1
    while (g <= ncol(floa)) {
      sed[f,g] <- sqrt(se_fscores[[f]]^2 + se_fscores[[g]]^2)
      g <- g+1
    }
    f <- f+1
  }
  #Bind all together
  f_char <- list()
  f_char[[1]] <- data.frame(cbind(av_rel_coef, nload, eigenvals, expl_var, reliability, se_fscores))
  row.names(f_char[[1]]) <- paste("f",1:ncol(loa), sep="")
  f_char[[2]] <- f_cor
  f_char[[3]] <- sed
  names(f_char) <- cbind("characteristics", "cor_zsc", "sd_dif")
  #cbind("Average reliability coefficient, Number of loading Q-sorts, Eigenvalues, Percentage of explained variance, Composite reliability, Standard error of factor scores", "Correlation coefficients between factors z-scores", "Standard errors of differences")
  return(f_char)
}