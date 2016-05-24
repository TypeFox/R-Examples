joinEQTL <- function(eqtlTemp){
  output <- list()
  # Join the ProbeLoc vectors
  temp <- c()
  for(i in 1:length(eqtlTemp))
  {
      temp <- c(temp,unlist(eqtlTemp[[i]][1]))
  }
  names(temp) <- c()
  output[[1]] <- temp

  # Join the testedGenome matrices
  temp <- matrix(NA,ncol=6,nrow=0)
  for(i in 1:length(eqtlTemp))
  {
      temp <- rbind(temp,eqtlTemp[[i]][2]$TestedSNP)
  }
  output[[2]] <- temp


  # Join the p.values vectors
  temp <- c()
  for(i in 1:length(eqtlTemp))
  {
      temp <- c(temp,unlist(eqtlTemp[[i]][3]))
  }
  names(temp) <- c()
  output[[3]] <- temp
  names(output) <- c("ProbeLoc","TestedSNP","p.values")
  output

}
