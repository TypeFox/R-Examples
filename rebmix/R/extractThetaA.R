.extractThetaA <- function(w, Theta)
{
  Names <- names(Theta)
  
  pdf <- unlist(Theta[grep("pdf", Names)])
  
  theta1 <- unlist(Theta[grep("theta1", Names)])
  
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(Theta[grep("theta2", Names)])
  
  theta2[is.na(theta2)] <- 0
  
  c <- length(w); d <- length(pdf) / c

  output <- array(data = list(NULL), dim = c(d, c), dimnames = NULL)

  for (i in 1:d) {
    for (j in 1:c) {
      output[[i, j]]$pdf <- as.character(pdf[(j - 1) * d + i])
      output[[i, j]]$theta1 <- as.numeric(theta1[(j - 1) * d + i])
      output[[i, j]]$theta2 <- as.numeric(theta2[(j - 1) * d + i])
    }	
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .extractThetaA
