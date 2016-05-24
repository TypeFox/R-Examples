.extractThetaB <- function(w, Theta)
{
  Names <- names(Theta)
  
  pdf <- unlist(Theta[grep("pdf", Names)])
  
  theta1 <- unlist(Theta[grep("theta1", Names)])
  
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(Theta[grep("theta2", Names)])
  
  theta2[is.na(theta2)] <- 0
  
  c <- length(w); d <- length(pdf) / c

  output <- array(data = list(NULL), dim = c, dimnames = NULL)

  for (j in 1:c) {
    output[[j]]$pdf <- as.character(pdf[((j - 1) * d + 1):(j * d)])
    output[[j]]$theta1 <- as.numeric(theta1[((j - 1) * d + 1):(j * d)])
    output[[j]]$theta2 <- matrix(theta2[((j - 1) * d * d + 1):(j * d * d)], nrow = d, byrow = TRUE)
  }	

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .extractThetaB
