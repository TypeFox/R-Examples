###############################################################################
# print.NPCD:                                                                 #
#                                                                             #
# Define functions to print the various classes of outputs generated from the #
# functions in this package, including AlphaNP, AlphaMLE, JMLE, and Qrefine.  #                                                          #
#                                                                             #
###############################################################################

#install.packages("R.oo")
#library(R.oo)

################################################################################

# Print outputs for Nonparametric Alpha estimation
# inputs: x: the output from AlphaNP
# outputs: print the estimated examinee attribute profiles

setMethodS3("print", "AlphaNP", function(x, ...)  
{
 out <- as.matrix(x$alpha.est)
 rowname.tmp <- colname.tmp <- NULL
 for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Examinee", i))
 for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
 rownames(out) <- rowname.tmp
 colnames(out) <- colname.tmp
 cat("The estimated examinee attribute profiles\n")
 cat(paste(paste("Method:", x$method), "\n"))
 print(out)
}
)

################################################################################

# Print outputs for MLE Alpha estimation
# inputs: x: the output from AlphaMLE
# outputs: print the estimated examinee attribute profiles

setMethodS3("print", "AlphaMLE", function(x, ...)  
{
  out <- as.matrix(x$alpha.est)
  rowname.tmp <- colname.tmp <- NULL
  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Examinee", i))
  for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
  rownames(out) <- rowname.tmp
  colnames(out) <- colname.tmp
  cat("The Estimated Examinee Attribute Profiles\n")
  cat("Method: conditional MLE\n")
  print(out)
}
)

################################################################################

# Print outputs for conditional MLE estimation of the item parameters
# inputs: x: the output from ParMLE
# outputs: print the estimated examinee attribute profiles

setMethodS3("print", "ParMLE", function(x, ...)  
{
  if (x$model %in% c('DINA', 'DINO', 'NIDA'))
  {
  	  out <- cbind(x$slip, x$se.slip, x$guess, x$se.guess)
	  rowname.tmp <- NULL
	  if (x$model %in% c('DINA','DINO'))
	  {
	    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
	  } else if (x$model == 'NIDA')
	  {
	    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Attribute", i))
	  }
	  rownames(out) <- rowname.tmp
	  colnames(out) <- c("slip","SE.slip","guess","SE.guess")
  }

  if (x$model == 'GNIDA')
  {
  	  x$slip[x$Q == 0] <- x$guess[x$Q == 0] <- x$se.slip[x$Q == 0] <- x$se.guess[x$Q == 0] <- NA
  	  out <- cbind(matrix(rbind(x$slip, x$se.slip), nrow=nrow(x$slip)), matrix(rbind(x$guess, x$se.guess), nrow=nrow(x$guess)))
	  rowname.tmp <- NULL
	  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
	  rownames(out) <- rowname.tmp
	  colname.tmp <- NULL
	  for (i in 1:ncol(x$slip)) colname.tmp <- c(colname.tmp, paste("slip", i, sep=''), paste("SE.slip", i, sep=''))
	  for (i in 1:ncol(x$guess)) colname.tmp <- c(colname.tmp, paste("guess", i, sep=''), paste("SE.guess", i, sep=''))
	  colnames(out) <- colname.tmp
  }
    
  if (x$model == 'RRUM')
  {
  	  x$r[x$Q == 0] <- x$se.r[x$Q == 0] <- NA
  	  out <- cbind(x$pi, x$se.pi, matrix(rbind(x$r, x$se.r), nrow=nrow(x$r)))
	  rowname.tmp <- NULL
	  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
	  rownames(out) <- rowname.tmp
	  colname.tmp <- c('pi', 'SE.pi')
	  for (i in 1:ncol(x$r)) colname.tmp <- c(colname.tmp, paste("r", i, sep=''), paste("SE.r", i, sep=''))
	  colnames(out) <- colname.tmp
  }
  
  cat("The Estimated Item Parameters\n")
  cat(paste(paste("Model:", x$model), "\n"))
  cat("Method: conditional MLE\n")
  print(out)
  
}
)

################################################################################

# Print outputs for JMLE estimation of the examinee attribute profiles and item parameters
# inputs: x: the output from JMLE
# outputs: print the estimated examinee attribute profiles and the item parameters

setMethodS3("print", "JMLE", function(x, ...)  
{
  # Item parameters
  if (x$model %in% c('DINA', 'DINO', 'NIDA'))
  {
  	  out <- cbind(x$par.est$slip, x$par.est$se.slip, x$par.est$guess, x$par.est$se.guess)
	  rowname.tmp <- NULL
	  if (x$model %in% c('DINA','DINO'))
	  {
	    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
	  } else if (x$model == 'NIDA')
	  {
	    for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Attribute", i))
	  }
	  rownames(out) <- rowname.tmp
	  colnames(out) <- c("slip","SE.slip","guess","SE.guess")
  }
  
  if (x$model == 'GNIDA')
  {
  	  x$par.est$slip[x$Q == 0] <- x$par.est$guess[x$Q == 0] <- x$par.est$se.slip[x$Q == 0] <- x$par.est$se.guess[x$Q == 0] <- NA
	  out <- cbind(matrix(rbind(x$par.est$slip, x$par.est$se.slip), nrow=nrow(x$par.est$slip)), matrix(rbind(x$par.est$guess, x$par.est$se.guess), nrow=nrow(x$par.est$guess)))
	  rowname.tmp <- NULL
	  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
	  rownames(out) <- rowname.tmp
	  colname.tmp <- NULL
	  for (i in 1:ncol(x$par.est$slip)) colname.tmp <- c(colname.tmp, paste("slip", i, sep=''), paste("SE.slip", i, sep=''))
	  for (i in 1:ncol(x$par.est$guess)) colname.tmp <- c(colname.tmp, paste("guess", i, sep=''), paste("SE.guess", i, sep=''))
	  colnames(out) <- colname.tmp
  }
  
  if (x$model == 'RRUM')
  {
  	  x$par.est$r[x$Q == 0] <- x$par.est$se.r[x$Q == 0] <- NA
  	  out <- cbind(x$par.est$pi, x$par.est$se.pi, matrix(rbind(x$par.est$r, x$par.est$se.r), nrow=nrow(x$par.est$r)))
	  rowname.tmp <- NULL
	  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
	  rownames(out) <- rowname.tmp
	  colname.tmp <- c('pi', 'SE.pi')
	  for (i in 1:ncol(x$par.est$r)) colname.tmp <- c(colname.tmp, paste("r", i, sep=''), paste("SE.r", i, sep=''))
	  colnames(out) <- colname.tmp
  }

  cat(paste(paste("Model:", x$model), "\n"))
  cat("Method: JMLE\n")
  cat(paste(paste("Convergence:", x$conv), "\n"))
  cat(paste(paste("Number of iterations:", x$n.ite), "\n"))
  cat("The Estimated Item Parameters\n")  
  print(out)
  
  # Examinee attribute profiles
  
  out <- as.matrix(x$alpha.est)
  rowname.tmp <- colname.tmp <- NULL
  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Examinee", i))
  for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
  rownames(out) <- rowname.tmp
  colnames(out) <- colname.tmp
  cat("\nThe Estimated Examinee Attribute Profiles\n")
  print(out)
}
)

################################################################################

# Print outputs for Q refinement
# inputs: x: the output from Qrefine
# outputs: print the refined Q-matrix and the modified entries

setMethodS3("print", "Qrefine", function(x, ...)  
{
  out <- as.matrix(x$modified.Q)
  rowname.tmp <- colname.tmp <- NULL
  for (i in 1:nrow(out)) rowname.tmp <- c(rowname.tmp, paste("Item", i))
  for (i in 1:ncol(out)) colname.tmp <- c(colname.tmp, paste("Attribute", i))
  rownames(out) <- rowname.tmp
  colnames(out) <- colname.tmp
  cat("The Modified Q-matrix\n")
  print(out)
  
  cat("\nThe Modified Entries\n")
  print(x$modified.entries) 
}
)

