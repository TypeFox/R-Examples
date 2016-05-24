CreatePssm <-
function(training.data, var)
{

  training.count <- nrow(training)
  amino <- var[[1]]
  register <- var[[2]]
  amino.count <- length(amino)
  register.count <- length(register)
  
  # Create trimer PSSM
  pssm <- array(c(rep(0,(21*7))), c(2, 21, 7), dimnames=list(c("DIMER", "TRIMER"), amino,letters[1:7]))
  cat("Computing PSSM from user-defined training data\n")
  for(i in 1:training.count)
  {
    if(training$type[i] == "TRIMER")
    {
      seq <- reg <- type <- c()
			type <- "TRIMER"
			seq <- as.matrix(strsplit(as.vector(training$sequence[i]),"")[[1]])
			reg <- as.matrix(strsplit(as.vector(training$register[i]),"")[[1]])
			
			ncol<-0
			for(n in 1:amino.count)
      {
        for(h in 1:register.count)
        {
          sect <- c()
          l.sect <- 0
          sect <- intersect(which(seq == amino[n]), which(reg == register[h]))
          l.sect <- length(sect)
          pssm[2,n,h] <- pssm[2,n,h] + l.sect
        }
      }
    }
    # Create dimer PSSM
    else if(training$type[i]=="DIMER")
    {
      seq <- reg <- type <- c()
			type <- "DIMER"
			seq <- as.matrix(strsplit(as.vector(training$sequence[i]),"")[[1]])
			reg <- as.matrix(strsplit(as.vector(training$register[i]),"")[[1]])
			
			ncol<-0
			for(n in 1:amino.count)
      {
        for(h in 1:register.count)
        {
          sect <- c()
          l.sect <- 0
          sect <- intersect(which(seq == amino[n]), which(reg == register[h]))
          l.sect <- length(sect)
          pssm[1,n,h] <- pssm[1,n,h] + l.sect
        }
      }
    }

  }
		
  return(pssm)
}
