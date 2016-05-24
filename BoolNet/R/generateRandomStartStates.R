
generateRandomStartStates <- function(network, n)
{
  mat <- matrix(nrow=n,ncol=length(network$genes))
  fixedPositions <- which(network$fixed != -1)
  nonFixedPositions <- which(network$fixed == -1)
  
  if (n > (2 ^ length(nonFixedPositions)))
    stop("The number of states to generate exceeds the total number of possible states!")

  if (length(fixedPositions) != 0)
    # fill fixed positions with the corresponding values
    mat[,fixedPositions] <- sapply(fixedPositions,function(x)
            rep(network$fixed[x],n))

  if (n != 2 ^ length(nonFixedPositions))
  {  
    # generate other positions randomly
    mat[,nonFixedPositions] <- round(runif(n=n*length(nonFixedPositions)))
  }
  else
  {
    mat[,nonFixedPositions] <- allcombn(2,length(nonFixedPositions)) - 1
  }
  
  # eliminate duplicates
  mat <- unique(mat)

  while (nrow(mat) != n)
  # if duplicates were removed, generate new states until the
  # desired number of states is reached
  {
    vec <- rep(0,length(network$genes))
    if (length(fixedPositions) != 0)
      # fill fixed positions
      vec[fixedPositions] <- sapply(fixedPositions,
                  function(x)network$fixed[x])
    # generate other positions randomly
    vec[nonFixedPositions] <- round(runif(n=length(nonFixedPositions)))
    mat <- unique(rbind(mat,vec))
  }
  #cat("Using initial states:\n")

  # print states and form a list
  res <- lapply(1:nrow(mat),function(i)
    {
      #cat(paste(mat[i,],collapse=""),"\n",sep="")
      mat[i,]
    })
  #cat("\n")
return(res);
}
