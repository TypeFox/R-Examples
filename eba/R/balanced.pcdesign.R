balanced.pcdesign <- function(nstimuli){
  # Create a (completely) balanced paired-comparison design
  # Last mod: Sep/30/2009, FW

  t <- nstimuli
  if(t < 3 | t > 26) stop("nstimuli must lie between 3 and 26")

  if(t %% 2){  # t odd, completely balanced design
    m1 <- matrix(letters[1:t], t, .5*(t + 1), TRUE)
    m2 <- matrix(letters[1:t], t, ncol(m1) - 1, TRUE)
    m2 <- as.matrix(m2[nrow(m2):1, ncol(m2):1])
    X <- cbind(m1[,1], matrix(paste(m2, m1[,-1], sep=""), nrow(m2), ncol(m2)))
    A <- as.character(t(X[,-1]))
  }else{       # t even, design is not completely balanced
    t <- t - 1
    m1 <- matrix(letters[1:t], t, .5*(t + 1), TRUE)
    m2 <- matrix(letters[1:t], t, ncol(m1) - 1, TRUE)
    m2 <- as.matrix(m2[nrow(m2):1, ncol(m2):1])
    fcol <- paste(letters[t + 1], m1[,1], sep="")  # first column
    idx <- seq(2, nrow(m1), 2)
    fcol[idx] <- paste(substr(fcol[idx],2,2), substr(fcol[idx],1,1), sep="")
    X <- cbind(fcol, matrix(paste(m2, m1[,-1], sep=""), nrow(m2), ncol(m2)))
    A <- as.character(t(X))
  }
  dimnames(X) <- NULL
  B <- paste(substr(A, 2, 2), substr(A, 1, 1), sep="")  # reversed list
  
  out <- list(pairs=X, listA=A, listB=B)
  out
}

