marginalRelevance <-
function(x, y)
  {
    e1 <- new.env()
    if(!is.factor(y))stop("error: y is not a factor")
    NF <- dim(x)[2]
    NK <- length(unique(y))
    u <- matrix(0, 1, NF)
    l <- matrix(0, 1, NF)
    xOrdered <- matrix(0, dim(x)[1], NF) 
    bestVars <- matrix(0, 1, NF)
    
    
    
    x_dot_f <- apply(x, 2, mean)
    for (k in 1 : NK){
      u <- u + sum(y == levels(y)[k]) * ((apply(x[y == levels(y)[k], ], 2, mean) - x_dot_f)^2)
      l <- l + apply((x[y == levels(y)[k], ] - kronecker(matrix(1, sum(y == levels(y)[k]), 1), t(as.matrix(apply(x[y == levels(y)[k], ], 2, mean)))))^2, 2, sum)
    }
    e1$score <- u/l
    e1$rank <- rank(-e1$score) #, ties.method = "random"
    for (i in 1 : NF){
      bestVars[i] <- which(e1$rank == i)
      xOrdered[ ,i] <- x[ ,which(e1$rank == i)]
    } 
    e1$bestVars <- bestVars
    e1$orderedData <- xOrdered
    e1 <- as.list(e1)
    class(e1) = "marginalRelevance"
    return(e1)
  }
