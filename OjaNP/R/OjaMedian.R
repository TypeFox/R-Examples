`ojaMedian` <-
function(X, alg="evolutionary", sp=1, na.action=na.fail, control=ojaMedianControl(...), ...){
   x <- y <- 1
   X<-na.action(X)
   if (!all(sapply(X, is.numeric))) 
        stop("'X' must be numeric")
   if (is.data.frame(X)) X <- as.matrix(X)
   
   alg <- match.arg(alg,c("evolutionary", "exact", "bounded_exact", "grid"))
   
   rows <- dim(X)[1]
   cols <- dim(X)[2]
   outvec <- c(1:cols)
   storage.mode(rows) <- "integer"
   storage.mode(cols) <- "integer"
   storage.mode(X) <- "double"
   storage.mode(outvec) <- "double"

  if(alg != "exact") SEED <- sample(1:5000, sp)

  if (alg=="evolutionary"){
    
    icsX <- ics(X, S1 = control$S1, S2 = control$S2, S1args = control$S1args, S2args = control$S2args)
    Z <- as.matrix(ics.components(icsX))
    B.inv <- solve(coef(icsX))
        
    output <- c(rep(0,cols))
    for(i in 1:sp){
    solution <-  .Call("ojaEvo", Z, as.numeric(control$sigmaInit), as.numeric(control$sigmaAda), control$adaFactor, as.numeric(control$iter), 
                        control$useAllSubsets, as.numeric(control$nSubsetsUsed), as.numeric(control$sigmaLog10Dec), 
                       control$storeSubDet);
    solution.X <- as.vector(solution$best %*% t(B.inv))
     output <- output + solution.X
     }
     RES <- output/sp
     }

  else if (alg=="grid"){
    
    icsX <- ics(X, S1 = control$S1, S2 = control$S2, S1args = control$S1args, S2args = control$S2args)
    Z <- as.matrix(ics.components(icsX))
    B.inv <- solve(coef(icsX))
    
    action <- 2
    param4 <- debug <- 0
    output <- c(rep(0,cols))
    for(i in 1:sp){
    solution<-.C("r_oja", rows, cols, Z, vec = outvec, y, as.integer(action), as.double(control$eps), as.double(control$chi2), as.integer(control$samples), as.integer(param4), as.integer(debug))
    solution.X <- as.vector(solution$vec %*% t(B.inv))
    output <- output + solution.X
    }
    RES <- output/sp
  }
  
  else if (alg=="exact"){
    action <- 1
    param2 <- param3 <- param4 <- debug <- 0
 #   debug = 1
    res<-.C("r_oja", rows, cols, X, vec = outvec, y, as.integer(action), as.double(control$maxlines), as.double(param2), as.integer(param3), as.integer(param4), as.integer(debug),1)
    RES <- res$vec
  }  
  else if (alg=="bounded_exact"){
    action <- 6
    param2 <- control$volume
    param3 <- control$boundedExact
    param4 <- debug <- 0
    #debug = 1
    res<-.C("r_oja", rows, cols, X, vec = outvec, y, as.integer(action), as.double(control$maxlines), as.double(param2), as.integer(param3), as.integer(param4), as.integer(debug),1)
    RES <- res$vec
  }
  names(RES)<-colnames(X)
  return(RES)
  }
