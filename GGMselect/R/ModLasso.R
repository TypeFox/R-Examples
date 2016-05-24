

calcModLasso <- function(X, Dmax, method, max.st,
                          max.iter, eps,
                         beta,
                         tau,  h,  T0, verbose) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Select a set of candidate edges, ranked by importance
  # INPUT
  #   X: n x p matrix
  #   Dmax: p dim. vector
  #   method: vector taking values in "EW", "LA", "C01"
  #   max.st: scalar ( min(n,p-1))
  #   max.iter : positive integer number (maximal number of
  #                iterations)
  #   eps  : positive real number (precision)
  # OUTPUT
  #   GrGlob: matrix with 3 columns and a variable number of rows
  #   GrGlob[,1] = nodes
  #   GrGlob[,2] = neighbour of the node
  #   GrGlob[,3] = regularisation parameter
  #      (ranked in decreasing order)
  # CALLED BY
  #   calcLarsNEW when family=LA or EW
  # ---------------------------------------------------------------
  
  GrGlob  <- NULL
  p <- dim(X)[2]
  n <- dim(X)[1]
  if (method == "LA")
    betainit <- rep(1,p-1)
  
  for (a in 1:p) {
    modAgarder <- max.st
    Y <- X[,a]
    Z <- X[,-a]

    if (method == "EW")
             betainit <- sapply(EW(x=Z,y=Y,beta,tau,h,T0, max.iter, eps),c)
    
    U <-  t(t(Z)*abs(betainit))

    ll <- lars(U, Y, normalize=FALSE, intercept=FALSE, max.steps=modAgarder, trace= verbose)
    action <- unlist(ll$action)
    sign.act <- sign(action)
    Vois <- ((1:p)[-a])[abs(action)]*sign.act
    switch(method,
           EW =
           {
             modAgarder <-min(modAgarder,length((1:length(action))[cumsum(sign.act) <= Dmax[a]]))
           },
           LA =
           {
             modAgarder <- min(modAgarder,length(action))
           }) #end switch
    
    result <- cbind(rep(a,modAgarder),
                    Vois[1:modAgarder],
                    ll$lambda[1:modAgarder])
    GrGlob <- rbind(GrGlob,result)
  } # end a

  dimnames(GrGlob)[[2]] <- c("a","Vois(a)","lambda")
  # the rows of GrGlob are ordered according to lambda
  ind <- order(GrGlob[,3],decreasing=TRUE)
  GrGlob <- GrGlob[ind,]
  return(GrGlob)
}

