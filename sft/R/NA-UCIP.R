estimateNAH <- function(RT, CR=NULL) {
  nt <- length(RT)

  if ( is.null(CR) | length(CR) != nt ) {
    CR <- rep(1, nt)
  }

  RTx <- sort(RT,index.return=TRUE)
  RT <- RTx$x
  CR <- as.logical(CR)[RTx$ix]

  Y <- rep(NA, nt)
  for (i in 1:nt ) { Y[i] <- sum(RT >= RT[i]) }

  H <- stepfun( RT[CR], c(0,cumsum(1/Y[CR]  )))
  H.v<-stepfun( RT[CR], c(0,cumsum(1/Y[CR]^2)))
  return(list(H=H, Var=H.v))
}
        

estimateNAK <- function(RT, CR=NULL) {
  nt <- length(RT)

  if ( is.null(CR) | length(CR) != nt ) {
    CR <- rep(1, nt)
  }

  RTx <- sort(RT,index.return=TRUE)
  RT <- RTx$x
  CR <- as.logical(CR)[RTx$ix]

  G <- rep(NA, nt)
  for (i in 1:nt ) { G[i] <- sum(RT <= RT[i]) }

  K <- stepfun( RT[CR], c(rev(- cumsum( rev(1/G[CR])   )),0), right=TRUE)
  K.v<-stepfun( RT[CR], c(rev(  cumsum( rev(1/G[CR]^2) )),0), right=TRUE)
  return(list(K=K, Var=K.v))
}


estimateUCIPor <- function(RT, CR=NULL) {
    allRT <- sort(c(RT, recursive=TRUE))
    ncond <- length(RT)
    nt <- length(allRT)

    if ( is.null(CR) | length(CR) != length(RT) ){
      CR <- vector("list", length(RT))
    }

    #Y <- matrix(NA, nt, ncond)
    #Ystepfun <- vector("list", ncond)
    #H <- vector("list", ncond)
    #H.v <- vector("list", ncond)

    #for( i in 1:ncond )  {
    #    sorted <- sort(RT[[i]], index.return=TRUE)
    #    RT[[i]] <- sorted$x
    #    CR[[i]] <- CR[[i]][sorted$ix]
    #    CR[[i]] <- as.logical(CR[[i]])
    #}

    Hucip <- rep(0, nt)
    Hucip.v <- rep(0, nt)

    for ( i in 1:ncond ) {
        Hi <- estimateNAH(RT[[i]], CR[[i]])
        Hucip <- Hucip + Hi$H(allRT)
        Hucip.v <- Hucip.v + Hi$Var(allRT)
    }

    Hucip <- stepfun(allRT, c(0, Hucip))
    Hucip.v <- stepfun(allRT, c(0, Hucip.v))

    return(list(H=Hucip, Var=Hucip.v))
}


estimateUCIPand <- function(RT, CR=NULL) {
    allRT <- sort(c(RT, recursive=TRUE))
    ncond <- length(RT)
    nt <- length(allRT)

    if ( is.null(CR) | length(CR) != length(RT) ) {
      CR <- vector("list", length(RT))
    }

    #Y <- matrix(NA, nt, ncond)
    #Ystepfun <- vector("list", ncond)
    #K <- vector("list", ncond)
    #K.v <- vector("list", ncond)

    #for( i in 1:ncond )  {
    #    sorted <- sort(RT[[i]], index.return=TRUE)
    #    RT[[i]] <- sorted$x
    #    CR[[i]] <- CR[[i]][sorted$ix]
    #    CR[[i]] <- as.logical(CR[[i]])
    #}

    Kucip <- rep(0, nt)
    Kucip.v <- rep(0, nt)

    for ( i in 1:ncond ) {
        Ki <- estimateNAK(RT[[i]], CR[[i]])
        Kucip <- Kucip + Ki$K(allRT)
        Kucip.v <- Kucip.v + Ki$Var(allRT)
    }

    Kucip <- stepfun(allRT, c(Kucip,0),right=TRUE)
    Kucip.v <- stepfun(allRT, c(Kucip.v,0),right=TRUE)

    return(list(K=Kucip, Var=Kucip.v))
}
