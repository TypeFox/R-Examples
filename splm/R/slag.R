## spatial lag of object according to listw or matrix

slag <- function(x, listw, maxlag=1, ...) {
    UseMethod("slag")
}

slag.default <- function(x, listw, maxlag=1, index, ...){
    ## needs a vector and a well-specified index
    if(length(x)!=length(index)) {
        stop("Argument and index lengths differ")
    }
    wx <- slagres(x=x, tind=index, listw=listw, maxlag=maxlag, ...)
    return(wx)
}

slag.pseries <- function(x, listw, maxlag=1, ...) {
    ## retrieve index attribute from pseries
    #ind <- attr(x, "index")[,1]
    tind <- attr(x, "index")[,2]

    wx <- slagres(x=x, tind=tind, listw=listw, maxlag=maxlag, ...)

    ## make it a regular pseries
    attr(wx, "index") <- attr(x, "index")
    class(wx) <- c("pseries", class(wx))

    return(wx)
}

slagres <- function(x, tind, listw, maxlag, ...) {
    ## all calculations done inside here
    ## check and if necessary transform
    if(class(listw)[1]=="matrix") {
        listw <- mat2listw(listw, ...)
    }
    ## if maxlag>1 then make higher-order W
    if(maxlag>1) {
        listw <- mat2listw(wlag(listw, maxlag))
    }

    ## unique values
    #unind <- unique(ind)
    tunind <- unique(tind)

    wx <- rep(NA, length(x))

    for(t. in 1:length(tunind)) {
        tpos <- tind==tunind[t.]
        xt <- x[tpos]
        wxt <- lag.listw(listw, xt)
        wx[tpos] <- wxt
    }
    return(wx)
}

wlag<-function(x, maxlag, std=TRUE) {
  ## accepts nb, listw or matrix
  ## returns the proximity matrix of all neighbours up to order=maxlag
  #require(spdep)

  ## convert in neighbours list
  cl1 <- class(x)[1]
  x <- switch(cl1,
              nb={x},
              matrix={mat2listw(x)$neighbours},
              listw={x$neighbours})

  n<-length(x)

  mynb<-nblag(x,maxlag=maxlag)

  mytot<-vector("list",n)

  for(i in 1:n) {
    mytot[[i]]<-mynb[[1]][[i]]
    for(j in 2:maxlag) mytot[[i]]<-c(mytot[[i]],mynb[[j]][[i]])
    ## reorder
    mytot[[i]]<-mytot[[i]][order(mytot[[i]])]
    }

  ## make lagged proximity matrix
  lagmat<-matrix(0,ncol=n,nrow=n)
  for(i in 1:n) lagmat[i,mytot[[i]]]<-1

  ## row-std. if requested
  if(std) lagmat<-lagmat/apply(lagmat,1,sum)

  return(lagmat)
  }





