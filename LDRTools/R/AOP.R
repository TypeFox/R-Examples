AOP <-
function(x, weights = "constant")
    {
    weights <- match.arg(weights, c("constant", "inverse", "sq.inverse"))
    
    RES <- switch(weights, constant = {AOP.constant(x)}
                         , inverse = {AOP.inverse(x)}
                         , sq.inverse = {AOP.sq.inverse(x)}
                         )
    RES
    }



# subfunction AOP.constant

AOP.constant <- function(x)
    {
    n.mat <- length(x)
    p <- dim(x[[1]])[1]

    # Individual k-estimates
    k <- sapply(x, matrix.trace)
    
    # Check for projectors with rank 0
    if(any(k==0)) { k0 <- sum(k==0)
                    MES <- paste("There are ", k0, " projectors with rank 0")
                    warning(MES)
                    if (k0==n.mat) {return(list(P = matrix(0,p,p), O = matrix(0,p,p), k = 0))}
                    x <- x[k!=0]
                    k <- k[k!=0]
                    }

    n.mat.new <- length(x)

    P.matbar <- Reduce("+",x) / n.mat.new
    E.P.matbar <- eigen(P.matbar)
    # Estimate of k 
    k.estim <- sum(E.P.matbar$values >= 0.5)
    O.mat <- E.P.matbar$vector[,1:k.estim]
    # AOP 
    P.mathat <- O2P(O.mat)
    list(P = P.mathat, O = O.mat, k = k.estim)
    }

# subfunction AOP.inverse
AOP.inverse <- function(x)
    { 
    n.mat <- length(x)
    p <- dim(x[[1]])[1]

    
    # Individual k-estimates
    k <- sapply(x, matrix.trace)
    
    # Check for projectors with rank 0
    if(any(k==0)) { k0 <- sum(k==0)
                    MES <- paste("There are ", k0, " projectors with rank 0")
                    warning(MES)
                    if (k0==n.mat) {return(list(P = matrix(0,p,p), O = matrix(0,p,p), k = 0))}
                    x <- x[k!=0]
                    k <- k[k!=0]
                    }

    n.mat.new <- length(x)
    x.weighted <-  mapply("/", x, k, SIMPLIFY=FALSE)
    P.star <- Reduce("+", x.weighted) / n.mat.new
    EVD.P.star <- eigen(P.star)
    f.criterion <- (1 - 2*cumsum(EVD.P.star$values))/ (2*(1:p))

    # Estimate of k 
    k.estim <- which.min(f.criterion)

    O.mat <- EVD.P.star$vector[,1:k.estim]
    # AOP
    P.mathat <- O2P(O.mat)
    list(P = P.mathat, O = O.mat, k = k.estim)
    }

# subfunction AOP.sq.inverse
AOP.sq.inverse <- function(x)
    { 
    n.mat <- length(x)
    p <- dim(x[[1]])[1]

    # Individual k-estimates
    k <- sapply(x, matrix.trace)
    
    # Check for projectors with rank 0.
    if(any(k==0)) { k0 <- sum(k==0)
                    MES <- paste("There are ", k0, " projectors with rank 0")
                    warning(MES)
                    if (k0==n.mat) {return(list(P = matrix(0,p,p), O = matrix(0,p,p), k = 0))}
                    x <- x[k!=0]
                    k <- k[k!=0]
                    }
        
    n.mat.new <- length(x)            
    x.weighted <-  mapply("/", x, sqrt(k), SIMPLIFY=FALSE)
    P.star <- Reduce("+", x.weighted) / n.mat.new
    EVD.P.star <- eigen(P.star)
    f.criterion <- cumsum(EVD.P.star$values) / sqrt(1:p)

    # Estimate of k
    k.estim <- which.max(f.criterion)
    
    O.mat <- EVD.P.star$vector[,1:k.estim]
    # AOP
    P.mathat <- O2P(O.mat)
    list(P = P.mathat, O = O.mat, k = k.estim)
    }
