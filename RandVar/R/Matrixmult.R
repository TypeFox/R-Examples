## "Matrixmultiplication"
setMethod("%*%", signature(x = "matrix", y = "EuclRandVariable"),
    function(x, y){
        if(!is.numeric(x))
            stop("x is no numeric matrix")
        if(ncol(x) != dimension(y))
            stop("the number of columns of x != dimension of y")

        fct1 <- NULL;  e <- NULL
        dimn <- y@Range@dimension
        map <- vector("list", nrow(x))
        for(i in 1:nrow(x)){
            if(length(y) == 1){
                if(!identical(all.equal(x[i,1:dimn], numeric(dimn), 
                                tolerance = .Machine$double.eps^0.75), TRUE)){
                    map[[i]] <- function(x){ f1 <- fct1; e %*% f1(x) }
                    body(map[[i]]) <- substitute({ f1 <- fct1; e %*% f1(x) },
                                            list(fct1 = y@Map[[1]], e = x[i,1:dimn]))
                }else{
                    map[[i]] <- function(x){ numeric(1) }
                }
            }else{
                if(!identical(all.equal(x[i, 1:dimn], numeric(dimn), 
                                tolerance = .Machine$double.eps^0.75), TRUE)){
                    fct <- function(x){ f1 <- fct1; e %*% f1(x) }
                    body(fct) <- substitute({ f1 <- fct1; e %*% f1(x) },
                                            list(fct1 = y@Map[[1]], e = x[i,1:dimn]))
                }else{
                    fct <- function(x) { numeric(1) }
                }
                for(j in 2:length(y)){
                    if(!identical(all.equal(x[i, ((j-1)*dimn + 1):(j*dimn)], numeric(dimn), 
                                tolerance = .Machine$double.eps^0.75), TRUE)){
                        body(fct) <- substitute({ f1 <- fct1; f2 <- fct2; f1(x) + e %*% f2(x)},
                                            list(fct1 = fct, fct2 = y@Map[[j]], 
                                                 e = x[i,((j-1)*dimn + 1):(j*dimn)]))
                    }
                }
                map[[i]] <- fct
            }
        }

        R <- new("EuclRandMatrix") 
        R@Map <- map
        R@Dim <- as.integer(c(nrow(x), 1))
        R@Domain <- y@Domain
        R@Range <- Reals()

        return(R)
    })
setMethod("%*%", signature(x = "matrix", y = "EuclRandMatrix"),
    function(x, y){
        dx <- dim(x)
        dy <- y@Dim
        dimn <- y@Range@dimension
        if(dx[2] != dy[1]*dimn)
            stop("non-conformable arguments")
        map <- vector("list", dx[1]*dy[2])

        for(i in 1:dy[2])
            for(j in 1:dx[1])
                map[[(i-1)*dx[1] + j]] <- Map(t(x[j,]) %*% y[,i])[[1]]
        
        y@Map <- map
        y@Dim <- c(dx[1], dy[2])
        y@Range <- Reals()

        return(y)
    })
setMethod("%*%", signature(x = "matrix", y = "EuclRandVarList"),
    function(x, y){
        if(!is.numeric(x))
            stop("x is no numeric matrix")
        if(ncol(x) != dimension(y))
            stop("the number of columns of x != dimension of y")
        
        fct1 <- NULL;  e <- NULL
        map <- vector("list", nrow(x))
        for(i in 1:nrow(x)){
            if(numberOfMaps(y) == 1){
                dimn <- y[[1]]@Range@dimension
                if(!identical(all.equal(x[i,1:dimn], numeric(dimn), 
                                tolerance = .Machine$double.eps^0.75), TRUE)){
                    map[[i]] <- function(x){ f1 <- fct1; e %*% f1(x) }
                    body(map[[i]]) <- substitute({ f1 <- fct1; e %*% f1(x) },
                                            list(fct1 = y[[1]]@Map[[1]], e = x[i,1:dimn]))
                }else{
                    map[[i]] <- function(x){ numeric(1) }
                }
            }else{
                comp <- 0
                ind2 <- 0
                for(k in 1:length(y)){
                    dimn <- y[[k]]@Range@dimension
                    for(j in 1:length(y[[k]])){
                        comp <- comp + 1
                        ind1 <- ind2 + 1
                        ind2 <- ind1 + dimn - 1
                        if(comp == 1){
                            if(!identical(all.equal(x[i, 1:dimn], numeric(dimn), 
                                            tolerance = .Machine$double.eps^0.75), TRUE)){
                                fct <- function(x){ f1 <- fct1; e %*% f1(x) }
                                body(fct) <- substitute({ f1 <- fct1; e %*% f1(x) },
                                                        list(fct1 = y[[k]]@Map[[1]], e = x[i,1:dimn]))
                            }else{
                                fct <- function(x) { numeric(1) }
                            }
                            next
                        }
                        if(!identical(all.equal(x[i, ind1:ind2], numeric(dimn), 
                                    tolerance = .Machine$double.eps^0.75), TRUE)){
                            body(fct) <- substitute({ f1 <- fct1; f2 <- fct2; f1(x) + e %*% f2(x)},
                                                list(fct1 = fct, fct2 = y[[k]]@Map[[j]], 
                                                     e = x[i, ind1:ind2]))
                        }
                    }
                }
                map[[i]] <- fct
            }
        }

        R <- new("EuclRandMatrix") 
        R@Map <- map
        R@Dim <- as.integer(c(nrow(x), 1))
        R@Domain <- y[[1]]@Domain
        R@Range <- Reals()

        return(R)
    })
setMethod("%*%", signature(x = "numeric", y = "EuclRandMatrix"),
    function(x, y){ t(x) %*% y })

setMethod("%*%", signature(x = "numeric", y = "EuclRandVariable"),
    function(x, y){ t(x) %*% as(y, "EuclRandMatrix") })

setMethod("%*%", signature(x = "EuclRandVariable", y = "matrix"),
    function(x, y){ t(y %*% x) })

setMethod("%*%", signature(x = "EuclRandMatrix", y = "matrix"),
    function(x, y){ t(t(y) %*% t(x)) })

setMethod("%*%", signature(x = "EuclRandVarList", y = "matrix"),
    function(x, y){ t(y %*% x) })

setMethod("%*%", signature(x = "EuclRandMatrix", y = "numeric"),
    function(x, y){ x %*% as.matrix(y) })

setMethod("%*%", signature(x = "EuclRandVariable", y = "numeric"),
    function(x, y){ t(x) %*% as.matrix(y) })

setMethod("%*%", signature(x = "EuclRandVariable", y = "EuclRandVariable"),
    function(x, y){
        nrvalues <- length(x)
        if(nrvalues != length(y))
            stop("non-conformable arguments")
        if(!compatibleDomains(x, y))
            stop("the domains of the two random variables are not compatible")

        fct1 <- NULL;  fct2 <- NULL
        fct <- function(x){ f1 <- fct1; f2 <- fct2; f1(x) %*% f2(x) }
        body(fct) <- substitute({ f1 <- fct1; f2 <- fct2; f1(x) %*% f2(x) },
                                list(fct1 = x@Map[[1]], fct2 = y@Map[[1]]))
        
        if(nrvalues > 1)
            for(i in 2:nrvalues)
                body(fct) <- substitute({ f <- fct; f1 <- fct1; f2 <- fct2; f(x) + f1(x) %*% f2(x) },
                                        list(fct = fct, fct1 = x@Map[[i]], fct2 = y@Map[[i]]))

        if(identical(x@Domain, y@Domain))
            Dom <- x@Domain
        else
            if(is(x@Domain, class(y@Domain)[1]))
                Dom <- x@Domain
            else
                Dom <- y@Domain

        R <- new("EuclRandMatrix") 
        R@Map <- list(fct)
        R@Dim <- as.integer(c(1, 1))
        R@Domain <- Dom
        R@Range <- Reals()

        return(R)
    })
setMethod("%*%", signature(x = "EuclRandMatrix", y = "EuclRandMatrix"),
    function(x, y){
        dx <- x@Dim
        dy <- y@Dim
        if(dx[2] != dy[1])
            stop("non-conformable arguments")
        if(!compatibleDomains(x, y))
            stop("the domains of the two random variables are not compatible")
        if(x@Range@dimension != y@Range@dimension)
            stop("the two random matrices have different ranges")
        map <- vector("list", dx[1]*dy[2])

        for(i in 1:dy[2])
            for(j in 1:dx[1])
                map[[(i-1)*dx[1] + j]] <- Map(x[j,] %*% y[,i])[[1]]
        
        x@Map <- map
        x@Dim <- c(dx[1], dy[2])
        x@Range <- EuclideanSpace(dimension = x@Range@dimension*y@Range@dimension)

        return(x)
    })
setMethod("%*%", signature(x = "EuclRandVariable", y = "EuclRandMatrix"),
    function(x, y){
        if(y@Dim[1] != 1)
            stop("non-conformable arguments")
        if(!compatibleDomains(x, y))
            stop("the domains of the two random variables are not compatible")

        R <- new("EuclRandMatrix") 
        R@Map <- x@Map
        R@Dim <- as.integer(c(length(x), 1))
        R@Domain <- x@Domain
        R@Range <- x@Range

        return(R %*% y) 
    })
setMethod("%*%", signature(x = "EuclRandMatrix", y = "EuclRandVariable"),
    function(x, y){
        nrvalues <- length(y)
        if(nrvalues != x@Dim[2])
            stop("non-conformable arguments")
        if(!compatibleDomains(x, y))
            stop("the domains of the two random variables are not compatible")

        R <- new("EuclRandMatrix") 
        R@Map <- y@Map
        R@Dim <- as.integer(c(nrvalues, 1))
        R@Domain <- y@Domain
        R@Range <- y@Range

        return(x %*% R) 
    })
