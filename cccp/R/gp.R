##
## R function for solving a geometric program
##
gp <- function(F0, g0, FList = list(), gList = list(),
               nno = NULL, A = NULL, b = NULL, optctrl = ctrl()){
    n <- ncol(F0)
    mnl <- length(FList)
    if(!identical(mnl, length(gList))){
        stop("Length of list objects 'FList' and 'gList' differ.\n")
    }
    mnl <- mnl + 1L
    FList <- c(list(F0), FList)
    gList <- c(list(g0), gList)
    gList <- lapply(gList, function(x) as.matrix(x))
    ## Checking for inequality constraints and
    ## creating objects for CONEC (epigraph form)
    G <- matrix(0, nrow = mnl, ncol = n)
    h <- matrix(0, nrow = mnl, ncol = 1)
    dims <- mnl
    cone <- "NLFC"
    K <- 1L
    sidx <- matrix(c(0, mnl - 1L), nrow = 1, ncol = 2)
    if(!is.null(nno)){
        cone <- c(cone, "NNOC")
        K <- 2L
        G <- rbind(G, nno$G)
        h <- rbind(h, nno$h)
        dims <- c(dims, nrow(nno$G))
        sidx <- rbind(sidx, c(mnl, nrow(G) -1L))
    } else {
        if(mnl == 1L){
            warning("No restrictions provided, trying solve().\n")
            ans <- try(solve(F0, -g0))
            if(class(ans) == "try-error"){
                stop("Solving unconstrained objective 'F0 * x + g0' failed.\n")
            } else {
                return(exp(ans))
            }
        }
    }
    G <- cbind(G, 0)
    G[1, ncol(G)] <- -1.0
    cList <- new(CONEC, cone, G, h, sidx, dims, K, n + 1L)
    ##
    ## Checking equality constraints
    ##
    ## checking whether x0 satisfies equality constraints
    if(is.null(A)){
        A <- matrix(0, nrow = 0, ncol = n + 1L)
    } else {
        A <- cbind(A, rep(0, nrow(A)))
    }
    if(is.null(b)){
        b <- matrix(0, nrow = 0, ncol = 1)
    }
    if(is.null(dim(b))){
        b <- matrix(b, ncol = 1)
    }
    ## Calling cpp-routine
    gpp(FList, gList, cList, A, b, optctrl)
}
