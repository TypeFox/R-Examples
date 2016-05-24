##
## Function for creating an object of reference-class 'DCP'
dcp <- function(x0, f0, g0, h0, cList = list(), nlfList = list(), nlgList = list(), nlhList = list(), A = NULL, b = NULL){
    x0 <- as.matrix(x0)
    n <- nrow(x0)
    K <- length(cList)
    mnl <- length(nlfList)
    ## Checking whether x0 is in the domain of nonlinear objective
    f0Dom <- is.nan(f0(x0))
    if(f0Dom){
        stop("Initial point 'x0' is not in the domain of nonlinear objective 'f0'.\n")
    }
    ##
    ## Checking provided non-linear constraints (if applicable)
    ##
    if(mnl > 0){
        if(!all(unlist(lapply(nlfList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlfList' are functions.\n")
        }
        if(!all(unlist(lapply(nlgList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlgList' are functions.\n")
        }
        if(!all(unlist(lapply(nlhList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlhList' are functions.\n")
        }
        fDom <- unlist(lapply(nlfList, function(fcc) fcc(x0)))
        idxnan <- which(is.nan(fDom))
        if(any(idxnan)){
            stop(paste("Initial point 'x0' is not in the domain of nonlinear convex constraint(s): ", idxnan, ".\n", sep = ""))
        }
        if(length(nlfList) != length(nlgList)){
            stop("Length of lists for nonlinear functions and gradient functions do differ.\n")
        }
        if(length(nlfList) != length(nlhList)){
            stop("Length of lists for nonlinear functions and Hessian functions do differ.\n")
        }
        ## Creating list-object of non-linear constraints, their Gradient and Hessian functions
        nList <- list(c(list(f0), nlfList), c(list(g0), nlgList), c(list(h0), nlhList))
        ## Creating objects related to NLFC
        mnl <- mnl + 1L
        Gnl <- matrix(0, nrow = mnl, ncol = n)
        hnl <- matrix(0, nrow = mnl, ncol = 1)
    } else {
        mnl <- 1L
        Gnl <- matrix(0, nrow = 1, ncol = n)
        hnl <- matrix(0, nrow = 1, ncol = 1)
        nList <- list(list(f0), list(g0), list(h0))
    }
    ##
    ## Checking/defining inequality constraints (epigraph form, right-adding 't')
    ##
    if(K > 0){
        cone <- unlist(lapply(cList, function(x) x[["conType"]]))
        if(!all(cone %in% c("NNOC", "SOCC", "PSDC"))){
            stop("List elements of cone constraints must be either created by calls to:\n'nnoc()', or 'socc()', or 'psdc()'.\n")
        }
        cone <- c("NLFC", cone)
        GList <- c(list(Gnl), lapply(cList, function(x) x[["G"]]))
        hList <- c(list(hnl), lapply(cList, function(x) x[["h"]]))
        dims <- c(mnl, as.integer(unlist(lapply(cList, function(x) x[["dims"]]))))
        K <- K + 1L
        G <- do.call("rbind", GList)
        G <- cbind(G, 0)
        G[1, ncol(G)] <- -1.0
        h <- do.call("rbind", hList)
        ridx <- cumsum(unlist(lapply(GList, nrow)))
        sidx <- cbind(c(0, ridx[-length(ridx)]), ridx - 1)
        cList <- new(CONEC, cone, G, h, sidx, dims, K, n + 1L)
    } else { ## case: no cone constraints, but nonlinear constraints, at least f0
        Gepi <- cbind(Gnl, 0)
        Gepi[1, ncol(Gepi)] <- -1.0
        sidx <- matrix(c(0, nrow(Gepi) - 1L), nrow = 1, ncol = 2) 
        nepi <- n + 1L
        cList <- new(CONEC, "NLFC", Gepi, hnl, sidx, mnl, 1L, nepi)
    } 
    ##
    ## Checking equality constraints
    ##
    ## checking whether x0 satisfies equality constraints
    if(!is.null(A)){
        A <- as.matrix(A)
        eq <- identical(as.vector(A %*% x0 - b), rep(0, nrow(A)))
        if(!eq){
            stop("Initial point 'x0' does not satisfy equality constraints.\n")
        }
    }
    if(is.null(A)){
        A <- matrix(0, nrow = 0, ncol = n + 1)
    } else {
        A <- cbind(A, rep(0, nrow(A)))
    }
    if(is.null(b)){
        b <- matrix(0, nrow = 0, ncol = 1)
    }
    if(is.null(dim(b))){
        b <- matrix(b, ncol = 1)
    } 
    
    ans <- new(DCP,
               x0 = rbind(x0, 0.0), ## set initial value of 't = 0.0'
               cList = cList,
               nList = nList,
               A = A,
               b = b
               )
    return(ans)
}
