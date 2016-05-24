##
## Function for creating an object of reference-class 'DNL'
dnl <- function(q, A = NULL, b = NULL, cList = list(),
                x0, nlfList = list(), nlgList = list(), nlhList = list()){
    if(is.matrix(q)){
        warning("Matrix provided for q, extracting first column for argument 'q'.\n")
        q <- q[, 1]
    }
    n <- length(q)
    if(is.matrix(x0)){
        warning("Matrix provided for x0, extracting first column for argument 'x0'.\n")
        x0 <- x0[, 1]
    }
    if(!identical(length(x0), length(q))){
        stop("Length of initial point 'x' is not equal to the dimension of the objective.\n")
    }
    if(is.null(A)){
        A <- matrix(0, nrow = 0, ncol = n)
    } 
    if(is.null(dim(A))){
        A <- matrix(A, nrow = 1)
    }
    if(is.null(b)){
        b <- numeric()
    }
    ##
    ## Checking provided non-linear constraints
    ##
    mnl <- length(nlfList)
    if(mnl < 1){
        warning("Empty list for non-linear convex constraints provided.\nReturning as DLP object.\n")
        cpd <- dlp(q = q, A = A, b = b, cList = cList)
        return(cpd)
    }
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
    nList <- list(nlfList, nlgList, nlhList)
    ## Creating objects related to NLFC
    Gnl <- matrix(0, nrow = mnl, ncol = n)
    hnl <- matrix(0, nrow = mnl, ncol = 1)
    K <- length(cList)
    if(K > 0){
        cone <- unlist(lapply(cList, function(x) x[["conType"]]))
        if(!all(cone %in% c("NNOC", "SOCC", "PSDC"))){
            stop("List elements of cone constraints must be either created by calls to:\n'nnoc()', or 'socc()', or 'psdc()'.\n")
        }
        cone <- c("NLFC", cone)
        GList <- lapply(cList, function(x) x[["G"]])
        GList <- c(list(Gnl), GList)
        G <- do.call("rbind", GList)
        hList <- lapply(cList, function(x) x[["h"]])
        hList <- c(list(hnl), hList)
        h <- do.call("rbind", hList)
        ridx <- cumsum(unlist(lapply(GList, nrow)))
        sidx <- cbind(c(0, ridx[-length(ridx)]), ridx - 1)
        dims <- c(mnl, as.integer(unlist(lapply(cList, function(x) x[["dims"]]))))
        K <- K + 1L
        cList <- new(CONEC, cone, G, h, sidx, dims, K, n)
    } else {
        cList <- new(CONEC, "NLFC", Gnl, hnl, mnl, 1L, n)
    }
    ans <- new(DNL,
               q = q,
               A = A,
               b = b,
               cList = cList,
               x0 = as.matrix(x0),
               nList = nList)
    return(ans)
}
