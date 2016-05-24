abouheif.moran <- function (x, W=NULL,
                            method=c("oriAbouheif","patristic","nNodes","Abouheif","sumDD"),
                            a=1, nrepet=999,alter=c("greater", "less", "two-sided")) {

    ## some checks
    ## if(!require(ade4)) stop("The ade4 package is not installed.")
    alter <- match.arg(alter)
    method <- match.arg(method)

    ## handle W
    if(!is.null(W)){ # W is provided
        if (any(W<0)) stop ("negative terms found in 'W'")
        if (nrow(W) != ncol(W)) stop ("'W' is not squared")
        W <- as.matrix(W)
    } else { # otherwise computed W from x, a phylo4d object
        if(!inherits(x, "phylo4d")) stop("if W is not provided, x has to be a phylo4d object")
        if (is.character(chk <- checkPhylo4(x))) stop("bad phylo4d object: ",chk)
        ##if (is.character(chk <- checkData(x))) stop("bad phylo4d object: ",chk) no longer needed
        W <- proxTips(x, method=method, a=a, normalize="row", symmetric=TRUE)
    }

    nobs <- ncol(W)
    ## W has to be symmetric
    W <- (W + t(W))/2

    ## take data from x if it is a phylo4d
    if(inherits(x, "phylo4d")){
        if (is.character(chk <- checkPhylo4(x))) stop("bad phylo4d object: ",chk)
        ## if (is.character(chk <- checkData(x))) stop("bad phylo4d object: ",chk) : no longer needed
        x <- tdata(x, type="tip")
    }

    ## main computations
    x <- data.frame(x)
    test.names <- names(x)
    x <- data.matrix(x) # convert all variables to numeric type

    if (nrow(x) != nobs) stop ("non convenient dimension")
    nvar <- ncol(x)
    res <- .C("gearymoran",
        param = as.integer(c(nobs,nvar,nrepet)),
        data = as.double(x),
        W = as.double(W),
        obs = double(nvar),
        result = double (nrepet*nvar),
        obstot = double(1),
        restot = double (nrepet),
        PACKAGE="adephylo"
    )
    res <- as.krandtest(obs=res$obs,sim=matrix(res$result,ncol=nvar, byrow=TRUE),
                        names=test.names,alter=alter)
    return(res)
} # end abouheif.moran
