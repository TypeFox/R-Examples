###################
# orthobasis.phylo
###################
orthobasis.phylo <- function(x=NULL, prox=NULL,
                             method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"), a=1){
    ## if(!require(phylobase)) stop("phylobase package is not installed")
    ## if(!require(ade4)) stop("ade4 package is not installed")

    ## handle arguments
    method <- match.arg(method)

    if(is.null(prox)){ # have to compute prox
        x <- as(x, "phylo4")
        if (is.character(checkval <- checkPhylo4(x))) stop(checkval)
        W <- proxTips(x, tips="all", method=method, a=a, normalize="row", symmetric=TRUE)
    } else { # prox is provided
        W <- as.matrix(prox)
        if(!is.matrix(W)) stop("W is not a matrix")
        if(ncol(W) != nrow(W)) stop("W is not a square matrix")
         diag(W) <- 0
        W <- 0.5 * (t(W) + W) # re-symmetrization
    }

    n <- nrow(W)


    ## main computation -> call to orthobasis.mat
    res <- orthobasis.mat(W, cnw=FALSE)

    ## build output
    row.names(res) <- rownames(W)
    names(res) <- paste("ME", 1:ncol(res))
    names(attr(res,"values")) <- names(res)
    attr(res,"call") <- match.call()
    attr(res,"class") <- c("orthobasis","data.frame")

    return(res)
} # end orthobasis.phylo





###########
# me.phylo
###########
me.phylo <- orthobasis.phylo
