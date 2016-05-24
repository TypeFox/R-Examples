## ###################################################################
##
## Maximum Cardinality Search (mcs)
## Maximum Cardinality Search for marked graphs (mcsmarked)
##
## Returns perfect ordering if it exists and character(0) otherwise
##
## ###################################################################

mcs <- function(object, root=NULL, index=FALSE){
  UseMethod("mcs")
}

## FIXME: mcs: returns character(0) if graph is not undirected. Should
## FIXME: mcs: instead signal an error??

mcs.default <- function(object, root=NULL, index=FALSE){
    cls <- match.arg(class( object ),
                     c("graphNEL","matrix","dgCMatrix","igraph"))
    mm <- coerceGraph( object, "matrix" )
    if ( !is.UGMAT(mm) )
        character(0) ##FIXME: mcs.default: Should perhaps be error...
    else
        mcsMAT( mm, root=root, index=index )
}


mcsMAT <- function (amat, vn = colnames(amat), root = NULL, index = FALSE)
{
    vn.old <- vn
    if (!is.null(root)){
        vn    <- c(root, setdiffPrim(vn, root))
        root2 <- match(vn, vn.old)-1
    } else {
        root2 <- 0:(ncol(amat)-1)
    }
    ##cat("mcsMAT:"); print(root2)

    a <- mcsMAT_( amat, root2 )

    if (index){
        if (a[1]<0){
            NA
        } else {
            a+1
        }
    } else {
        if (a[1]<0){
            character(0)
        } else {
            vn.old[a+1]
        }
    }
}



mcsmarked <- function (object, discrete=NULL, index = FALSE){
  UseMethod("mcsmarked")
}

mcsmarked.default <- function (object, discrete=NULL, index = FALSE){
    cls <- match.arg(class( object ),
                     c("graphNEL","igraph","matrix","dgCMatrix"))
    switch(cls,
           "graphNEL" ={
               if (is.null(discrete))
                   mcsMAT(graphNEL2dgCMatrix(object), index=index)
               else
                   mcsmarkedMAT(graphNEL2dgCMatrix(object), discrete=discrete, index = index)
           },
           "igraph"   ={
               if (is.null(discrete))
                   mcsMAT(igraph::get.adjacency(object), index=index)
               else
                   mcsmarkedMAT(igraph::get.adjacency(object), discrete=discrete, index = index)
           },
           "dgCMatrix"=,
           "matrix"   ={
               if (is.null(discrete))
                   mcsMAT(object, index=index)
               else
                   mcsmarkedMAT(object, discrete=discrete, index = index)
           })
}



## FIXME: mcsmarkedMAT: candidate for C++ implementation.
mcsmarkedMAT <- function(amat, vn = colnames(amat), discrete = NULL, index = FALSE) {

    nv   <- length(vn)

    if (is.null(discrete)){
        return(mcsMAT(amat, vn=vn, index=index))
    } else {
        if (is.logical(discrete)){
            discrete <- as.numeric(discrete)
        }
        if (is.numeric(discrete)){
            if ( length(discrete) != nv ){
                stop("'discrete' is numeric or logical but does not have the correct length")
            } else {
                vn.ext <- c(".", vn)
                idx <- c(1, discrete)
            }
        } else {
            if (is.character(discrete)){
                vn.ext <- c(".", vn)
                idx <- match(discrete, vn.ext)
            } else {
                stop ("'discrete' is not a character")
            }
        }
    }

    amat.ext <- as(Matrix(0, nrow=nv+1L, ncol=nv+1L), "dgCMatrix")
    amat.ext[2:(nv+1),2:(nv+1)] <- amat
    amat.ext[idx, 1L] <- 1L
    amat.ext[1L, idx] <- 1L
    ans <- mcsMAT(amat.ext, vn=vn.ext, index=index)
    if (length(ans)>0)
        ans <- ans[-1L]
    ans
}


