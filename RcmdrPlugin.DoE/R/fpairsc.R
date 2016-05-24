fpairsc <- function(factors, G1, class=1){
    nfactors <- length(factors)
    ## create pairs for contour plots
    if (!class %in% c(1,2,3,4)) stop("class must be an integer from 1 to 4")
    if (!is.numeric(G1)) stop("G1 must be numeric")
    if (length(G1)==0) stop("At least one effect must be in group 1")
    if (length(G1)==1 & class==1) stop("no plots required!")
    if (!is.character(factors)) stop("factors must be character")
    if (!all(G1%%1==0)) stop("G1 must consist of integer values")
    if (!length(unique(G1))==length(G1)) stop("non-unique values in G1")
    if (!nfactors>=max(G1)) stop("G1 must not contain numbers larger than nfactors")
    if (!min(G1)>=1) stop("G1 must not contain numbers smaller than 1")
    G2 <- setdiff(1:nfactors,G1)
    
    ## maximum number of factors for resolution V
    if (length(G1)>1)
    requirement <- apply(matrix(factors[G1[combn(length(G1),2)]],nrow=2),2,"paste",collapse="*")
    else requirement <- character(0)
    if (class==3)
    requirement <- c(requirement, outer(factors[G1],factors[G2],FUN=function(X,Y) paste(pmin(X,Y),pmax(X,Y),sep="*")))
    if (class==2 & length(G2)>=2){
    requirement <- c(requirement, apply(matrix(factors[G2[combn(length(G2),2)]],nrow=2),
                            2,"paste",collapse="*"))
                            }
    if (class==4)
    requirement <- c(outer(factors[G1],factors[G2],FUN=function(X,Y) paste(pmin(X,Y),pmax(X,Y),sep="*")))
    requirement
}
