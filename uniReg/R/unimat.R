unimat <- function(p,m){
#######################################################################################
# Generates the matrix A of unimodality-constraints with mode m (for p coefficients). #
#######################################################################################
    if(!(p%%1==0 && p>=1 && is.finite(p))){stop("p should be a finite whole number >=1.")}
    if(!(m%%1==0 && m>=1 && p>=m)){stop("m should be a whole number between 1 and p.")}
    
    Amon <- function(q){return(cbind(diag(-1,q-1),rep(0,q-1)) + cbind(rep(0,q-1),diag(1,q-1)))}
    return(rbind(cbind(Amon(m),matrix(0,nrow=m-1,ncol=p-m)),
                 cbind(matrix(0,nrow=p-m,ncol=m-1),-Amon(p-m+1))))
}
