`TrenchInverse` <-
function( G )
{
    EPS<-.Machine$double.eps # 1+EPS==1, machine epsilon
    if (!is.matrix(G)) 
        stop("Error: argument G must be symmetric Toeplitz matrix")
    if ((!is.numeric(G))||any(is.na(G)))
        stop("Error: invalid matrix - nonnumerics or na's present")
    a<-G[1,]
    G2<-toeplitz(a)
    if (sum(abs(G2-G))>EPS)
        stop("Error: argument G must be symmetric Toeplitz matrix")
    if (length(a)==1)
        return(1/G)
#order of matrix is >1
    scz<-a[1]
    if (scz < EPS) stop("error in (1,1)-entry of G, variance <= 0")
    a<-a/scz
    out<-.C("trenchR", as.double(a), as.integer(length(a)), EPS,
      tr = array(0, dim=c(length(a),length(a))),  fault = as.integer(1),
      PACKAGE="ltsa" )
    fault<-out$fault
    if (fault == 1) 
            stop("error: matrix singular")
    (out$tr)/scz
}

