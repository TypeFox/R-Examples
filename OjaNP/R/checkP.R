`checkP` <-
function(p,n,k,silent,type){
    if (!is.null(p)) if(length(p)!=1 | p < 0 | p > 1) stop("'p' must be NULL or between 0 and 1")
    constant <- 6e+06
    if (type=="sign")
        N <- choose(n,k-1)
    else if (type=="rank")
        N <- choose(n,k)
    else if (type=="signedrank")
        N <- 2^k*choose(n,k)
    else
        N <- choose(2*n,k)
    ti <- N*k^3
    if (is.null(p)){    
        if (ti <= constant)   
            p <- 1
        else
            p <- 0.25*constant/ti
    }else{
        if (p > 1) p <- 1
        if (p*N < 1) p <- 1/(N-1)
        if (p*ti > constant){
            longRunTimeMessage(p=p,k=k,N=N,silent=silent)
        }
    } 
    return(p)
}
