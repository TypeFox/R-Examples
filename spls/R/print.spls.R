
# print fit

"print.spls" <-
function( x, ... )
{
    xmat <- x$x
    p <- ncol(xmat)
    A <- x$A
    xAnames <- colnames(xmat)[A]
    q <- ncol(x$y)    
    eta <- x$eta
    K <- x$K
    kappa <- x$kappa
    select <- x$select
    fit <- x$fit
    
    if ( q == 1 )
    {
        cat( "\nSparse Partial Least Squares for an univariate response\n" )
        cat( "----\n")
        cat( paste("Parameters: eta = ",eta,", K = ",K,"\n",sep="") )
    }
    if ( q > 1 )
    {
        cat( "\nSparse Partial Least Squares for multivariate responses\n" )
        cat( "----\n")
        cat( paste("Parameters: eta = ",eta,", K = ",K,", kappa = ",kappa,"\n",sep="") )
    }
    cat( paste("PLS algorithm:\n",
        select," for variable selection, ",fit," for model fitting\n\n",sep="") )
    
    cat( paste("SPLS chose ",length(A)," variables among ",p," variables\n\n",sep='') )
    cat( "Selected variables: \n" )
    if ( !is.null(xAnames) )
    {
        for (i in 1:length(A))
        {
            cat( paste(xAnames[i],'\t',sep='') )
            if ( i%%5==0 ) { cat('\n') }
        }
    } else
    {
        for (i in 1:length(A))
        {
            cat( paste(A[i],'\t',sep='') )
            if ( i%%5==0 ) { cat('\n') }
        }
    }
    cat('\n')
}
