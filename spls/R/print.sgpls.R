
# print fit

"print.sgpls" <-
function( x, ... )
{
    xmat <- x$x
    p <- ncol(xmat)
    y.u <- unique(x$y)
    q <- length(y.u)
    eta <- x$eta
    K <- x$K
    
    if ( q == 2 )
    {
        cat( "\nSparse Generalized Partial Least Squares for binary classification\n" )
    }
    if ( q > 2 )
    {
        cat( "\nSparse Generalized Partial Least Squares for multicategory classification\n" )
    }
    cat( "----\n")
    cat( paste("Parameters: eta = ",eta,", K = ",K,"\n",sep="") )
    
    if ( q == 2 )
    {
        A <- x$A
        xAnames <- colnames(xmat)[A]
        
        cat( paste("SGPLS chose ",length(A)," variables among ",p," variables\n\n",sep='') )
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
    if ( q > 2 )
    {
        # total active set
        
        A.all <- x$A.all
        xAnames <- colnames(xmat)[A.all]
        
        cat( paste("SGPLS chose ",length(A.all)," variables among ",p," variables\n\n",sep='') )
        cat( "Selected variables: \n" )
        if ( !is.null(xAnames) )
        {
            for (i in 1:length(A.all))
            {
                cat( paste(xAnames[i],'\t',sep='') )
                if ( i%%5==0 ) { cat('\n') }
            }
        } else
        {
            for (i in 1:length(A.all))
            {
                cat( paste(A.all[i],'\t',sep='') )
                if ( i%%5==0 ) { cat('\n') }
            }
        }
        cat('\n\n')
        
        
        # class-comparison-specific active set
        
        A.bygroup <- x$A.bygroup
        
        for ( i in 1:(q-1) )
        {
            A.bygroup.i <- A.bygroup[[i]]
            xAnames <- colnames(xmat)[ A.bygroup.i ]
            
            cat( paste("For class 0 vs. ",i," comparison,\n",sep="") )
            cat( paste("SGPLS chose ",length(A.bygroup.i)," variables among ",p," variables\n\n",sep='') )
            cat( "Selected variables: \n" )
            if ( !is.null(xAnames) )
            {
                for (i in 1:length(A.bygroup.i))
                {
                    cat( paste(xAnames[i],'\t',sep='') )
                    if ( i%%5==0 ) { cat('\n') }
                }
            } else
            {
                for (i in 1:length(A.bygroup.i))
                {
                    cat( paste(A.bygroup.i[i],'\t',sep='') )
                    if ( i%%5==0 ) { cat('\n') }
                }
            }
            cat('\n\n')
        }
    }
}
