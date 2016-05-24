
# correct wrong specification of parameters

"correctp" <-
function( x, y, eta, K, kappa, select, fit )
{    
    # eta
    
    if ( min(eta)<0 | max(eta)>=1 )
    {
        if ( max(eta)==1 )
        {
            stop('eta should be strictly less than 1!')
        }    
        if ( length(eta)==1 )
        {
            stop('eta should be between 0 and 1!')
        } else
        {
            stop('eta should be between 0 and 1! \n  Choose appropriate range of eta!')
        }
    }
    
    # K
    
    if ( max(K) > ncol(x) )
    {
        stop('K cannot exceed the number of predictors! Pick up smaller K!')
    }
    if ( max(K) >= nrow(x) )
    {
        stop('K cannot exceed the sample size! Pick up smaller K!')
    }
    if ( min(K)<=0 | !all(K%%1==0) )
    {
        if ( length(K)==1 )
        {
            stop('K should be a positive integer!')
        } else
        {
            stop('K should be a positive integer! \n  Choose appropriate range of K!')
        }    
    }
    
    # kappa
    
    if ( kappa>0.5 | kappa<0 )
    {
        cat('kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n')
        kappa <- 0.5
    }   
    
    # select
    
    if ( select!="pls2" & select!="simpls" )
    {
        cat('Invalid PLS algorithm for variable selection.\n')
        cat('pls2 algorithm is used. \n\n')
        select <- 'pls2'
    }
    
    # fit
    
    fits <- c("simpls", "kernelpls", "widekernelpls", "oscorespls")
    if ( !any(fit==fits) )
    {
        cat('Invalid PLS algorithm for model fitting\n')
        cat('simpls algorithm is used. \n\n')
        fit <- 'simpls'
    }
    
    list( K=K, eta=eta, kappa=kappa, select=select, fit=fit )
}
