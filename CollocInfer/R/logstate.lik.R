### This transforms a likelihood for data given a state x given in the
### functions in "more" into a likelihood given exp(x) and associated
### derivatives

make.logstate.lik <- function()
{

logstate.lik.fun <- function(data,times,y,p,more)
{
    y = exp(y)
    x = more$fn(data,times,y,p,more$more)
    
    return(x)
}

logstate.lik.dfdx <- function(data,times,y,p,more)
{
     y = exp(y)
     x = more$dfdx(data,times,y,p,more$more)

     for(i in 1:dim(x)[2]){
        x[,i] = x[,i]*y[,i]
     }
     return(x)
}

logstate.lik.dfdp <- function(data,times,y,p,more)
{
     y = exp(y)
     x = more$dfdp(data,times,y,p,more$more)
     return(x)
}

logstate.lik.d2fdx2 <- function(data,times,y,p,more)
{
     x1 = logstate.lik.dfdx(data,times,y,p,more)
     
     y = exp(y)
     x = more$d2fdx2(data,times,y,p,more$more)

     for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
              x[,i,j] = x[,i,j]*y[,i]*y[,j]
        }
        x[,i,i] = x[,i,i] + x1[,i]
     }
     return(x)
}

logstate.lik.d2fdxdp <- function(data,times,y,p,more)
{
     y = exp(y)
     x = more$d2fdxdp(data,times,y,p,more$more)

     for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
           x[,i,j] = x[,i,j]*y[,i]
        }
     }
     return(x)
}



    return( 
        list(
            fn = logstate.lik.fun,
            dfdx = logstate.lik.dfdx,
            dfdp = logstate.lik.dfdp,
            d2fdx2 = logstate.lik.d2fdx2,
            d2fdxdp = logstate.lik.d2fdxdp
        )
    )
}
