######################################################
# Finite differencing for ODE RHS
######################################################


make.findif.ode <- function()
{

findif.ode<-list()
findif.ode$fn <- function(times,y,p,more)
{
    x = more$fn(times,y,p,more$more)
}

findif.ode$dfdx <- function(times,y,p,more)
{
    x1 = more$fn(times,y,p,more$more)
    x = array(0,c(dim(x1),ncol(y)))

    for(i in 1:ncol(y)){
        ty = y
        ty[,i] = y[,i] + more$eps
        x[,,i] = (more$fn(times,ty,p,more$more)-x1)/more$eps
    }
    return(x)
}

findif.ode$dfdp <- function(times,y,p,more)
{
    x1 = more$fn(times,y,p,more$more)
    x = array(0,c(dim(x1),length(p)))

    for(i in 1:length(p)){
        tp = p
        tp[i] = p[i] + more$eps
        x[,,i] = (more$fn(times,y,tp,more$more)-x1)/more$eps        
    }
    return(x)
}

findif.ode$d2fdx2 <- function(times,y,p,more)
{
    x1 = findif.ode$dfdx(times,y,p,more)
    x = array(0,c(dim(x1),ncol(y)))
    
    for(i in 1:ncol(y)){
        ty = y
        ty[,i] = y[,i] + more$eps
        x[,,,i] = (findif.ode$dfdx(times,ty,p,more)-x1)/more$eps
    }
    return(x)
}

findif.ode$d2fdxdp <- function(times,y,p,more)
{
    x1 = findif.ode$dfdx(times,y,p,more)
    x = array(0,c(dim(x1),length(p)))
    
    for(i in 1:length(p)){
        tp = p
        tp[i] = p[i] + more$eps
        x[,,,i] = (findif.ode$dfdx(times,y,tp,more)-x1)/more$eps
    }

    return(x)

}


    return(findif.ode)
}
