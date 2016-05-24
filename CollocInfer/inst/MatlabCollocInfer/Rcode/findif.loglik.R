
make.findif.loglik = function()
{
findif.loglik.fun <- function(data,times,y,p,more)
{
    x = more$fn(data,times,y,p,more$more)
}

findif.loglik.dfdx <- function(data,times,y,p,more)  # deriv wrt state
{
    x1 = more$fn(data,times,y,p,more$more)
    x = array(0,dim(y))

    for(i in 1:ncol(y)){
        ty = y
        ty[,i] = y[,i] + more$eps
        x[,i] = (more$fn(data,times,ty,p,more$more)-x1)/more$eps
    }
    return(x)
}


findif.loglik.dfdy <- function(data,times,y,p,more)  # deriv wrt response
{
    x1 = more$fn(data,times,y,p,more$more)
    x = array(0,dim(data))

    for(i in 1:ncol(data)){
        tdata = data
        tdata[,i] = data[,i] + more$eps
        x[,i] = (more$fn(tdata,times,y,p,more$more)-x1)/more$eps
    }
    return(x)
}


findif.loglik.dfdp <- function(data,times,y,p,more)
{

    x1 = more$fn(data,times,y,p,more$more)
    x = array(0,c(length(x1),length(p)))

    for(i in 1:length(p)){
        tp = p
        tp[i] = p[i] + more$eps
        x[,i] = (more$fn(data,times,y,tp,more$more)-x1)/more$eps        
    }
    return(x)
}

findif.loglik.d2fdx2 <- function(data,times,y,p,more)
{
    x1 = findif.loglik.dfdx(data,times,y,p,more)
    x = array(0,c(dim(x1),ncol(y)))

    for(i in 1:ncol(y)){
        ty = y
        ty[,i] = y[,i] + more$eps
        x[,,i] = (findif.loglik.dfdx(data,times,ty,p,more)-x1)/more$eps
    }
    return(x)
}

findif.loglik.d2fdy2 <- function(data,times,y,p,more)
{
    x1 = findif.loglik.dfdy(data,times,y,p,more)
    x = array(0,c(dim(x1),ncol(data)))

    for(i in 1:ncol(data)){
        tdata = data
        tdata[,i] = data[,i] + more$eps
        x[,,i] = (findif.loglik.dfdy(tdata,times,y,p,more)-x1)/more$eps
    }
    return(x)
}

findif.loglik.d2fdxdy <- function(data,times,y,p,more)
{
    x1 = findif.loglik.dfdx(data,times,y,p,more)
    x = array(0,c(dim(x1),ncol(data)))

    for(i in 1:ncol(data)){
        tdata = data
        tdata[,i] = data[,i] + more$eps
        x[,,i] = (findif.loglik.dfdx(tdata,times,y,p,more)-x1)/more$eps
    }
    return(x)
}


findif.loglik.d2fdxdp <- function(data,times,y,p,more)
{
    x1 = findif.loglik.dfdx(data,times,y,p,more)
    x = array(0,c(dim(x1),length(p)))
    
    for(i in 1:length(p)){
        tp = p
        tp[i] = p[i] + more$eps
        x[,,i] = (findif.loglik.dfdx(data,times,y,tp,more)-x1)/more$eps
    }

    return(x)
}


findif.loglik.d2fdydp <- function(data,times,y,p,more)
{
    x1 = findif.loglik.dfdy(data,times,y,p,more)
    x = array(0,c(dim(x1),length(p)))
    
    for(i in 1:length(p)){
        tp = p
        tp[i] = p[i] + more$eps
        x[,,i] = (findif.loglik.dfdy(data,times,y,tp,more)-x1)/more$eps
    }

    return(x)
}

#make.findif.loglik <- function()
#{
    return(
        list(
            fn = findif.loglik.fun,
            dfdx = findif.loglik.dfdx,
            dfdy = findif.loglik.dfdy,
            dfdp = findif.loglik.dfdp,
            d2fdx2 = findif.loglik.d2fdx2,
            d2fdy2 = findif.loglik.d2fdy2,
            d2fdxdy = findif.loglik.d2fdxdy,
            d2fdxdp = findif.loglik.d2fdxdp,
            d2fdydp = findif.loglik.d2fdydp
        )
    )
}
