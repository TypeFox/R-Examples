######################################################
# Finite differencing for process variance
######################################################

make.findif.var = function()
{

findif.var.fun <- function(times,y,p,more)
{
    x = more$var.fn(times,y,p,more$more)
}

findif.var.dfdx <- function(times,y,p,more)
{
    x1 = more$var.fn(times,y,p,more$more)
    x = array(0,c(dim(x1),ncol(y)))

    for(i in 1:ncol(y)){
        ty = y
        ty[,i] = y[,i] + more$eps
        x[,,,i] = (more$var.fn(times,ty,p,more$more)-x1)/more$eps
    }
    return(x)
}

findif.var.dfdp <- function(times,y,p,more)
{
    x1 = more$var.fn(times,y,p,more$more)
    x = array(0,c(dim(x1),length(p)))

    for(i in 1:length(p)){
        tp = p
        tp[i] = p[i] + more$eps
        x[,,,i] = (more$var.fn(times,y,tp,more$more)-x1)/more$eps        
    }
    return(x)
}

findif.var.d2fdx2 <- function(times,y,p,more)
{
    x1 = findif.var.dfdx(times,y,p,more)
    x = array(0,c(dim(x1),ncol(y)))
    
    for(i in 1:ncol(y)){
        ty = y
        ty[,i] = y[,i] + more$eps
        x[,,,,i] = (findif.var.dfdx(times,ty,p,more)-x1)/more$eps
    }
    return(x)
}

findif.var.d2fdxdp <- function(times,y,p,more)
{
    x1 = findif.var.dfdx(times,y,p,more)
    x = array(0,c(dim(x1),length(p)))
    
    for(i in 1:length(p)){
        tp = p
        tp[i] = p[i] + more$eps
        x[,,,,i] = (findif.var.dfdx(times,y,tp,more)-x1)/more$eps
    }

    return(x)

}


#make.findif.var <- function()
#{
    return(
        list(
            fn = findif.var.fun,
            dfdx = findif.var.dfdx,
            dfdp = findif.var.dfdp,
            d2fdx2 = findif.var.d2fdx2,
            d2fdxdp = findif.var.d2fdxdp
        )
    )
}
