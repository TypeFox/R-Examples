### CONSTANT VARIANCE FUNCTIONS

    


make.cvar <- function()
{

cvar<-list()
cvar$var.fn<- function(t,y,p,more)  
{
    n = ncol(y)
    
    more = checkmore.cvar(more,n)
    pmat = 0*more$mat
    
    pmat[more$sub[,1:2]] = pmat[more$sub[,1:2]] + p[more$sub[,3]]
    
    if(nrow(pmat) > 1){
        pmat = pmat - diag(diag(pmat))
        pmat[more$sub[,c(2,1)]] = pmat[more$sub[,c(2,1)]] + p[more$sub[,3]]    
    }
    
    pmat = pmat + more$mat    
    
    return( array(pmat,c(1,dim(pmat))) )
}



cvar$var.dfdx <- function(t,y,p,more)
{
   if( !is.null(more$mat) ){ ny = nrow(more$mat) }
   else{ ny = ncol(y) }
   nx = ncol(y)
   
    return( array(0,c(1,ny,ny,nx)) )
}


cvar$var.dfdp <- function(t,y,p,more)
{
    n = ncol(y)

    more = checkmore.cvar(more,n)

    r = array(0,c(1,dim(more$mat),length(p)))
    
    ind = cbind(rep(1,2*dim(more$sub)[1]), rbind( more$sub, more$sub[,c(2,1,3)]) )

    r[ind] = 1
    
    return(r)
}


cvar$var.d2fdxdp <- function(t,y,p,more)
{
   if( !is.null(more$mat) ){ ny = nrow(more$mat) }
   else{ ny = ncol(y) }  
   nx = ncol(y)
   
    return( array(0,c(1,rep(ny,2),nx,length(p))) )
}

cvar$var.d2fdx2 <- function(t,y,p,more)
{
    if( !is.null(more$mat) ){ ny = nrow(more$mat) }
    else{ ny = ncol(y) }
    nx = ncol(y)
    
    return( array(0,c(1,ny,ny,nx,nx)) )
}


checkmore.cvar <- function(more,n)    # checks additional arguments to cvar
{
    if(is.null(more)){ more$sub = matrix(1:n,n,3) }    # default to diagonal covariance
    else{ 
        if(is.null(more$sub)){ more$sub = matrix(0,0,3) }
    }
    if(is.null(more$mat)){ 
        if(dim(more$sub)[1] == 0){ more$mat = diag(rep(1,n)) }
        else{ more$mat = matrix(0,n,n) }
    }
    
    return(more)
}

    
return(cvar) 
} 
