make.exptrans <- function()
{

exptrans.fun <- function(times,y,p,more)
{ 

    y = exp(y)
    x = more$fn(times,y,p,more$more)

    return(x)
}

exptrans.dfdx <- function(times,y,p,more)
{                  
     y  = exp(y)
     x  = more$dfdx(times,y,p,more$more)

     for(i in 1:dim(x)[2]){
	     for(j in 1:dim(x)[3]){
		      x[,i,j] = x[,i,j]*y[,j]
        }
     }
     return(x)
}

exptrans.dfdp <- function(times,y,p,more)
{
     y = exp(y)
     x = more$dfdp(times,y,p,more$more)

     return(x)
}

exptrans.d2fdx2 <- function(times,y,p,more)
{  
     x1 = exptrans.dfdx(times,y,p,more)
     y = exp(y)
     x = more$d2fdx2(times,y,p,more$more)
      
     for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
           for(k in 1:dim(x)[4]){
              x[,i,j,k] = x[,i,j,k]*y[,j]*y[,k]
           }
           x[,i,j,j] = x[,i,j,j] + x1[,i,j]
        }
     }
     
     return(x)
}

exptrans.d2fdxdp <- function(times,y,p,more)
{   
     y = exp(y)
     x = more$d2fdxdp(times,y,p,more$more)

      for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
           for(k in 1:dim(x)[4]){
		          x[,i,j,k] = x[,i,j,k]*y[,j]
           }
        }
     }
     return(x)
}



  more = list(fn = exptrans.fun,
    dfdx = exptrans.dfdx,
    dfdp = exptrans.dfdp,
    d2fdx2 = exptrans.d2fdx2,
    d2fdxdp = exptrans.d2fdxdp,
    extras = NULL)
}
