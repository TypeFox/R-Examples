make.logtrans <- function()
{

logtrans.fun <- function(times,y,p,more)
{ 
	  y = exp(y)
    x = more$fn(times,y,p,more$more)
    x = x/y

    return(x)
}

logtrans.dfdx <- function(times,y,p,more)
{                  
     x1 = logtrans.fun(times,y,p,more)
     y = exp(y)
     x = more$dfdx(times,y,p,more$more)
     

     for(i in 1:dim(x)[2]){
	     for(j in 1:dim(x)[3]){
		      x[,i,j] = x[,i,j]*y[,j]/y[,i]
        }
        x[,i,i] = x[,i,i] - x1[,i]
     }
     return(x)
}

logtrans.dfdp <- function(times,y,p,more)
{
     y = exp(y)
     x = more$dfdp(times,y,p,more$more)
    
     for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
            x[,i,j] = x[,i,j]/y[,i]
        }
     }
     return(x)
}

logtrans.d2fdx2 <- function(times,y,p,more)
{  
     x1 = logtrans.dfdx(times,y,p,more)
     y = exp(y)
     x = more$d2fdx2(times,y,p,more$more)
      
     for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
           for(k in 1:dim(x)[4]){
              x[,i,j,k] = x[,i,j,k]*y[,j]*y[,k]/y[,i]
           }
        }
     }
     
     for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
           
           x[,i,i,j] = x[,i,i,j] - x1[,i,j]
           
           x[,i,j,i] = x[,i,j,i] - x1[,i,j]
           x[,i,j,j] = x[,i,j,j] + x1[,i,j]
          
        }
     }

     return(x)
}

logtrans.d2fdxdp <- function(times,y,p,more)
{   
     x1 = logtrans.dfdp(times,y,p,more)
     y = exp(y)
     x = more$d2fdxdp(times,y,p,more$more)
     
     
      for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[3]){
           for(k in 1:dim(x)[4]){
		          x[,i,j,k] = x[,i,j,k]*y[,j]/y[,i]
           }
        }
     }
     
     for(i in 1:dim(x)[2]){
        for(j in 1:dim(x)[4]){
           x[,i,i,j] = x[,i,i,j] - x1[,i,j]
        }
     }

     return(x)
}



  more = list(fn = logtrans.fun,
    dfdx = logtrans.dfdx,
    dfdp = logtrans.dfdp,
    d2fdx2 = logtrans.d2fdx2,
    d2fdxdp = logtrans.d2fdxdp,
    extras = NULL)
}
