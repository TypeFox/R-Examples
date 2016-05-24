make.Henon <- function()
{

Henon.ode <- function(times,y,p,more)
{
    r = y
    r[1] = 1- p[1]*y[1]^2 + y[2]
    r[2] = p[2]*y[1]
    
    return(r)
}

Henon.fun <- function(times,y,p,more)
{
    r = y
    r[,1] = 1- p[1]*y[,1]^2 + y[,2]
    r[,2] = p[2]*y[,1]
    
    return(r)
}

Henon.dfdx <- function(times,y,p,more)
{
    r = array(0,c(length(times),ncol(y),ncol(y)))
      
    r[,1,1] = -2*p[1]*y[,1]
    r[,1,2] = 1
    r[,2,1] = p[2]
    
    return(r)
}

Henon.dfdp <- function(times,y,p,more)
{
    r = array(0,c(length(times),ncol(y),length(p)))
    
    r[,1,1] = -y[,1]^2
    r[,2,2] = y[,1]
    
    return(r)
}


Henon.d2fdx2 <- function(times,y,p,more)
{
    r = array(0,c(length(times),rep(ncol(y),3)))
    
    r[,1,1,1] = -2*p[1]
    
    return(r)
}


Henon.d2fdxdp <- function(times,y,p,more)
{
    r = array(0,c(length(times),ncol(y),ncol(y),length(p)))
    
    r[,1,1,1] = -2*y[,1]
    r[,2,1,2] = 1
    
    return(r)
}



  return( 
    list(
        fn = Henon.fun,
        ode = Henon.ode,
        dfdx = Henon.dfdx,
        dfdp = Henon.dfdp,
        d2fdx2 = Henon.d2fdx2,
        d2fdxdp = Henon.d2fdxdp
    )
  )
}
