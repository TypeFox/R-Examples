Circle.fun <- function(times,y,p,more)
{
  r = y;
  r[,1] = p[1]*y[2];
  r[,2] = -p[2]*y[1];
  return(r)
}



Circle.dfdx <- function(times,y,p,more)
{
  r = array(0,c(dim(y),ncol(y)))
  r[,1,2] = p[1];
  r[,2,1] = -p[2]; 
  return(r)
}


Circle.dfdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),length(p)))
  r[,1,1] = y[2];
  r[,2,2] = -y[1];
  return(r)
}


Circle.d2fdx2 <- function(times,y,p,more)
{
  r = array(0,c(dim(y),ncol(y),ncol(y)))
  return(r)
}


Circle.d2fdxdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2,length(p)))
  r[,1,2,1] = 1;
  r[,2,1,2] = -1;
  return(r)
}


make.Circle <- function()
{
  return( 
    list(
        fn = Circle.fun,
        dfdx = Circle.dfdx,
        dfdp = Circle.dfdp,
        d2fdx2 = Circle.d2fdx2,
        d2fdxdp = Circle.d2fdxdp
    )
  )
}

