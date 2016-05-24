make.fhn.input <- function()
{

fhn.fun <- function(times,y,p,more)
{
  I = more$infn(times)
  r = y;
  r[,'V'] = (y[,'V'] - y[,'V']^3/3 - y[,'R'])   + I 
  r[,'R'] = (y[,'V'] + p['a'] - p['b']*y[,'R'])*p['c'] 
  return(r)
}


fhn.dfdx <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2))
  dimnames(r) = list(NULL,colnames(y),colnames(y))

  r[,'V','V'] = 1 - y[,'V']^2
  r[,'V','R'] =  -1
  r[,'R','V'] = p['c']
  r[,'R','R'] = -p['b']*p['c']
  return(r)
}


fhn.dfdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),length(p)))
  dimnames(r) = list(NULL,colnames(y),names(p))

  r[,'R','a'] = p['c'];
  r[,'R','b'] = -y[,'R']*p['c'];
  r[,'R','c'] = y[,'V']+p['a']-p['b']*y[,'R'];
  return(r)
}


fhn.d2fdx2 <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2,2))
  dimnames(r) = list(NULL,colnames(y),colnames(y),colnames(y))

  r[,'V','V','V'] = -2*y[,'V']
  return(r)
}


fhn.d2fdxdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2,length(p)))
  dimnames(r) = list(NULL,colnames(y),colnames(y),names(p))

  r[,'R','V','c'] = 1
  r[,'R','R','b'] = -p['c']
  r[,'R','R','c'] = -p['b']

  return(r)
}



  return( 
    list(
        fn = fhn.fun,
        dfdx = fhn.dfdx,
        dfdp = fhn.dfdp,
        d2fdx2 = fhn.d2fdx2,
        d2fdxdp = fhn.d2fdxdp
    )
  )
}

