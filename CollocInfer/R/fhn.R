make.fhn <- function()
{

fhn.fun <- function(times,y,p,more)
{
  r = y;
  r[,'V'] = p['c']*(y[,'V'] - y[,'V']^3/3 + y[,'R'])
  r[,'R'] = -(y[,'V'] -p['a'] + p['b']*y[,'R'])/p['c'] 
  return(r)
}


fhn.fun.ode <- function(times,y,p)
{
  r = y;
  dimnames(r) = dimnames(y);
  r['V'] = p['c']*(y['V'] - y['V']^3/3 + y['R'])
  r['R'] = -(y['V'] -p['a'] + p['b']*y['R'])/p['c'] 

  return(list(r))
}


fhn.dfdx <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2))
  dimnames(r) = list(NULL,colnames(y),colnames(y))

  r[,'V','V'] = p['c'] - p['c']*y[,'V']^2
  r[,'V','R'] =  p['c']
  r[,'R','V'] = (-1/p['c'])
  r[,'R','R'] = (-p['b']/p['c'])
  return(r)
}


fhn.dfdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),length(p)))
  dimnames(r) = list(NULL,colnames(y),names(p))

  r[,'V','c'] =  (y[,'V']-y[,'V']^3/3+y[,'R']);
  r[,'R','a'] = 1/p['c'];
  r[,'R','b'] = (-y[,'R']/p['c']);
  r[,'R','c'] = ((y[,'V']-p['a']+p['b']*y[,'R'])/(p['c']^2));
  return(r)
}


fhn.d2fdx2 <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2,2))
  dimnames(r) = list(NULL,colnames(y),colnames(y),colnames(y))

  r[,'V','V','V'] = -2*p['c']*y[,'V']
  return(r)
}


fhn.d2fdxdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2,length(p)))
  dimnames(r) = list(NULL,colnames(y),colnames(y),names(p))

  r[,'V','V','c'] = 1 - y[,'V']^2
  r[,'V','R','c'] = 1
  r[,'R','V','c'] = 1/p['c']^2
  r[,'R','R','b'] = -1/p['c']
  r[,'R','R','c'] = p['b']/p['c']^2

  return(r)
}



  return( 
    list(
        fn = fhn.fun,
        fn.ode = fhn.fun.ode,
        dfdx = fhn.dfdx,
        dfdp = fhn.dfdp,
        d2fdx2 = fhn.d2fdx2,
        d2fdxdp = fhn.d2fdxdp
    )
  )
}

