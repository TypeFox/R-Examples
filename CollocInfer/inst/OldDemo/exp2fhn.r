exp2fhn.fun <- function(times,y,p,more)
{
  r = y;
  r[,'V'] = y[,'V']* p['c']*(log(y[,'V']) - (log(y[,'V']))^3/3 + log(y[,'R']))
  r[,'R'] = -y[,'R']*(log(y[,'V']) -p['a'] + p['b']*log(y[,'R']))/p['c'] 
  return(r)
}


exp2fhn.fun.ode <- function(times,y,p)
{
  r = y;
  dimnames(r) = dimnames(y);
  r['V'] = y['V']* p['c']*(log(y['V']) - (log(y['V']))^3/3 + log(y['R']))
  r['R'] = -y['R']*(log(y['V']) -p['a'] + p['b']*log(y['R']))/p['c'] 
  return(list(r))
}


exp2fhn.dfdx <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2))
  dimnames(r) = list(NULL,colnames(y),colnames(y))

  r[,'V','V'] = p['c']*(log(y[,'V']) - (log(y[,'V']))^3/3 + log(y[,'R']))+p['c']*(1-(log(y[,'V']))^2)
  r[,'V','R'] =  p['c']*y[,'V']/y[,'R']
  r[,'R','V'] = -y[,'R']/(y[,'V']*p['c'])
  r[,'R','R'] = -p['b']/p['c']*log(exp(1)*y[,'R'])-1/p['c']*log(y[,'V']*exp(-p['a']))
  return(r)
}


exp2fhn.dfdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),length(p)))
  dimnames(r) = list(NULL,colnames(y),names(p))

  r[,'V','c'] =  (log(y[,'V'])-(log(y[,'V']))^3/3+log(y[,'R']))*y[,'V'];
  r[,'R','a'] = y[,'R']/p['c'];
  r[,'R','b'] = (-y[,'R']*log(y[,'R']))/p['c'];
  r[,'R','c'] = (log(y[,'V'])-p['a']+p['b']*log(y[,'R']))*y[,'R']/p['c']^2;
  return(r)
}


exp2fhn.d2fdx2 <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2,2))
  dimnames(r) = list(NULL,colnames(y),colnames(y),colnames(y))

  r[,'V','V','V'] = p['c']/y[,'V']*(1-(log(y[,'V']))^2-2*log(y[,'V']))
  r[,'V','V','R'] = p['c']/y[,'R']
  r[,'V','R','V'] = p['c']/y[,'R']
  r[,'V','R','R'] = - p['c']*y[,'V']/y[,'R']^2
  r[,'R','V','V'] = y[,'R']/(y[,'V']^2*p['c'])
  r[,'R','V','R'] = -1/(p['c']*y[,'V'])
  r[,'R','R','V'] = -1/(p['c']*y[,'V'])
  r[,'R','R','R'] = -p['b']/(p['c']*y[,'R'])
  return(r)
}


exp2fhn.d2fdxdp <- function(times,y,p,more)
{
  r = array(0,c(dim(y),2,length(p)))
  dimnames(r) = list(NULL,colnames(y),colnames(y),names(p))

  r[,'V','V','c'] = log(y[,'V'])-(log(y[,'V']))^3/3+log(y[,'R']) + 1- log(y[,'V'])^2
  r[,'V','R','c'] = y[,'V']/y[,'R']
  r[,'R','V','c'] = 1/p['c']^2*y[,'R']/y[,'V']
  r[,'R','R','a'] = 1/p['c']
  r[,'R','R','b'] = -1/p['c']*(1+ log(y[,'R']))
  r[,'R','R','c'] = p['b']/p['c']^2*log(exp(1)*y[,'R'])+1/p['c']^2*log(y[,'V']*exp(-p['a']))

  return(r)
}


make.exp2fhn <- function()
{
  return( 
    list(
        fn = exp2fhn.fun,
        fn.ode = exp2fhn.fun.ode,
        dfdx = exp2fhn.dfdx,
        dfdp = exp2fhn.dfdp,
        d2fdx2 = exp2fhn.d2fdx2,
        d2fdxdp = exp2fhn.d2fdxdp
    )
  )
}

