
make.SEIR <- function()
{
SEIR.ode <- function(t,y,parms)
{

    p    = parms$p
    more = parms$more

    beta = more$beta.fun(t,p,more)


    tmpvec = beta*y[,'S']*(p['i']+y[,'I'])
    r = y
    r[,'S'] = -tmpvec + p['mu']                        - p['nu']*y[,'S']
    r[,'E'] =  tmpvec - p['sigma']*y[,'E']             - p['nu']*y[,'E']
    r[,'I'] =  p['sigma']*y[,'E'] - p['gamma']*y[,'I'] - p['nu']*y[,'I']

    return(list(r))
}


SEIR.fun <- function(t,y,p,more)
{

    beta = more$beta.fun(t,p,more)

    tmpvec = beta*y[,'S']*(p['i']+y[,'I'])
    r = y
    r[,'S'] = -tmpvec + p['mu']                        - p['nu']*y[,'S']
    r[,'E'] =  tmpvec - p['sigma']*y[,'E']             - p['nu']*y[,'E']
    r[,'I'] =  p['sigma']*y[,'E'] - p['gamma']*y[,'I'] - p['nu']*y[,'I']

    return(r)
}


SEIR.dfdx <- function(t,y,p,more)
{
    beta = more$beta.fun(t,p,more)

    r = array(0,c(length(t),ncol(y),ncol(y)))

    dimnames(r) = list(NULL,colnames(y),colnames(y))

    betaYI = beta*(p['i']+y[,'I'])
    betaYS = beta*y[,'S']
    r[,'S','S'] = -betaYI - p['nu']
    r[,'S','I'] = -betaYS
    r[,'E','S'] =  betaYI
    r[,'E','E'] = -p['sigma'] - p['nu']
    r[,'E','I'] =  betaYS
    r[,'I','E'] =  p['sigma']
    r[,'I','I'] = -p['gamma'] - p['nu']
 
    return(r)
}

SEIR.dfdp <- function(t,y,p,more)
{
    beta    = more$beta.fun(t,p,more)
    dbetadp = more$beta.dfdp(t,p,more)

    r = array(0,c(length(t),ncol(y),length(p)))

    dimnames(r) = list(NULL,colnames(y),names(p))

    tmpmat = diag(y[,'S']*(p['i']+y[,'I']))%*%dbetadp
    betaYS = beta*y[,'S']
    r[,'S',more$beta.ind] = -tmpmat
    r[,'E',more$beta.ind] =  tmpmat
    r[,'S','mu']    = 1
    r[,'S','i']     = -betaYS
    r[,'E','i']     =  betaYS
    r[,'S','nu']    = -y[,'S']
    r[,'E','nu']    = -y[,'E']
    r[,'I','nu']    = -y[,'I']
    r[,'E','sigma'] = -y[,'E']
    r[,'I','sigma'] =  y[,'E']   
    r[,'I','gamma'] = -y[,'I']
    return(r)
}


SEIR.d2fdx2 <- function(t,y,p,more)
{
    beta = more$beta.fun(t,p,more)

    r = array(0,c(length(t),ncol(y),ncol(y),ncol(y)))

    dimnames(r) = 
         list(NULL,colnames(y),colnames(y),colnames(y))

    r[,'S','S','I'] = -beta
    r[,'S','I','S'] = -beta
    r[,'E','S','I'] =  beta
    r[,'E','I','S'] =  beta

    return(r)
}


SEIR.d2fdxdp <- function(t,y,p,more)
{
    beta    = more$beta.fun(t,p,more)
    dbetadp = more$beta.dfdp(t,p,more)

    r = array(0,c(length(t),ncol(y),ncol(y),length(p)))
    dimnames(r) = list(NULL,colnames(y),colnames(y),names(p))

    tmpdiag1 = diag(p['i']+y[,'I']) %*% dbetadp
    tmpdiag2 = diag(y[,'S'])        %*% dbetadp

    r[,'S','S',more$beta.ind] = -tmpdiag1
    r[,'S','I',more$beta.ind] = -tmpdiag2
    r[,'E','S',more$beta.ind] =  tmpdiag1
    r[,'E','I',more$beta.ind] =  tmpdiag2
    r[,'S','S','i']     = -beta
    r[,'E','S','i']     =  beta    
    r[,'S','S','nu']    = -1
    r[,'E','E','nu']    = -1
    r[,'I','I','nu']    = -1   
    r[,'E','E','sigma'] = -1
    r[,'I','E','sigma'] =  1    
    r[,'I','I','gamma'] = -1

    return(r)
}



  return(
    list(
        ode.fn = SEIR.ode,
        fn = SEIR.fun,
        dfdx = SEIR.dfdx,
        dfdp = SEIR.dfdp,
        d2fdx2 = SEIR.d2fdx2,
        d2fdxdp = SEIR.d2fdxdp
    )
  )
}


make.var.SEIR <- function()
{

SEIR.var.fun <- function(t,y,p,more)
{
    beta = more$beta.fun(t,p,more)

    r = array(0,c(length(t),ncol(y),ncol(y)))
    dimnames(r) = list(NULL,colnames(y),colnames(y))
    
    r[,'S','S'] = p['mu'] + beta*y[,'S']*(p['i'] + 
                  y[,'I']) + p['nu']*y[,'S']
    r[,'E','E'] = beta*y[,'S']*(p['i'] +y [,'I']) + 
                  p['sigma']*y[,'E'] + p['nu']*y[,'E']
    r[,'I','I'] = p['sigma']*y[,'E'] + p['gamma']*y[,'I'] + 
                  p['nu']*y[,'I']

    r[,'S','E'] = -beta*y[,'S']*(p['i']+y[,'I'])
    r[,'E','S'] = r[,'S','E']

    r[,'E','I'] = p['sigma']*y[,'E']
    r[,'I','E'] = r[,'E','I']

    return(r)
}

SEIR.var.dfdx <- function(t,y,p,more)
{
    p = more$p.fun(t,more$pdef)
    beta = more$beta.fun(t,p,more)

    r = array(0,c(length(t),ncol(y),ncol(y),ncol(y)))
    dimnames(r) = list(NULL,colnames(y),colnames(y),
                       colnames(y))


    r[,'S','S','S'] = beta*(p['i']+y[,'I']) + p['nu']
    r[,'S','S','I'] = beta*y[,'S']

    r[,'E','E','S'] = beta*(p['i']+y[,'I'])
    r[,'E','E','E'] = p['sigma'] + p['nu']
    r[,'E','E','I'] = beta*y[,'S']

    r[,'I','I','E'] = p['sigma']
    r[,'I','I','I'] = p['gamma'] + p['nu']

    r[,'S','E','S'] = -beta*(p['i']+y[,'I'])
    r[,'S','E','I'] = -beta*y[,'S']
    r[,'E','S',] = r[,'S','E',]

    r[,'E','I','E'] = p['sigma']
    r[,'I','E',] = r[,'E','I',]

    return(r)
}

SEIR.var.dfdp <- function(t,y,p,more)
{
    beta    = more$beta.fun(t,p,more)
    dbetadp = more$beta.dfdp(t,p,more)

    r = array(0,c(length(t),ncol(y),ncol(y),length(p)))
    dimnames(r) = list(NULL,colnames(y),colnames(y),names(p))

    tmpdiag1 = diag(y[,'S']*(p['i'] + y[,'I']) )%*% dbetadp
    tmpdiag2 = diag(y[,'S']*(p['i']+y[,'I']))   %*% dbetadp
    r[,'S','S',more$beta.ind] =  tmpdiag1   
    r[,'E','E',more$beta.ind] =  tmpdiag2    
    r[,'S','E',more$beta.ind] = -tmpdiag1
    r[,'E','S',more$beta.ind] = r[,'S','E',more$beta.ind]

    r[,'S','S','mu'] = 1

    r[,'S','S','i'] =  beta*y[,'S']
    r[,'S','E','i'] = -beta*y[,'S']
    r[,'E','S','i'] = r[,'S','E','i']
    r[,'E','E','i'] =  beta*y[,'S']
    
    r[,'S','S','nu'] = y[,'S']
    r[,'E','E','nu'] = y[,'E']
    r[,'I','I','nu'] = y[,'I']
    
    r[,'E','E','sigma'] = y[,'E']
    r[,'E','I','sigma'] = y[,'E']
    r[,'I','E','sigma'] = y[,'E']
    
    r[,'I','I','gamma'] = y[,'I']

    return(r)
}


SEIR.var.d2fdx2 <- function(t,y,p,more)
{
    beta = more$beta.fun(t,p,more)

    r = array(0,c(length(t),rep(ncol(y),4)))
    dimnames(r) = list(NULL,colnames(y),colnames(y),colnames(y),colnames(y))

    r[,'S','S','S','I'] = beta
    r[,'S','S','I','S'] = beta

    r[,'E','E','S','I'] = beta
    r[,'E','E','I','S'] = beta

    r[,'S','E','S','I'] = -beta
    r[,'S','E','I','S'] = -beta
    r[,'E','S',,] = r[,'S','E',,]

    return(r)
}

SEIR.var.d2fdxdp <- function(t,y,p,more)
{
    beta    = more$beta.fun(t,p,more)
    dbetadp = more$beta.dfdp(t,p,more)

    r = array(0,c(length(t),rep(ncol(y),3),length(p)))
    dimnames(r) = list(NULL,colnames(y),colnames(y),colnames(y),names(p))

    r[,'S','S','S',more$beta.ind] = diag(p['i']+y[,'I'])%*%dbetadp
    r[,'S','S','I',more$beta.ind] = diag(y[,'S'])%*%dbetadp    

    r[,'E','E','S',more$beta.ind] = diag(p['i']+y[,'I'])%*%dbetadp
    r[,'E','E','I',more$beta.ind] = diag(y[,'S'])%*%dbetadp
    
    r[,'S','E','S',more$beta.ind] = -diag(p['i']+y[,'I'])%*%dbetadp
    r[,'S','E','I',more$beta.ind] = -diag(y[,'S'])%*%dbetadp

    r[,'S','S','S','i'] = beta
    r[,'E','E','S','i'] = beta
    r[,'S','E','S','i'] = -beta
    
    r[,'S','S','S','nu'] = 1
    r[,'E','E','E','nu'] = 1
    r[,'I','I','I','nu'] = 1
    
    r[,'E','E','E','sigma'] = 1
    r[,'E','I','E','sigma'] = 1
    r[,'I','E','I','sigma'] = 1

    r[,'I','I','I','gamma'] = 1
    
    r[,'E','S',,] = r[,'S','E',,]

    return(r)
}



  return(
    list(
        var.fn = SEIR.var.fun,
        var.dfdx = SEIR.var.dfdx,
        var.dfdp = SEIR.var.dfdp,
        var.d2fdx2 = SEIR.var.d2fdx2,
        var.d2fdxdp = SEIR.var.d2fdxdp
    )
  )
}
