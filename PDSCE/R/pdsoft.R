pdsoft <-
function(s, lam,  tau=1e-4, init=c("soft", "diag", "dense", "user"),
s0=NULL, i0=NULL, standard=TRUE,tolin=1e-8, tolout=1e-8, maxitin=1e4, 
                  maxitout=1e3, quiet=TRUE)
{
  p=dim(s)[1]
  init=match.arg(init)

  if(sum(lam)==0) init="dense"
    

  if(standard)
  {
    dhat = sqrt(diag(s))
    dhat.inv = 1/dhat
    S=dhat.inv * s * rep(dhat.inv, each = p)

    if(init=="diag")
    {
      s0=diag(p)
      i0=diag(p)
    } else if( init=="soft")
    {
      tmp=abs(S) - lam
      tmp=tmp*(tmp > 0)
      S.soft=sign(S)*tmp
      diag(S.soft)=diag(S)
      oe=eigen(S.soft, symmetric=TRUE)
      evs=oe$val/2 + sqrt( oe$val^2 + 4*tau) /2
      s0=tcrossprod(oe$vec*rep(evs, each=p), oe$vec)
      i0=tcrossprod(oe$vec*rep(1/evs, each=p), oe$vec)
    } else if( init=="dense")
    {
      oe=eigen(S, symmetric=TRUE)
      evs=oe$val/2 + sqrt( oe$val^2 + 4*tau) /2
      s0=tcrossprod(oe$vec*rep(evs, each=p), oe$vec)
      i0=tcrossprod(oe$vec*rep(1/evs, each=p), oe$vec)
    } else
    {
      if( is.null(s0) & (!is.null(i0)) )
        s0=qr.solve(i0)
      if( (!is.null(s0)) & is.null(i0) )
        i0=qr.solve(s0)
      if(is.null(s0) & is.null(i0) )
        stop("Error: must specify either s0 or i0 with user initialization") 
    }
  }else ## not standard
  {
    S=s
    if(init=="diag")
    {
      s0=diag(diag(S))
      i0=diag(1/diag(S))
    } else if( init=="soft")
    {
      tmp=abs(S) - lam
      tmp=tmp*(tmp > 0)
      S.soft=sign(S)*tmp
      diag(S.soft)=diag(S)
      oe=eigen(S.soft, symmetric=TRUE)
      evs=oe$val/2 + sqrt( oe$val^2 + 4*tau) /2
      s0=tcrossprod(oe$vec*rep(evs, each=p), oe$vec)
      i0=tcrossprod(oe$vec*rep(1/evs, each=p), oe$vec)      
    } else if( init=="dense")
    {
      oe=eigen(S, symmetric=TRUE)
      evs=oe$val/2 + sqrt( oe$val^2 + 4*tau) /2
      s0=tcrossprod(oe$vec*rep(evs, each=p), oe$vec)
      i0=tcrossprod(oe$vec*rep(1/evs, each=p), oe$vec)
    } else
    {
      if( is.null(s0) & (!is.null(i0)) )
        s0=qr.solve(i0)
      if( (!is.null(s0)) & is.null(i0) )
        i0=qr.solve(s0)
      if(is.null(s0) & is.null(i0) )
        stop("Error: must specify either s0 or i0 with user initialization") 
    }   
  }  

  if(is.null(dim(lam)[1]))
  {
    lam = matrix(lam, nrow=p, ncol=p)
  }

  if( sum(lam) > 0 )
  {
    Soff=S
    diag(Soff)=0
    tolmult=sum(abs(Soff))/2
    S = as.double(S)
    i0 = as.double(i0)
    s0 = as.double(s0)
    lam = as.double(lam)
    b = as.double(tau)
    tolin = as.double(tolin*(tolmult/p))
    tolout = as.double(tolout*tolmult)
    totalout=1
    mode(totalout) = "integer"
    mode(p) = "integer"
    mode(maxitin) = "integer"
    mode(maxitout) = "integer"
    mode(quiet) = "integer"
    tmp=matrix(0, nrow=p, ncol=p)
    tmp=as.double(tmp)
    tmp2=matrix(0, nrow=p, ncol=p)
    tmp2=as.double(tmp2)

    tosoft=as.double(rep(0,p))
    
    coutput=.C("pdsc",S=S, Sigma=s0, Omega=i0, tosoft=tosoft, pin=p, lam=lam, tauin=b, tolin=tolin,
            maxitin=maxitin, tolout=tolout, maxitout=maxitout, totalout=totalout)

    sigma = matrix(coutput$Sigma, nrow=p, ncol=p)
    omega = matrix(coutput$Omega, nrow=p, ncol=p)
  
    if(!quiet)
    {
      cat("Total outer iterations = ", coutput$totalout, "\n")
    }
  } else ## lam == 0
  {
    sigma=s0
    omega=i0
  }
 
  if(standard)
  {
    theta=sigma
    theta.inv=omega
    sigma = dhat*theta*rep(dhat, each=p)
    omega = dhat.inv *theta.inv * rep(dhat.inv, each=p)
  } else
  {
    theta=NULL
    theta.inv=NULL
  }
  return(list(omega=omega, sigma=sigma, theta=theta, theta.inv=theta.inv))
}

