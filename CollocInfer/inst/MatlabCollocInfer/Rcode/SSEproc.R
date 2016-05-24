make.SSEproc <- function()
{

  ##############################################################################

  SSEproc <- function(coefs,bvals,pars,more)
  {
    devals  = as.matrix(bvals$bvals%*%coefs)
    ddevals = as.matrix(bvals$dbvals%*%coefs)
    colnames(devals)  = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    f = make.SSElik()$fn(ddevals,more$qpts,devals,pars,more)
    return( sum(f) )
  }

  ##############################################################################

  dSSEproc.dc <- function(coefs,bvals,pars,more)
  {
    devals  = as.matrix(bvals$bvals%*%coefs)
    ddevals = as.matrix(bvals$dbvals%*%coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames    
    #  g1 = dSSE.dx(ddevals,more$qpts,devals,pars,more)
    g1 = make.SSElik()$dfdx(ddevals,more$qpts,devals,pars,more)
    weights = checkweights(more$weights,more$whichobs,g1)
    g2 = weights*(ddevals-more$fn(more$qpts,devals,pars,more$more))
    g = as.vector( as.matrix(t(bvals$bvals)%*%g1 + 2*t(bvals$dbvals)%*%g2) )
    return(g)
  }

  ##############################################################################

  dSSEproc.dp <- function(coefs,bvals,pars,more)
  {
    devals  = as.matrix(bvals$bvals%*%coefs)
    ddevals = as.matrix(bvals$dbvals%*%coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames    
    #  g = dSSE.dp(ddevals,more$qpts,devals,pars,more)
    g = make.SSElik()$dfdp(ddevals,more$qpts,devals,pars,more)
    g = apply(g,2,sum)
    return(g)
  }

  ##############################################################################

  d2SSEproc.dc2 <- function(coefs,bvals,pars,more)
  {
    devals = as.matrix(bvals$bvals%*%coefs)
    ddevals =  as.matrix(bvals$dbvals%*%coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames    
    #  H1 = d2SSE.dx2(ddevals,more$qpts,devals,pars,more)
    H1 = make.SSElik()$d2fdx2(ddevals,more$qpts,devals,pars,more)
    H2 = more$dfdx(more$qpts,devals,pars,more$more)
    #  H = array(0,c(rep(dim(bvals$bvals)[2],2),rep(dim(devals)[2],2)))
    weights = checkweights(more$weights,more$whichobs,H1[,,1])
    H = list(len=dim(bvals$bvals)[2])
    for(i in 1:dim(devals)[2]){
      H[[i]] = list(len=dim(devals))
      for(j in 1:dim(devals)[2]){
        H[[i]][[j]] = t(bvals$bvals)%*%diag(H1[,i,j])%*%bvals$bvals - 
            2*t(bvals$dbvals)%*%diag(H2[,i,j]*weights[,i])%*%bvals$bvals -
            2*t(bvals$bvals)%*%diag(H2[,j,i]*weights[,j])%*%bvals$dbvals
      }
      H[[i]][[i]] = H[[i]][[i]] + 2*t(bvals$dbvals)%*%diag(weights[,i])%*%bvals$dbvals 
    }
    H = blocks2mat(H)
    return(H)
  }

  ##############################################################################

  d2SSEproc.dcdp <- function(coefs,bvals,pars,more)
  {
    devals  =  as.matrix(bvals$bvals%*%coefs)
    ddevals =  as.matrix(bvals$dbvals%*%coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames    
    #  H1 = d2SSE.dxdp(ddevals,more$qpts,devals,pars,more)
    H1 = make.SSElik()$d2fdxdp(ddevals,more$qpts,devals,pars,more)
    H2 = 2*more$dfdp(more$qpts,devals,pars,more$more)
    weights = checkweights(more$weights,more$whichobs,H1[,,1])
    H = c()
    for(i in 1:length(pars)){
      H = cbind(H,as.vector(as.matrix(t(bvals$bvals)%*%H1[,,i] -
            t(bvals$dbvals)%*%(weights*H2[,,i]))))
    }
    return(H)
  }

  ##############################################################################

  return(
        list(
            fn      = SSEproc,
            dfdc    = dSSEproc.dc,
            dfdp    = dSSEproc.dp,
            d2fdc2  = d2SSEproc.dc2,
            d2fdcdp = d2SSEproc.dcdp
       )
  )
}






