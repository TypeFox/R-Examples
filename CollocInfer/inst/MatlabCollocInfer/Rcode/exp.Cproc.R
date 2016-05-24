########### Likelihood for Rates of Change ################ 
# 
# Here we note that the derivative is modeled as a respose
# dependent on the state, and represented as "y" in 
# in calls to gradient evalutations
########################################################### 


make.exp.Cproc <- function()
{

exp.Cproc <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals$bvals%*%coefs))
   ddevals = (as.matrix(bvals$dbvals%*%coefs))*devals

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    
   f = more$fn(ddevals,more$qpts,devals,pars,more$more)

   return( sum(f) )
}


exp.dCproc.dc <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals$bvals%*%coefs))
   ddevals = (as.matrix(bvals$dbvals%*%coefs))*devals

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    
   g1 = more$dfdx(ddevals,more$qpts,devals,pars,more$more)
   g2 = more$dfdy(ddevals,more$qpts,devals,pars,more$more)

   g = as.vector( t(bvals$bvals)%*%(g1*devals) + 
                 t(bvals$dbvals)%*%(g2*devals) +
                 t(bvals$bvals)%*%(g2*ddevals) )

   return(g)
}


exp.dCproc.dp <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals$bvals%*%coefs))
   ddevals = (as.matrix(bvals$dbvals%*%coefs))*devals

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

   g = more$dfdp(ddevals,more$qpts,devals,pars,more$more)

   g = apply(g,2,sum)

   return(g)
}



exp.d2Cproc.dc2 <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals$bvals%*%coefs))
   ddevals = (as.matrix(bvals$dbvals%*%coefs))*devals

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

  g1 = more$dfdx(ddevals,more$qpts,devals,pars,more$more)
  g2 = more$dfdy(ddevals,more$qpts,devals,pars,more$more)

  H1 = more$d2fdx2(ddevals,more$qpts,devals,pars,more$more)
  H2 = more$d2fdxdy(ddevals,more$qpts,devals,pars,more$more)
  H3 = more$d2fdy2(ddevals,more$qpts,devals,pars,more$more)

 # H = array(0,c(rep(dim(bvals$bvals)[2],2),rep(dim(devals)[2],2)))
   H = list(len=dim(bvals$bvals)[2])
   
  for(i in 1:dim(devals)[2]){
    ibvals = diag(devals[,i])%*%bvals$bvals
    idbvals = diag(ddevals[,i])%*%bvals$bvals+diag(devals[,i])%*%bvals$dbvals
    H[[i]] = list(len=dim(devals))
    for(j in 1:dim(devals)[2]){
        jbvals = diag(devals[,j])%*%bvals$bvals
        jdbvals = diag(ddevals[,j])%*%bvals$bvals+diag(devals[,j])%*%bvals$dbvals
        H[[i]][[j]] = t(ibvals)%*%diag(H1[,i,j])%*%jbvals +
            t(ibvals)%*%diag(H2[,i,j])%*%jdbvals + 
            t(idbvals)%*%diag(H2[,j,i])%*%jbvals + 
            t(idbvals)%*%diag(H3[,i,j])%*%jdbvals
    }
    H[[i]][[i]] = H[[i]][[i]] + t(ibvals)%*%diag(g1[,i])%*%bvals$bvals +
                t(idbvals)%*%diag(g2[,i])%*%bvals$bvals +
                t(ibvals)%*%diag(g2[,i])%*%bvals$dbvals    
    
#    t(bvals$dbvals)%*%diag(devals[,i]*g2[,i])%*%bvals$bvals +
#                t(bvals$bvals)%*%diag(devals[,i]*g2[,i]*devals[,i])%*%bvals$dbvals
#                t(bvals$bvals)%*%diag(ddevals[,i]*g2[,i]*devals[,i]+devals[,i]*g1[,i]*devals[,i])%*%bvals$bvals
  }

  H = blocks2mat(H)

  return(H)
}



exp.d2Cproc.dcdp <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals$bvals%*%coefs))
   ddevals = (as.matrix(bvals$dbvals%*%coefs))*devals

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

   H1 = more$d2fdxdp(ddevals,more$qpts,devals,pars,more$more)
   H2 = more$d2fdydp(ddevals,more$qpts,devals,pars,more$more)

   H = c()

   for(i in 1:length(pars)){
      H = cbind(H,as.vector(t(bvals$bvals)%*%(devals*H1[,,i] + ddevals*H2[,,i]) +
            t(bvals$dbvals)%*%(devals*H2[,,i])))
   }

   return(H)
}



    return(
        list(
            fn = exp.Cproc,
            dfdc = exp.dCproc.dc,
            dfdp = exp.dCproc.dp,
            d2fdc2 = exp.d2Cproc.dc2,
            d2fdcdp = exp.d2Cproc.dcdp
        )
    ) 
} 






