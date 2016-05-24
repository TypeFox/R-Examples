########### Likelihood for Rates of Change ################ 
# 
# Here we note that the derivative is modeled as a respose
# dependent on the state, and represented as "y" in 
# in calls to gradient evalutations
########################################################### 

make.Cproc <- function()
{

Cproc <- function(coefs,bvals,pars,more) 
{ 
   devals = as.matrix(bvals$bvals%*%coefs)
   ddevals = as.matrix(bvals$dbvals%*%coefs)

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
   
   f = more$fn(ddevals,more$qpts,devals,pars,more$more)

   return( sum(f) )
}


dCproc.dc <- function(coefs,bvals,pars,more) 
{ 
   devals = as.matrix(bvals$bvals%*%coefs)
   ddevals =as.matrix( bvals$dbvals%*%coefs)

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

   g1 = more$dfdx(ddevals,more$qpts,devals,pars,more$more)
   g2 = more$dfdy(ddevals,more$qpts,devals,pars,more$more)

  g = as.vector( t(bvals$bvals)%*%g1 + t(bvals$dbvals)%*%g2 )

  return(g)
}


dCproc.dp <- function(coefs,bvals,pars,more) 
{ 
  devals = as.matrix(bvals$bvals%*%coefs)
  ddevals = as.matrix(bvals$dbvals%*%coefs)

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

  g = more$dfdp(ddevals,more$qpts,devals,pars,more$more)

  g = apply(g,2,sum)

  return(g)
}



d2Cproc.dc2 <- function(coefs,bvals,pars,more) 
{ 
  devals = as.matrix(bvals$bvals%*%coefs)
  ddevals = as.matrix(bvals$dbvals%*%coefs)

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

  H1 = more$d2fdx2(ddevals,more$qpts,devals,pars,more$more)
  H2 = more$d2fdxdy(ddevals,more$qpts,devals,pars,more$more)
  H3 = more$d2fdy2(ddevals,more$qpts,devals,pars,more$more)

#  H = array(0,c(rep(dim(bvals$bvals)[2],2),rep(dim(devals)[2],2)))

#  for(i in 1:dim(devals)[2]){
#	for(j in 1:dim(devals)[2]){
#		H[,,i,j] = t(bvals$bvals)%*%diag(H1[,i,j])%*%bvals$bvals +
#            t(bvals$bvals)%*%diag(H2[,i,j])%*%bvals$dbvals + 
#            t(bvals$dbvals)%*%diag(H2[,j,i])%*%bvals$bvals + 
#            t(bvals$dbvals)%*%diag(H3[,i,j])%*%bvals$dbvals 
#	}
#  }

#  H = array(0,c(rep(dim(bvals$bvals)[2],2),rep(dim(devals)[2],2)))

 H = list(len=dim(bvals$bvals)[2])
  for(i in 1:dim(devals)[2]){
  H[[i]] = list(len=dim(devals))
    for(j in 1:dim(devals)[2]){
        H[[i]][[j]] = t(bvals$bvals)%*%diag(H1[,i,j])%*%bvals$bvals +
            t(bvals$bvals)%*%diag(H2[,i,j])%*%bvals$dbvals + 
            t(bvals$dbvals)%*%diag(H2[,j,i])%*%bvals$bvals + 
            t(bvals$dbvals)%*%diag(H3[,i,j])%*%bvals$dbvals             
    }
  }


  H = blocks2mat(H)

  return(H)
}



d2Cproc.dcdp <- function(coefs,bvals,pars,more) 
{ 
  devals = as.matrix(bvals$bvals%*%coefs)
  ddevals = as.matrix(bvals$dbvals%*%coefs)

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

  H1 = more$d2fdxdp(ddevals,more$qpts,devals,pars,more$more)
  H2 = more$d2fdydp(ddevals,more$qpts,devals,pars,more$more)

  H = c()

  for(i in 1:length(pars)){
	  H = cbind(H,as.vector(t(bvals$bvals)%*%H1[,,i] +
			t(bvals$dbvals)%*%H2[,,i]))
  }

  return(H)
}



    return(
        list(
            fn = Cproc,
            dfdc = dCproc.dc,
            dfdp = dCproc.dp,
            d2fdc2 = d2Cproc.dc2,
            d2fdcdp = d2Cproc.dcdp
        )
    ) 
} 






