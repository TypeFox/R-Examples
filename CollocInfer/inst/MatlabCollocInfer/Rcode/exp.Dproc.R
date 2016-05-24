########### Likelihood for Discrete-Time Dynamics ################ 

make.exp.Dproc <- function()
{

exp.Dproc <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals[1:(nrow(bvals)-1),]%*%coefs))
   ddevals = exp(as.matrix(bvals[2:nrow(bvals),]%*%coefs))

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    
   f = more$fn(ddevals,more$qpts,devals,pars,more$more)

   return( sum(f) )
}


exp.dDproc.dc <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals[1:(nrow(bvals)-1),]%*%coefs))
   ddevals = exp(as.matrix(bvals[2:nrow(bvals),]%*%coefs))

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

   g1 = more$dfdx(ddevals,more$qpts,devals,pars,more$more)*devals
   g2 = more$dfdy(ddevals,more$qpts,devals,pars,more$more)*ddevals

  g = as.vector( t(bvals[1:(nrow(bvals)-1),])%*%g1 + t(bvals[2:nrow(bvals),])%*%g2 )

  return(g)
}


exp.dDproc.dp <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals[1:(nrow(bvals)-1),]%*%coefs))
   ddevals = exp(as.matrix(bvals[2:nrow(bvals),]%*%coefs))

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

  g = more$dfdp(ddevals,more$qpts,devals,pars,more$more)

  g = apply(g,2,sum)

  return(g)
}



exp.d2Dproc.dc2 <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals[1:(nrow(bvals)-1),]%*%coefs))
   ddevals = exp(as.matrix(bvals[2:nrow(bvals),]%*%coefs))

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
   
   g1 = more$dfdx(ddevals,more$qpts,devals,pars,more$more)*devals
   g2 = more$dfdy(ddevals,more$qpts,devals,pars,more$more)*ddevals   

  H1 = more$d2fdx2(ddevals,more$qpts,devals,pars,more$more)
  H2 = more$d2fdxdy(ddevals,more$qpts,devals,pars,more$more)
  H3 = more$d2fdy2(ddevals,more$qpts,devals,pars,more$more)

  #H = array(0,c(rep(dim(bvals[1:(nrow(bvals)-1),])[2],2),rep(dim(devals)[2],2)))
   H = list(len=dim(bvals)[2])
   
  for(i in 1:dim(devals)[2]){
    H[[i]] = list(len=dim(devals))
    for(j in 1:dim(devals)[2]){
        H[[i]][[j]] = t(bvals[1:(nrow(bvals)-1),])%*%diag(devals[,i]*H1[,i,j]*devals[,j])%*%bvals[1:(nrow(bvals)-1),] +
            t(bvals[1:(nrow(bvals)-1),])%*%diag(devals[,i]*H2[,i,j]*ddevals[,j])%*%bvals[2:nrow(bvals),] + 
            t(bvals[2:nrow(bvals),])%*%diag(ddevals[,i]*H2[,j,i]*devals[,j])%*%bvals[1:(nrow(bvals)-1),] + 
            t(bvals[2:nrow(bvals),])%*%diag(ddevals[,i]*H3[,i,j]*ddevals[,j])%*%bvals[2:nrow(bvals),] 
    }
    H[[i]][[i]] = H[[i]][[i]] + t(bvals[1:(nrow(bvals)-1),])%*%diag(g1[,i])%*%bvals[1:(nrow(bvals)-1),] +
                t(bvals[2:nrow(bvals),])%*%diag(g2[,i])%*%bvals[2:nrow(bvals),]
  }

  H = blocks2mat(H)

  return(H)
}



exp.d2Dproc.dcdp <- function(coefs,bvals,pars,more) 
{ 
   devals = exp(as.matrix(bvals[1:(nrow(bvals)-1),]%*%coefs))
   ddevals = exp(as.matrix(bvals[2:nrow(bvals),]%*%coefs))

    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames

  H1 = more$d2fdxdp(ddevals,more$qpts,devals,pars,more$more)
  H2 = more$d2fdydp(ddevals,more$qpts,devals,pars,more$more)

  H = c()

  for(i in 1:length(pars)){
      H = cbind(H,as.vector(t(bvals[1:(nrow(bvals)-1),])%*%(devals*H1[,,i]) +
            t(bvals[2:nrow(bvals),])%*%(ddevals*H2[,,i])))
  }

  return(H)
}



    return(
        list(
            fn = exp.Dproc,
            dfdc = exp.dDproc.dc,
            dfdp = exp.dDproc.dp,
            d2fdc2 = exp.d2Dproc.dc2,
            d2fdcdp = exp.d2Dproc.dcdp
        )
    ) 
} 






