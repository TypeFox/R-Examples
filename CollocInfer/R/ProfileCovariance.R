################################################################################
# This file contains definitions for the functions to obtain Newey-West based
# estimates of parameter estimates resulting for the generalized profile
# proceedure. 
#
# The main function is 'Profile.covariance'. 
#  Utility functions blocks2mat, trimr, Newey.West and NeweyWest.Var are
#    also defined.
################################################################################

Profile.covariance <- function(pars,active=NULL,times,data,coefs,lik,proc,
                               in.meth='nlminb',control.in=NULL,eps=1e-6,
                               GN=FALSE)
{
    # First, we'll need to allow for repeated measurements
    if(length(dim(data)) == 3){
       index = matrix(1:(dim(data)[1]*dim(data)[2]),dim(data)[1],dim(data)[2],
                      byrow=FALSE)
       times = rep(x=times,times=dim(data)[3])
       data = matrix(data,dim(data)[1]*dim(data)[2],dim(data)[3])
    }
    else{ index = matrix(1:dim(data)[1],dim(data)[1],1) }
    if(length(dim(coefs)) == 3){
       coefs = matrix(coefs,dim(coefs)[1]*dim(coefs)[2],dim(coefs)[3])
    }
 
    check.lik.proc.data.coefs(lik,proc,data,times,coefs)

    if(is.null(active)){ active = 1:length(pars) }

    apars = pars[active]
   
    g = ProfileDP(pars=apars,allpars=pars,times=times,data=data,coefs=coefs,
                  lik=lik,proc=proc,active=active,sumlik=FALSE)
    if(!is.matrix(g)){ g = matrix(g,length(g),length(active)) }

    if(!GN){

      H = matrix(0,length(apars),length(apars))
      #  see file OutOptimization for functions ProfileDP and ProfileDP.AllPar
      gg = ProfileDP(pars=apars,allpars=pars,times=times,data=data,
                     coefs=coefs,lik=lik,proc=proc,active=active,sumlik=TRUE)
      for(i in 1:length(apars)){
        if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
        if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
        if(file.exists('counter.tmp')){file.remove('counter.tmp')}

          tpars = apars
          tpars[i] = tpars[i] + eps

          tf = ProfileErr(tpars,pars,times,data,coefs,lik,proc,in.meth=in.meth,
                          control.in=control.in,active=active)
          tg = ProfileDP(tpars,pars,times,data,coefs,lik,proc,active=active,
                         sumlik=TRUE)

          H[,i] = (tg-gg)/eps

        if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
        if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
        if(file.exists('counter.tmp')){file.remove('counter.tmp')}
      }
    }
    else{
        H = t(g[,active])%*%g[,active]
    }
    
    Covar = 0
    for(i in 1:ncol(index)){       
      Covar = Covar + NeweyWest.Var( 0.5*(t(H)+H) ,g[index[,i],],5)
    }        
    return( Covar )
}



################################################################################
#
# Some Utilities
#
################################################################################

blocks2mat = function(H)   # List of matrices -> large matrix
{
  if(is.spam(H[[1]][[1]])){
      if(length(H[[1]])>1){ out = cbind(H[[1]][[1]],H[[1]][[2]]) }
      else{ out = H[[1]][[1]] }
     if(length(H[[1]]) > 2){for(j in 3:length(H[[1]])){
        out = cbind(out,H[[1]][[j]])
      }}  
  
     if(length(H) > 1){for(i in 2:length(H)){
      if(length(H[[i]])>1){ tout = cbind(H[[i]][[1]],H[[i]][[2]]) }
      else{ tout = H[[i]][[1]] }
      if(length(H[[i]])>2){for(j in 3:length(H[[i]])){
        tout = cbind(tout,H[[i]][[j]])
      }}
      out = rbind(out,tout)
    }}
  }
  else{
    rowdims = rep(0,length(H))
    coldims = rep(0,length(H[[1]]))
    for(i in 1:length(H[[1]])){ coldims[i] = ncol(H[[1]][[i]]) }
    for(i in 1:length(H)){ rowdims[i] = ncol(H[[i]][[1]]) }
    
    if(is.matrix(H[[1]][[1]])){ out = matrix(0,sum(rowdims),sum(coldims)) }
    else{ out = Matrix(0,sum(rowdims),sum(coldims),sparse=TRUE) }

    rowdims = cumsum(c(0,rowdims))
    coldims = cumsum(c(0,coldims))


    for(i in 1:length(H)){
      for(j in 1:length(H[[i]])){
        out[(rowdims[i]+1):rowdims[i+1],(coldims[j]+1):coldims[j+1]] = 
               H[[i]][[j]]
      }
    }
  }
  return(out)
}

################################################################################

#blocks2mat = function(H)  # List of matrices -> large matrix
#{
#
#    out = c()
#
#
#
#    for(i in 1:length(H)){
#      tout = H[[i]][[1]]
#      if(length(H[[i]])>1){
#        for(j in 2:length(H[[i]])){
#        print(c(i,j))
#            if(inherits(H[[i]][[j]],'dgCMatrix')|inherits(H[[i]][[j]],
#                         'dgeMatrix')){ tout=cBind(tout,H[[i]][[j]]) }
#            else{ tout = cbind(tout,H[[i]][[j]]) }
#        }
#      }
#      if(i > 1){
#        if(inherits(tout,'dgCMatrix')|inherits(tout,'dgeMatrix')){
#          print('still sparse')
#          print(c(dim(out),dim(tout)))
#          out=rBind(out,tout)
#          print(is.matrix(out))
#        }
#        else{ out = rbind(out,tout) }
#      }
#      else{ out = tout }
#    }
#
#    print('hello')
#    return(out)
#}

## Newey West Calculations, with thanks to Steve Ellner


## GAUSS trimr function: trims n1 rows from the start and n2 rows from the end
## of a matrix or vector 

################################################################################

trimr <- function (a,n1,n2) {
        da<-dim(a); 
        if(is.null(da)) {a[(n1+1):(length(a)-n2)]}
        else {a[(n1+1):(da[1]-n2),]}
}

################################################################################

Newey.West <-function(x,y,maxlag) {
        w=1-(1:maxlag)/(maxlag+1); w=w/length(x); 
        out=mean(x*y); 
        for(i in 1:maxlag) {
            out=out + w[i]*sum(trimr(x,i,0)*trimr(y,0,i)) +
                      w[i]*sum(trimr(y,i,0)*trimr(x,0,i))
        }
        return(out)     
} 

################################################################################

NeweyWest.Var = function(H,g,maxlag)       
{
    V = solve(H)
    I = 0*H    
    if(is.null(maxlag)){ 
        n = nrow(g)
        maxlag = max(5,n^(0.25))
    }              
    if(maxlag > 0){
        for(i in 1:ncol(g)){
            for(j in i:ncol(g)){
                I[i,j] = Newey.West(g[,i],g[,j],maxlag)
                I[j,i] = I[i,j]
            }    
        }
    }
    return( V%*%(I+ t(g)%*%g)%*%V  )
}

