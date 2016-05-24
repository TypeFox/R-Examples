plot.pda.fd = function(x, whichdim=1,npts=501,...)
{
  # This basically plots the elements of bwtlist, allowing the user
  # to specify how the functions are collected.

#  rangval = pdaList$resfdlist[[1]]$basis$rangeval
  rangval = x$resfdlist[[1]]$basis$rangeval

#  m = length(pdaList$resfdlist)
  m = length(x$resfdlist)
  tfine = seq(rangval[1],rangval[2],length.out=npts)

  whichdim=unique(sort(whichdim))

#  bwtlist = pdaList$bwtlist
  bwtlist = x$bwtlist

  # Firstly the one-variable case, do we plot all the functions
  # on one plot or not?
  if(m == 1){
    d = length(bwtlist)

    if(whichdim == 3){
      par(mfrow=c(d,1))
      for(i in 1:d){
        titlestr = paste('Coefficient for Derivative',i-1)
        plot(bwtlist[[i]]$fd,main=titlestr,...)
        }
    }
    else{
      betamat = matrix(0,npts,d)
      legendstr = c()

      for(i in 1:d){
        betamat[,i] = eval.fd(tfine,bwtlist[[i]]$fd)
        legendstr = c(legendstr,paste('Deriv',i))
      }
      xlabstr = names(bwtlist[[1]]$fd$fdnames)[[1]]
      ylabstr = names(bwtlist[[1]]$fd$fdnames)[[3]]

      matplot(tfine,betamat,type='l',lty=c(1:d),xlab=xlabstr,ylab=ylabstr,...)
      legend(x='topleft',legend=legendstr,lty=c(1:d),...)
    }
  }

  # Otherwise, we can plot by any combination of variables,
  # equations and derivatives.

  else{
    d = length(bwtlist[[1]][[1]])

    xlabstr = names(bwtlist[[1]][[1]][[1]]$fd$fdnames)[[1]]
    ylabstr = names(bwtlist[[1]][[1]][[1]]$fd$fdnames)[[3]]

    betamat = array(0,c(npts,m,m,d))
    legendstr = array('',c(m,m,d))

    for(i in 1:m){
      for(j in 1:m){
        for(k in 1:d){
                betamat[,i,j,k] = eval.fd(tfine,bwtlist[[i]][[j]][[k]]$fd)
                legendstr[i,j,k] = paste('var',i,'eq',j,'deriv',k)
        }
      }
    }

    if(length(whichdim)==1){
      if(whichdim==1){
        par(mfrow=c(m,1))
        for(i in 1:m){
          tbetamat = matrix(betamat[,i,,],npts,m*d,byrow=FALSE)
          tlegendstr = as.vector(legendstr[i,,])
          matplot(tfine,tbetamat,type='l',lty=c(1:(d*m)),col=c(1:(d*m)),xlab=xlabstr,ylab=ylabstr,...)
          legend(x='topleft',legend=tlegendstr,lty=c(1:(d*m)),col=c(1:(d*m)),...)
        }
      }
      if(whichdim==2){
        par(mfrow=c(m,1))
        for(j in 1:m){
          tbetamat = matrix(betamat[,,j,],npts,m*d,byrow=FALSE)
          tlegendstr = as.vector(legendstr[,j,])
          matplot(tfine,tbetamat,type='l',lty=c(1:(d*m)),col=c(1:(d*m)),xlab=xlabstr,ylab=ylabstr,...)
          legend(x='topleft',legend=tlegendstr,lty=c(1:(d*m)),col=c(1:(d*m)),...)
        }
      }
      if(whichdim==3){
        par(mfrow=c(d,1))
        for(k in 1:d){
          tbetamat = matrix(betamat[,,,k],npts,m*m,byrow=FALSE)
          tlegendstr = as.vector(legendstr[,,k])
          matplot(tfine,tbetamat,type='l',lty=c(1:(m*m)),col=c(1:(m*m)),xlab=xlabstr,ylab=ylabstr,...)
          legend(x='topleft',legend=tlegendstr,lty=c(1:(m*m)),col=c(1:(m*m)),...)
        }
      }
    }
    else if(length(whichdim)==2){
      if(whichdim[1]==1){
        if(whichdim[2]==2){
          par(mfrow=c(m,m))
          for(i in 1:m){
            for(j in 1:m){
              matplot(tfine,betamat[,i,j,],type='l',lty=c(1:d),col=c(1:d),xlab=xlabstr,ylab=ylabstr,...)
              legend(x='topleft',legend=legendstr[i,j,],lty=c(1:d),col=c(1:d),...)
            }
          }
        }
        if(whichdim[2]==3){
          par(mfrow=c(m,d))
          for(i in 1:m){
            for(k in 1:d){
              matplot(tfine,betamat[,i,,k],type='l',lty=c(1:m),col=c(1:m),xlab=xlabstr,ylab=ylabstr,...)
              legend(x='topleft',legend=legendstr[i,,k],lty=c(1:m),col=c(1:m),...)
            }
          }
        }
      }
      else{
        par(mfrow=c(m,d))
        for(j in 1:m){
          for(k in 1:d){
            matplot(tfine,betamat[,,j,k],type='l',lty=c(1:m),col=c(1:m),xlab=xlabstr,ylab=ylabstr,...)
            legend(x='topleft',legend=legendstr[,j,k],lty=c(1:m),col=c(1:m),...)
          }
        }
      }
    }
    else{
      for(j in 1:m){
        #X11()
        dev.new()
        par(mfrow=c(m,d))
        for(i in 1:m){
          for(k in 1:d){
            plot(bwtlist[[i]][[j]][[k]]$fd,main=legendstr[i,j,k],...)
          }
        }
      }
    }
  }
}
