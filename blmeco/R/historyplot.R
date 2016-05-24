historyplot <- function(fit, parameter){
  
  if(class(fit)=="stanfit"){
    nchains <- fit@sim$chains
    matsim <- as.matrix(fit)
    namespars <- sapply(strsplit(colnames(matsim), "[", fixed=TRUE),function(x){x[[1]]})
    indexpar <- namespars==parameter
    matsim <- matsim[,indexpar]
    colkey <- rainbow(nchains)
    if(is.null(dim(matsim))){ 
      nsim <- length(matsim)/nchains
      par(mar=c(3,5,0.5, 0.5))
      plot(matsim[1:nsim], type="l", 
             ylab=parameter)
        for(i in 2:nchains) 
          lines(matsim[(nsim*(i-1)+1):(nsim*i)], 
                col=colkey[i])
      }
    if(!is.null(dim(matsim))){
      nsim <- nrow(matsim)/nchains
      npars <- ncol(matsim)
      nplots <- ceiling(npars/5)
      for(k in 1:nplots){
        par(mfrow=c(5,1), mar=c(3,5,0.5, 0.5), ask=TRUE)
          for(j in 1:min(npars-5*(k-1), 5)){
          plot(matsim[1:nsim,j+5*(k-1)], type="l", 
             ylab=colnames(matsim)[j+5*(k-1)])
          for(i in 2:nchains) 
            lines(matsim[(nsim*(i-1)+1):(nsim*i),j+5*(k-1)], 
                                col=colkey[i])
        }
      }
    }  
  }# close stanfit
  
  if(class(fit)=="bugs"){
    nchains <- fit$n.chains
    matsim <- fit$sims.array
    parnames <- dimnames(fit$sims.array)[[3]]
    namespars <- sapply(strsplit(parnames, "[", fixed=TRUE),function(x){x[[1]]})
    indexpar <- namespars==parameter
    matsim <- matsim[,,indexpar]
    colkey <- rainbow(nchains)
    if(length(dim(matsim))==2){ 
      nsim <- dim(matsim)[1]
      par(mar=c(3,5,0.5, 0.5))
      plot(matsim[,1], type="l", 
           ylab=parameter)
      for(i in 2:nchains) 
        lines(matsim[,i], 
              col=colkey[i])
    }
    if(length(dim(matsim))==3){
      nsim <- dim(matsim)[1]
      npars <- dim(matsim)[3]
      nplots <- ceiling(npars/5)
      for(k in 1:nplots){
        par(mfrow=c(5,1), mar=c(3,5,0.5, 0.5), ask=TRUE)
        for(j in 1:min(npars-5*(k-1), 5)){
          plot(matsim[,1,j+5*(k-1)], type="l", 
               ylab=dimnames(matsim)[[3]][j+5*(k-1)])
          for(i in 2:nchains) 
            lines(matsim[,i,j+5*(k-1)], 
                  col=colkey[i])
        }
      }
    }  
  }# close bugs  
}
  
  