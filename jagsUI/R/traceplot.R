
#Handle masking issues

if (!isGeneric("traceplot")) {
  setGeneric("traceplot",
             function(x, ...)
               standardGeneric("traceplot"),
             useAsDefault = function(x, ...) coda::traceplot(x, ...))
}  

#Set traceplot method
setClass("jagsUI")

setMethod("traceplot", signature(x = "jagsUI"),
  function (x, parameters=NULL, ...) {
    if(class(x)!="jagsUI"){stop('Requires jagsUI object as input')}
    samples <- x$samples    
    if(is.null(parameters)){params <- names(as.data.frame(samples[[1]]))
    } else {params <- translate.params(x,parameters)}

    nparams <- length(params)
    nchains <- x$mcmc.info[[1]]
    rhat <- x$Rhat
            
    sep <- par()$ask
    par(ask=TRUE)
    
    xmax <- x$mcmc.info$n.samples / x$mcmc.info$n.chains
            
    col=c('red','blue','green','yellow','orange','violet')
    if(nchains>6){col=rainbow(nchains)}
            
    for (i in 1:nparams){
      if(nchains>1){
        rhat <- gelman.diag(samples[,params[i]],autoburnin=FALSE)$psrf[1]
        title <- paste('Trace of ',params[i],', Rhat = ',round(rhat,2),sep="")
      } else {title <- paste('Trace of ',params[i],sep="") }
        plot(x = 1:xmax, y = samples[,params[i]][[1]], main = title, xlab="Iterations", ylab="Value",type="l", col=col[1],
            ylim=range(samples[,params[i]]))
        if(nchains>1){
          for (j in 2:nchains){
            lines(x = 1:xmax, y = samples[,params[i]][[j]],type="l", col=col[j])
        }}
      }
            
    on.exit(par(ask=sep))            
  }
)