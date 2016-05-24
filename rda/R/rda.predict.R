predict.rda <- function(object, x, y, xnew, prior, alpha, delta,
                        type=c("class", "posterior", "nonzero"), 
                        trace=FALSE, ...){
  this.call <- match.call()

  if ( missing(object) ) {
    stop("A rda fit object must be supplied.")
  }

  if ( missing(xnew) ) {
    stop("A new data to predict must be supplied.")
  }

  if ( missing(alpha) ) {
    alpha <- object$alpha
  }

  if ( missing(delta) ) {
    delta <- object$delta
  } 

  if ( missing(prior) ) {
    prior <- object$prior
  }
  
  type <- match.arg(type)
  switch(type, 
         class={
           tmp <- rda(x=x, y=y, xnew=xnew, ynew=NULL,
                      prior=prior, alpha=alpha, delta=delta,
                      regularization=object$reg, genelist=FALSE, trace=trace)
           tmp$call <- this.call
           drop(tmp$yhat.new)
         }, 
         posterior={
           tmp <- rda(x=x, y=y, xnew=xnew, ynew=NULL,
                      prior=prior, alpha=alpha, delta=delta,
                      regularization=object$reg, genelist=FALSE, trace=trace)
           tmp$call <- this.call
           tmp1 <- apply(exp(tmp$posterior), c(1, 2, 3), sum)
           tmp2 <- exp(tmp$posterior)
           for(i in 1:dim(tmp2)[1]){
             for(j in 1:dim(tmp2)[2]){
               for(k in 1:dim(tmp2)[3]){
                 tmp2[i, j, k, ] <- tmp2[i, j, k, ]/tmp1[i, j, k]
               }
             }
           }
           drop(tmp2)
         }, 
         nonzero={
           tmp <- rda(x=x, y=y, xnew=xnew, ynew=NULL,
                      prior=prior, alpha=alpha, delta=delta,
                      regularization=object$reg, genelist=TRUE, trace=trace)
           tmp$call <- this.call
           drop(tmp$gene.list)
         })
}

  
