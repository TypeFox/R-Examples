genelist.rda <- function(x, y, alpha, delta, prior=table(y)/length(y),
                         gnames=NULL, regularization="S"){
  this.call <- match.call()
  
  #if ( missing(fit) ) {
  #  stop("A rda fit object must be supplied.")
  #}
  
  if ( missing(x) || missing(y) ) {
    stop("You forgot to supply the data.")
  }

  if ( missing(alpha) || missing(delta) ) {
    stop("You must supply both parameters alpha and delta.")
  }

  #if ( missing(prior) ) {
  #  prior <- fit$prior
  #}

  if ( is.null(gnames) || missing(gnames) ) {
    if ( is.null(dimnames(x)[[1]]) ) {
      gnames <- as.character(1:nrow(x))
    }
    else {
      gnames <- dimnames(x)[[1]]
    }
  }  

  if ( length(gnames) != nrow(x) ) {
    stop("The gene names you provided is not of the same length as the original gene list.")
  }

  tmp <- rda(x=x, y=y, alpha=alpha, delta=delta, prior=prior,
             regularization=regularization, genelist=TRUE, trace=FALSE)
  return(gnames[as.logical(tmp$gene.list[1, 1, ])])
}

  
