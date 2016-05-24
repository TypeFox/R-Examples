`summary.obsSens` <-
function(object, digits=3, ...){
  d <- dim(object$beta)
  nr <- d[1]
  ne <- length(object$beta)
  d[1] <- 2*d[1]
  dn <- dimnames(object$beta)
  dn[[1]] <- c( rbind(dn[[1]], '') )
  out <- array('',d, dn)

  out[ seq(1, length=ne, by=2) ] <- paste('   ',format( object$beta, digits=digits))
  out[ seq(2, length=ne, by=2) ] <- paste( '(', format(object$lcl, digits=digits),
                           ',', format(object$ucl, digits=digits), ')', sep='')

  attr(out,'log') <- object$log
  attr(out,'xname') <- object$xname
  attr(out,'type') <- object$type
  class(out) <- 'summary.obsSens'

  out
}

