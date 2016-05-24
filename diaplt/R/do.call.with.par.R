# do.call.with.par
# http://code.google.com/p/cowares-excel-hello/source/browse/trunk/util_r/
#
# Copyright (C) 2013 Tomizono
# Fortitudinous, Free, Fair, http://cowares.nobody.jp
#                            http://paidforeveryone.wordpress.com/

# extract graphic parameters from dots argument
# and call some graphical functions
#
#  do.call.with.par('plot.default', list(...), 
#                   xlim=lim[,1], ylim=lim[,2], 
#                   x=NA, xlab='', ylab='', xaxt=axt[1], yaxt=axt[2])
#  do.call.with.par('text', list(...), dots.win=T,
#                   x=pt[1,1], y=pt[3,2], col=col.label,
#                   labels=text)

only.graphic.pars <- function(pars, what='plot.default') {
  if(length(pars) == 0) return(list())

  gnames <- names(par(no.readonly=T))
  pnames <- names(formals(args(what)))
  gpnames <- unique(c(pnames[pnames != '...'], gnames))

  selectgname <- pars[gpnames]
  to.rm.na <- names(selectgname)

  if(is.null(to.rm.na)) list() else 
    selectgname[!is.na(to.rm.na)]
}

do.call.with.par <- function(what, args, dots.win=FALSE, ...) {
  pars <- 
    if(dots.win) {
      modifyList(only.graphic.pars(args, what), list(...))
    } else {
      modifyList(list(...), only.graphic.pars(args, what))
    }
  do.call(what, pars)
}

extract.from.dots <- function(which, default=NULL, ...) {
  x <- list(...)[[which]]
  if(is.null(x)) default else x
}

