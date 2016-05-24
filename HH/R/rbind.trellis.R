## rbind.trellis <- function(..., deparse.level=1,
##                           combineLimits=TRUE, useOuterStrips=TRUE) {
##   dddots <- list(...)
##   dimdddots <- sapply(dddots, dim)
##   if (is.list(dimdddots) || length(dim(dimdddots)) > 1)
##     stop("Only one-dimensional trellis objects can be used in rbind.trellis or cbind.trellis.",
##          call.=FALSE)
##   if (any(dimdddots != dimdddots[1]))
##     stop("All one-dimensional trellis objects in rbind.trellis or cbind.trellis must have the same dim value.",
##          call.=FALSE)
##   cdddots <- do.call(c, c(dddots, list(layout=c(dim(dddots[[1]]), length(dddots)))))
##   dddnames <- names(dddots)
##   if (is.null(dddnames)) dddnames <- LETTERS[1:length(dddots)]
##   cdddots$condlevels <- list(dddots[[1]]$condlevels[[1]], dddnames)
##   cdddots$index.cond <- list(1:length(dddots[[1]]$condlevels[[1]]), 1:length(dddots))
##   cdddots$perm.cond <- 1:2
##   dim(cdddots$packet.sizes) <- dim(cdddots)
##   dimnames(cdddots$packet.sizes) <- dimnames(cdddots)
##   if (useOuterStrips) cdddots <- useOuterStrips(cdddots)
##   if (combineLimits) cdddots <- combineLimits(cdddots)
##   cdddots
## }


rbind.trellis <- function(..., deparse.level=1,
                          combineLimits=TRUE, useOuterStrips=TRUE) {
  dddots <- list(...)
  dim.dddots <- lapply(dddots, dim)
  dim2.dddots <- sapply(dddots,
                      function(x) {
                        dim.x <- dim(x)
                        if (length(dim.x) > 2)
                          stop("Only one-or two-dimensional trellis objects can be used in rbind.trellis or cbind.trellis.",
                               call.=FALSE)
                        ## if (is.null(length(dim.x))) return(length(x)) ## not for trellis objects
                        if (length(dim.x) == 1)
                          c(dim.x, 1)
                        else
                          dim.x
                      })
  dimnames.dddots <- lapply(seq(along=dddots),
                            function(i) {
                              x <- dddots[[i]]
                              dn <- dimnames(x)
                              if (length(dn) == 1) {
                                name.i <- names(dddots)[i]
                                if (is.null(name.i) || nchar(name.i)==0)
                                  name.i <- LETTERS[i]
                                dn <- c(dn, list(name.i))
                              }
                              dn
                            })
  ddd.layout=c(sum(dim2.dddots[2,]), dim2.dddots[1,1])
  cdddots <- c(...,
               layout=ddd.layout,
               x.same=combineLimits,
               y.same=combineLimits)

  dimnames.start <- list(
    unlist(sapply(dimnames.dddots, `[[`, 2)),
    dimnames.dddots[[1]][[1]])

  mdddots <- matrix.trellis(cdddots, nrow=ddd.layout[1], ncol=ddd.layout[2], byrow=TRUE,
                            dimnames=list(
                              dimnames.start[[1]][1:ddd.layout[1]],
                              dimnames.start[[2]][1:ddd.layout[2]]))

  if (useOuterStrips)
    useOuterStrips(mdddots)
  else
    mdddots
}

cbind.trellis <- function(..., deparse.level=1,
                           combineLimits=TRUE, useOuterStrips=TRUE) {
  dddots <- list(...)
  t.dddots <- lapply(dddots, transpose)
  ddd.names <- names(dddots)
  if (is.null(ddd.names))
    ddd.names <- LETTERS[seq(along=dddots)]
  missing.names <- is.na(ddd.names) | nchar(ddd.names)==0
  ddd.names[missing.names] <- LETTERS[missing.names]

  for (ii in seq(along=t.dddots)) {
    if (length(dimnames(dddots[[ii]]))==1)
      dimnames(t.dddots[[ii]])[[2]] <- ddd.names[ii]
  }

  result <- transpose(do.call(rbind.trellis,
                              c(t.dddots,
                                list(deparse.level=deparse.level,
                                     combineLimits=combineLimits, useOuterStrips=useOuterStrips))))
  if (useOuterStrips)
    useOuterStrips(result)
  else
    result
}
