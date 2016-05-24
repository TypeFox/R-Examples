morph.metrop <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, morph,
    ...)
UseMethod("morph.metrop")

morph.metrop.morph.metropolis <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, morph, ...) {
  if (missing(morph)) {
    morph <- obj$morph
    obj$final <- obj$morph.final
  } else {
    # if the transformation was changed, transform the last state from the
    # original space to be the initial state.
    obj$final <- morph$transform(obj$final)
  }

  if (missing(outfun)) outfun <- obj$outfun
  if (missing(blen)) blen <- obj$blen
  if (missing(nspac)) nspac <- obj$nspac
  if (missing(debug)) debug <- obj$debug
  if (missing(scale)) scale <- obj$scale
  
  morphed.obj <- metrop.metropolis(obj,
                                   nbatch=nbatch,
                                   blen=blen,
                                   nspac=nspac,
                                   scale=scale,
                                   outfun=morph$outfun(outfun),
                                   debug=debug,
                                   ...)
  
  unmorphed.obj <- .morph.unmorph(morphed.obj, morph, outfun)
  return(unmorphed.obj)
}

morph.metrop.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, morph, ...) {

  if (missing(morph)) morph <- morph.identity()
  if (missing(outfun)) outfun <- NULL
  
  morphed.obj <- metrop.function(morph$lud(obj),
                                 initial=morph$transform(initial),
                                 nbatch=nbatch,
                                 blen=blen,
                                 scale=scale,
                                 outfun=morph$outfun(outfun),
                                 debug=debug,
                                 ...)
  
  unmorphed.obj <- .morph.unmorph(morphed.obj, morph, outfun)
  return(unmorphed.obj)
}

.morph.unmorph <- function(obj, morph, outfun) {
  obj$morph       <- morph
  obj$morph.final <- obj$final
  obj$final       <- morph$inverse(obj$final)
  obj$outfun      <- outfun
  class(obj) <- c("mcmc", "morph.metropolis")
  return(obj)
}
