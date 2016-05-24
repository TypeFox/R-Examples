##' (Noisy) BBOB 2009 test function generator.
##'
##' @param dimensions Size of parameter space.
##' @param fid Function ID, valid values are 1 to 24.
##' @param iid Instance ID, defaults to \code{1L}.
##' @param noiseSeed Seed for the noise random number generator,
##'   defaults to \code{1L}.
##' @return A \code{soo_function}.
##'
##' @note Due to the way the instances are generated, the function
##' value at the optimal parameter settings (as returned by
##' \code{\link{global_minimum}}) may differ slightly from the optimal
##' function value. These differences are in the order of
##' \eqn{10^{-16}}. 
##'
##' Also note that the random number generator used for the noisy test
##' functions is shared by all instanciated test functions. This means
##' that if you run multiple trials in parallel  within the same
##' interpreter, your results will not necessarily be repeatable.
##' 
##' @author \R interface by Olaf Mersmann. Original C code graciously
##' provided by the BBOB team (Anne Auger, Hans-Georg Beyer, Nikolaus
##' Hansen, Steffen Finck, Raymond Ros, Marc Schoenauer, Darrell
##' Whitley)
##' 
##' @references For a complete description of the 24 test functions
##' and much more background information please see the
##' \href{http://coco.gforge.inria.fr/doku.php?id=bbob-2009-downloads}{BBOB homepage}
##' 
##' @export
##' @useDynLib soobench do_bbob_opt do_bbob_eval do_set_bbob_noise_seed
##' @rdname bbob2009_function.Rd
bbob2009_function <- function(dimensions, fid, iid) {
  fid <- as.integer(fid)
  stopifnot(fid < 25)
  iid <- as.integer(iid)
  opt <- .Call(do_bbob_opt, fid, iid, as.integer(dimensions))
  f <- function(x) {}
  body(f) <- substitute(.Call(do_bbob_eval, fid, iid, x),
                        list(fid=fid, iid=iid))
  
  soo_function(name=sprintf("BBOB 2009 problem %i (instance %i)", fid, iid),
               id=sprintf("bbob2009-%02i-%02i-%id", fid, iid, dimensions),
               fun=f,
               dimensions=dimensions,
               lower_bounds=rep(-5, dimensions),
               upper_bounds=rep(5, dimensions),
               best_par=opt[[1]],
               best_value=opt[[2]])
}

##' @rdname bbob2009_function.Rd
##' @export
noisy_bbob2009_function <- function(dimensions, fid, iid, noiseSeed=1L) {
  fid <- as.integer(fid)
  realfid <- 100L + fid
  stopifnot(100 < realfid, realfid < 125)
  iid <- as.integer(iid)
  .Call(do_set_bbob_noise_seed, as.integer(noiseSeed))
  opt <- .Call(do_bbob_opt, realfid, iid, as.integer(dimensions))
  f <- function(x) {}
  body(f) <- substitute(.Call(do_bbob_eval, fid, iid, x),
                        list(fid=realfid, iid=iid))

  if ((104 <= realfid && realfid <= 106) ||
      (110 <= realfid && realfid <= 112)) {
    scale <- max(1, sqrt(dimensions) / 8)
    opt[[1]] <- scale * 0.75 * opt[[1]]
  }

  soo_function(name=sprintf("Noisy BBOB 2009 problem %i (instance %i)", fid, iid),
               id=sprintf("noisy-bbob2009-%02i-%02i-%id", fid, iid, dimensions),
               fun=f,
               dimensions=dimensions,
               lower_bounds=rep(-5, dimensions),
               upper_bounds=rep(5, dimensions),
               best_par=opt[[1]],
               best_value=opt[[2]])
}
