##' simple sim.
##'
##' Deprecated. Should just use sample or many other means
##' @param no.samples number of samples
##' @param f function to simulate from
##' @param ... passed to f
##' @return a sample
##' @export
simple.sim <- function(no.samples,f,...) {
  sample <-1:no.samples
  for (i in 1:no.samples) {
    sample[i]<-f(...)
  }
  sample
}
