plot.fwsim <-
function(x, which = 1L, ...) {
  if (!is(x, "fwsim")) stop("x must be a fwsim object")

  plot(1L:x$pars$G, x$expected_pop_sizes, type = "b", lty = 2, xaxt = "n", xlab = "Generation", ylab = "Population size")
  axis(1L:x$pars$G, 1L:x$pars$G)
  points(1L:x$pars$G, x$pop_sizes, type = "b", lty = 1)
  legend("topleft", c("Actual", "Expected"), lty = 1:2)
    
  return(invisible(NULL))
}

