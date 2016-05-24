##################################################################
## spacings

spacing_equal <- function(sp = unit(0.3, "lines")) {
  if (!is.unit(sp)) sp <- unit(sp, "lines")
  function(d, condvars = NULL) lapply(d, function(x) if(x > 1) rep(sp, x - 1) else NA)
}
class(spacing_equal) <- "grapcon_generator"

spacing_dimequal <- function(sp) {
  if (!is.unit(sp)) sp <- unit(sp, "lines")
  function(d, condvars = NULL)
    lapply(seq_along(d), function(i) if(d[i] > 1) rep(sp[i], d[i] - 1) else NA)
}
class(spacing_dimequal) <- "grapcon_generator"

spacing_increase <- function(start = unit(0.3, "lines"), rate = 1.5) {
  if (!is.unit(start)) start <- unit(start, "lines")
  function(d, condvars = NULL) {
    sp <- start * rev(cumprod(c(1, rep.int(rate, length(d) - 1))))
    spacing_dimequal(sp)(d = d, condvars = condvars)
  }
}
class(spacing_increase) <- "grapcon_generator"

spacing_highlighting <- function(start = unit(0.2, "lines"), rate = 1.5) {
 if (!is.unit(start)) start <- unit(start, "lines")
  function(d, condvars = NULL)
    c(spacing_increase(start, rate)(d, condvars)[-length(d)],
      list(unit(rep(0, d[length(d)]), "lines")))
}
class(spacing_highlighting) <- "grapcon_generator"

spacing_conditional <- function(sp = unit(0.3, "lines"),
                                start = unit(2, "lines"), rate = 1.8) {
  condfun <- spacing_increase(start, rate)
  equalfun <- spacing_equal(sp)
  equalfun2 <- spacing_equal(start)
  function(d, condvars) {
    if (length(d) < 3)
      return(spacing_equal(sp)(d, condvars))
    condvars <- seq(condvars)
    ret <- vector("list", length(d))
    ret[condvars] <- if (length(condvars) < 3)
      equalfun2(d[condvars])
    else
      condfun(d[condvars])
    ret[-condvars] <- equalfun(d[-condvars])
    ret
  }
}
class(spacing_conditional) <- "grapcon_generator"
