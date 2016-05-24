

.onLoad <- function(libname, pkgname) {
  if (is.null(rvpar())) {
    rvpar(
          oorm = TRUE,
          rvlwd = 2.0,
          rvcol = "default",
          rvlex = function (lwd) 1.5,
          rvpoint = c("median", "50%", "95%"),
          point.sample = 400,
          line.sample = 20,
          summary.dimnames = TRUE,
          summary.quantiles.numeric = c(0.01, 0.025, 0.25, 0.50, 0.75, 0.975, 0.99),
          summary.quantiles.integer = c(0, 0.01, 0.025, 0.25, 0.50, 0.75, 0.975, 0.99, 1),
          print.digits = 2
          )
  }
  if (is.null(rvpar("n.sims"))) {
    setnsims(4000)
  }
}

##.onAttach <- function(libname, pkgname) {
##  packageStartupMessage("Package rv attached.\n")
##}

