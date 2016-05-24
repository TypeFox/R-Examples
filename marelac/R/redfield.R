## -----------------------------------------------------------------------------
## the Redfield ratio
## -----------------------------------------------------------------------------
redfield <- function(q, species, method = c("mol", "mass"),
                     ratio = c(C = 106, H = 263, O = 110, N = 16, P = 1)) {
  if (!is.numeric(q)) stop("q must be numeric")
  method   <- match.arg(method)
  speciess <- names(ratio)
  species  <- match.arg(species, speciess)

  if (method == "mass") {
    ratio <- with(marelac::atomicweight,{
      ratio * c(C, H, O, N, P)
    })
  }
  p <- match(species, speciess)
  as.data.frame(q %o% ratio / ratio[p])
}

