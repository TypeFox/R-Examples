
## =============================================================================
## Perturb     : adds the numerical differencing value to a value
## =============================================================================

# internal function #

perturb <- function (value , pert=1e-8) {
    # A small, positive value, not too small
  pmax(abs(value) * pert,pert)
}

