# ========================================================================
# rvRhat; rvneff; rvnchains - convenience functions for some attributes
# ========================================================================
#

rvRhat <- function(x) {
  unlist(rvattr(x, "Rhat"))
}

rvneff <- function(x) {
  unlist(rvattr(x, "n.eff"))
}

rvnchains <- function(x) {
  unlist(rvattr(x, "n.chains"))
}

