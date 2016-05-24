make.branches.mkn.expokit <- function(cache, control) {
  control <- check.control.expokit(control)
  tol <- control$tol
  m <- control$m
  
  function(y, len, pars, t0, idx) {
    .Call("r_branches_mkn_expokit",
          pars[["Q"]], pars[["iq"]], pars[["jq"]], pars[["qnorm"]],
          len, y, m, tol)
  }
}

make.all.branches.mkn.expokit <- function(cache, control) {
  ## This message is important, as this doesn't always seem to work
  ## (with code -42, which suggests that I've missed something
  ## somewhere).
  message("Using experimental expokit code")
  branches <- make.branches.mkn.expokit(cache, control)
  function(pars, intermediates, preset=NULL) {
    pars.sparse <- expm.expokit.sparse.pars(pars)
    all.branches.matrix(pars.sparse, cache,
                        initial.conditions.mkn,
                        branches, preset)
  }
}

check.control.expokit <- function(control) {
  control <- modifyList(list(tol=1e-8, m=5), control)
  control$m <- check.integer(control$m)
  control
}
