ror.addPreference <- function(ror, a, b) {
  a = a-1
  b = b-1
  .jcall(ror$model, "V", method="addPreference", as.integer(a), as.integer(b))
}

combineConstraints <- function(...) {
  allConst = list(...)

  lhs <- c()
  dir <- c()
  rhs <- c()

  for (const in allConst) {
    lhs <- rbind(lhs, const$lhs)
    dir <- c(dir, const$dir)
    rhs <- c(rhs, const$rhs)
  }

  return(list(lhs=lhs, dir=dir, rhs=rhs))
}
