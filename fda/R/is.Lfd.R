is.Lfd <- function(Lfdobj) {
#  check whether LFDOBJ is a linear differential operator
  if (inherits(Lfdobj, "Lfd")) return(TRUE) else return(FALSE)
}
