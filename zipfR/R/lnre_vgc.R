lnre.vgc <- function (model, N, m.max=0, variances=FALSE)
{
  if (!inherits(model, "lnre")) stop("first argument must belong to a subclass of 'lnre'")
  if (! (is.numeric(N) && all(N >= 0)))
    stop("'N' argument must be a vector of non-negative numbers")
  if (!missing(m.max) && !(length(m.max) == 1 && is.numeric(m.max) && 0 <= m.max && m.max <= 9))
    stop("'m.max' must be a single integer in the range 1 ... 9")

  E.V <- EV(model, N)
  E.Vm.list <- list()
  if (m.max > 0) E.Vm.list <- lapply(1:m.max, function (.m) EVm(model, .m, N))

  if (!variances) {
    vgc(N=N, V=E.V, Vm=E.Vm.list, expected=TRUE)
  }
  else {
    V.V <- VV(model, N)
    V.Vm.list <- list()
    if (m.max > 0) V.Vm.list <- lapply(1:m.max, function (.m) VVm(model, .m, N))
    vgc(N=N, V=E.V, VV=V.V, Vm=E.Vm.list, VVm=V.Vm.list, expected=TRUE)
  }
}
