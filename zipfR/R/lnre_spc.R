lnre.spc <- function (model, N=NULL, variances=FALSE, m.max=100)
{
  if (! inherits(model, "lnre")) stop("first argument must belong to a subclass of 'lnre'")
  if (missing(N)) N <- N(model)  # if specified as default value, R complains about "recursive reference"
  if (! (is.numeric(N) && length(N) == 1 && all(N > 0)))
    stop("'N' argument must be a single positive number")
  if (variances && missing(m.max)) m.max <- 50

  m.vec <- 1:m.max
  E.Vm <- EVm(model, m.vec, N)
  E.V <- EV(model, N)

  if (variances) {
    V.Vm <- VVm(model, m.vec, N)
    V.V <- VV(model, N)
    spc(Vm=E.Vm, m=m.vec, VVm=V.Vm, N=N, V=E.V, VV=V.V, m.max=m.max, expected=TRUE)
  }
  else {
    spc(Vm=E.Vm, m=m.vec, N=N, V=E.V, m.max=m.max, expected=TRUE)
  }  
}
