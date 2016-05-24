spc.vector <- function (obj, m.min=1, m.max=15, all=FALSE)
{
  if (! inherits(obj, "spc")) stop("argument must be object of class 'spc'")
  if (all) m.max <- max(obj$m)
  if (m.min < 1) stop("'m.min' must be >= 1")
  if (m.max < m.min) stop("'m.max' must be >= 'm.min'")

  Vm.spc(obj, m.min:m.max)
}
