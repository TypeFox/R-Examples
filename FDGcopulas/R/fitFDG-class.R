## FDGcopula object 
setClass(
  Class="fitFDG",
  representation=representation(
    estimate = "numeric",
    var.est = "matrix",
    optimalvalues = "numeric",
    convergence = "list",
    FDGcopula = "FDGcopula")
  )
