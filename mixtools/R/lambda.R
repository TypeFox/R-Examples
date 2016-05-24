lambda=function (z, x, xi, h = NULL, kernel = c("Gaussian", "Beta", 
    "Triangle", "Cosinus", "Optcosinus"), g = 0) 
{
  if (is.null(h)) {
    cat("WARNING! BANDWIDTH MUST BE SPECIFIED!", "\n")
  }
  X.mat=cbind(1,(xi-x))
  kernel <- match.arg(kernel)
  inwindow <- (abs((xi - x)/h) <= 1)
  if (kernel == "Gaussian") {
    W=(kern.G(x, xi, h) * inwindow)
  }
  else if (kernel == "Beta") {
    W=(kern.B(x, xi, h, g) * inwindow)
  }
  else if (kernel == "Triangle") {
    W=(kern.T(x, xi, h) * inwindow)
  }
  else if (kernel == "Cosinus") {
    W=(kern.C(x, xi, h) * inwindow)
  }
  else if (kernel == "Optcosinus") {
    W=(kern.O(x, xi, h) * inwindow)
  }
  A=try(solve(t(X.mat)%*%(W*X.mat)), silent=TRUE)
  if(class(A)=="try-error") {
    A=ginv(t(X.mat)%*%(W*X.mat))
  }
  beta.x=A%*%t(X.mat)%*%(W*z)
  beta.x
}