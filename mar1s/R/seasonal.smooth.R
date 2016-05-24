`seasonal.smooth` <-
function(x, basis = create.fourier.basis(nbasis = 3), lambda = 0, ...)
{
  result <- ts(start = c(0, 1), end = c(1, 0),
               frequency = frequency(x))  
  fd <- smooth.basisPar(time(x)[is.finite(x)], x[is.finite(x)],
                        basis, lambda = lambda)$fd  
  result[] <- predict(fd, as.vector(time(result)))
  attr(result, 'fd') <- fd
  return(result)
}
