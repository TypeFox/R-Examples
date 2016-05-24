`eq` <-
function (x, y, tol=.Machine$double.eps^.5) identical(all.equal(x, y, tolerance=tol), TRUE)

