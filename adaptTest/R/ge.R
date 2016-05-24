`ge` <-
function (x, y, tol=.Machine$double.eps^.5) identical(x>y || eq(x,y,tol=tol), TRUE)

