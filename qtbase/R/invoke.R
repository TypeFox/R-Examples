qinvoke <- function(x, method, ...) {
  .Call("qt_qinvoke", x, method, FALSE, list(...), PACKAGE="qtbase")
}

qinvokeStatic <- function(x, method, ...) {
  .Call("qt_qinvokeStatic", x, method, list(...), PACKAGE="qtbase")
}

qinvokeSuper <- function(x, method, ...) {
  .Call("qt_qinvoke", x, method, TRUE, list(...), PACKAGE="qtbase")
}
