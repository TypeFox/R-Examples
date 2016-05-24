values <-
function(x, data=mydata, ...) {


  # get actual variable name before potential call of data$x
  x.name <- deparse(substitute(x)) 

  # get data frame name
  dname <- deparse(substitute(data))

  # get conditions and check for data existing
  xs <- .xstatus(x.name, dname)
  in.global <- xs$ig 

  # see if variable exists in the data frame, if x not in Global Env or function call 
  if (!missing(x) && !in.global) .xcheck(x.name, dname, data)

  if (!in.global) x.call <- eval(substitute(data$x))
  else {  # vars that are function names get assigned to global
    x.call <- x
    if (is.function(x.call)) x.call <- eval(substitute(data$x))
  }

  print(x.call, ...)

}
