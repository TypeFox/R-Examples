 ########## R function: aspmPenRead ##########

# For extracting information from a 
# penalized spline term.

# Last changed: 16 JUN 2006


"aspmPenRead" <-
  function (term) 
{
  arg.list <- substring(term, 3, (nchar(term) - 1))
  var.name <- break.string(arg.list, ",")[1]
  var.val <- eval(parse(text = var.name))
  out <- arg.search(arg.list, "knots=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    knots <- default.knots(var.val)
  if (present == TRUE) 
    knots <- spmArgRead(arg)$val
  out <- arg.search(arg.list, "var.knot=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    var.knot <- NULL
  if (present == TRUE) 
    var.knot <- spmArgRead(arg)$val
  out <- arg.search(arg.list, "adap=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    adap <- TRUE
  if (present == TRUE) 
    adap <- spmArgRead(arg)$val
  out <- arg.search(arg.list, "spar=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    spar <- NULL
  if (present == TRUE) 
    spar <- spmArgRead(arg)$val
  out <- arg.search(arg.list, "adf=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    adf <- "miss"
  if (present == TRUE) 
    adf <- spmArgRead(arg)$val
  out <- arg.search(arg.list, "basis=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    basis <- "tps"
  if (present == TRUE) 
    basis <- spmArgRead(arg)$val
  out <- arg.search(arg.list, "degree=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    if (basis == "trunc.poly") 
      degree <- 1
    else degree <- 3
  if (present == TRUE) 
    degree <- spmArgRead(arg)$val
  return(list(name = var.name, var = var.val, knots = knots, 
              adap = adap, var.knot = var.knot, spar = spar, adf = adf, 
              degree = degree, basis = basis))
}
