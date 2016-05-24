########## R function: aspmKrigeRead ##########

# For extracting information from a bivariate
# penalized spline term.

# Last changed: 16 JUN 2006


"aspmKrigeRead" <-
  function (term) 
{
  arg.list <- substring(term, 3, (nchar(term) - 1))
  var.name <- break.string(arg.list, ",")[1:2]
  var.val <- cbind(eval(parse(text = var.name[1])), eval(parse(text = var.name[2])))
  out <- arg.search(arg.list, "knots=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) {
    out <- arg.search(arg.list, "k=")
    arg <- out$arg
    if (!is.null(arg)) {
      num.knots <- spmArgRead(out$arg)$val
      knots <- default.knots.2D(var.val[, 1], var.val[, 2], num.knots)
    }
    else {
      knots <- default.knots.2D(var.val[, 1], var.val[,2])
      num.knots <- nrow(knots)
    }
  }
  if (present == TRUE) {
    knots <- spmArgRead(arg)$val
    num.knots <- nrow(knots)
  }
  out <- arg.search(arg.list, "var.knot=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    var.knot <- NULL
  if (present == TRUE) {
    var.knot <- spmArgRead(arg)$val
  }
  out <- arg.search(arg.list, "adap=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    adap <- TRUE
  if (present == TRUE) 
    adap <- spmArgRead(arg)$val
  out <- arg.search(arg.list, "degree=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    degree <- 2
  if (present == TRUE) 
    degree <- spmArgRead(arg)$val
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
  out <- arg.search(arg.list, "bdry=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) 
    bdry <- NULL
  if (present == TRUE) 
    bdry <- spmArgRead(arg)$val
  return(list(name = var.name, var = var.val, knots = knots, 
              adap = adap, var.knot = var.knot, num.knots = num.knots, 
              spar = spar, adf = adf, bdry = bdry, degree = degree))
}
