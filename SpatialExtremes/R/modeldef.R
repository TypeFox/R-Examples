modeldef <- function(data, formula){
  ##This function computes the design matrix and the penalization
  ##matrix if penalized smoothing splines are used

  if (substr(formula[3], 1, 2) == "rb")
    ans <- modeldef.rb(data, formula)

  else{
    cov.names <- colnames(data)

    if (is.null(cov.names))
      stop("''coord'' and/or ``marg.cov'' must have named columns")

    cov.names.uniq <- unique(cov.names)
    
    if (length(cov.names.uniq) != length(cov.names)){
      warning("''coord'' and/or ``marg.cov'' have duplicates named columns. Omiting.")
      data <- data[,cov.names.uniq]
    }

    ans <- modeldef.lm(data, formula)
  }

  return(ans)
}

modeldef.lm <- function(data, formula){
  ##This function computes the design matrix from any valid R formula

  if (class(formula) != "formula")
    stop("''loc.form'', ''scale.form'' and ''shape.form'' must be R formulas")

  init.fun <- function(y)
    lm(formula, data = as.data.frame(cbind(y = y, data)))$coeff
  
  formula.tmp <- formula[-2]
  dsgn.mat <- model.matrix(formula.tmp, data = as.data.frame(data))

  ##The number of ``purely parametric'' parameters to be estimated
  n.ppar <- ncol(dsgn.mat)
  
  return(list(dsgn.mat = dsgn.mat, pen.mat = 0, degree = 0, knots = 0,
              type = "lm", penalty.tot = 0, formula = formula, data = data,
              init.fun = init.fun, n.ppar = n.ppar))
}

modeldef.rb <- function(data, formula){
  ##This function computes the design matrix X for using radial basis
  ##as well as the penalty matrix K

  model <- eval(formula[[3]], envir = as.data.frame(data))

  dsgn.mat <- model$dsgn.mat
  pen.mat <- model$pen.mat
  degree <- model$degree
  knots <- model$knots
  penalty <- model$penalty
  data <- model$data
  n.ppar <- model$n.ppar
  
  init.fun <- function(y)
    rbpspline(y, data, knots, degree, penalty)$beta

  return(list(dsgn.mat = dsgn.mat, pen.mat = pen.mat, degree = degree,
              knots = knots, type = "rb", penalty.tot = penalty^degree,
              init.fun = init.fun, penalty = penalty, formula = formula,
              data = data, n.ppar = n.ppar))
}
