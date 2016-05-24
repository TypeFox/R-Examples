ssfa <- function(formula, data = NULL, data_w = NULL, intercept = TRUE, pars = NULL, par_rho = TRUE, 
                 form = "cost")
  {  
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- mf[-1]
  w <- as.matrix(data_w)
  tm <- attr(mf, "terms")
  intercept <- attr(tm, "intercept") == 1
  
  # Error messages
  if(!is.vector(y)) stop("y is not a vector")
  if(!is.data.frame(x)) stop("x is not a data frame")
  if(length(y) != nrow(x)) stop("x and y lengths differ")
  
  ssfa.fit(y = y, x = x, w=w, intercept = intercept, pars = pars, form=form, par_rho=par_rho)
}
