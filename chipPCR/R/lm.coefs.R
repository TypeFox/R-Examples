lm.coefs <- function(x, y, method.reg) {
  
  # Define the method for the linear regression
  # lmrob (robustbase) uses MM-type estimators, rfit (Rfit)
  # uses a rank-based estimation model for linear regression,
  # and lm (stats) an ordinary least squares 
  
  method.reg <- check.method(c("lmrob", "rfit", "least", "rq"), method.reg)
  
  lm.fit <- function(y, x, method) {
    switch(method,
           lmrob = do.call(function(x, y) lmrob(y ~ x), c(list(x = x, y = y))),
           #TO DO: check if Rfit can be called other way
           rfit = do.call(function(x, y) Rfit::rfit.default(y ~ x), c(list(x = x, y = y))),
           least = do.call(function(x, y) lm(y ~ x), c(list(x = x, y = y))),
           rq = do.call(function(x, y) rq(y ~ x), c(list(x = x, y = y)))
           )  
  }
  
  if(method.reg == "rfit") {
    op.warn <- options("warn")[[1]]
    options(warn=2)
  }
  tried.fit <- try(lm.fit(y, x, method = method.reg), 
                   silent = TRUE)
  
  if(method.reg == "rfit") {
    options(warn=op.warn)
  }
  
  
  if(class(tried.fit) != "try-error" && is.null(tried.fit[["converged"]]))
    tried.fit[["converged"]] <- TRUE
  if (class(tried.fit) != "try-error" && tried.fit[["converged"]] == TRUE) { 
    coefficients <- data.frame(tried.fit[1])
  } else { 
    warning(paste0("Chosen method ", method.reg, " failed to converge. 
                   Performed linear regression. "))
    coefficients <- data.frame(lm.fit(y, x, method = "least")[1]) 
  } 
  coefficients
}