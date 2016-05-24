######################################################################
# Formula interface to ccif
######################################################################

ccif.formula <-  function(formula, 
                          data,
                          ...) {
                                
  if (!inherits(formula, "formula"))
    stop("uplift: method is only for formula objects")
  
  
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data"), 
                names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  special <- "trt"
  mt <- if(missing(data)) terms(formula, special) else terms(formula, special, data = data)
  mf$formula <- mt  
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  attr(Terms, "intercept") <- 0
  trt.var <- attr(Terms, "specials")$trt
  if (length(trt.var) > 0)
    ct <- mf[, trt.var] else
      stop("upliftRF: formula does not include a control/treatment variable using trt().")    
  y <- model.response(mf, "numeric")
  var_names <- attributes(Terms)$term.labels[-(trt.var-1)]
  x <- model.frame(terms(reformulate(var_names)),
                    mf)                  
  res <- ccif(x = x, y = y, ct = ct, ...)
  cl <- match.call()
  cl[[1]] <- as.name("ccif")
  res$call <- cl
  return(res)
  
}  

### END FUN
