
design <- function (formula, data)
{
  call <- match.call()

  
  if (missing(data)) 
    data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data"), names(mf), 0L) 
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  
  mf <- eval(mf, parent.frame())   
  mt <- attr(mf, "terms")

  
  Y <- model.response(mf, "any")
     if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) 
      names(Y) <- nm
    } 
  
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  
  return(X)
}

