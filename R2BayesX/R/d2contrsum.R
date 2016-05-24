d2contrsum <- function(x)
{
  if(is.null(x) || !is.data.frame(x))
    stop("x must be a data.frame!")
  nx <- names(x)
  for(i in 1L:length(x)) {
    if(is.factor(x[[i]])) { 
      var <- x[[i]]
      var <- C(var, contr.sum)
      contr <- attr(var, "contrasts")
      colnames(contr) <- rownames(contr)[1:(nrow(contr) - 1)]
      attr(var, "contrasts") <- contr
      x[[i]] <- var
    }
  }

  return(x)
}
