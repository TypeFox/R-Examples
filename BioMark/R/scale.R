scalefun <- function(sc.p = c("none", "log", "sqrt", "pareto", "auto"))
{
  if (is.null(sc.p)) sc.p <- "none"
  sc.p <- match.arg(sc.p)
  
  function(X) {
    switch(sc.p,
           "none" = return(scale(X, scale = FALSE)),
           "log" = return(scale(log(X), scale = FALSE)),
           "sqrt" = return(scale(sqrt(X), scale = FALSE)),
           "pareto" = return(scale(X, scale = sqrt(apply(X, 2, sd)))),
           "auto" = return(scale(X, scale = TRUE)))
  }
}
