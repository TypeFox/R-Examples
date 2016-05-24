"chanames" <-
function (x, allow.null = TRUE) 
{
  if (is.mcmc.list(x)) {
    if (is.null(names(x))) 
      if (allow.null) 
        NULL
      else paste("chain", 1:length(x), sep = "")
    else names(x)
  }
  else NULL
}

"chanames<-" <-
function (x, value) 
{
  if (is.mcmc.list(x)) 
      names(x) <- value
    else stop("Not an mcmc.list object")
    x
}

"varnames" <-
function (x, allow.null = TRUE) 
{
  if (!is.mcmc(x) && !is.mcmc.list(x)) 
    return(NULL)
  y <- if (is.mcmc(x)) 
    dimnames(x)[[2]]
  else if (is.mcmc.list(x)) 
    dimnames(x[[1]])[[2]]
  if (is.null(y) && !allow.null) 
    y <- paste("var", 1:nvar(x), sep = "")
  return(y)
}

"varnames<-" <-
function (x, value) 
{
    if (is.mcmc(x)) {
        if (length(dim(x)) < 2) {
            dim(x) <- c(length(x), 1)
        }
        colnames(x) <- value
    }
    else if (is.mcmc.list(x)) {
        for (i in 1:nchain(x)) varnames(x[[i]]) <- value
    }
    else stop("Not an mcmc or mcmc.list object")
    x
}

"nchain" <-
function (x) 
{
    if (is.mcmc(x)) 
        1
    else if (is.mcmc.list(x)) 
        length(x)
    else NULL
}

"nvar" <-
function (x) 
{
  
  if (is.mcmc(x)) {
    if (is.matrix(x)) ncol(x) else 1
  }
  else if (is.mcmc.list(x)) {
    if (is.matrix(x[[1]])) ncol(x[[1]]) else 1
  }
  else NULL
}

"niter" <-
function (x) 
{
  if (is.mcmc(x)) {
    if (is.matrix(x)) nrow(x) else length(x)
  }
  else if (is.mcmc.list(x)) {
    if (is.matrix(x[[1]])) nrow(x[[1]]) else length(x[[1]])
  }
  else NULL
}












