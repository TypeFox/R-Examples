#  Return subsets of 'pdb' or 'atoms' objects which meet conditions.

subset.atoms <- function(x, subset, drop = FALSE, reindex.all = TRUE, ...)
{
  if (missing(subset)) 
    r <- TRUE
  else
  {
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must evaluate to logical")
    r <- r & !is.na(r)
  }
  x <- x[r, , drop = drop]
  
  if(nrow(x) == 0) x <- NULL
  
  if(reindex.all) x <- reindex.atoms(x)
  return(x)
}

subset.pdb <- function(x, subset, drop = FALSE, reindex.all = TRUE, ...)
{
  if (missing(subset)) 
    r <- TRUE
  else
  {
    e <- substitute(subset)
    r <- eval(e, x$atoms, parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must evaluate to logical")
    r <- r & !is.na(r)
  }
  if(any(r)) x$atoms <- x$atoms[r, , drop = drop]
  else x["atoms"] <- list(NULL)
  r <-     x$conect$eleid.1 %in% x$atoms$eleid
  r <- r & x$conect$eleid.2 %in% x$atoms$eleid
  if(any(r)) x$conect <- x$conect[r,]
  else x["conect"] <- list(NULL)
  
  if(reindex.all) x <- reindex.pdb(x)
  
  return(x)
}
