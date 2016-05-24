#  Merge by row two data frames of class 'coords'.

merge.coords <- function(x, y, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.coords(y)) stop("'y' must be an object of class 'coords'")
  
  if(basis(x) != basis(y)) stop("'x' and 'y' basis differ")
  
  to.return <- rbind(x,y)
  return(to.return)
}

merge.atoms <- function(x, y, reindex = TRUE, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  if(!is.atoms(y)) stop("'y' must be an object of class 'atoms'")
  
  if(basis(x) != basis(y)) stop("'x' and 'y' basis differ")
  
  y.resid <- as.integer(y$resid + max(x$resid))
  y.eleid <- as.integer(y$eleid + max(x$eleid))
  
  y$resid <- y.resid
  y$eleid <- y.eleid
  
  to.return <- rbind(x,y)
  if(reindex) to.return <- reindex.atoms(to.return)
  
  return(to.return)
}

# Merge two objects of class 'pdb'.

merge.pdb <- function(x, y, reindex = TRUE, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  if(!is.pdb(y)) stop("'y' must be an object of class 'pdb'")
  
  if(basis(x) != basis(y)) stop("'x' and 'y' basis differ")
  
  if(any(x$cryst1$abc != y$cryst1$abc) | any(x$cryst1$abg != y$cryst1$abg) | any(x$cryst1$sgroup != y$cryst1$asgroup))
    warning("Different periodical boundary conditions are defined for 'x' and 'y'. 'x$cryst1' has been kept.")
  
  title  <- unique(c(x$title , y$title ))
  remark <- unique(c(x$remark, y$remark))
  cryst1 <- x$cryst1
  
  y$conect$eleid.1 <- y$conect$eleid.1 + max(x$atoms$eleid)
  y$conect$eleid.2 <- y$conect$eleid.2 + max(x$atoms$eleid)

  eleid.1 <- c(x$conect$eleid.1, y$conect$eleid.1)
  eleid.2 <- c(x$conect$eleid.2, y$conect$eleid.2)
  conect      <- conect(eleid.1, eleid.2)

  atoms <- merge.atoms(x$atoms, y$atoms, reindex = FALSE)
  
  to.return   <- pdb(atoms, cryst1, conect, remark, title)
  if(reindex) to.return <- reindex.pdb(to.return)
  
  return(to.return)
}

