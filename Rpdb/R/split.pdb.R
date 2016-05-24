#  split and unsplit objects of class 'pdb'.

split.pdb <- function(x, f, drop = FALSE, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  atoms <- split(x$atoms, f, drop)

  to.return <- lapply(atoms,pdb, x$cryst1, x$conect, x$remark, x$title)
  to.return <- lapply(to.return,
                   function(x){
                     r <-     x$conect$eleid.1 %in% x$atoms$eleid
                     r <- r & x$conect$eleid.2 %in% x$atoms$eleid
                     if(any(r)) x$conect <- x$conect[r,]
                     else x["conect"] <- list(NULL)
                     return(x)
                   }
  )
  
  return(to.return)
}

unsplit.pdb <- function(value, f, drop = FALSE, ...)
{
  if(!all(unlist(lapply(value, is.pdb))))
    stop("'value' must be a list containing only 'pdb' objects")
  
  title  <- value[[1]]$title
  remark <- value[[1]]$remark
  cryst1 <- value[[1]]$cryst1
  
  atoms <- lapply(value, function(x) return(x$atoms))
  atoms <- unsplit(atoms, f, drop)
  
  conect <- lapply(value, function(x) return(x$conect))
  conect <- do.call(rbind, conect)
  rownames(conect) <- 1:nrow(conect)
  
  to.return <- pdb(atoms, cryst1, conect, remark, title)
  
  return(to.return)
}