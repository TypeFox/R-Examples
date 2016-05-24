#  Reinitialize the indexation of residues and elements IDs.

reindex <- function(...)
  UseMethod("reindex")

reindex.atoms <- function(x, eleid = TRUE, resid = TRUE, ...)
{
  if(eleid)
  {
    eleid <- as.factor(x$eleid)
    levels(eleid) <- 1:nlevels(eleid)
    x$eleid <- as.integer(as.character(eleid))
  }
  if(resid)
  {
    resid <- as.factor(x$resid)
    levels(resid) <- 1:nlevels(resid)
    x$resid <- as.integer(as.character(resid))
  }
  rownames(x) <- 1:nrow(x)
  return(x)
}

reindex.pdb <- function(x, eleid = TRUE, resid = TRUE, ...)
{
  if(eleid)
  {
    if(anyDuplicated(x$atoms$eleid)){
      x$atoms$eleid <- 1:natom(x)
      if(!is.null(x$conect)){
        cat("Recalculating connectivity\n")
        x$conect <- conect(x)      
      }
    }
    else{
      if(!is.null(x$conect)){
        eleid.1 <- factor(x$conect$eleid.1, levels = x$atoms$eleid)
        eleid.2 <- factor(x$conect$eleid.2, levels = x$atoms$eleid)
        levels(eleid.1) <- 1:nlevels(eleid.1)
        levels(eleid.2) <- 1:nlevels(eleid.2)
        eleid.1 <- as.integer(as.character(eleid.1))
        eleid.2 <- as.integer(as.character(eleid.2))
        x$conect <- conect(eleid.1, eleid.2)
      }
      eleid <- factor(x$atoms$eleid, levels = x$atoms$eleid)
      neleid <- nlevels(eleid)
      levels(eleid) <- 1:neleid
      x$atoms$eleid <- as.integer(as.character(eleid))
    }
  }
  if(resid)
  {
    resid <- as.factor(x$atoms$resid)
    levels(resid) <- 1:nlevels(resid)
    x$atoms$resid <- as.integer(as.character(resid))
  }
  rownames(x$atoms) <- 1:nrow(x$atoms)
  return(x)
}
