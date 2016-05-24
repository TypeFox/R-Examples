# Automatically generated from all.nw using noweb

as.data.frame.pedigree <- function(x, ...) {

  dadid <- momid <- rep(0, length(x$id))
  dadid[x$findex>0] <- x$id[x$findex]
  momid[x$mindex>0] <- x$id[x$mindex]
  df <- data.frame(id=x$id, dadid=dadid, momid=momid, sex=x$sex)
  
  if(!is.null(x$affected))
    df$affected = x$affected
  
  if(!is.null(x$status))
    df$status = x$status
  return(df)
}
