no.intercept <- function(mm) 
{
  if(attributes(mm)$assign[1] == 0) {
    remove.int <- which(attributes(mm)$assign == 0)
    mm <- mm[, -remove.int, drop = FALSE]
    attributes(mm) <- list(dim = dim(mm), dimnames = dimnames(mm), 
                           assign = attributes(mm)$assign[-remove.int])
  } else return(mm)
  mm
}