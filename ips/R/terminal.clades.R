## package: ips
## called by: phylo2mafft, {megaptera}stepF
## author: Christoph Heibl (at gmx.net)
## last update: 2014-07-30

terminal.clades <- function(phy){
  obj <- lapply(1:Ntip(phy), sister, phy = phy)
  for ( i in seq_along(obj)) 
    obj[[i]] <- sort(c(obj[[i]], i))
  obj <- unique(obj)
  is.nested <- function(x, y){
    identical(sort(union(x, y)), x)
  }
  id <- vector(length = length(obj))
  for ( i in seq_along(obj) ){
    id[i] <- !any(sapply(obj[-i], is.nested, x = obj[[i]]))
  }
  obj[id]
#   obj <- unlist(obj[sapply(obj, length) == 1])
#   matrix(obj, ncol = 2, byrow = TRUE)
}