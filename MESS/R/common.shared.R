common.shared <- function (id, ...) 
{
  UseMethod("common.shared")
}


common.shared.pedigreeList <- function (id, ...) 
{
#    if (!require("Matrix")) { 
#        stop("suggested package Matrix not installed") 
#    } 
    
    famlist <- unique(id$famid)
    nfam <- length(famlist)
    matlist <- vector("list", nfam)
    for (i in 1:length(famlist)) {
      family <- id[i]
      temp <- common.shared(family)
      matlist[[i]] <- as(forceSymmetric(temp), "dsCMatrix")
    }
    result <- bdiag(matlist)
    
    if (any(duplicated(id$id))) {
      dimnames(result) <- list(NULL, paste(id$famid, id$id, 
                                           sep = "/"))
    } else { dimnames(result) <- list(id$id, id$id) }
    result
}

common.shared.pedigree <- function (id, ...) 
{
  n <- length(id$id)

  if (n == 1) 
    return(matrix(1, 1, 1, dimnames = list(id$id, id$id)))
  if (any(duplicated(id$id))) 
    stop("All id values must be unique")
  
  temp <- matrix(rep(1, n*n), n)
  dimnames(temp) <- list(id$id, id$id)
  temp
}
