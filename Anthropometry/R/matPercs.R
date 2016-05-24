matPercs <- function(archoids,data){
  percs <- list()
  for(j in 1 : length(archoids)){
    percs[[j]] <- sapply(1 : dim(data)[2], percentilsArchetypoid,
                         archoids[j], data, 0)
  }
  mat <- matrix(unlist(percs), nrow = ncol(data),
                  ncol = length(percs), byrow = FALSE)
  return(mat)
}