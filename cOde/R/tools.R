#' reduceSensitivities
#' @param sens Named character, the sensitivity equations
#' @param vanishing Character, names of the vanishing sensitivities
#' @details Given the set \code{vanishing} of vanishing sensitivities, the algorithm
#' determins sensitivities that vanish as a consequence of the first set.
#' @return Named character, the sensitivity equations with zero entries for vanishing
#' sensitivities.
#' @export
reduceSensitivities <- function(sens, vanishing) {

  sensvar <- names(sens)
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) v[2]))
  ini.zero <- which(senssplit.1 != senssplit.2)
  ini.nonzero <- which(senssplit.1 == senssplit.2) 
  
  exit <- FALSE
  sensvar.zero <- sensvar[ini.zero]
  sensvar.nonzero <- sensvar[ini.nonzero]
  while(!exit){
    find_nonzero <- unlist(lapply(sensvar.zero, function(s){
      allSyb <- union(s,sensvar[sensvar %in%getSymbols(sens[s])])
      nDpf <- setdiff(allSyb,vanishing)
      nIni <- intersect(sensvar.nonzero,allSyb)
      nonzero <- (length(nDpf)+length(nIni) > 0)
    }))
    if(any(find_nonzero)){
      sensvar.nonzero <- c(sensvar.nonzero,sensvar.zero[find_nonzero])
      sensvar.zero <- sensvar.zero[!find_nonzero]
    }else{
      exit <- TRUE
    }
  }
  sens[sensvar.zero] <- "0"
  sens <- replaceSymbols(sensvar.zero, "0", sens)
  attr(sens,"is.zero") <- sensvar.zero
  #return(list(sens, sensvar.zero,sensvar.nonzero))
  return(sens)
}
