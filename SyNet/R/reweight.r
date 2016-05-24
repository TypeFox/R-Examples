reweight <- function(wm, similarity = TRUE, t1t2 = if(similarity) 
            c(0, 1) else c(Inf,0), normalized = TRUE){
  stopifnot(is.matrix(wm))
  stopifnot(isSymmetric(wm))
  stopifnot(any(!is.na(wm)))
  if(similarity) {
    noaffthr <- min(t1t2)
    closeaffthr <- max(t1t2)
    if(noaffthr < 0 | closeaffthr > 1) 
      stop("Affinity thresholds should be inside [0, 1] for similarity/sympatry matrix")
    fun <- ">=" 
    nocredit <- "<"
    mpond <- ifelse(wm >= closeaffthr, 1, wm) 
    diag(mpond) <- 1  
  }  else { 
    noaffthr <- max(t1t2)
    closeaffthr <- min(t1t2)
    fun <- "<="
    nocredit <- ">"
    mpond <- ifelse(wm <= closeaffthr, 0, wm)
    diag(mpond) <- 0  
  }
  rew <- matrix(0, nrow(wm), nrow(wm))
  for(i in 1:nrow(mpond))
    for (j in i:nrow(mpond)) {
       if(do.call(nocredit, list(mpond[i, j], noaffthr))) next
       aux <- do.call(fun, list(mpond[i,j], mpond[i,])) + do.call(fun, list(mpond[i,j], mpond[j,])) 
       rew[i, j] <- rew[j, i] <- sum(aux==2) + 1/(2 + sum(aux == 0))
  }
  if(normalized) return(rew/max(rew)) else return(rew)
}

