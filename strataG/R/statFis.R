#' @rdname popStructStat
#' @export
#' 
statFis <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  if(ploidy(g) < 2 | nStrata(g) == 1) {
    return(list(
      stat.name = "Fis", 
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statFis_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g)
  )
  
  if(length(result) == 1) keep.null <- FALSE
  p.val <- if(length(result) == 1) NA else mean(result >= result[1], na.rm = TRUE) 
  return(list(
    stat.name = "Fis", 
    result = c(estimate = result[1], p.val = p.val),
    null.dist = if(keep.null) result[-1] else NULL
  ))
}