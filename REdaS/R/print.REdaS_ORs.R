print.REdaS_ORs <- function(x, ...){
  
  if(class(x) != "REdaS_ORs") stop('"x" must be an object of class "REdaS_ORs"', call. = FALSE)
  
  if(all(is.na(unlist(x$ORs)))) stop('no results due to frequencies of 0 in all tables.\n\tsee ?odds_ratios', call. = FALSE)
  
  res <- matrix(unlist(x$ORs), ncol=5L, byrow=TRUE)
  
  dimnames(res) <- list(
    unlist(lapply(x$comps, function(txt){ paste0(txt[1], " vs. ", txt[2]) })),
    c("OR", "log(OR)", "Std. Err.", "z-value", "p-value")
  )
  
  if(interactive()) writeLines("")
  cat("Odds Ratios\n")
  printCoefmat(res, P.values = TRUE, has.Pvalue = TRUE)
  if(interactive()) writeLines("")
  invisible(res)
  
}
