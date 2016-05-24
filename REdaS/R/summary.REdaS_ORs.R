summary.REdaS_ORs <- function(object, ...){
  
  if(class(object) != "REdaS_ORs") stop('"object" must be an object of class "REdaS_ORs"', call.=FALSE)
  
  if(all(is.na(unlist(object$ORs)))) stop('no results due to frequencies of 0 in all tables.\n\tsee ?odds_ratios', call.=FALSE)

  n_tabs <- length(object$tables)
  
  if(interactive()) writeLines("")
  writeLines(paste0(ifelse(n_tabs == 1L, "", "\tPairwise "), "Odds Ratios"))
  
  for(i in seq_len(n_tabs)){
    cat("\n")
    print(object$tables[[i]])
    writeLines(paste0("\nOdds-Ratio = ", round(object$ORs[[i]]$or, 3)))
    writeLines(paste0("Log(Odds-Ratio) = ", round(object$ORs[[i]]$lor, 3), ", Standard Error = ", round(object$ORs[[i]]$se, 3)))
    pv <- my.format.pval(object$ORs[[i]]$p, 5L)
    writeLines(paste0("z-value = ", round(object$ORs[[i]]$z, 3), ", p-value", ifelse(pure_all.equal(0, round(object$ORs[[i]]$p, 5L)), " ", " = "), pv))
    if(i < n_tabs) cat("\n")
  }
  if(interactive()) writeLines("")
  
}
