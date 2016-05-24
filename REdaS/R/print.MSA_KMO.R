print.MSA_KMO <- function(x, stats = c("both", "MSA", "KMO"), vars = "all", sort = FALSE, show = "all", digits = getOption("digits"), ...){
  stats <- match.arg(stats)
  pMSA <- ifelse(stats %in% c("both", "MSA"), TRUE, FALSE)
  pKMO <- ifelse(stats %in% c("both", "KMO"), TRUE, FALSE)
  
  if(is.character(vars) && length(vars) == 1L && vars == "all"){ var_idx <- seq_along(x$MSA) } else {
    if(!is.numeric(vars)){ stop('"vars" must be "all" or a numeric vector of indices for the variables.') }
    if(!all(vars %in% seq_along(x$MSA))){ stop('"vars" contains invalid indices.') }
    var_idx <- vars
  }
  MSAres <- round(x$MSA[var_idx], digits)
  if(sort){
    MSAres <- sort(MSAres, decreasing = FALSE)
  }
  if(!(is.character(show) && length(show) == 1L && show == "all")) MSAres <- MSAres[show]
  
  if(pMSA){ writeLines("\nKaiser-Meyer-Olkin Statistics\n") } else
          { writeLines("\nKaiser-Meyer-Olkin Statistic") }
  writeLines(paste0("Call: ", deparse(x$call)))

  if(pMSA){
    writeLines("\nMeasures of Sampling Adequacy (MSA):")
    print.default(MSAres)
  }
  if(pKMO){
    writeLines(paste0("\nKMO-Criterion: ", round(x$KMO, digits)))
  }
  writeLines("")
  invisible(list(MSA=MSAres, KMO=round(x$KMO, digits)))
}
