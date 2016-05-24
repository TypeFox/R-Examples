#@ei.object -- ei object
#@... options
#Extracts information from ei object

eiread <- function(ei.object, ...){
  function.list <- list("betab" = .betaB, "betaw" = .betaW,
                        "phi" = .phi, "sbetab" = .sbetab,
                        "sbetaw" = .sbetaw, "psisims" = .psisims,
                        "bounds" = .bounds, "CI80b" =.CI80b,
                        "CI80w" = .CI80w, "abounds" = .abounds,
                        "aggs" = .aggs, "maggs" = .maggs,
                        "VCaggs" = .VCaggs, "eaggbias" = .eaggbias,
                        "goodman" = .goodman)
  dec <- ei.object$precision
  arguments <- list(...)
  results <- list()
  for (arg in arguments) {
    if (arg %in% names(function.list))
      results[[arg]] <-
        floor(function.list[[arg]](ei.object)*10^dec)/10^dec
    else
      results[[arg]] <- NA
  }
  if (length(results)==1)
    results <- results[[1]]
  if (length(results) <1)
    warning("qi results object is empty")
  #class(results) <- "eiread"
  results
}
