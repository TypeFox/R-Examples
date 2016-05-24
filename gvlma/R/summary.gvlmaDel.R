"summary.gvlmaDel" <-
function(object, allstats = TRUE, ...)
{
  gvlmaDelobj <- object
  nobs <- nrow(gvlmaDelobj)
  oddcols <- seq(1, ncol(gvlmaDelobj), by = 2)
  Stats <- gvlmaDelobj[,oddcols]
  nstats <- ncol(Stats)
  pvals <- gvlmaDelobj[,-oddcols]
  AllQs <- lapply(as.list(gvlmaDelobj), summary)
  StatQs <- lapply(as.list(Stats), summary)
  pvalQs <- lapply(as.list(pvals), summary)
  AllQsMat <- matrix(unlist(AllQs), nrow = ncol(gvlmaDelobj), byrow= TRUE)
  StatQsMat <- matrix(unlist(StatQs), nrow = length(oddcols), byrow = TRUE)
  pvalQsMat <- matrix(unlist(pvalQs), nrow = length(oddcols), byrow = TRUE)
  rownames(AllQsMat) <- names(gvlmaDelobj)
  rownames(StatQsMat) <- names(StatQs)
  rownames(pvalQsMat) <- names(pvalQs)
  colnames(AllQsMat) <- names(AllQs[[1]])
  colnames(StatQsMat) <- names(StatQs[[1]])
  colnames(pvalQsMat) <- names(pvalQs[[1]])

  # unusual observations
  upperfencestat <- StatQsMat[,5] + 3*(StatQsMat[,5] - StatQsMat[,2])
  lowerfencestat <- StatQsMat[,2] - 3*(StatQsMat[,5] - StatQsMat[,2])
  upperfencepval <- pvalQsMat[,5] + 3*(pvalQsMat[,5] - pvalQsMat[,2])
  lowerfencepval <- pvalQsMat[,2] - 3*(pvalQsMat[,5] - pvalQsMat[,2])
  hiStats <- Stats > matrix(upperfencestat,
                            nrow = nobs, ncol = nstats, byrow = TRUE)
  loStats <- Stats < matrix(lowerfencestat,
                            nrow = nobs, ncol = nstats, byrow = TRUE)
  hipvals <- pvals > matrix(upperfencepval,
                            nrow = nobs, ncol = nstats, byrow = TRUE)
  lopvals <- pvals < matrix(lowerfencepval,
                            nrow = nobs, ncol = nstats, byrow = TRUE)
  unusualStats <- hiStats | loStats
  unusualpvals <- hipvals | lopvals
  unusualobs <- unusualStats | unusualpvals
######################################################################  
  cat("\nGlobal test deletion statistics.\n")
  cat("\nLinear Model:\n",
      deparse(attr(gvlmaDelobj, "lmcall")), "\n")
  cat("\ngvlma call:\n",
      deparse(attr(gvlmaDelobj, "gvlmacall")), "\n")
  cat("\n\nSummary values:\n")
  print(AllQsMat)
  if (allstats)
    statloop <- c("Global Stat", "Directional Stat1",
               "Directional Stat2", "Directional Stat3",
               "Directional Stat4")
  else statloop <- "Global Stat"
  statindx <- 0
  for (nm in statloop)
    {
      statindx <- statindx + 1
      unusualobsStat <- unusualobs[,statindx]
      unusualindx <- which(unusualobsStat)
      unusualobsprint <- data.frame(Stats[unusualindx, statindx],
                                    pvals[unusualindx, statindx],
                                    row.names = unusualindx)
      names(unusualobsprint) <- paste(c("Delta", ""),
                                      nm, c("(%)", "p-value"), sep = " ") 
      cat("\n\nUnusual observations for ", nm, ":\n", sep = "")
      print(unusualobsprint)
    }
  invisible(as.data.frame(unusualobs))
}

