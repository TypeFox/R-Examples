stats.acrGenes <-
function (diff, alpha, tail) {
  # Descr:    calculates statistics across genes
  # Deps:     statshelpers.qntls
  #           statshelpers.cv
  # I/p:      diff
  #           alpha
  #           tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.acrGenes", fg="red"),
        sep="")
  }

  # Perform calculations per row (i.e per MCMC generation)
  acrossG_Sum = rowSums(diff)
  acrossG_Mean = rowMeans(diff)
  acrossG_Median = rowMedians(diff)
  acrossG_Mode = rowModes(diff)
  acrossG_CV = statshelpers.cv(diff)

  # "acrossG" is a matrix with four columns hereafter
  acrossG = cbind(acrossG_Sum, 
                  acrossG_Mean, 
                  acrossG_Median,
                  acrossG_Mode,
                  acrossG_CV)

  qntls = statshelpers.qntls(acrossG, alpha, tail)
  sigSgns = statshelpers.sigsgn(qntls, tail)

  subCol1 = c(apply(acrossG, MARGIN=2, mean))
  subCol2 = c(apply(acrossG, MARGIN=2, sd))
  subCol3 = c(sigSgns)
  outData = paste(round(subCol1, 2), 
               paste("(", "\u00B1", round(subCol2, 2), ")", sep=""), 
               subCol3)

  return(outData)
}
