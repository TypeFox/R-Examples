stats.perGene <-
function (diff, alpha, tail) {
  # Descr:    calculates statistics per gene
  # Deps:     statshelpers.qntls
  # I/p:      diff
  #           alpha
  #           tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.perGene", fg="red"), 
        sep="")
  }

  # Structure of "diff" at this point:
  #        [gene1]   [gene2]  [gene3]
  # [1,]  -2.72395  -9.48870 45.79565
  # [2,]  14.97560  38.06380 88.60285
  # [3,]   8.76280  21.09895 50.51950
  qntls = statshelpers.qntls(diff, alpha, tail)
  sigSgns = statshelpers.sigsgn(qntls, tail)
  # Calculation of means
  # Note: "MARGIN=2" means "across a column"
  subCol1 = c(apply(diff, MARGIN=2, mean))
  # Calculation of stdv
  subCol2 = c(apply(diff, MARGIN=2, sd))
  # Assignment of signif. values
  subCol3 = c(sigSgns)
  # "\u00B1" is unicode sign of "plusminus"
  # Output format:   mean(+-1 SD)[:space:]signif.level
  #                  Mean and SD are rounded to two decimal places.
  outData = paste(round(subCol1, 2), 
               paste("(", "\u00B1", round(subCol2, 2), ")", sep=""), 
               subCol3)

  return(outData)
}
