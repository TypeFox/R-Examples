statshelpers.qntls <-
function (ind, alpha, tail) {
  # Descr:  generates a distribution of quantiles;
  #         quantile levels are hardcoded in variable "qntlLevels"
  # Deps:   (none)
  #         ind = a data frame
  #         alpha
  #         tail

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.qntls",
        fg="red"), sep="")
  }

  #################################
  # Set quantile levels via alpha #
  #################################
  if (tail=="1l" | tail=="1r") {
    qntlLevels = c(alpha, 1-alpha)
  }
  # Adjusting alpha value for two-tailed test
  if (tail=="2") {
    qntlLevels = c((alpha/2), 1-(alpha/2))
  }

  qntls = t(apply(ind, MARGIN=2, quantile, qntlLevels, na.rm=T))
  #qntls = rbind(qntls, quantile(acrossGenes, qntlLevels, na.rm=T))

  return(qntls)
}
