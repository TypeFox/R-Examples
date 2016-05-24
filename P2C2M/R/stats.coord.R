stats.coord <-
function(ind, loci, alpha) {
  # Descr:  generate significance tables
  # Deps:   stats.outmaker
  # I/p:    ind
  #         loci
  #         alpha
  # Note:   CV = coefficient of variance
  #         avgs = arithmetic means
  #         qntls = quantiles

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  slctn = get("P2C2M_flg_dscrStat", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.coord", fg="red"), "\n",
        sep="")
  }

##############################
# 1. Setting number of tails #
##############################

  # T-tests for all descriptive statistics are two-tailed, because there is no
  # a priori reason in which direction they should differ.
  tailL = list()
  tailL$LCWT = "2"
  tailL$COAL = "2"
  tailL$NDC = "2"
  tailL$GSI = "2"
  ## T-tests for NDC should be left one-tailed, because trees not compliant 
  ## with the coalescent model have a higher number of deep coalescences.
  #NDCtail = "1l"
  ## T-tests for GSI should be right one-tailed, because trees not compliant 
  ## with the coalescent model have lower values.
  #GSItail = "1r"

#############################
# 2. Inferring significance #
#############################
  perGene = acrGenes = list()
  for (s in slctn) {
    perGene[[s]] = stats.perGene(ind[[s]]$dif, alpha, tailL[[s]])
    acrGenes[[s]] = stats.acrGenes(ind[[s]]$dif, alpha, tailL[[s]])
  }

#######################
# 3. Combining output #
#######################
  outList = list()
  # perGene output
  outList$perGene = sapply(perGene, cbind)
  rownames(outList$perGene) = c(loci)
  # acrossGene output
  outList$acrGenes = sapply(acrGenes, cbind)
  rownames(outList$acrGenes) = c("Sum", "Mean", "Median", "Mode", "CV")
  # Naming rows of output
  names = c()
  for (stat in slctn) {
    names = c(names, paste(stat, "[", tailL[[stat]], "]", sep=""))
  }
  colnames(outList$perGene) = colnames(outList$acrGenes) = names
  
  return(outList)
}
