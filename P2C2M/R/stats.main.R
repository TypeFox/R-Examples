stats.main <-
function (path, xml.file, loci, resultData, prmFile) {
  # Descr:    coordinates executing of modules
  # Deps:     (various)
  # I/p:      path = absolute path to Working directory
  #           xml.file = name of infile
  #           loci = names of loci
  #           resultData
  #           prmFile

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> stats.main", fg="red"), "\n",
        sep="")
  }

##########################
# 1. Summarizing results #
##########################
  loghelpers.prntmngr("Summarizing results", uprFlg=T)
                                               
  # Setting outdata
  outData = list()
  # Setting alpha values
  alphaValues = c(0.1, 0.05, 0.01)
  for (val in alphaValues) {
    # Coordinating the calculation of the statistics
    results = stats.coord(resultData, loci, val)
    valStr = paste("alpha", as.character(val), sep="")
    outData[[valStr]] = results
  }

#####################
# 2. Writing legend #
#####################
  legend = "Differences between the posterior and the posterior predictive distributions. Each cell contains the following information in said order: mean, standard deviation, significance level. Codes in square brackets indicate the number of tails. Alpha values are automatically adjusted for the number of tails."
  outData$legend = legend
return(outData)
}
