readhelpers.makeCFtable <-
function(ind) {
  # Descr:  generating a constituent-frequency table
  # Deps: -
  # I/p:  ind
  #
  # The command will perform the following conversion:
  # FROM:
  #         Var1            Freq
  #    [1,] "alpinus"        "4"
  # TO:
  #        spec              V2 
  #   Var1 "alpinus"         "1"
  #   Var1 "alpinus"         "2"
  #   Var1 "alpinus"         "3"
  #   Var1 "alpinus"         "4"

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readhelpers.makeCFtable", fg="red"), 
        sep="")
  }

  ind = as.matrix(as.data.frame(table(ind[,1])))
  aList = c()
  for (i in 1:length(ind[,1])) {
    aList = c(aList, rep(ind[i,1], times=ind[i,2]))
  }
  outd = as.matrix(as.data.frame(cbind(aList, 1:length(aList))))

  return(outd)
}
