calchelpers.dmvparse <-
function(sTree, nSp) {
  # Descr:    parsing demographic data (dmv and dmt) of a branch
  # Deps:     -
  # I/p:      sTree
  #           nSp

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.dmvparse", fg="red"),
        sep="")
  }
  
    # DEBUGLINES:
    #cat("\nsTree\n"); print(sTree)
    #cat("\nsTree$edge\n"); print(sTree$edge)
    #cat("\nsTree$edge.length\n"); print(sTree$edge.length)
    #cat("\nsTree$dmv\n"); print(sTree$dmv)

    dmvD = cbind(sTree$edge[,2], sTree$edge.length)                     # generating rows of node number - branch length pairs
    dmvD = rbind(dmvD, c((nSp+1), Inf))                                 # adding another branch of length Inf
    dmvD = cbind(dmvD, sTree$dmv)                                       # adding dmv values
    dmvD = dmvD[order(dmvD[,1]),]                                       # order the matrix by the first column
    stBt = ape::branching.times(sTree)                                  # calc branching times of the species tree (via ape-function)

    # TFL may not be necessary, as stBt may already be sorted
    stBt = stBt[order(as.numeric(names(stBt)))]                         # sort the branching times by their names (which are numbers)

    pre = structure(rep(0, nSp), .Names=c(1:nSp))                       # add x zeros to the beginning of branching times list, where x = number of species in sTree
    stBt = c(pre, stBt)

    dmvD = cbind(dmvD, stBt)                                            # add the column "stBt" to the matrix "dmvD"
    colnames(dmvD) = c("node", "length", "dmv", "sbt")
    rownames(dmvD) = c(1:length(stBt))

return(dmvD)
}
