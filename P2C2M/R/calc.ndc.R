calc.ndc <-
function(gTree, sTree, assoc) {
  # Descr:  returns the number of deep coalescences for an entire tree
  # Deps:   calc.parse
  #         calchelpers.gtreeparse
  # I/p:    sTree
  #         gTree
  #         assoc
  # Note:   ndc = "number of deep coalescences"

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.ndc",fg="red"),sep="")
  }

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD
  ndc = c()
  for (node in nodes) {
    tempData = calchelpers.gtreeparse(sTree, gTree, assoc, dmvD, node)
    ndc = c(ndc, tempData$lNd-1)
  }

  return(sum(ndc))
}
