calc.parse <-
function(sTree, assoc) {
  # Descr:  parses species tree nodes for metric calculation
  # Deps:   calchelpers.dmvparse
  # I/p:    sTree
  #         assoc

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.parse",fg="red"),sep="")
  }

  # DEBUGLINES:
  #cat("\nassoc\n"); print(assoc)
  #cat("\nsTree$tip.label\n"); print(sTree$tip.label)
  
  spNames = sTree$tip.label
  n_sp = length(spNames)
  #nBr = (2*n_sp)-1
  tiplist = list()
  for(i in 1:n_sp) {
    tiplist[[spNames[i]]] = assoc[which(assoc[,1]==spNames[i]),2]
  }
  dmvD = calchelpers.dmvparse(sTree, n_sp)                              # returns demographic info of the species tree
  n_tips_per_sp = lapply(tiplist, length)                               # calculate number of tips per species

  if(any(n_tips_per_sp==1)) {                                           # evaluate number of tips per species tree
    tmp = dmvD[-which(n_tips_per_sp==1),]                               # Remove dmv values for terminals (i.e., where n_tips_per_sp==1)
  } else {
    tmp = dmvD
  }
  nodes = tmp[,"node"]
  names(nodes) = NULL

  outD = list()
  outD$nodes = nodes
  outD$dmvD = dmvD

  return(outD)
}
