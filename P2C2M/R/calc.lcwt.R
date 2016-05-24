calc.lcwt <-
function(gTree, sTree, assoc, ploidy) {
  # Descr:    calculates the probability of a whole gene tree
  #           by calculating the probability of subtrees defined 
  #           by their MRCA (i.e. their node); done for all 
  #           nodes in a tree, values are summed up thereafter
  # Deps:     calc.parse
  #           calchelpers.brprob
  # I/p:      sTree
  #           gTree
  #           assoc
  #           ploidy
  # Note:     gtp = "gene tree probability" (a form of coalescent likelihood)

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> calc.lcwt",fg="red"),sep="")
  }

  handle = calc.parse(sTree, assoc)
  nodes = handle$nodes
  dmvD = handle$dmvD
  
  # DEBUGLINES:
  #cat("\nnodes\n"); print(nodes)
  #cat("\ndmvD\n"); print(dmvD)

  lnP = c()
  for(node in nodes) {
    lnP = c(lnP, log(calchelpers.brprob(sTree, gTree, assoc, 
                                        ploidy, dmvD, node)))
  }

  return(sum(lnP))
}
