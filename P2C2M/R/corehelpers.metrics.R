corehelpers.metrics <-
function(gTree, pTree, pTreeNames, sTree, 
                               assoc, ploidy, descrStats, singleAllele) {
  # Descr:    corehelpers.metrics
  # Deps:     calc.lcwt
  #           calc.ndc
  #           calc.coal
  # I/p:      gTree
  #           pTree
  #           pTreeNames
  #           sTree
  #           assoc
  #           ploidy
  #           descrStats
  #           singleAllele

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> corehelpers.metrics", fg="red"), sep="")
  }

  outL = list()

  if ("GSI" %in% descrStats) {
    gsi = calc.gsi(gTree, assoc, singleAllele)
    outL$GSI = frmtMntse(gsi, 4)                                        # Setting a specific number of significands
  }

  if ("LCWT" %in% descrStats) {
    gtp = calc.lcwt(gTree, sTree, assoc, ploidy)
    outL$LCWT = frmtMntse(gtp, 4)
  }

  if ("NDC" %in% descrStats) {
    ndc = calc.ndc(gTree, sTree, assoc)
    outL$NDC = frmtMntse(ndc, 4)
  }

  if ("COAL" %in% descrStats) {
      
    # DEBUGLINES:
    #cat("\ngTree\n"); print(gTree)
    #cat("\npTree\n"); print(pTree)
    #cat("\npTreeNames\n"); print(pTreeNames)
    #cat("\nassoc\n"); print(assoc)
    
    ray = calc.coal(gTree, pTree, pTreeNames, assoc)
    outL$COAL = frmtMntse(ray, 4)
  }

  outD = unlist(outL)
  names(outD) = NULL
  outD = toString(outD)

  return(outD)
}
