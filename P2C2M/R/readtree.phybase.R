readtree.phybase <-
function(inFn) {
  # Descr:  read tree of format phybase
  # Deps:   (none)
  # I/p:    inFn

  beastVers = get("P2C2M_flg_beastV", envir=P2C2M_globalVars)
  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)

  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readtree.phybase", fg="red"), 
        sep="")
  }

  pySc = system.file("exec", paste("BEAST2phybase_", beastVers, ".py", sep=""), package="P2C2M")  # Specify name of parsing script
  system(paste("python2", pySc, inFn))
  pTrees = phybase::read.tree.string(paste(rmext(inFn), ".P2C2M.phyb", sep=""))  # Read tree via PHYBASE
  pTrees$tree = gsub("^ .", "", pTrees$tree)                            # Remove leading white space from tree specs

  return(pTrees)
}
