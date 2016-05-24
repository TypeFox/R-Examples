readhelpers.parsetrees <-
function(inFn) {
  # Descr:  reading in trees
  # Deps: -
  # I/p:  inFn

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n",xtermStyle::style("DEBUG> readhelpers.parsetrees",fg="red"),sep="")
  }

  outD = list()
  tree = ape::read.nexus(inFn)
  outD[["full"]] = tree                                                 # Pass on the full tree data

  lines = scan(file=inFn, what="", sep="\n", quiet=T)                   # Load inFn line by line
  outD[["treestrings"]] = lines[grep("tree STATE", lines)]              # Extract all those lines that contain keyword "tree STATE"

  return(outD)
}
