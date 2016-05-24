readtree.gtree <-
function(inFn, locus) {
  # Descr:  read gene tree    
  # Deps:   ?     
  # I/p:    inFn
  #         locus

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readtree.gtree", fg="red"), sep="")
  }

#################
# 1. Load trees #
#################
  treeD = readhelpers.parsetrees(inFn)                                  # Reading in the trees

  # DEBUGLINES:
  #cat("\ntreeD\n"); print(treeD)

##################################
# 2. Extract branch and metadata #
##################################
  treestrings = gsub("tree STATE_.*\\[&R\\] ", "", treeD$treestrings)   # Remove everything in lines except tree definition
  br_and_meta_D = lapply(treestrings, readhelpers.extract_br_and_meta, locus, sTreeFlg=F)  # Extracting branch- and metadata

  # DEBUGLINES:
  #cat("\nbr_and_meta_D\n"); print(br_and_meta_D)

################################################
# 3. Add branch- and metadata to phylo-objects #
################################################
  #LEGACY: tree = treeD$info
  outD = list()

  #LEGACY: if(class(tree)=="multiPhylo") {
  #LEGACY:   for(i in 1:length(tree)) {
    for (i in 1:length(treeD$full)) {
      #mode(br_and_meta_D[[i]]) = "numeric"
      #tree[[i]][["rate"]] = as.matrix(br_and_meta_D[[i]][,c(-1,-2)])    # Append branch rates to trees
      #outD[[i]] = list()
      #outD[[i]][["rate"]] = as.matrix(br_and_meta_D[[i]][,c(-1,-2)])   # Append branch rates to trees
      treeD_singleTree = treeD$full[[i]]                                # Additions possible only to phylo objects, not to multiphylo objects
      treeD_singleTree[["rate"]] = as.numeric(br_and_meta_D[[i]][,c(-1,-2)])   # Append branch rates to trees
      outD[[i]] = treeD_singleTree
    }
  #LEGACY: }
  #LEGACY: if(class(tree)=="phylo") {
  #LEGACY:   mode(br_and_meta_D[[1]]) = "numeric"
  #LEGACY:   tree[["rate"]] = as.numeric(br_and_meta_D[[1]][,c(-1,-2)])
  #LEGACY: }

  # DEBUGLINES:
  #cat("\noutD\n"); print(outD)

###########################
# 4. Summarize and return #
###########################
  if (verboseBool) {
    #LEGACY: cat("\tN of tree specs loaded: ", length(treeD$data), "\n", sep="")
    cat("\tN of tree specs loaded: ", length(outD), "\n", sep="")
  }
  return(outD)
}
