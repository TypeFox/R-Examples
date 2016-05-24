readtree.stree <-
function(inFn) {
  # Descr:  read species tree
  # Deps:   ?
  # I/p:    inFn
  # Note:   little change compared to Noah's scripts

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readtree.stree", fg="red"), sep="")
  }

#################
# 1. Load trees #
#################
  treeD = readhelpers.parsetrees(inFn)

  # DEBUGLINES:
  #cat("\ntreeD\n"); print(treeD)

##################################
# 2. Extract branch and metadata #
##################################
  treestrings = gsub("tree STATE_[[:digit:]]+[[:space:]]\\[[^]]*\\] = \\[&R\\] ", "", treeD$treestrings)  # extract only the tree specification
  br_and_meta_D = lapply(treestrings, readhelpers.extract_br_and_meta, locus=NA, sTreeFlg=T)

  # DEBUGLINES:
  #cat("\nbr_and_meta_D\n"); print(br_and_meta_D)

################################################
# 3. Add branch- and metadata to phylo-objects #
################################################
  #LEGACY: tree = treeD$info
  outD = list()

  #LEGACY: # Notes by Noah: Appending branch width stats to trees; vector "check" is 
  #LEGACY: # included as a verification and should exactly match vector "edge.length"
  if(ncol(br_and_meta_D[[1]])==6) {
    #LEGACY: if(class(tree)=="multiPhylo") {
      #LEGACY: for(i in 1:length(tree)) {
      for(i in 1:length(treeD$full)) {
        #LEGACY: #tree[[i]][["check"]]<-as.numeric(br_and_meta_D[[i]][,2])
        #outD[[i]] = list()
        treeD_singleTree = treeD$full[[i]]                              # Additions possible only to phylo objects, not to multiphylo objects
        #treeD_singleTree[["rate"]] = as.numeric(br_and_meta_D[[i]][,c(-1,-2)])  # Append branch rates to trees
        treeD_singleTree[["dmt"]] = as.numeric(as.matrix(br_and_meta_D[[i]][,3]))
        treeD_singleTree[["dmv_start"]] = as.numeric(br_and_meta_D[[i]][,4])
        treeD_singleTree[["dmv_end"]] = as.numeric(br_and_meta_D[[i]][,5])
        treeD_singleTree[["dmv"]] = as.numeric(br_and_meta_D[[i]][,6])
        outD[[i]] = treeD_singleTree
      #LEGACY: }
    }
    #LEGACY: if(class(tree)=="phylo") {
      #LEGACY: #	tree[["check"]]<-as.numeric(br_and_meta_D[[1]][,2])
      #LEGACY: tree[["dmt"]] = as.numeric(br_and_meta_D[[1]][,3])
      #LEGACY: tree[["dmv_start"]] = as.numeric(br_and_meta_D[[1]][,4])
      #LEGACY: tree[["dmv_end"]] = as.numeric(br_and_meta_D[[1]][,5])
      #LEGACY: tree[[i]][["dmv"]] = as.numeric(br_and_meta_D[[i]][,6])
    #LEGACY: }
  }

  if(ncol(br_and_meta_D[[1]])==3) {
    #LEGACY: if(class(tree)=="multiPhylo") {
      for(i in 1:length(treeD$full)) {
        #LEGACY: # tree[[i]][["check"]]<-as.numeric(br_and_meta_D[[i]][,2])
        treeD_singleTree = treeD$full[[i]]                              # Additions possible only to phylo objects, not to multiphylo objects
        #treeD_singleTree[["rate"]] = as.numeric(br_and_meta_D[[i]][,c(-1,-2)])  # Append branch rates to trees
        treeD_singleTree[[i]][["dmv"]] = as.numeric(br_and_meta_D[[i]][,3])
        outD[[i]] = treeD_singleTree
      }
    #LEGACY: }
    #LEGACY: if(class(tree)=="phylo") {
      #LEGACY: #	tree[["check"]]<-as.numeric(br_and_meta_D[[1]][,2])
      #LEGACY: tree[["dmv"]] = as.numeric(br_and_meta_D[[1]][,3])
    #LEGACY: }
  }

###########################
# 4. Summarize and return #
###########################
  if (verboseBool) {
    cat("\tN of tree specs loaded: ", length(outD), "\n", sep="")
  }
  return(outD)
}
