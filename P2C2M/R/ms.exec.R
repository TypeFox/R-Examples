ms.exec <-
function(assoc, stree, nTips, ploidy, prmFHndl) {
  ######################
  # Function "ms.exec" #
  ######################
  # Descr:  conduction simulations via Hudson's ms   
  # Deps:   calchelpers.nodetips
  # I/p:    stree = the species tree, inferred via the multispecies 
  #                coalescent model
  #         nTips = the total number of tips to be simulated
  #         assoc = a dataframe with two columns and n rows
  #         prmFHndle
  #
  # Note:   "assoc": The first column contains the tip labels of the 
  #                  population tree and the second contains the number of 
  #                  alleles to be sampled from that population.
  #
  #         "ploidy": refers to *BEAST's "ploidy" in the xml files
  #                   and modifies the DMV values. When all loci have the same 
  #                   ploidy, it should be left as 1. When ploidy varies, it 
  #                   should be 0.5 for mitoch. and 2 for diploid nuclear.
  #
  #         "stree": all trees must be species trees with associated 
  #                 'dmv' values.
  #
  #         dmv-values: "stree$dmv" = Nu haploid alleles, or 2Nu individuals; 
  #                     hence, multiply by 2 to get the standard 4Nu 
  #                     population scaled mutation rate. BEAST's.dmvparse is 
  #                     equal to Nu, where N is the pop size in alleles. 
  #                     Hence: 2Nu in diploid individuals.

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  nReps = get("P2C2M_flg_nReps", envir=P2C2M_globalVars)
  
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> ms.exec", fg="red"), sep="")
  }

###########################################
# 1. Preparatory actions, adjusting theta #
###########################################
  n_sp = length(stree$tip.label)

  theta = stree$dmv*2*ploidy
  # relTheta-values are theta-values devided by the root dmv (which is the 
  # last one in the list, hence: tail(x,1)
  relTheta = theta/tail(theta, 1)

  # Branch lengths of the species tree also must be divided by root theta
  stree$edge.length = stree$edge.length/tail(theta, 1)
  # Set up a data frame of branching times, where the row names (which equal
  # the node number) are their own column
  brTms = as.data.frame(sort(ape::branching.times(stree)))
  brTms[,2] = rownames(brTms)
  # nodeNames is a list of node names plus the root
  nodeNames = c(stree$edge[,2], (n_sp+1))

##############################
# 2. Converting associations #
##############################
  # FROM: alpinus,m206395a
  #       alpinus,m206395b
  #       alpinus,m206400a
  #       alpinus,m206400b
  # TO:   alpinus, 4
  assoc = as.matrix(as.data.frame(table(assoc[,1])))

###################################
# 3. Set initial pop sizes for ms #
###################################
  popSizes = ""
  for(i in 1:n_sp) {
    spName = which(stree$tip.label==assoc[i,1])
    # stree$dmv and stree$edge are in the same order
    size_subpop = relTheta[which(stree$edge[,2]==spName)]
   popSizes = paste(popSizes, "-n", i, size_subpop)
  }

##############################
# 4. Set island model for ms #
##############################
  # "length(stree$tip.label)" provides the number of tips in the species tree
  islModel = paste("-I", length(stree$tip.label), paste(assoc[,2], collapse=' '))

############################################################
# 5. Represent the species tree as past demographic events #
############################################################
  demogrEvents = ""
  for (r in 1:nrow(brTms)) {

    # "brTms" are the ages of all nodes, sorted from smallest to largest
    # "nodeNames" is a list of node names plus the root
    brTmName = brTms[r, 2]
    brTmVal = brTms[r, 1]

    # Calculate a scaling factor
    scaledBrTme = relTheta[which(nodeNames==brTmName)]

    # Take those branches that have the shortes length
    children = sort(stree$edge[stree$edge[,1] == brTmName, 2])
    child1 = sort(calchelpers.nodetips(stree, children[1]))[1]
    child2 = sort(calchelpers.nodetips(stree, children[2]))[1]
    tips = sort(c(child1, child2))

    # Check which of the elements in 
    # "stree$tip.label[tips[2]]" are also in "assoc[,1]"
    popI = which(assoc[,1]==stree$tip.label[tips[2]])
    popJ = which(assoc[,1]==stree$tip.label[tips[1]])

    # -ej t i j: Move all lineages in subpopulation i to subpopulation j at time t
    # -en t i x: Set subpop i size to x*N0 at time t and growth rate to zero
    demogrEvents = paste(demogrEvents, "-ej", brTmVal, popI, popJ, 
                                       "-en", brTmVal, popJ, scaledBrTme)
  }

###################################################
# 6. Combining specifications for full ms command #
###################################################
  # Path to Hudson's ms (Hudson 2002)
  pathToMS = system.file("msdir", "ms", package="P2C2M")
  
  # "nTips" specifies the number of tips in the gene tree
  # "islModel" specifies the number of tips in the species tree and the composition of the number of tips in the gene tree
  # "popSizes" specifies the population sizes in ech subpopulation
  # "demogrEvents" specifies the species tree topology as past demographic events
  fullCmd = paste(pathToMS, nTips, nReps, islModel, popSizes,
                  demogrEvents, "-T")

#################
# 7. Execute ms #
#################
  # The grep command catches the tree spec via the trailing semicolon;
  # the first backslash necessary due to R, the second due to bash shell
  msOutp = system(paste(fullCmd, "| grep \\;"), intern=T)

##############################
# 8. Read in simulated trees #
##############################
  if (nReps == 1) {
    trees = list()
    trees[[1]] = ape::read.tree(text=msOutp)
  }
  else {
    trees = ape::read.tree(text=msOutp)
  }

  # Undoing the adjustment of species tree brlens 
  # as was necessary for ms-simulation
  if (length(trees) > 1) {
    for (t in 1:length(trees)) {
      trees[[t]]$edge.length = trees[[t]]$edge.length*tail(theta, 1)
    }
  }
  else {
    trees$edge.length = trees$edge.length*tail(theta, 1)
  }

  return(trees)
}
