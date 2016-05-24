calchelpers.gtreeparse <-
function(sTree, gTree, assoc, dmvD, node) {
# Descr:    parsing branches for likelihood calculations
# Deps:     calchelpers.descend
# I/p:      sTree
#           gTree
#           assoc
#           dmvD
#           node
# Note:     branching times = distance from each node to the tips, 
#           under the assumption that the tree is ultrametric
#           gBt = gene tree branching times (plural!)
#           fBt = first branching time
#           lBt = last branching time
#           fNd = first node
#           lNd = last node

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.gtreeparse", fg="red"),
        sep="")
  }

    # returns the subtree that starts at node <node>
    subtree = calchelpers.descend(sTree, gTree, assoc, node)
    # infers the branching times for a subtree (here the gene tree)
    gBt = sort(ape::branching.times(subtree))
    # generate a matrix with branching times in first column and node 
    # names in second column
    gBt = c(0, gBt)
    gBt = cbind(gBt, length(gBt):1)
    # get the branching times for a node in the species tree 
    fBt = dmvD[node,"sbt"]
    # get the branching times for a node in the species tree + the branch 
    # lengths for that node
    lBt = (dmvD[node,"sbt"] + dmvD[node,"length"])
    # get those node IDs (column 2), whose gene tree branching times (column 1) 
    # are equal the largest of the species tree branching times
    fBtMax = max(gBt[gBt[,1]<=fBt,1])
    fNd = gBt[gBt[,1]==fBtMax,2]
    # get those node IDs (column 2), whose  gene tree branching times (column 1) 
    # that equal the largest of the species tree branching times + branch lengths
    lBtMax = max(gBt[gBt[,1]<=lBt,1])
    lNd = gBt[gBt[,1]==lBtMax,2]

    outd = list("gBt"=gBt, "fBt"=fBt, "lBt"=lBt, "fNd"=fNd, "lNd"=lNd)

return(outd)
}
