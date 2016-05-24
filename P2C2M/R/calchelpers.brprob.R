calchelpers.brprob <-
function(sTree, gTree, assoc, ploidy, dmvD, node) {
# Descr:    organising likelihood calculations for particular branches
#           Calculate prob of gene tree part, which is a subset 
#           of the species tree.
#           More general: Calculate the prob of subtrees defined 
#           by their MRCA (i.e. their node). 
# Deps:     calchelpers.gtreeparse
#           calchelpers.probcalc
# I/p:      sTree
#           gTree
#           assoc
#           ploidy
#           dmvD = demographic data
#           node
# Note:     branching times = distance from each node to the tips, 
#            under the assumption that the tree is ultrametric

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.brprob", fg="red"),
        sep="")
  }

    ind = calchelpers.gtreeparse(sTree, gTree, assoc, dmvD, node)

    gBt = ind$gBt
    fBt = ind$fBt
    lBt = ind$lBt
    fNd = ind$fNd
    lNd = ind$lNds

    if(fBt!=0) {prob=1}
    else
        {
        # Retain those gTree branching times (i.e. column 1) whose node names 
        # (i.e. column 2) fall within "fBt" and "lBt"
        gBt = gBt[gBt[,1]>=fBt & gBt[,1]<=lBt,]
        # Add lNd row at bottom of gBt-matrix
        gBt = rbind(gBt, c(lBt, lNd))
        # Perform actual prob calculation
        prob = calchelpers.probcalc(gBt, dmvD, ploidy, node, fBt, lBt, fNd, lNd)
        }

return(prob)
}
