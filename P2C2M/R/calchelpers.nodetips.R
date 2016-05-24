calchelpers.nodetips <-
function(tree, node) {
  # Descr:    returns the tip numbers for a given node
  # Deps:     calchelpers.nodetips
  # I/p:      tree = a phylog. tree with numbered internal nodes and numbered tips
  #           node = a node in that tree, specified as an integer

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.nodetips", fg="red"),
        sep="")
  }

    # DEBUGLINES:
    #cat("\ntree$edge\n"); print(tree$edge)
    #cat("\nnode\n"); print(node)

    n_tips = length(tree$tip.label)
    if (node <= n_tips) {node}
    #if (node <= n_tips) {                                              # Potential improvements than line above
    #    outD = node
    #    outD
    #}
    else 
        {
        outD = numeric()
        k = tree$edge[which(tree$edge[,1] == node),2]                   # CURRENT PROBLEM UNDER LCWT: When diff. allele numbers, node nicht in tree$edge[,1]
        for (j in k)
            {
            if (j <= n_tips) {outD = c(outD, j)}
            else {outD = c(outD, calchelpers.nodetips(tree, j))}
            }
        outD        
        }
}
