calchelpers.probcalc <-
function(gBt, dmvD, ploidy, node, fBt, lBt, fNd, lNd) {
# Descr:    calculating likelihoods for particular branches
# Deps:     -
# I/p:      gBt
#           dmvD
#           ploidy
#           node
#           fBt
#           lBt
#           fNd
#           lNd
# Note:     The input variable "gBt" contains two columns: 
#           branching time differences and node ID
#           branching times = distance from each node to the tips, 
#           under the assumption that the tree is ultrametric

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> calchelpers.probcalc", fg="red"),
        sep="")
  }

# 1. Calculating differences in branching times
    value = length(gBt[,1])
    # waitTms = branching time differences between two branches
    waitTms = gBt[2:value,1] - gBt[1:value-1,1]
    # append node IDs again (so that "waitTms" mimicks "gBt"), 
    # with the exception of last node ID in list
    waitTms = cbind(waitTms, gBt[1:value-1,2])
# 2. Calculating lists of values
    dmv = (2*dmvD[node,"dmv"])*ploidy
    # lambda is a list of values, each calculated according to the following formula
    lambda = (waitTms[,2]*(waitTms[,2]-1))/dmv
    # exponent is a list of values; Euler-constant to the power of (-lambda*waitTms[,1])
    exponent = exp(-lambda*waitTms[,1])
# 3. Calculating products for each list
    # Calculate product of exponent list, given that lambda is never 0
    exponent = prod(exponent[lambda!=0])
    # Calculate product of lambda list, however disregarding last list element
    lambda = prod(lambda[1:(length(lambda)-1)])
# 4. Calculating overall probability
    prob = lambda * exponent

return(prob)
}
