statshelpers.sigsgn <-
function (qntls, tail) {
  ##################################
  # Function "statshelpers.sigsgn" #
  ##################################
  # Descr:    applies significance signs
  # Deps:     -
  # I/p:      qntls = a data frame of two columns

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.sigsgn", 
        fg="red"), sep="")
  }

  sigSgns = rep(0, length(qntls[,1]))

  # Schemes for one-tailed test
  if (tail=="1l") {
    # Read: "Whichever elements of quants are smaller than zero,
    #        are considered significant and, hence, receive a star 
    #        in the first column."
    sigSgns[which(qntls[,1] > 0)] = "*"
    sigSgns[which(qntls[,1] < 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,1] == 0)] = "n.s."
  }
  if (tail=="1r") {
    sigSgns[which(qntls[,2] < 0)] = "*"
    sigSgns[which(qntls[,2] > 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,2] == 0)] = "n.s."
  }

  # Scheme for two-tailed test
  if (tail=="2") {
    sigSgns[which(qntls[,1] > 0 | qntls[,2] < 0)] = "*"
    sigSgns[which(qntls[,1] < 0 & qntls[,2] > 0)] = "n.s."
    # rare cases
    sigSgns[which(qntls[,1] == 0 & qntls[,2] > 0)] = "n.s."
    sigSgns[which(qntls[,1] < 0 & qntls[,2] == 0)] = "n.s."
  }

  return(sigSgns)
}
