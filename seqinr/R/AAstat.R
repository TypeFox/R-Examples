data(sysdata, envir=environment())
AAstat <- function(seq, plot = TRUE){
  #
  # seq is a protein sequence as a vector of (upper case) chars.
  #
  AAP <- SEQINR.UTIL$AA.PROPERTY
  tutu <- lapply(names(AAP), function(x) which(seq %in% AAP[[x]]))
  names(tutu) <- names(AAP)

  n.items <- length(tutu) # Number of physoco-chemical properties
  n.res <- length(seq)    # Number of residues in the protein
  
  if(plot == TRUE){
    coul <- rainbow(n.items)
    plot(c(0, n.res), c(0, n.items + 1), type = "n", axes = FALSE, 
      ann = FALSE, xlim = c(0, n.res + 1))
    title(xlab = "Position of the residues along the sequence")
    axis(2, at = seq(1.5, 10, 1), labels = names(tutu), col.lab = "blue", las = 1, cex.axis = 0.8)
    axis(1, at = seq(0, n.res, 15), labels = seq(0, n.res, 15), col.axis = "blue")
    lapply(seq_len(n.items), function(x){
      segments(tutu[[x]], x, tutu[[x]], x + 1, col = coul[x], lwd = 2)
      rect(0, x, n.res, 1, lwd = 2)
      })
    rect(0, n.items + 1, n.res, 1, lwd = 2)
  }

  res1 <- lapply(tutu,function(x){length(x)/n.res})
  res2 <- table(factor(seq, levels = levels(SEQINR.UTIL$CODON.AA$L)))
  res3 <- computePI(seq)
  return(list(Compo = res2, Prop = res1, Pi = res3))
}
