Bon_EV <- function(pvalue, alpha)
{
  new_MTP_adjp <- rep(length(pvalue))
  BH <- p.adjust(pvalue, "BH")
  qobj <- qvalue(p = pvalue)
  qvalues <- qobj$qvalues
  ngene <- length(pvalue)
  for (i in 1: ngene)
  { 
    adjpv <- ngene*(pi0est(pvalue)$pi0)/sum(BH <= alpha, na.rm=TRUE)*pvalue
    new_MTP_adjp[i] <- min(adjpv[i], 1)
  }
  mylist <- list(raw_P_value = pvalue, BH_adjp = BH, Storey_adjp = qvalues, Bon_EV_adjp = new_MTP_adjp)
  return(mylist)
}



