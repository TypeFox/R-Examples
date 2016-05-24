plotnae <- function(d, mpl, return.nae = FALSE, ...){
  tra <- TraFromMpl(mpl)
  d$to[which(d$delta == 0)] <- "cens"
  nae <- mvna(data = d, state.names = 1:ncol(tra), tra = tra, cens.name = "cens")
  plot(nae, bty = "n", ...)
  if(return.nae){
    return(nae)
  }
}