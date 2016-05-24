gather <- function(P, i, G.Adj) {
  ind <- c()
  vari <- P$children[[i]]$var
  for (k in (1:length(P$children))[-i]) {
    if (vari %in% P$children[[k]]$cond) {
      if (!dSep(G.Adj, P$children[[k]]$var, vari, setdiff(P$children[[k]]$cond, vari))) {
        ind <- c(ind, k)
      } else {
        P$children[[k]]$cond <- setdiff(P$children[[k]]$cond, vari)
      }
    }
  }
  return(list(ind, P))
}