join <- function(J, D, vari, cond, G.Adj) {
  J.new <- character()
  D.new <- character()
  if (length(J) == 0) {
    J.new <- vari
    D.new <- cond
    return(list(J.new, D.new))
  }

  if (length(intersect(J, cond)) > 0 && vari %in% D) {
    return(list(J, D))
  }

  cond.uni <- union(D, cond)

  if (length(cond.uni) > 0) {
    ds <- powerset(cond.uni, nonempty = FALSE)
    n <- length(ds)
    for (i in 1:n) {
      a.set <- union(setdiff(ds[[i]], setdiff(D, vari)), setdiff(setdiff(D, vari), ds[[i]]))
      b.set <- union(setdiff(ds[[i]], setdiff(cond, J)), setdiff(setdiff(cond, J), ds[[i]]))
      if (wrap.dSep(G.Adj, J, a.set, setdiff(D, a.set)) && 
          wrap.dSep(G.Adj, vari, b.set, setdiff(cond, b.set))) {
        J.new <- union(J, vari)
        D.new <- ds[[i]]
        return(list(J.new, D.new))
      }
    } 
  } else {
    J.new <- union(J, vari)
    D.new <- cond
    return(list(J.new, D.new))
  }

  return(list(J, D))

}
