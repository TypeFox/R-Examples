zero.children <-
function(P) {
  if (P$recursive) {
    if (length(P$children) == 0) return(TRUE)
    zeros <- logical(length(P$children))
    for (i in 1:length(P$children)) {
      zeros[i] <- zero.children(P$children[[i]])
    }
    if (all(zeros)) return(TRUE)
    else return(FALSE)
  }
  if (length(P$var) == 0) return(TRUE)
  else return(FALSE)
}