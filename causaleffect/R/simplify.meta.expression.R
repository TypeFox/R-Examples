simplify.meta.expression <-
function(P, D, to) {
  if (P$recursive) {
    for (i in 1:length(P$children)) {
      P$children[[i]] <- simplify.meta.expression(P$children[[i]], D, to)
    }
  } else {
    domain <- which(letters == P$domain)
    if (length(domain) > 0) return(simplify.expression(P, D[[domain]], to[[domain]]))
    else return(simplify.expression(P, D[[1]], to[[1]]))
  } 
  return(P)
}
