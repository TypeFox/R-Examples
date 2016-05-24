simplify <- 
function(P, G.Adj) {
  factorials <- sapply(1:length(P$children), factorial)
  j <- 0
  while (j < length(P$sumset)) {
    j <- j + 1
    vars <- unlist(lapply(P$children, FUN = function(x) x$var))
    i <- which(vars == P$sumset[j])
    k <- 1
    n <- 1
    terms <- gather(P, i, G.Adj)
    I <- terms[[1]]
    P <- terms[[2]]
    q <- length(I) + 1
    omega <- perm(I, n, i)
    fail.terms <- list()
    J <- character()
    D <- character()    
    while (k <= q) {  
      joint <- join(J, D, P$children[[omega[k]]]$var, P$children[[omega[k]]]$cond, G.Adj)
      if (length(joint[[1]]) <= length(J)) {
        #fail.terms <- c(fail.terms, list(omega[1:(k-1)]))
        #n <- n + 1
        #while (n < factorials[q-1]) {
        #  omega <- perm(I, n, q)
        #  if (is.prefix(fail.terms, omega)) n <- n + 1
        #  else break
        #}
        #if (n > factorials[q-1]) break
        J <- character()
        D <- character()
        k <- 1
        break
      } else {
        k <- k + 1
        J <- joint[[1]]
        D <- joint[[2]]
      }
    }
    if (k == q + 1) {
      P <- factorize(J, D, P, omega, q)     
      P$children[[i]] <- NULL
      P$sumset <- P$sumset[-j]
      j <- 0 
    }
  }
  return(P)       
}