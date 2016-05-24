mst <- function(x) {
  numpoints <- nrow(x)
  if(numpoints == 1) return (list(path = 0, wght = 1, xy0 = 1, xy1 = 1))
  caso <- order(x[lower.tri(x)])
  anc <- vector("integer", numpoints)
  path <- 0
  nanc <- 0
  inmst <- vector("logical", numpoints)
  cl <- ceiling(-0.5 + numpoints - sqrt(0.25 + numpoints*(numpoints-1) - 2*caso))
  rw <- numpoints*(1 - cl) + 0.5*cl*(cl + 1) + caso
  flag <- 0
  i <- 0
  pt1 <- pt2 <- c()
  nlz <- suma <- array(0, numpoints)
  while (flag < (numpoints - 1))
  {
    i <- i + 1
    aux <- 2*inmst[cl[i]] + inmst[rw[i]]
    if(aux == 3 & anc[cl[i]] == anc[rw[i]]) next
    if(aux == 0) {inmst[c(cl[i], rw[i])] <- TRUE; nanc <- nanc + 1; anc[c(cl[i], rw[i])] <- nanc}
    if(aux == 1) {inmst[cl[i]] <- TRUE; anc[cl[i]] <- anc[rw[i]]}
    if(aux == 2) {inmst[rw[i]] <- TRUE; anc[rw[i]] <- anc[cl[i]]}
    if(anc[cl[i]] != anc[rw[i]]) anc[anc == anc[rw[i]]] <- anc[cl[i]]
    path <- path + x[rw[i], cl[i]]
    flag <- flag + 1
    suma[c(rw[i], cl[i])] <- suma[c(rw[i], cl[i])] + x[rw[i], cl[i]]
    nlz[c(rw[i], cl[i])] <- nlz[c(rw[i], cl[i])] + 1
    pt1 <- c(pt1, rw[i])
    pt2 <- c(pt2, cl[i])
  }
  suma <- suma/nlz
  return (list(path = path, wght = suma/sum(suma), xy0 = pt1, xy1 = pt2))
}

