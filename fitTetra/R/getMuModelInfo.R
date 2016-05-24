getMuModelInfo <-
function(mutype, ng=1) {
#Note: order must match gEMMax.mu.ratiodist !
#Note: list items 1:7 correspond to mutype 0:6 !
  modelinfo <- list()
  modelinfo[[1]] <- list(
    name = paste("free, ng=",ng,sep=""),
    npar = ng) # 0=mu.free: ng means
  modelinfo[[2]] <- list(
    name = "b2",
    npar = 3) # 1 : c1, c2, f
  modelinfo[[3]] <- list(
    name = "b1",
    npar = 2) # c, f
  modelinfo[[4]] <- list(
    name = "b2,q2",
    npar = 5) # c1, c2, d1, d2, f
  modelinfo[[5]] <- list(
    name = "b1,q2",
    npar = 4) # c, d1, d2, f
  modelinfo[[6]] <- list(
    name = "b2,q",
    npar = 4) # c1, c2, d, f
  modelinfo[[7]] <- list(
    name = "b1,q",
    npar = 3) # c, d, f

  modelinfo[[mutype+1]]
}
