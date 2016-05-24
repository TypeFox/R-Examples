test.best <- function() {
    data(williams)
    checkTrue(all(getDesign(searchCrossOverDesign(s=10, v=5, p=5, start.designs=list(williams5t), n=c(200,10)))==williams5t))
}

test.evaluation <- function() {
  design1 <- t(rbind(c(1,1,2,2),
                     c(2,2,1,1),
                     c(1,1,2,2),
                     c(2,2,1,1),
                     c(1,2,2,1),
                     c(2,1,1,2)))
  design2 <- t(rbind(c(1,1,2,1),
                     c(2,2,1,2),
                     c(1,1,2,1),
                     c(2,2,1,2),
                     c(1,2,2,1),
                     c(2,1,1,2)))
  design3 <- t(rbind(c(1,1,2,2),
                     c(2,2,1,1),
                     c(1,1,2,2),
                     c(2,2,1,1),
                     c(1,1,2,2),
                     c(2,2,1,1)))
  
  general.carryover(design1, model=1)
  general.carryover(design2, model=1)
  general.carryover(design3, model=1)
  design.efficiency(design1)
  design.efficiency(design2)
  design.efficiency(design3)
  
  model <- 1
  v <- 2
  H <- Crossover:::linkMatrix(model=1, v=2)
  for (design in list(design1, design2, design3)) {
    rcD1 <- Crossover:::rcd_R(design, v, model=model) # R-Code
    rcD2 <- Crossover:::rcd(design, v, model=model) # C-Code
    all(rcD1==rcD2)
    X1 <- rcD1
    #X1 <- Crossover:::getRCDesignMatrix(rcD1, v*v+v) # R-Code
    A1 <- Crossover:::infMatrix_R(X1, v, model=model) # R-Code
    A1b <- Crossover:::infMatrix_R(X1, v, model=model, method=2) # R-Code
    all(A1==A1b)
    Csub <- contrMat(n = rep(1, v), type = "Tukey")
    C <- as.matrix(cbind(Csub, matrix(0, dim(Csub)[1], v)))
    general.carryover(design, model=model)
    diag(C %*% ginv(t(H) %*% A1 %*% H) %*% t(C))
    diag(ginv(t(H) %*% A1 %*% H) %*% t(C) %*% C)
    A2 <- Crossover:::infMatrix(X1, v, model=model) # C-Code
    A2
    rcD1
    X1
    A1
    t(H)%*%A1%*%H
  }
  general.carryover(design1, model=model)
}