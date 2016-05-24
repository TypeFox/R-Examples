Gdist = function(gen,method=1){
 
  # matches with dist.genpop from adegenet

  nloc = ncol(gen)
  recode = function(gen) cbind(gen,(gen-2)*-1)
  X = recode(gen)
  
  # NEI
  if (method == 1) {
    cat("Nei distance\n")
    d = tcrossprod(X)
    vec = sqrt(diag(d))
    d = d/vec[col(d)]
    d = d/vec[row(d)]
    d = -log(d)
    d = as.dist(d)
  }
  # Edwards
  else if (method == 2) {
    cat("Edwards distance\n")
    X2 = sqrt(X)
    d = tcrossprod(X2)
    d = 1 - d/(nloc*2)
    diag(d) = 0
    d = sqrt(d)
    d = as.dist(d)
  }
  # Reynolds
  else if (method == 3) {
    cat("Reynolds distance\n")
    X = X/2
    denomi = tcrossprod(X)
    vec = rowSums((X)^2)
    d = -2 * denomi + vec[col(denomi)] + vec[row(denomi)]
    diag(d) = 0
    denomi = 2 * (nloc - denomi)
    diag(denomi) = 1
    d = d/denomi
    d = sqrt(d)
    d = as.dist(d)
  }
  # Rogers
  else if (method == 4) {
    cat("Rogers distance\n")
    nlig = nrow(X)
    loc.fac = colnames(gen)
    kX = lapply(split(X, loc.fac[col(X)]), matrix, nrow = nlig)
    dcano = function(mat) {
      daux = tcrossprod(mat)
      vec = diag(daux)
      daux = -2*daux+vec[col(daux)]+vec[row(daux)]
      diag(daux) = 0
      daux = sqrt(0.5*daux)
      return(daux)
    }
    d = matrix(0, nlig, nlig)
    for (i in 1:length(kX)) {
      d = d + dcano(kX[[i]])
    }
    d = d/length(kX)/sqrt(2)
    d = as.dist(d)
  }  
return(d)
}