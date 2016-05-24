
MCreml = function(y,K,X=NULL,MC=300,samp=300){
  anyNA = function(x) any(is.na(x))
  if(samp>=length(y)){stop("Sample size has to be smaller than sample space")}
  if(nrow(K)!=length(y)){stop("Kinship and response variable have incompatible dimensions")}
  if(ncol(K)!=length(y)){stop("Kinship and response variable have incompatible dimensions")}
  n = MC; t = samp
  moda=function (x){
    it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
      y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
      diffs=sup-inf; i=min(which(diffs==min(diffs)))
      M=median(y[i:(i+k)]); it=it-1}; return(M)}
  h2 = c(); Vg = c(); Ve = c()
  for(i in 1:n){
    R = sample(1:length(y),t)
    if(any(is.na(y[R]))){mis = which(is.na(y[R])); R=R[-mis] }
    fit = reml(y[R],X=X,K=K[R,R])
    Vg = c( fit$VC[1], Vg )
    Ve = c( fit$VC[2], Ve )
    h2 = c( fit$VC[3], h2 )
  }
  H = unlist(h2); G = unlist(Vg); E = unlist(Ve)
  samples = cbind(G,E,H); rownames(samples) = 1:MC
  mode.Vg = moda(G); mode.Ve = moda(E); mode.h2 = moda(H)
  VC = c(mode.Vg,mode.Ve,mode.h2); names(VC) = c("Vg","Ve","h2")
  result = list("samples"=samples,"modes"=VC)
  return(result)
}

Gdist = function(gen,method=1){
  # from adegenet{dist.genpop}
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
    d = as.dist(d)}
  # Edwards
  else if (method == 2) {
    cat("Edwards distance\n")
    X2 = sqrt(X)
    d = tcrossprod(X2)
    d = 1 - d/(nloc*2)
    diag(d) = 0
    d = sqrt(d)
    d = as.dist(d)}
  # Reynolds
  else if (method == 3){
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
    d = as.dist(d)}
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
      return(daux)}
    d = matrix(0, nlig, nlig)
    for (i in 1:length(kX)) {
      d = d + dcano(kX[[i]])}
    d = d/length(kX)/sqrt(2)
    d = as.dist(d)}  
  return(d)
}