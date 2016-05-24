"profiles.drm" <-
  function(n.categories,n.repetitions, structure="exchangeable"){
    match.arg(structure,c("exchangeable","M","M2"))
    if(structure=="exchangeable")
      dep <- "I"
    else(dep <- structure)
    v <- expand.grid(rep(list(1:n.categories), n.repetitions))
    dims <- dim(kronecker.drm(v[1, ], nclass = n.categories, dep = dep))
    w <- apply(v, 1, function(i, nclass, dep) {
      kronecker.drm(i, nclass = nclass, dep = dep)
    }, nclass = n.categories, dep = dep)
    if (length(grep("M", dep)) > 0) {
      w <- array(w, dim = c(dims, nrow(v)))
    }
    w
  }







