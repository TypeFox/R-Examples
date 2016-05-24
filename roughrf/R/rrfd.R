rrfd = function(dat, yvar = ncol(dat), tr, te, mispct, number.trees, m) {
  rdms = function(dt, pct, kpc) {
    nr = nrow(dt)
    nc = ncol(dt)
    nd = dt
    nm = floor(nr * pct)
    for (z in 1:nc) {
      if (!(z %in% kpc)) 
        nd[sample(nr, nm, replace = F), z] = NA
    }
    mfix(nd, 1)
  }
  
  fin = matrix(0, length(te), number.trees)
  for (i in 1:number.trees) {
    pmr = rdms(dat[tr, ], mispct, yvar)
    rf = randomForest(pmr[, -yvar], pmr[, yvar], ntree = 1, mtry = m)
    fin[, i] = as.numeric(predict(rf, dat[te, -yvar])) - 1
  }
  list(pred = fin)
}
