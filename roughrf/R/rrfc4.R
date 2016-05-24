rrfc4 = function(dat, yvar = ncol(dat), tr, te, mispct, number.trees) {
  hotdeckimp = function(dt, mispct, yvar = ncol(dt)) {
    nr = nrow(dt)
    nc = ncol(dt)
    nd = dt
    nm = floor(nr * mispct)
    for (z in 1:nc) {
      if (!(z %in% yvar)) 
        nd[sample(nr, nm, replace = F), z] = 
        nd[sample(nr, nm, replace = F), z]
    }
    nd
  }
  
  fin = matrix(0, length(te), number.trees)
  for (i in 1:number.trees) {
    pmr = hotdeckimp(dat[tr, ], mispct, yvar)
    rf = randomForest(pmr[, -yvar], pmr[, yvar], ntree = 1)
    fin[, i] = as.numeric(predict(rf, dat[te, -yvar])) - 1
  }
  list(pred = fin)
}