rrfc7 = function(dat, yvar = ncol(dat), tr, te, mispct, number.trees) {
  rdms = function(dt, pct, kpc) {
    nr = nrow(dt)
    nc = ncol(dt)
    nd = dt
    nm = floor(nr * pct)
    for (z in 1:nc) {
      if (!(z %in% kpc)) 
        nd[sample(nr, nm, replace = F), z] = NA
    }
    nd
  }
  
  fin = matrix(0, length(te), number.trees)
  for (i in 1:number.trees) {
    misdat = rdms(dat[tr, ], mispct, yvar)
    pmr = rfImpute(misdat[, -yvar], misdat[, yvar])
    rf = randomForest(pmr[, -1], pmr[, 1], ntree = 1)
    fin[, i] = as.numeric(predict(rf, dat[te, -yvar])) - 1
  }
  list(pred = fin)
}