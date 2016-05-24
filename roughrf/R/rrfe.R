rrfe = function(dat, yvar = ncol(dat), tr, te, mispct, number.trees, k) {
  rf = randomForest(dat[tr, -yvar], dat[tr, yvar])
  vp = varImpPlot(rf)
  
  
  rdms = function(dt, pct, kpc, k) {
    colmis = c(1:length(vp))[rbinom(length(vp), 1, (vp/max(vp))^(k)) == 
                               1]
    nr = nrow(dt)
    nc = ncol(dt)
    nd = dt
    nm = floor(nr * 0.01)
    for (z in 1:nc) {
      if ((z %in% colmis)) 
        nd[sample(nr, nm, replace = F), z] = NA
    }
    mfix(nd, 1)
  }
  
  
  fin = matrix(0, length(te), number.trees)
  for (i in 1:number.trees) {
    pmr = rdms(dat[tr, ], mispct, yvar, k)
    rf = randomForest(pmr[, -yvar], pmr[, yvar], ntree = 1)
    fin[, i] = as.numeric(predict(rf, dat[te, -yvar])) - 1
  }
  list(pred = fin)
} 