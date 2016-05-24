rrfc5 = function(dat, yvar = ncol(dat), tr, te, mispct, number.trees) {
  regdat = function(dat) {
    datmis = dat
    nc = ncol(dat)
    for (i in 1:nc) {
      if (is.factor(dat[, i])) 
        datmis[, i] = predict(multinom(dat[, i] ~ ., data = dat[, 
                                                                -i]))
      if (!is.factor(dat[, i])) 
        datmis[, i] = predict(glm(dat[, i] ~ ., data = dat[, -i]))
    }
    datmis
  }
  regdat = regdat(dat[tr, ])
  
  regimp = function(dt, mispct, yvar = ncol(dt)) {
    nr = nrow(dt)
    nc = ncol(dt)
    nd = dt
    nm = floor(nr * mispct)
    for (z in 1:nc) {
      smpcs = sample(nr, nm, replace = F)
      if (!(z %in% yvar)) 
        nd[smpcs, z] = regdat[smpcs, z]
    }
    nd
  }
  
  fin = matrix(0, length(te), number.trees)
  for (i in 1:number.trees) {
    pmr = regimp(dat[tr, ], mispct, yvar)
    rf = randomForest(pmr[, -yvar], pmr[, yvar], ntree = 1)
    fin[, i] = as.numeric(predict(rf, dat[te, -yvar])) - 1
  }
  list(pred = fin)
}