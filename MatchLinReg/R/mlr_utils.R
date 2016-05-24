mlr.add.int <- function(Xcore, X) { # TODO: this is buggy, seems to include binary columns twice
  Xret <- array(NA, dim = c(nrow(X), ncol(Xcore)*ncol(X)))
  colnames(Xret) <- rep(NA, ncol(Xret))
  colnames.xcore <- colnames(Xcore)
  colnames.x <- colnames(X)
  c <- 0
  colnames.standard <- c()
  for (i in 1:ncol(Xcore)) {
    for (j in 1:ncol(X)) {
      colname.tmp <- colname.standardize(paste(colnames.xcore[i], colnames.x[j], sep = ":"))
      
      if (!(colname.tmp %in% colnames.standard) && 
            !(length(unique(Xcore[, i]))==2 && colnames.xcore[i] == colnames.x[j])) {
        c <- c + 1
        colnames(Xret)[c] <- colname.tmp
        Xret[, c] <- Xcore[, i] * X[, j]
        colnames.standard <- c(colnames.standard, colname.tmp)
      }
    }
  }
  return (Xret[, 1:c])
}

mlr.add.step <- function(X, thresh = 20, ncuts = 3) {
  # 1) determine which columns are numeric
  is.numeric <- which(apply(X, 2, function(x) {
    length(unique(x)) >= thresh
  }))
  # 2) define binary variables along quantile values
  ret <- c()
  for (i in 1:length(is.numeric)) {
    qtmp <- quantile(X[, is.numeric[i]], probs = seq(from = 0.0, to = 1.0, length.out = ncuts + 2))
    for (j in 1:ncuts) {
      xtmp <- 1*(X[, is.numeric[i]] < qtmp[j + 1])
      ret <- cbind(ret, xtmp)
      colnames(ret)[ncol(ret)] <- paste(colnames(X)[is.numeric[i]], j, sep = "")
    }
  }
  return (ret)
}

mlr.generate.Z.o <- function(X, interaction.order = 3, step.funcs = TRUE, step.thresh = 20, step.ncuts = 3) {
  X.core <- X
  for (i in 1:(interaction.order - 1)) {
    X.core <- mlr.add.int(X.core, X)
    if (i == 1) {
      X.all <- X.core
    } else {
      X.all <- cbind(X.all, X.core)
    }
  }
  if (step.funcs) {
    X.all <- cbind(X.all, mlr.add.step(X, thresh = step.thresh, ncuts = step.ncuts))
  }
  return (X.all)
}

mlr.smd <- function(tr, X) {
  idx.1 <- which(tr == 1)
  idx.0 <- which(tr == 0)
  
  smd.vec <- apply(X, 2, function(x) {
    (mean(x[idx.1]) - mean(x[idx.0])) / sd(x)
  })
  
  return (smd.vec)
}
