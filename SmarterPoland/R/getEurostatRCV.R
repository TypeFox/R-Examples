getEurostatRCV <-
function(kod = "educ_iste", ...) {
  dat <- getEurostatRaw(kod, ...)
  dat2 <- t(as.data.frame(strsplit(as.character(dat[,1]), split=",")))
  cnames <- strsplit(colnames(dat)[1], split="[,\\\\]")[[1]]
  colnames(dat2) <- cnames[-length(cnames)]

  measures <- colnames(dat)[-1]
  df <- do.call("rbind", lapply(measures, function(x) {
    data.frame(dat2, x, dat[, x])
  }))
  rownames(df) <- apply(df[,-ncol(df),drop=FALSE], 1, paste, collapse="_")
  colnames(df) <- c(cnames, "value")
  df
}

