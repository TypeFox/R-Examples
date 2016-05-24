imputeQs <- function(data) {
  dat <- data
  A <- dat[(sapply(dat, function(x) is.numeric(x)) == T)]
  Missing.Vars <- colnames(A)[apply(is.na(A), 2, any)]
  New.Vars <- sapply(Missing.Vars, function(this.var) {
    cut(dat[, this.var], include.lowest = T,
        breaks = unique(quantile(dat[, this.var], na.rm = T)), 
        ordered_result = T)
  })
  New.Vars <- apply(New.Vars, 2, function(x) ifelse(x == "" | is.na(x), "Missing", x))
  colnames(New.Vars) <- Missing.Vars
  New.Vars <- data.frame(New.Vars)
  dat[, Missing.Vars] <- New.Vars[, Missing.Vars]
  dat.factors <- dat[(sapply(dat, function(x) is.factor(x)) == T)]
  dat2 <- apply(dat.factors, 2, function(x) ifelse(x == "" | is.na(x), "Missing", x))
  dat3 <- as.data.frame(dat2)
  names(dat3) <- names(dat.factors)
  dat[, names(dat3)] <- dat3
  return(dat)
}
