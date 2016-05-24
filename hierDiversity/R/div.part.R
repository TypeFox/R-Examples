div.part <-
function(dat, group, q = 1) {
  if (is.vector(dat)) {
    dat <- matrix(dat, nrow = 1)
  }
  ass.level <- unique(as.vector(group))
  res.table <- array(dim = c(length(ass.level), 5), 
    dimnames = list(NULL, c("alpha", "beta", "gamma", 
      "turnover", "homogeneity")))
  for (i in 1:length(ass.level)) {
    if (is.matrix(group)) {
      x <- dat[apply(group == ass.level[i], 1, 
        sum) > 0, ]
    }
    else {
      x <- dat[group == ass.level[i], ]
    }
    alpha <- dz(x, lev = "alpha", q = q)
    if (is.vector(x) == FALSE) {
      beta <- dz(x, lev = "beta", q = q)
      gamma <- dz(x, lev = "gamma", q = q)
      turnover <- (beta - 1)/(dim(x)[1] - 1)
      homogeneity <- ((1/beta) - (1/dim(x)[1]))/(1 - 
        1/dim(x)[1])
    }
    else if (is.vector(x) == TRUE) {
      beta <- 1
      gamma <- alpha
      turnover <- 0
      homogeneity <- 1
    }
    res.table[i, 1] <- round(alpha, digits = 3)
    res.table[i, 2] <- round(beta, digits = 3)
    res.table[i, 3] <- round(gamma, digits = 3)
    res.table[i, 4] <- round(turnover, digits = 3)
    res.table[i, 5] <- round(homogeneity, digits = 3)
  }
  return(data.frame(ass.level, res.table))
}
