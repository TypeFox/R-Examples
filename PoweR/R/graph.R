graph <- function(matrix.pval, xi = c(seq(0.001,0.009,by=0.001),seq(0.010,0.985,by=0.005),seq(0.990,0.999,by=0.001)),
			      type = c("pvalue.plot", "pvalue.discrepancy", "size.power"), center = FALSE, scale = FALSE) {
	
  type <- match.arg(type)
		
  # Retrieve informations from matrix.pval
  pj <- matrix.pval$pvals
  N <- nrow(pj)
  null.dist <- matrix.pval$null.dist
  stat.indices <- matrix.pval$stat.indices
  n <- matrix.pval$n
  M <- matrix.pval$M
  alter <- matrix.pval$alter
  parstats <- matrix.pval$parstats
  method <- matrix.pval$method

  # Using these values of pj, we calculate Fxi
  Fxi <- calcFx(pj,xi)

### P value plots ###
  if (type == "pvalue.plot") plot.pvalue(Fxi)

### P value discrepancy plots ###
  if (type == "pvalue.discrepancy") plot.discrepancy(Fxi)

### Size-power curves
  if (type == "size.power") {
    pnull <- many.pval(stat.indices, law.index=null.dist, n, M, N, alter, law.pars=NULL, parstats, null.dist, method, center = center, scale = scale)$pvals
    Fxnull <- calcFx(pnull, xi)
    plot.sizepower(Fxi,Fxnull)
  }
	
}


