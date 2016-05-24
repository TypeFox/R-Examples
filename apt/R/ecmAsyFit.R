ecmAsyFit <- function(y, x, lag = 1, split = TRUE, 
    model = c("linear", "tar", "mtar"), thresh = 0) 
{
  if(lag < 1) stop("\nLag for ECM must be 1 or larger.\n") 
  if(!is.ts(y) || !is.ts(x)) stop("Please provide time series data.\n")
  if (!identical(tsp(y), tsp(x))) {
    stop("Properties of y and x are different.\n")
  }
  model <- match.arg(model)

  name.y <- deparse(substitute(y))
  name.x <- deparse(substitute(x))
  A <- start(y); B <- end(y); Q <- tsp(y)[3]
  z <- ts(data = residuals(lm(y~x)), start = A, end = B, frequency = Q)  
  lz <- lag(z, k = -1);  dz <- diff(z); ldz <- lag(dz, k = -1)

  if(model == "linear") {ind <- ifelse(test = lz >= 0,      yes = 1, no=0)} 
  if(model == "tar"   ) {ind <- ifelse(test = lz >= thresh, yes = 1, no=0)}
  if(model == "mtar"  ) {ind <- ifelse(test =ldz >= thresh, yes = 1, no=0)}
  pos <- lz * ind; neg <- lz * (1 - ind)
  
  dx <- diff(x); dy <- diff(y)
  xx <- bsLag(ts.union(dx, dy), lag = lag, prefix = "diff.", 
      var.name = c(name.x, name.y), suffix = ".t_", include.orig = TRUE)
  xp <- ifelse(test = xx >= 0, yes = xx, no = 0)
  xn <- ifelse(test = xx < 0,  yes = xx, no = 0)
  if (tsp(xx)[1] > tsp(pos)[1]) {aa <- start(xx)} else {aa <- start(pos)} 
  if (tsp(xx)[2] < tsp(pos)[2]) {bb <- end(xx)  } else {bb <- end(pos)  }   
  da <- window(cbind(xx, pos, neg), start = aa, end = bb, frequency = Q)
  db <- window(cbind(xp,xn,pos,neg), start = aa, end = bb, frequency = Q)
  colnames(da) <- c(colnames(xx), "ECT.t_1.pos", "ECT.t_1.neg")
  colnames(db) <- c(paste(colnames(xx), ".pos", sep = ""), 
      paste(colnames(xx), ".neg", sep = ""), "ECT.t_1.pos", "ECT.t_1.neg")
  db <- db[, c(((1 + lag) * 0 + 1):((1 + lag) * 1), 
               ((1 + lag) * 2 + 1):((1 + lag) * 3),
               ((1 + lag) * 1 + 1):((1 + lag) * 2), 
               ((1 + lag) * 3 + 1):((1 + lag) * 4 + 2))]

  DepVar.x <- da[, 1]; DepVar.y <- da[, lag+2]
  if (split) { 
    X. <- db[, c(-1, -(lag + 2), -(2 * lag + 3), -(3 * lag + 4))]
  } else {
    X. <- da[, c(-1, -(lag + 2))]
  }
  ecm.x <- lm(DepVar.x ~ 1 + X.)
  ecm.y <- lm(DepVar.y ~ 1 + X.)
  
  result <- listn(y, x, lag, split, model, IndVar=X., name.x, name.y,
    ecm.y, ecm.x)
  if(split) {result$data = db} else { result$data = da}   
  if(model != "linear") {result$thresh = thresh}
  class(result) <- c("ecm", "ecmAsyFit")
  return(result)
} 