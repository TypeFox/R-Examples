ransy <-
function(scaling, driver.peak) {
  pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, "/ProcessedTable.csv", sep="")
  x <- read.csv(pwd.n, sep=",", header=TRUE)
  x.x <- x[,2:ncol(x)]
  rownames(x.x) <- x[,1]
  x.t <- x.x
  one = matrix(rep(1, ncol(x.t)), nrow=1)
  driver = x.t[,driver.peak]%*%one
  D = x.t/driver
  m = matrix(colMeans(D), nrow=1)
  sd = matrix(apply(D, 2, sd), nrow=1)
  R = m/sd
  R[,driver.peak] = 1
  R[,driver.peak] = max(R)
  Rt = t(R)
  library(gplots)
  plot(Rt, type="h", main = paste("RANSY (", rownames(x.x)[driver.peak], ")", sep=""), ylab = paste("Mean/sd of ratio with ", rownames(x.x)[driver.peak], sep=""), xlab = "Variables")
  text(Rt, labels=colnames(x.x), cex=0.6)
  dirout = paste(getwd(), "/RANSY/", sep="")
  dir.create(dirout)
  out = paste(dirout, "ransy_", colnames(x.x)[driver.peak], ".pdf", sep="")
  dev.copy2pdf(file=out)
}
