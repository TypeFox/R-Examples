volcano <-
function(file, plot.vol) {
 pwdfile=paste(getwd(), "/Univariate/DataTable.csv", sep="")
 file=pwdfile
 x <- read.csv(file, sep=",", header=TRUE)
 x.x = x[,3:ncol(x)]
 rownames(x.x) = x[,2]
 k = matrix(x[,1], ncol=1)
 x.n = cbind(k, x.x)
 sorted = x.n[order(x.n[,1]),]
 sorted.x = as.matrix(sorted[,-1], ncol=ncol(sorted)-1)
 g = c()
 for (i in 1:nrow(sorted)) {
  if (any(g == sorted[i,1])) {g=g}
  else {g=matrix(c(g,sorted[i,1]), ncol=1)}
 }
NoF=nrow(g)
dirout.fc = paste(getwd(), "/Univariate/Fold_Changes/", sep="")
 dir.create(dirout.fc)
dirout.vol = paste(getwd(), "/Univariate/Volcano_Plots/", sep="")
 dir.create(dirout.vol)
for (i in 1:NoF) {
  for (j in 1:NoF) {
   if (i < j) {
    ni=paste("r.",i,".csv",sep="")
    nj=paste("r.",j,".csv",sep="")
    pwdi = paste(getwd(), "/Univariate/Groups/", ni, sep="")
    pwdj = paste(getwd(), "/Univariate/Groups/", nj, sep="")
    pv = paste(getwd(), "/Univariate/Pvalues/Pvalues_", i, "vs", j, ".csv", sep="")
    I=read.csv(pwdi, header=TRUE)
    J=read.csv(pwdj, header=TRUE)
    I = I[,-1]
    J = J[,-1]
meanI = matrix(colMeans(I), ncol=ncol(I))
meanJ = matrix(colMeans(J), ncol=ncol(J))
MeanI = matrix(rep(NA, ncol(I)), nrow=1)
MeanJ = matrix(rep(NA, ncol(I)), nrow=1)
  for (m in 1:ncol(I)) {
    if (meanI[,m] < 0 | meanJ[,m] < 0) {
      MeanI[,m] = 1
      MeanJ[,m] = 1
    } else {
	MeanI[,m] = meanI[,m]
	MeanJ[,m] = meanJ[,m]
      }
  }
	
FC = matrix(MeanI/MeanJ, nrow=ncol(I))
rownames(FC) = colnames(I)
fc.csvfile = paste("Fold_Change_", i, "vs", j, ".csv", sep="")
write.csv(FC, paste(dirout.fc, fc.csvfile, sep=""))
PV = read.csv(pv, header=TRUE)
PV = matrix(PV[,-1], ncol=1)
logfc = log2(FC)
logpv = -log10(PV)
colpv=matrix(rep(NA, nrow(PV)), ncol=1)
  for (p in 1:nrow(PV)) {
    if (logfc[p,] < -0.3219281 | logfc[p,] > 0.2630344) {
      if (logpv[p,] > 1.30103) {
      colpv[p,] = "navy"
      } else {
      colpv[p,] = "dark grey"
	}
    } else { 
      colpv[p,] = "dark grey"
      }
  }
max.fc = 1.3*(max(abs(logfc)))
V = paste(dirout.vol, "VolcanoPlot_", i, "vs", j, ".pdf", sep="")
pospv=matrix(rep(NA, nrow(PV)), ncol=1)
    for (p in 1:nrow(PV)) {
        if (logfc[p,] < 0) {
            pospv[p,] = 2
        } else { 
        pospv[p,] = 4
    }
   }
pdf(V)
plot(logfc, logpv, col=colpv, pch = 19, xlim=c(-max.fc,max.fc), xlab = "Log2 (Fold Change)", ylab = "Log10 (Pvalue)", main = paste("Volcano Plot ", i, " vs ", j, sep=""), sub = "(Variables in Blue are significant (Pvalue<0.05) and showed Fold Changes >1.2 or <0.8)")
text(logfc, logpv, labels=colnames(sorted.x), cex=0.8, pos=pospv, col=colpv)
axis(2, at = c(-1,150), pos=c(-0.3219281,0), col="blue", lwd=0.3)
axis(2, at = c(-1,150), pos=c(0.2630344,0), col="blue", lwd=0.3)
axis(1, at = c(-150,150), pos=c(1.30103,0), col="blue", lwd=0.3)
dev.off()
if (plot.vol) {
	dev.new()
plot(logfc, logpv, col=colpv, pch = 19, xlim=c(-max.fc,max.fc), xlab = "Log2 (Fold Change)", ylab = "Log10 (Pvalue)", main = paste("Volcano Plot ", i, " vs ", j, sep=""), sub = "(Variables in Blue are significant (Pvalue<0.05) and showed Fold Changes >1.2 or <0.8)")
text(logfc, logpv, labels=colnames(sorted.x), cex=0.8, pos=pospv, col=colpv)
axis(2, at = c(-1,150), pos=c(-0.3219281,0), col="blue", lwd=0.3)
axis(2, at = c(-1,150), pos=c(0.2630344,0), col="blue", lwd=0.3)
axis(1, at = c(-150,150), pos=c(1.30103,0), col="blue", lwd=0.3)
}
}}}
}
