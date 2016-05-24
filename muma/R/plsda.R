plsda <-
function (scaling) {
pwd.x = paste(getwd(), "/Preprocessing_Data_", scaling, "/ProcessedTable.csv", sep="")
x = read.csv(pwd.x, header=TRUE)
 x.x = x[,2:ncol(x)]
 rownames(x.x) = x[,1]
 pwdK = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.csv", sep="")
 k = read.csv(pwdK, header=TRUE)
 k.x = matrix(k[,-1], ncol=1)
 x.n = cbind(k.x, x.x)
 sorted = x.n[order(x.n[,1]),]
 g = c()
 for (i in 1:nrow(sorted)) { 
  if (any(g == sorted[i,1])) {g=g} 
  else {g=matrix(c(g,sorted[i,1]), ncol=1)}
 }
dimB = nrow(g)*nrow(sorted)
B = matrix(rep(NA, dimB), ncol=nrow(g))
for (i in 1:nrow(sorted)) {
  for (j in 1:nrow(g)) {
    if (sorted[i,1] == j) { 
      B[i,j] = 1}
    else {
      B[i,j] = 0
    }
  }
}
library(pls)
sorted.x = sorted[,-1]
sorted.un = matrix(unlist(sorted.x), ncol=ncol(sorted.x))
P = plsr(B ~ sorted.un, ncomp = nrow(g)-1, method = c("kernelpls"), validation = "CV")
rownames(P$scores) =  rownames(sorted.x)
rownames(P$loadings) =  colnames(sorted.x)
dirout = paste(getwd(), "/PLS-DA_", scaling, "/", sep="")
dir.create(dirout)
out.score = paste(dirout, "PLSDA_Scores_", scaling, ".csv", sep="")
write.csv(P$scores, out.score)
out.load = paste(dirout, "PLSDA_Loadings_", scaling, ".csv", sep="")
write.csv(P$loadings, out.load)
k = matrix(sorted[,1], ncol=1)
tutticolors=matrix(c(1,2,3,4,5,6,7,8,"rosybrown4", "green4", "navy", "purple2", "orange", "pink", "chocolate2", "coral3", "khaki3","thistle","turquoise3","palegreen1","moccasin","olivedrab3","azure4","gold3","deeppink"), ncol=1)
     col=c()
     for(i in 1:nrow(k)) {
      col=c(col, tutticolors[k[i,],])
     }
if (ncol(P$scores) == 1) {
  xlab = "Samples"
  ylab = "Score values Component 1"
  plot(P$scores[,1], col=col, pch=19, xlab=c(xlab), ylab=c(ylab), main = paste("PLS-DA Score Plot (", scaling, ")", sep=""))
  lim1 = nrow(P$scores)*2
  axis(1, at=c(-lim1,lim1), col="grey", pos=c(0,0), labels=FALSE, lwd=1)
  text(P$scores[,1], col=col, cex=0.5, pos=1, labels=rownames(P$scores))
  pwdout=paste(dirout, "ScorePlot_PLSDA_1Component_", scaling, ".pdf", sep="")
  dev.copy2pdf(file=pwdout)
  Max.pc2 = 1.1*(max(P$loadings[,1]))
  Min.pc2 = 1.1*(min(P$loadings[,1]))
  Mpc2=c(Min.pc2,Max.pc2)
  dev.new()
  plot(P$loadings[,1], ylim=Mpc2, main = paste("PLS-DA Loading Plot (", scaling, ")", sep=""), xlab="Variables", ylab="W*c values (Component1)")
  text(P$loadings[,1], cex=0.7, pos=1, labels=rownames(P$loadings))
  pwdout1=paste(dirout, "W*cPlot_PLSDA_1Component_", scaling, ".pdf", sep="")
  dev.copy2pdf(file=pwdout1)
} else {
pairs = c()
     if (ncol(P$scores) >= 10) {pairs = c(10)} else {pairs = c(ncol(P$scores))}

     pairs(P$scores[,1:pairs],col=col)
     pairs = paste(dirout, "Pairs_PLSDA_", scaling, ".pdf", sep="")
     dev.copy2pdf(file=pairs)
}
p.v = matrix(P$Xvar, ncol=1)
	p.v.csv = paste(dirout, "PLSDA_P_", scaling, ".csv", sep="")
	write.csv(p.v, file=p.v.csv)
p.vtot = matrix(P$Xtotvar, ncol=1)
	p.vtot.csv = paste(dirout, "PLSDA_Ptot_", scaling, ".csv", sep="")
	write.csv(p.vtot, file=p.vtot.csv, row.names = FALSE)
}
