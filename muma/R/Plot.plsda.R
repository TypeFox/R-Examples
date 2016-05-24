Plot.plsda <-
function (pcx, pcy, scaling) {
pwd.score = paste(getwd(), "/PLS-DA_", scaling, "/PLSDA_Scores_", scaling, ".csv", sep="")
score = read.csv(pwd.score, header=TRUE)
score.x = score[,-1]
rownames(score.x) = score[,1]
pwd.load = paste(getwd(), "/PLS-DA_", scaling, "/PLSDA_Loadings_", scaling, ".csv", sep="")
loading = read.csv(pwd.load, header=TRUE)
loading.x = loading[,-1]
rownames(loading.x) = loading[,1]
pwd.p = paste(getwd(), "/PLS-DA_", scaling, "/PLSDA_P_", scaling, ".csv", sep="")
p = read.csv(pwd.p, header=TRUE)
p.x = matrix(p[,-1], ncol=1)
pwd.ptot = paste(getwd(), "/PLS-DA_", scaling, "/PLSDA_Ptot_", scaling, ".csv", sep="")
p = read.csv(pwd.ptot, header=TRUE)
pvar.a = p.x[pcx,]/p
pvar.b = p.x[pcy,]/p
pvar.ai = round(pvar.a*100,1)
pvar.bi = round(pvar.b*100,1)
cum = pvar.ai + pvar.bi
xlab = paste("Component", pcx, "(", pvar.ai, "%)", sep="")
ylab = paste("Component", pcy, "(", pvar.bi, "%)", sep="")
max.pc1 = 1.3*(max(abs(score.x[,pcx])))
 max.pc2 = 1.3*(max(abs(score.x[,pcy])))
 lim = c()
  if (max.pc1 > max.pc2) {lim = c(-max.pc1,max.pc1)} else {lim = c(-max.pc2,max.pc2)}
pwdK = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.csv", sep="")
 k = read.csv(pwdK, header=TRUE)
k.s = k[order(k[,2]),]
tutticolors=matrix(c(1,2,3,4,5,6,7,8,"rosybrown4", "green4", "navy", "purple2", "orange", "pink", "chocolate2", "coral3", "khaki3","thistle","turquoise3","palegreen1","moccasin","olivedrab3","azure4","gold3","deeppink"), ncol=1)
     col=c()
     for(i in 1:nrow(k.s)) {
      col=c(col, tutticolors[k.s[i,2],])
     }
plot(score.x[,pcx], score.x[,pcy], col=col, pch=19, xlab = c(xlab), ylab = c(ylab), xlim = lim, ylim = lim, sub = paste("Cumulative Proportion of Variance Explained = ", cum, "%", sep=""), main = paste("PLS-DA Score Plot (", scaling, ")", sep=""))
text(score.x[,pcx], score.x[,pcy], col=col, cex=0.5, labels=rownames(score.x), pos=1)
axis(1, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
axis(2, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
library(car)
dataEllipse(score.x[,pcx], score.x[,pcy], levels = c(0.95), add=TRUE, col = "black", lwd = 0.4, plot.points=FALSE, center.cex=0.2)
dirout = paste(getwd(), "/PLS-DA_", scaling, "/", sep="")
scor = paste(dirout, "ScorePlot_PLS-DA_", scaling, ".pdf", sep="")
dev.copy2pdf(file=scor)
Max.pc1 = 1.1*(max(loading.x[,pcx]))
 Min.pc1 = 1.1*(min(loading.x[,pcx]))
 Mpc1=c(Min.pc1,Max.pc1)
 Max.pc2 = 1.1*(max(loading.x[,pcy]))
 Min.pc2 = 1.1*(min(loading.x[,pcy]))
 Mpc2=c(Min.pc2,Max.pc2)
dev.new()
plot(loading.x[,pcx], loading.x[,pcy], xlim = Mpc1, ylim = Mpc2, xlab = paste("w*c values ", pcx, sep=""), ylab = paste("w*c values ", pcy, sep=""), main = paste("PLS-DA Loading Plot (", scaling, ")", sep="")) 
text(loading.x[,pcx], loading.x[,pcy], labels=rownames(loading.x), cex=0.7, pos=1)
axis(1, at=Mpc1*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
axis(2, at=Mpc2*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
load = paste(dirout, "W*cPlot_PLS-DA_", scaling, ".pdf", sep="")
dev.copy2pdf(file=load)
}
