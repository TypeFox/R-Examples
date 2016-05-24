oplsda <-
function(scaling) {
  	pwd.x = paste(getwd(), "/Preprocessing_Data_", scaling, "/ProcessedTable.csv", sep="")
	x = read.csv(pwd.x, header=TRUE)
 	x.x = x[,2:ncol(x)]
 	rownames(x.x) = x[,1]
 	pwdK = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.csv", sep="")
	k = read.csv(pwdK, header=TRUE)
 	k.x = matrix(k[,-1], ncol=1)
 	x.n = cbind(k.x, x.x)
 	sorted = x.n[order(x.n[,1]),]
	k = matrix(sorted[,1], ncol=1)
 	g = c()
 	for (i in 1:nrow(sorted)) { 
  		if (any(g == sorted[i,1])) {
  			g=g
  		} else {
  			g=matrix(c(g,sorted[i,1]), ncol=1)
  			}
 	}
	Y = matrix(rep(NA, nrow(sorted)), ncol=1)
	for (i in 1:nrow(sorted)) {
  		for (l in 1:2) {
    		if (sorted[i,1] == l) { 
      			Y[i,] = 0}
    		else {
      			Y[i,] = 1
    		}
  		}
	}
  X = as.matrix(sorted[,-1], ncol=ncol(sorted)-1)
  nf = 1
  T=c()
  P=c()
  C=c()
  W=c()
  Tortho=c()
  Portho=c()
  Wortho=c()
  Cortho=c()
for (j in 1:nf) {  
  w = (t(X)%*%Y)%*%solve(t(Y)%*%Y)
  w1 = t(w)%*%w
  w2 = abs(sqrt(w1))
  w = w%*%solve(w2)
  t = (X%*%w)%*%solve(t(w)%*%w)
  t1 = t(t)%*%t
  c = t(Y)%*%t%*%solve(t1)
  c1 = t(c)%*%c
  u = Y%*%c%*%solve(c1)
  u1 = t(u)%*%u
  u2 = abs(sqrt(u1))
  #conv = abs((u2-o2)%*%solve(o2))
  p = (t(X)%*%t)%*%solve(t1)
  wortho = p - w
  wortho1 = t(wortho)%*%wortho
  wortho2 = abs(sqrt(abs(wortho1)))
  wortho = wortho%*%solve(wortho2)
  tortho = X%*%wortho%*%solve(t(wortho)%*%wortho)
  tortho1 = t(tortho)%*%tortho
  portho = t(X)%*%tortho%*%solve(tortho1)
  cortho = t(Y)%*%tortho%*%solve(tortho1)
  X = X - tortho%*%t(portho)
 
  T=matrix(c(T,t))
  P=matrix(c(P,p))
  C=matrix(c(C,c))
  W=matrix(c(W,w))
  Tortho=matrix(c(Tortho,tortho))
  Portho=matrix(c(Portho,portho))
  Wortho=matrix(c(Wortho,wortho))
  Cortho=matrix(c(Cortho,cortho))
}
T = matrix(T, ncol=nf)
T = scale(T, scale=FALSE, center=TRUE)
P = matrix(P, ncol=nf)
C = matrix(C, ncol=nf)
W = matrix(W, ncol=nf)
Tortho = matrix(Tortho, ncol=nf)
Portho = matrix(Portho, ncol=nf)
Cortho = matrix(Cortho, ncol=nf)
Wortho = matrix(Wortho, ncol=nf)
Xortho = Tortho%*%t(Portho)
max.pc1 = 1.3*(max(abs(T[,nf])))
 max.pc2 = 1.3*(max(abs(Tortho[,nf])))
 lim = c()
  if (max.pc1 > max.pc2) {lim = c(-max.pc1,max.pc1)} else {lim = c(-max.pc2,max.pc2)}
tutticolors=matrix(c(1,2,3,4,5,6,7,8,"rosybrown4", "green4", "navy", "purple2", "orange", "pink", "chocolate2", "coral3", "khaki3","thistle","turquoise3","palegreen1","moccasin","olivedrab3","azure4","gold3","deeppink"), ncol=1)
     col=c()
     for(i in 1:nrow(k)) {
      col=c(col, tutticolors[k[i,],])
     }
plot(T[,nf], Tortho[,1], col=col, pch=19, xlim = lim, ylim = lim, xlab="T score [1]", ylab = "Orthogonal T score [1]", main = "OPLS-DA score scatter plot")
text(T[,nf], Tortho[,1], col=col, labels=rownames(sorted), cex=0.5, pos=1)
axis(1, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
axis(2, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
library(car)
dataEllipse(T[,nf], Tortho[,1], levels = c(0.95), add=TRUE, col = "black", lwd = 0.4, plot.points=FALSE, center.cex=0.2)
dirout = paste(getwd(), "/OPLS-DA", scaling, "/", sep="")
dir.create(dirout)
pwdxdef = paste(dirout, "X_deflated.csv", sep="")
write.csv(X, pwdxdef)
scor = paste(dirout, "ScorePlot_OPLS-DA_", scaling, ".pdf", sep="")
dev.copy2pdf(file=scor)
pwdT = paste(dirout, "TScore_Matrix.csv", sep="")
write.csv(T, pwdT)
pwdTortho = paste(dirout, "TorthoScore_Matrix.csv", sep="")
write.csv(T, pwdTortho)
#S-plot
s = as.matrix(sorted[,-1], ncol=ncol(sorted)-1)
p1 = c()
for (i in 1:ncol(s)) {
	scov = cov(s[,i], T)
	p1 = matrix(c(p1, scov), ncol=1)
}
pcorr1 = c()
for (i in 1:nrow(p1)) {
	den = apply(T, 2, sd)*sd(s[,i])
	corr1 = p1[i,]/den
	pcorr1 = matrix(c(pcorr1, corr1), ncol=1)
}
pwdp1 = paste(dirout, "p1_Matrix.csv", sep="")
write.csv(p1, pwdp1)
pwdpcorr1 = paste(dirout, "pcorr1_Matrix.csv", sep="")
write.csv(pcorr1, pwdpcorr1)
dev.new()
plot(p1, pcorr1, xlab="p[1]", ylab ="p(corr)[1]", main = paste("S-plot (OPLS-DA) ", scaling, sep=""))
text(p1, pcorr1, labels=colnames(s), cex=0.5, pos=1)
splot = paste(dirout, "SPlot_OPLS-DA_", scaling, ".pdf", sep="")
dev.copy2pdf(file=splot)
#PCA_OPLS
pc.all <- prcomp(X, center=FALSE, scale=FALSE)
     p.v <- matrix(((pc.all$sdev^2)/(sum(pc.all$sdev^2))), ncol = 1)
     p.i <- round(p.v*100,1)
     p.z <- matrix(1,nrow(p.i),1)
     p.f <- cbind(p.i,p.z)
     dirout.pca = paste(dirout, "PCA_OPLS/", sep="")
     dir.create(dirout.pca)
     write.csv(p.f, paste(dirout.pca, "PCA_P_OPLS", sep=""))
     write.csv(pc.all$x, paste(dirout.pca, "PCA_OPLS_ScoreMatrix.csv", sep=""))
     write.csv(pc.all$rotation, paste(dirout.pca,"PCA_OPLS_LoadingsMatrix.csv", sep=""))
     cum = p.f[1,1] + p.f[2,1]
     lim = c()
	 max.pc1 = 1.3*(max(abs(pc.all$x[,1])))
 	 max.pc2 = 1.3*(max(abs(pc.all$x[,2])))
 	 if (max.pc1 > max.pc2) {lim = c(-max.pc1,max.pc1)} else {lim = c(-max.pc2,max.pc2)}
	 pca <- paste("PC1 (",p.f[1,1],") %")
 	 pcb <- paste("PC2 (",p.f[2,1],") %")
	 xlab = c(pca)
 	 ylab = c(pcb)
	 D <- paste(dirout.pca, "PCA_OPLS_ScorePlot.pdf", sep="")
	 pdf(D)
	 plot(pc.all$x[,1], pc.all$x[,2], col=col, xlab = xlab, ylab = ylab, xlim = lim, ylim = lim, pch=19, sub = paste("Cumulative Proportion of Variance Explained = ", cum, "%", sep=""), main = "PCA Score Plot on orthogonal-deflated X")
 axis(1, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 axis(2, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 library(car)
 dataEllipse(pc.all$x[,1], pc.all$x[,2], levels = c(0.95), add=TRUE, col = "black", lwd = 0.4, plot.points=FALSE, center.cex=0.2)
 text(pc.all$x[,1], pc.all$x[,2], col=col, cex=0.5, labels=rownames(sorted), pos=1)
 dev.off()
 pca.load <- paste("Loading PC1 (",p.f[1,1],") %")
 pcb.load <- paste("Loading PC2 (",p.f[2,1],") %")
 Max.pc1 = 1.1*(max(pc.all$rotation[,1]))
 Min.pc1 = 1.1*(min(pc.all$rotation[,1]))
 Mpc1=c(Min.pc1*2,Max.pc1*2)
 Max.pc2 = 1.1*(max(pc.all$rotation[,2]))
 Min.pc2 = 1.1*(min(pc.all$rotation[,2]))
 Mpc2=c(Min.pc2*2,Max.pc2*2)
 E = paste(dirout.pca, "PCA_OPLS_LoadingPlot.pdf", sep="")
 pdf(E)
 plot(pc.all$rotation[,1], pc.all$rotation[,2], xlab = pca.load, ylab = pcb.load, xlim = c(Min.pc1,Max.pc1), ylim = c(Min.pc2, Max.pc2), main = "PCA Loading Plot on orthogonal-deflated X", sub = paste("Cumulative Proportion of Variance Explained = ", cum, "%", sep=""))
 axis(1, at=Mpc1, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 axis(2, at=Mpc2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 text(pc.all$rotation[,1], pc.all$rotation[,2], labels=rownames(pc.all$rotation), cex=0.6, pos=1)
 dev.off()
 }
