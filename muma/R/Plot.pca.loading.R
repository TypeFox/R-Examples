Plot.pca.loading <-
function(pcx,pcy,scaling) {
 loadi = paste(getwd(), "/PCA_Data_", scaling, "/PCA_LoadingsMatrix.csv", sep="")
 Loading <- read.csv(loadi, sep=",", header=TRUE)
 Loading.x <- Loading[,2:ncol(Loading)]
 rownames(Loading.x) <- Loading[,1]
 ppppp = paste(getwd(), "/PCA_Data_", scaling, "/PCA_P", sep="")
 Pvar <- read.csv(ppppp, sep=",", header=TRUE)
 Pvar.x <- Pvar[,2:ncol(Pvar)]
 rownames(Pvar.x) <- Pvar[,1]
 cum = Pvar[pcx,2] + Pvar[pcy,2]
 pca <- paste("Loadings PC",pcx," (",Pvar[pcx,2],") %")
 pcb <- paste("Loadings PC",pcy," (",Pvar[pcy,2],")%")
 lim.load = c()
 Max.pc1 = 1.1*(max(Loading.x[,pcx]))
 Min.pc1 = 1.1*(min(Loading.x[,pcx]))
 Mpc1=c(Min.pc1*2,Max.pc1*2)
 Max.pc2 = 1.1*(max(Loading.x[,pcy]))
 Min.pc2 = 1.1*(min(Loading.x[,pcy]))
 Mpc2=c(Min.pc2*2,Max.pc2*2)
 dev.new()
 plot(Loading.x[,pcx], Loading.x[,pcy], xlab = pca, ylab = pcb, xlim = c(Min.pc1,Max.pc1), ylim = c(Min.pc2, Max.pc2), main = paste("PCA Loading Plot (", scaling, ")", sep=""), sub = paste("Cumulative Proportion of Variance Explained = ", cum, "%", sep=""))
 axis(1, at=Mpc1, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 axis(2, at=Mpc2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 text(Loading.x[,pcx], Loading.x[,pcy], labels=rownames(Loading.x), cex=0.6, pos=1)
 E = paste(getwd(), "/PCA_Data_", scaling, "/LoadingPlot_PC",pcx,"vsPC",pcy,".pdf", sep="")
 dev.copy2pdf(file=E)
}
