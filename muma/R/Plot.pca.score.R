Plot.pca.score <-
function(pcx,pcy,scaling) {
 score = paste(getwd(), "/PCA_Data_", scaling, "/PCA_ScoreMatrix.csv", sep="")
 ppppp = paste(getwd(), "/PCA_Data_", scaling, "/PCA_P", sep="")
 Score <- read.csv(score, sep=",", header=TRUE)
 Score.x <- Score[,2:ncol(Score)]
 rownames(Score.x) <- Score[,1]
 pwdK = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.csv", sep="")
 k = read.csv(pwdK)
 Pvar <- read.csv(ppppp, sep=",", header=TRUE)
 Pvar.x <- Pvar[,2:ncol(Pvar)]
 rownames(Pvar.x) <- Pvar[,1]
 pca <- paste("PC",pcx," (",Pvar[pcx,2],") %")
 pcb <- paste("PC",pcy," (",Pvar[pcy,2],")%")
 cum = Pvar[pcx,2] + Pvar[pcy,2]
 xlab = c(pca)
 ylab = c(pcb)
 lim = c()
 max.pc1 = 1.3*(max(abs(Score.x[,pcx])))
 max.pc2 = 1.3*(max(abs(Score.x[,pcy])))
 if (max.pc1 > max.pc2) {lim = c(-max.pc1,max.pc1)} else {lim = c(-max.pc2,max.pc2)}
 tutticolors=matrix(c(1,2,3,4,5,6,7,8,"rosybrown4", "green4", "navy", "purple2", "orange", "pink", "chocolate2", "coral3", "khaki3","thistle","turquoise3","palegreen1","moccasin","olivedrab3","azure4","gold3","deeppink"), ncol=1)
     col=c()
     for(i in 1:nrow(k)) {
      col=c(col, tutticolors[k[i,2],])
     }
 dev.new()
 plot(Score.x[,pcx], Score.x[,pcy], col=col, xlab = xlab, ylab = ylab, xlim = lim, ylim = lim, pch=19, sub = paste("Cumulative Proportion of Variance Explained = ", cum, "%", sep=""), main = paste("PCA Score Plot (", scaling, ")", sep=""))
 axis(1, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 axis(2, at=lim*2, pos=c(0,0), labels=FALSE, col="grey", lwd=0.7)
 library(car)
 dataEllipse(Score.x[,pcx], Score.x[,pcy], levels = c(0.95), add=TRUE, col = "black", lwd = 0.4, plot.points=FALSE, center.cex=0.2)
 text(Score.x[,pcx], Score.x[,pcy], col=col, cex=0.5, labels=rownames(Score.x), pos=1)
 D <- paste(getwd(), "/PCA_Data_", scaling, "/ScorePlot_PC",pcx,"vsPC",pcy,".pdf", sep="")
 dev.copy2pdf(file=D)
 }
