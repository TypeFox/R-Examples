### R code from vignette source 'Correlplot.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Correlplot.Rnw:73-75
###################################################
#install.packages("Correlplot")
library("Correlplot")


###################################################
### code chunk number 2: scatterplot
###################################################
data(PearsonLee)
Mheight <- PearsonLee$Mheight
Dheight <- PearsonLee$Dheight
plot(Mheight,Dheight,asp=1,xlab="Mother's height (cm)",
     ylab="Daughter's height (cm)",pch=19,cex=0.05)


###################################################
### code chunk number 3: biplot
###################################################
X <- cbind(Mheight,Dheight)
n <- nrow(X)
Xt <- scale(X)/sqrt(n-1)
res.svd <- svd(Xt)
Fs <- sqrt(n-1)*res.svd$u
Gp <- res.svd$v%*%diag(res.svd$d)
plot(Fs[,1],Fs[,2],asp=1,pch=19,cex=0.05,xlab="First principal component",
     ylab="Second principal component")
arrows(0,0,3*Gp[,1],3*Gp[,2],col="red",lwd=2)
textxy(3*Gp[,1],3*Gp[,2],colnames(X),cex=1)


###################################################
### code chunk number 4: Correlplot.Rnw:133-136
###################################################
M <- Gp%*%t(Gp)
alpha <- acos(M[1,2])
ang <- alpha*180/pi


###################################################
### code chunk number 5: linangplot
###################################################
lin.out <- linangplot(Mheight,Dheight,cex=0.05)


###################################################
### code chunk number 6: Correlplot.Rnw:167-171
###################################################
library(xtable)
data(students)
R <- cor(students)
xtable(R,digits=3,caption="Correlation matrix for student grades on 5 subjects (Mec=Mecanics,Vec=Vectors,Alg=Algebra,Ana=Analysis,Sta=Statistics).",label="tab:sampleR")


###################################################
### code chunk number 7: pcaplot
###################################################
data(students)
R <- cor(students)
out.eigen <- eigen(R)
V <- out.eigen$vectors
D <- diag(out.eigen$values)
F <- V%*%sqrt(D)
plot(F[,1],F[,2],pch=19,asp=1,xlim=c(-1,1),ylim=c(-1,1))
origin()
arrows(0,0,F[,1],F[,2])
textxy(F[,1],F[,2],colnames(R),cex=1)


###################################################
### code chunk number 8: Correlplot.Rnw:225-227
###################################################
out.pfa <- pfa(students)
L <- out.pfa$La


###################################################
### code chunk number 9: pfaplot
###################################################
plot(L[,1],L[,2],pch=19,asp=1,xlim=c(-1,1),ylim=c(-1,1))
origin()
arrows(0,0,L[,1],L[,2])
text(L[,1],L[,2],colnames(students),cex=1,pos=c(1,2,2,2,3))


###################################################
### code chunk number 10: Correlplot.Rnw:245-246
###################################################
Rhatpfa <- L[,1:2]%*%t(L[,1:2])


###################################################
### code chunk number 11: Correlplot.Rnw:252-257
###################################################
rownames(Rhatpfa) <- rownames(R)
colnames(Rhatpfa) <- colnames(R)
xtable(Rhatpfa,digits=3,
       caption="Least squares approximation, using scalar products, to the correlation matrix obtained by PFA.",
       label="tab:pfahat")


###################################################
### code chunk number 12: Correlplot.Rnw:267-269
###################################################
correlogram(R,labs=colnames(R),main="",
            xlim=c(-1.3,1.3),ylim=c(-1.3,1.3))


###################################################
### code chunk number 13: Correlplot.Rnw:276-277
###################################################
correlogram(R,labs=colnames(R),main="",xlim=c(-1.3,1.3),ylim=c(-1.3,1.3))


###################################################
### code chunk number 14: Correlplot.Rnw:285-287
###################################################
angles <- fit_angles(R)
Rhatcor <- angleToR(angles)


###################################################
### code chunk number 15: Correlplot.Rnw:292-297
###################################################
rownames(Rhatcor) <- rownames(R)
colnames(Rhatcor) <- colnames(R)
xtable(Rhatcor,digits=3,
       caption="Approximation to the correlation matrix obtained by a correlogram.",
       label="tab:corcos")


###################################################
### code chunk number 16: Correlplot.Rnw:306-308
###################################################
correlogram(R,ifun="lincos",labs=colnames(R),main="",
            xlim=c(-1.3,1.3),ylim=c(-1.3,1.3))


###################################################
### code chunk number 17: Correlplot.Rnw:317-319
###################################################
theta_lin <- fit_angles(R)
Rhatcorlin <- angleToR(theta_lin,ifun="lincos")


###################################################
### code chunk number 18: Correlplot.Rnw:324-329
###################################################
rownames(Rhatcorlin) <- rownames(R)
colnames(Rhatcorlin) <- colnames(R)
xtable(Rhatcorlin,digits=3,
       caption="Approximation to the correlation matrix obtained by a linear correlogram.",
       label="tab:corcoslin")


