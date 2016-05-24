### R code from vignette source 'BHH2.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=76)
library(BHH2)


###################################################
### code chunk number 2: PermutationTest
###################################################
# Permutation test for Tomato Data
#cat("Tomato Data (not paired):\n")
data(tomato.data)
attach(tomato.data)
a <- pounds[fertilizer=="A"]
b <- pounds[fertilizer=="B"]
permtest(b,a)
detach()


###################################################
### code chunk number 3: PermutationTest
###################################################
# Permutation test for Boy's Shoes Example
#cat("Shoes Data (paired):\n")
data(shoes.data)
permtest(shoes.data$matB-shoes.data$matA)


###################################################
### code chunk number 4: BHH2.Rnw:142-146
###################################################
par(mar=c(4,1,2,1),mgp=c(2,1,0),cex=0.7)
data(tab03B1)
#stem(tab03B1$yield)
dotPlot(tab03B1$yield,main="Dot plot: Industrial Process Example",xlab="yield")


###################################################
### code chunk number 5: BHH2.Rnw:155-161
###################################################
par(mar=c(4,1,2,1),mgp=c(2,1,0),cex=0.7)
data(tab03B2)
plt <- dotPlot(tab03B2$diff10,xlim=2.55*c(-1,+1),xlab="differences",
    main="Dot plot: Reference Distribution of Differences")
segments(1.3,0,1.3,max(plt$y),lty=2)  #vertical line at x=1.3
text(1.3,max(plt$y),labels=" 1.30",adj=0)


###################################################
### code chunk number 6: BHH2.Rnw:193-198
###################################################
par(mfrow=c(1,1),mar=c(3,1,2,1),cex=0.7)
data(penicillin.data)
penicillin.aov <- aov(yield~blend+treat,data=penicillin.data)
anovaPlot(penicillin.aov,main="Anova plot: Penicillin Manufacturing Example",
    labels=TRUE,cex.lab=0.6)


###################################################
### code chunk number 7: BHH2.Rnw:209-213
###################################################
par(mfrow=c(1,1),mar=c(3,1,2,1),cex=0.7)
data(poison.data)
poison.lm <- lm(y~poison*treat,data=poison.data)
anovaPlot(poison.lm,main="Anova plot: Toxic Agents Example",cex.lab=0.6)


###################################################
### code chunk number 8: BHH2.Rnw:233-238
###################################################
par(mfrow=c(1,1),mar=c(3,1,2,1),cex=0.7)
data(corrosion.data)
corrosion.aov <- aov(resistance~heats+run+coating+heats:coating,data=corrosion.data)
anovaPlot(corrosion.aov,main="Anova plot: Corrosion Resistance Example",
    cex.lab=0.6)


###################################################
### code chunk number 9: BHH2.Rnw:260-265
###################################################
par(mar=c(4,3,2,1),mgp=c(2,1,0),cex=0.7)
# Lambda Plot tracing F values.
data(poison.data)
lambdaPlot(poison.lm,lambda=seq(-2,1,by=.1),stat="F",global=FALSE,cex=0.6,
    main="Lambda Plot: Toxic Agents Example")


###################################################
### code chunk number 10: BHH2.Rnw:280-286
###################################################
# Lambda Plot tracing t values.
par(mar=c(4,3,2,1),mgp=c(2,1,0),cex=0.7)
data(woolen.data)
woolen.lm <- lm(y~x1+x2+x3+I(x1^2)+I(x2^2)+I(x3^2)+I(x1*x2)+I(x1*x3)+I(x2*x3)+I(x1*x2*x3),data=woolen.data)
lambdaPlot(woolen.lm,main="Lambda plot: Woolen Thread Example (2nd order model)",
    stat="t",cex=0.6)


###################################################
### code chunk number 11: BHH2.Rnw:294-298
###################################################
par(mar=c(4,3,2,1),mgp=c(2,1,0),cex=0.7)
# Lambda Plot tracing F values.
lambdaPlot(lm(y~x1+x2+x3,data=woolen.data),lambda=seq(-1,1,length=31),
    main="Lambda plot: Woolen Thread Example (1st order model)",stat="F",global=TRUE)


###################################################
### code chunk number 12: DesignMatrix
###################################################
print(X <- ffDesMatrix(5,gen=list(c(-5,1,2,3,4))))


###################################################
### code chunk number 13: DesignMatrix
###################################################
ffFullMatrix(X,x=c(1,2,3,4),maxInt=2,blk=X[,5])$Xa


###################################################
### code chunk number 14: Subsets
###################################################
subsets(n=5,r=3,v=c("x","y","z","A","B"))


