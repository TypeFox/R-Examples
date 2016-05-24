lmerAll <-
function(z){
data <- z
modF <- lmer(y~x + (x|Series) + (x|Ind), data=data, REML=TRUE)
modIDS <- lmer(y~x + (x|Series) + (1|Ind), data=data, REML=TRUE)
modIDInt <- lmer(y~x + (x|Series), data=data, REML=TRUE)

d1<- (-2*logLik(modIDInt)) + (2*logLik(modIDS))
d2<- (-2*logLik(modIDS)) + (2*logLik(modF))

pwI<- if (0.5*pchisq(d1[1], 0, lower.tail=FALSE) + 0.5*pchisq(d1[1], 1, lower.tail=FALSE) < 0.05) 1 else 0
pwS<- if (0.5*pchisq(d2[1], 1, lower.tail=FALSE) + 0.5*pchisq(d2[1], 2, lower.tail=FALSE) < 0.05) 1 else 0

vcs<- VarCorr(modF)
VCVID <- as.vector(vcs$Ind[,])
VCVSeries <- as.vector(vcs$Series[,])
Residuals <- summary(modF)$sigma^2
Intercept <-  summary(modF)$coef[1,1]
Slope <- summary(modF)$coef[2,1]
RInt <- vcs$Ind[1,1]/(vcs$Series[1,1]+vcs$Ind[1,1])
RSlope <- vcs$Ind[2,2]/(vcs$Series[2,2]+vcs$Ind[2,2])
Rrn <- (vcs$Ind[1,1] + vcs$Ind[2,2] + (2*vcs$Ind[2,1]))/((vcs$Series[1,1] + vcs$Series[2,2] + (2*vcs$Series[2,1])) + (vcs$Ind[1,1] + vcs$Ind[2,2] + (2*vcs$Ind[2,1])))
list(VCVSeries, VCVID, pwI, pwS, RInt, RSlope, Rrn, Residuals, Intercept, Slope)
}
