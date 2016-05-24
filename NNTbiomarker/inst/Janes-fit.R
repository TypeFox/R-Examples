Fig1C = png::readPNG("/Users/Roger/Dropbox/_HOME/Paper_on_Biomarker_Studies_and_NNT/Janes-et-al-Figure-clean.png")
Fig1C = Fig1C[ , , 1, drop=T]
dim(Fig1C)
samplesize = 20000
somepoints = data.frame(
  y=sample(x = 1:nrow(Fig1C), size = samplesize, replace=T),
  x=sample(x = 1:ncol(Fig1C), size = samplesize, replace=T)
)


somepoints$
  values = apply(as.matrix(somepoints), 1,
  function(r) Fig1C[nrow(Fig1C)+1-r[1], r[2] ]
  )

summary(somepoints$values)

plot(somepoints$x, somepoints$y,
     cex = 1.1 - somepoints$values,
     pch=".",
     col=1 - (somepoints$values > 0.9832),
     xlab = "cutoff", ylab="DFS% at 5 years")
# To center, we mark the 4 corners, clockwise from lower left:
locationsCorners = locator(4)
xMin =  mean(locationsCorners$x[c(1,2)])
xMax = mean(locationsCorners$x[c(3,4)])
yMin =  mean(locationsCorners$y[c(1,4)])
yMax = mean(locationsCorners$y[c(2,3)])
somepoints$yCentered = (somepoints$y - yMin) / (yMax - yMin)
somepoints$xCentered = (somepoints$x - xMin) / (xMax - xMin)
plot(somepoints$xCentered, somepoints$yCentered,
     cex = 1.1 - somepoints$values,
     pch=".",
     col=1 - (somepoints$values > 0.9832),
     xlim=0:1, ylim=0:1,
  xlab = "cutoff", ylab="DFS% at 5 years")

# Now we mark the two LHS, the center, and the two RHS points.
# (top to bottom)

locationsOnCurves = locator(5)
locationsOnCurves = as.data.frame(locationsOnCurves)

points(locationsOnCurves, pch=as.character(1:5))

xdf = data.frame(x=seq(0,1,length.out = samplesize))
xvalues = xdf[[1]]

lmTreat = lm(y ~ poly(x,2), data = locationsOnCurves[c(2,3,4), ])
predictTreat = predict.lm(lmTreat, newdata = xdf)
lines(xvalues, predictTreat, col="green")

lmWait = lm(y ~ poly(x,2), data = locationsOnCurves[c(1,3,5), ])
predictWait = predict.lm(lmWait, newdata = xdf)
lines(xvalues, predictWait, col="red")

benefit = predictTreat - predictWait


# The following mean benefits are not appropriate:
# the NNT view requires that the patients in the group be indistinguishable;
# but we cannot discard the actual value of the test result.
# For each possible cutoff, what is the mean benefit of treating the positives?
# marginalBenefitPos = sapply(xvalues,
#                             function(cutoff)mean(benefit[xvalues > cutoff]))
# ## The last one is NaN because the condition is all F.
# nntPos = 1/marginalBenefitPos
#
# # For each possible cutoff, what is the mean benefit of treating the negatives?
# marginalBenefitNeg = sapply(xvalues,
#                             function(cutoff)mean(benefit[xvalues <= cutoff]))
# nntNeg = 1/marginalBenefitNeg
# nntNeg = pmin(nntNeg, 1000)
# nntNeg[nntNeg <= 0] = 1000
# xvaluesNeg = xvalues[nntNeg > 0]
# #nntNeg = nntNeg[nntNeg > 0]

NNTlower = 10
NNTupper = 50
plot(xvalues, nntNeg, xlim = c(min(xvaluesNeg), 1), log="y", ylim=c(2,1000))
abline(h=NNTlower)

plot(xvalues, nntPos,
     pch="P",
     xlab="cutoff")
abline(h=c(NNTlower,NNTupper))

plot(c(xvalues, xvalues), c(nntPos, nntNeg),
     pch=c(rep("P", length(xvalues)), rep("N", length(xvalues))),
     xlab="cutoff", log="y",
     xlim = c(min(xvaluesNeg), 1), ylim=c(2,1000))
abline(h=c(NNTlower,NNTupper))

### Where predictWait > predictTreat,
###  we assume that Pr(Treat harms)
#                  = predictWait - predictTreat
#                  = Pr(OK | W but not | T)
#       Pr(Treat helps) = 0
#       Pr(OK regardless) = predictTreat
#       Pr(D  regardless) = 1 - predictWait
### Where predictWait < predictTreat,
###  we assume that Pr(Treat helps)
#                  = predictTreat - predictWait
#                  = Pr(OK | T but not | W)
#       Pr(Treat harms) = 0
#       Pr(OK regardless) = predictWait
#       Pr(D  regardless) = 1 - predictTreat
###
pW = predictWait
pT = predictTreat
probVec = data.frame(OKregardless=pmin(pW, pT),
            TreatHelps=pmax(pT-pW, 0),
            WaitHelps=pmax(pW-pT, 0),
            Dregardless=1 - pmax(pW, pT)
)
plot(xvalues, probVec[,1], ylim=0:1, type ="n")
for(col in 1:4) lines(xvalues, probVec[,col], col=col+1)
legend(0,.7, legend=colnames(probVec), text.col = 1+1:4)

sampleSize = 1e4

janesData = data.frame(
  X = jitter(sample(xvalues, sampleSize, replace = T)))
pW = predict.lm(lmWait, newdata = data.frame(x=janesData$X))
pT = predict.lm(lmTreat, newdata = data.frame(x=janesData$X))
probVec = data.frame(OKregardless = pmin(pW, pT),
                     TreatHelps = pmax(pT-pW, 0),
                     WaitHelps = pmax(pW-pT, 0),
                     Dregardless = 1 - pmax(pW, pT) )
janesData$group = sapply(1:nrow(janesData), function(J) {
  names(which(rmultinom(1, 1, probVec[J, ])[,1]==1))
}
)
janesData$class = (janesData$group=="TreatHelps")
janesData = janesData[order(janesData$X), ]
#lmJonesData = lm(data=janesDataTreat, )

ROCplots(data= janesData, NNTlower=10, NNTupper=40)


RSquantiles = c(0, 11, 20, 28, 41, 112) # From Janes Fig 1 panel C.
quintiles = seq(0, 1, by = 0.2)
RSlm.out = lm(RSquantiles ~ quintiles)
plot(quintiles, RSquantiles)

janesDataSmooth= data.frame(
  X = rep(xvalues, times=2) ,
  class = rep(0:1, each=samplesize),
  weights = c(1-probVec$TreatHelps, probVec$TreatHelps)
)

ROCplots(       whichPlots = "nntRange",
  data= janesDataSmooth,
        NNTlower=10, NNTupper=50)

# "80% overall DFS;   44% marker negative".  Marker negative means marker <
mean(janesDataSmooth$X < 42)

## Conditional NNT's
nnt = 1/benefit
nnt[nnt <= 0] = 1000
xvalues[argmin(nnt,50)]
xvalues[argmin(nnt,10)]
abline(h=xvalues[argmin(nnt,50)], lwd=2)
abline(h=xvalues[argmin(nnt,10)], lwd=2)

# epxloring why this line is in theNNTpos too big region.
nnt[argmin(nnt,10)]


lines(nnt, xvalues, lwd=2)
plot.default(xvalues, benefit)
abline(h=0)
plot(xvalues, nnt, log="y", ylim=c(1,1000), xlab="Recurrence Score quantile", ylab="NNT")
abline(h=10)  ### NNTlower?
abline(h=20)  ### In discomfort zone?
abline(h=50)  ### NNTupper?
