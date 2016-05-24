library(metaplus)

# speed execution this uses just 99 bootstrap replications
# for final results this should be changed to 999
nosims <- 99

mag1 <- metaplus(yi,sei,slab=study,data=mag)

summary(mag1)

plot(mag1)

# can modify plot to obtain Odds Ratios
plot(mag1,atransf=exp, at=log(c(.01,.1,1,10,100)),xlab="Odds Ratio")

mag2 <- metaplus(yi,sei,slab=study,random="t-dist",data=mag)

summary(mag2)

summary(testOutliers(mag2,R=nosims))

mag3 <- metaplus(yi,sei,slab=study,random="mixture",data=mag)
 
summary(mag3)

summary(testOutliers(mag3,R=nosims))

cdp1 <- metaplus(yi,sei,plotci=TRUE,slab=study,data=cdp)
summary(cdp1)

plot(cdp1)

cdp2 <- metaplus(yi,sei,plotci=TRUE,slab=study,random="t-dist",data=cdp)
summary(cdp2)

summary(testOutliers(cdp2,R=nosims))

cdp3 <- metaplus(yi,sei,plotci=TRUE,slab=study,random="mixture",data=cdp)
summary(cdp3)

cdp3.outlierProbs <- outlierProbs(cdp3)
plot(cdp3.outlierProbs)

summary(testOutliers(cdp3,R=nosims))

plot(cdp1,extrameta=list(cdp2,cdp3))

marinho1 <- metaplus(meaneffect,seeffect,plotci=TRUE,slab=study,data=marinho)
summary(marinho1)

marinho2 <- metaplus(meaneffect,seeffect,plotci=TRUE,slab=study,random="t-dist",data=marinho)
summary(marinho2)

marinho3 <- metaplus(meaneffect,seeffect,plotci=TRUE,slab=study,random="mixture",data=marinho)
summary(marinho3)

plot(marinho1,extrameta=list(marinho2,marinho3))

marinho3.outlierProbs <- outlierProbs(marinho3)
plot(marinho3.outlierProbs)

summary(testOutliers(marinho3,R=nosims))

exercise1 <- metaplus(smd,sqrt(varsmd),mods=duration,plotci=TRUE,slab=study,data=exercise)
summary(exercise1)

exercise2 <- metaplus(smd,sqrt(varsmd),mods=duration,plotci=TRUE,random="t-dist",data=exercise)
summary(exercise2)

summary(testOutliers(exercise2,R=nosims))

exercise3 <- metaplus(smd,sqrt(varsmd),mods=eduration,plotci=TRUE,random="mixture",data=exercise)
summary(exercise3)

summary(testOutliers(exercise3,R=nosims))

exercise3.outlierProbs <- outlierProbs(exercise3)
plot(exercise3.outlierProbs)

exercise$duration4 <- exercise$duration-4
exercise$duration8 <- exercise$duration-8
exercise$duration12 <- exercise$duration-12

exercise.wk4 <- metaplus(smd,sqrt(varsmd),mods=duration4,plotci=TRUE,label="Random Mixture (Week 4)",
                         slab=study,random="mixture",data=exercise)
exercise.wk8 <- metaplus(smd,sqrt(varsmd),mods=duration8,plotci=TRUE,label="Random Mixture (Week 8)",
                         slab=study,random="mixture",data=exercise)
exercise.wk12 <- metaplus(smd,sqrt(varsmd),mods=duration12,plotci=TRUE,label="Random Mixture (Week 12)",
                          slab=study,random="mixture",data=exercise)

exercise.nodurn <- metaplus(smd,sqrt(varsmd),plotci=TRUE,label="Random Mixture (No Duration)",
                            slab=study,random="mixture",data=exercise)

plot(exercise.nodurn,extrameta=list(exercise.wk4,exercise.wk8,exercise.wk12))
