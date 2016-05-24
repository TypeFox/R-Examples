surface <- ice[ice$Location=='surface',c("Treatment","T1930")]
intra <- ice[ice$Location=='intramuscular',"T1930"]
newdata <- cbind(surface,intra)
names(newdata) <- c('Treatment','SurfTemp','IntraTemp')
anova(lm(SurfTemp - IntraTemp ~ Treatment, newdata))
