library(RadOnc); data(RadOnc)

## generate plots for a single DVH
stomach <- janedoe[["STOMACH"]]
par(mfrow=c(2,1))
plot(stomach,main="Stomach (Cumulative DVH)", type="cumulative")
points(1:5*1000, stomach[paste("V", 1:5*1000, "cGy", sep="")])
text(1:5*1000, stomach[paste("V", 1:5*1000, "cGy", sep="")],stomach[paste("V", 1:5*1000, "cGy", sep="")],pos=4)
plot(stomach,main="Stomach (Differential DVH)", type="differential", volume="absolute")

## generate plots for multiple DVHs
par(mfrow=c(1,1))
plot(janedoe,plot.type="individual",col=rainbow(length(janedoe)),main="Individual DVHs")
plot(janedoe[c("CTV","PTV")],johndoe[c("LIVER","STOMACH","DUODENUM")],plot.type="grouped",col=c("red","blue"),density=c(NA,20),angle=c(NA,90),width="sd",multiplier=0.5,main="Grouped DVHs")

## generate plots for statistical comparison among two groups of DVHs
plot(janedoe[c("CTV","PTV","DUODENUM")],johndoe[c("LIVER","STOMACH","SMALL_BOWEL","BODY",".*KIDNEY")],plot.type="wilcox",col=c("red","blue"),panel.lower="grouped",width="mad",multiplier=0.5,main="Wilcox Comparison")
