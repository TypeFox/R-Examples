## ----message=FALSE-------------------------------------------------------
library(lctools)
data(GR.Municipalities)
names(GR.Municipalities@data)

## ------------------------------------------------------------------------
myDF<-cbind(GR.Municipalities@data$Income01,GR.Municipalities@data$X, GR.Municipalities@data$Y)
myDF[!complete.cases(myDF),]
myDF.new<-na.omit(myDF)
bw<-12
wt<-'Binary'
sp.G<-spGini(myDF.new[,2:3],bw,myDF.new[,1],wt)
knitr::kable(sp.G, format = "pandoc", digits = 5)

## ----warning = FALSE-----------------------------------------------------
spGini.Sim20<-mc.spGini(Nsim=19,bw,myDF.new[,1],myDF.new[,2],myDF.new[,3],wt)
spGini.Sim20$pseudo.p

