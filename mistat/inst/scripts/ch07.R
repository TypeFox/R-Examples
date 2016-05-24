###################################################
### Chap07Start
###################################################
library(mistat)


###################################################
### AcceptanceSampling01
###################################################
library(AcceptanceSampling)

as.data.frame(
  find.plan(PRP=c(0.01, 0.05), 
            CRP=c(0.08, 0.05), 
            type="hypergeom", N=100))


###################################################
### TableAcceptanceSamplingN100AB05
###################################################
p0 <- rep(c(0.01, 0.03), each=10)
pt <- rep(seq(0.05, to=0.32, by=0.03), 2)

res <-as.data.frame(
  find.plan(PRP=c(p0[1], 0.05), 
          CRP=c(pt[1], 0.05), 
          type="hypergeom", N=100))
for(i in 2:20){
  res <- rbind(res,
               find.plan(PRP=c(p0[i], 0.05), 
                         CRP=c(pt[i], 0.05), 
                         type="hypergeom", N=100))
}
res
rm(res)


###################################################
### TableAcceptanceSamplingN100A10B20
###################################################
p0 <- rep(c(0.01, 0.03), each=10)
pt <- rep(seq(0.05, to=0.32, by=0.03), 2)

res <-as.data.frame(
  find.plan(PRP=c(p0[1], 0.20), 
          CRP=c(pt[1], 0.10), 
          type="hypergeom", N=100))
for(i in 2:20){
  res <- rbind(res,
               find.plan(PRP=c(p0[i], 0.1), 
                         CRP=c(pt[i], 0.2), 
                         type="hypergeom", N=100))
}
res
rm(res, i, p0, pt)


###################################################
### PlotOCpN100n50c1
###################################################
OC <- OC2c(n=50, 
           c=1, 
           type="hypergeom", 
           N=100, 
           pd=seq(0, 0.15, 
                  length.out=100))

plot(OC, type="l")

rm(OC)


###################################################
### TableOCpN100n50c1
###################################################
OC <- OC2c(n=50, 
           c=1, 
           type="hypergeom", 
           N=100, 
           pd=seq(0, 0.15, 
                  length.out=100))


data.frame(p=round(OC@pd[seq(1, 100, length.out=20)], 3),
  OCp=round(OC@paccept[seq(1, 100, length.out=20)], 4))

rm(OC)


###################################################
### PlotOCsequential
###################################################
library(Dodge)

PlanDesign <- SeqDesignBinomial(AQL=0.01, 
                                alpha=0.05, 
                                LQL=0.05, 
                                beta=0.05, 
                                Plots=FALSE)

invisible(
  SequentialBinomial(
    PlanDesign, 
    Plots=TRUE))

layout(1)
rm(PlanDesign)


###################################################
### PlotAOQN1000n250c5AOQ
###################################################
invisible(
  SSPlanBinomial(N=1000, 
                 n=250, 
                 Ac=5, 
                 p=seq(0.005, 0.035, 
                       length.out=200), 
                 Plots=TRUE))

layout(1)


###################################################
### PlotAOQN1000n250c5ATI
###################################################
invisible(
  SSPlanBinomial(N=1000, 
                 n=250, 
                 Ac=5, 
                 p=seq(0.005, 0.035, 
                       length.out=200), 
                 Plots=TRUE))

layout(1)


###################################################
### Chap07End
###################################################
detach(package:AcceptanceSampling)
detach(package:Dodge)
detach(package:mistat)
