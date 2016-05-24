###################################################
### Chap11Start
###################################################
library(mistat)


###################################################
### Dpvalue01
###################################################
X <- c(1.1, 0.3, -0.7, -0.1)

M <- 200

set.seed(123)

Di <- matrix(sample(x=c(-1,1), 
                    size=length(X)*M, 
                    replace=TRUE), 
             nrow=M)

Xi <- matrix(X, 
             nrow=M, 
             ncol=length(X), 
             byrow=TRUE)

sum(rowMeans(Di*Xi) >= mean(X))/M

rm(X, M, Di, Xi)


###################################################
### Dpvalue02
###################################################
X <- c(0.8, 0.6, 0.3, -0.1, 1.1, -0.2, 0.3, 0.5, 0.5, 0.3)

M <- 200

set.seed(123)

Di <- matrix(sample(x=c(-1,1), 
                    size=length(X)*M, 
                    replace=TRUE), 
             nrow=M)

Xi <- matrix(X, 
             nrow=M, 
             ncol=length(X), 
             byrow=TRUE)

Means <- rowMeans(Di*Xi)

sum(rowMeans(Di*Xi) >= mean(X))/M

stem(Means)

rm(Di, Xi, M, X, Means)


###################################################
### AnovaHADPAS
###################################################
data(HADPAS)                         

HADPAS$diska <- as.factor(HADPAS$diska)
HADPAS$hyb <- as.factor(HADPAS$hyb)  

AovH <- aov(res3 ~ diska + hyb,      
            data=HADPAS)                

summary(AovH)                        

tail(confint(AovH), 5)               

rm(AovH)


###################################################
### PlotEffectKeyboard
###################################################
data(KEYBOARDS)

boxplot(errors ~ keyboard, data=KEYBOARDS, ylab="Errors")


###################################################
### PlotEffectKeyboard
###################################################
boxplot(errors ~ job, data=KEYBOARDS, ylab="Errors")


###################################################
### PlotEffectKeyboard
###################################################
boxplot(errors ~ typist, data=KEYBOARDS, ylab="Errors")


###################################################
### FullFactorialPistonSim
###################################################
library(DoE.base)

Factors <- list(
  m=c(30, 45, 60), 
  k=c(1500, 3000, 4500))

FacDesign <- fac.design(
  factor.names=Factors,
  randomize=TRUE,
  replications=5,
  repeat.only=TRUE)

Levels <- data.frame(
  lapply(
    lapply(FacDesign, 
           as.character), 
    as.numeric),
  s=0.01, 
  v0=0.005,
  p0=95000,
  t=293,
  t0=350)

Ps <- pistonSimulation(m=Levels$m,
                       s=Levels$s,
                       v0=Levels$v0,
                       k=Levels$k,
                       p0=Levels$p0,
                       t=Levels$t,
                       t0=Levels$t0,
                       each=1,
                       seed=123)

FacDesign <- add.response(
  design=FacDesign, 
  response=Ps$seconds)

summary(
  aov(Ps.seconds ~ m*k, 
      data=FacDesign))

rm(Levels, Factors)


###################################################
### PlotEffectSpringCoeffOnCycleTime
###################################################
boxplot(seconds ~ k, 
        data=Ps)

rm(Ps)


###################################################
### PlotInteractionPlotOnCycleTime
###################################################
with(FacDesign,
     interaction.plot(
       x.factor=m, 
       trace.factor=k, 
       response=Ps.seconds,
       type="b",
       pch=15:18))

rm(FacDesign)


###################################################
### FractionalFactorialFrF2
###################################################
library(FrF2)
FrF2(nfactors=5, resolution=5)


###################################################
### FullFactorialDesign
###################################################
Design <- fac.design(nlevels=2, 
                     nfactors=5)

head(Design, 3)

tail(Design, 3)

rm(Design)


###################################################
### FactorialDesign5Factors
###################################################
Factors <- list(
  m=c(30, 60), 
  s=c(0.005, 0.02), 
  v0=c(0.002, 0.01),
  k=c(1000, 5000),
  t=c(290, 296))

FacDesign <- fac.design(
  factor.names=Factors,
  randomize=TRUE,
  replications=5,
  repeat.only=TRUE)

Levels <- data.frame(
  lapply(
    lapply(FacDesign, as.character), 
    as.numeric),
  p0=90000,
  t0=340, stringsAsFactors=F)

Ps <- pistonSimulation(m=Levels$m,
                       s=Levels$s,
                       v0=Levels$v0,
                       k=Levels$k,
                       p0=Levels$p0,
                       t=Levels$t,
                       t0=Levels$t0,
                       each=1,
                       seed=123)

FacDesign <- add.response(
  design=FacDesign, 
  response=Ps$seconds)

summary(
  aov(Ps.seconds ~ (m+s+v0+k+t)^2, 
      data=FacDesign))

rm(Levels, Ps, Factors)


###################################################
### PlotMainEffectsOnCycleTime5Factors
###################################################
MEPlot(obj=FacDesign, 
       main="")


###################################################
### PlotInteractionOnCycleTime5Factors
###################################################
IAPlot(obj=FacDesign, 
       main="")

rm(FacDesign)


###################################################
### FullFactorial3n
###################################################
data(STRESS)                                   

summary(                                       
  lm(stress ~ (A+B+C+I(A^2)+I(B^2)+I(C^2))^3,  
     data=STRESS))                             

summary(                                       
  aov(stress ~ (A+B+C)^3 +I(A^2)+I(B^2)+I(C^2),
      data=STRESS))                            


###################################################
### PlotMainEffectsOnStress
###################################################
Stress2 <- data.frame(
  lapply(STRESS[,!names(STRESS) %in% "stress"], 
         as.factor), 
  Stress=STRESS$stress)                          

layout(matrix(1:4, 2, byrow=TRUE))

with(Stress2,
     interaction.plot(
       x.factor=A, 
       trace.factor=rep(0, length(A)), 
       response=Stress, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(160, 280)))

with(Stress2,
     interaction.plot(
       x.factor=B, 
       trace.factor=rep(0, length(A)), 
       response=Stress, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(160, 280)))

with(Stress2,
     interaction.plot(
       x.factor=C, 
       trace.factor=rep(0, length(A)), 
       response=Stress, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(160, 280)))

layout(1)


###################################################
### PlotInteractionOnStress
###################################################
layout(matrix(1:4, 2, byrow=TRUE))

with(Stress2,
     interaction.plot(
       x.factor=A, 
       trace.factor=B, 
       response=Stress, 
       legend=TRUE,
       type="b",
       pch=15:18,
       ylim=c(160, 340)))

with(Stress2,
     interaction.plot(
       x.factor=A, 
       trace.factor=C, 
       response=Stress, 
       legend=TRUE,
       type="b",
       pch=15:18,
       ylim=c(160, 340)))

with(Stress2,
     interaction.plot(
       x.factor=B, 
       trace.factor=C, 
       response=Stress, 
       legend=TRUE,
       type="b",
       pch=15:18,
       ylim=c(160, 340)))

layout(1)

rm(Stress2)


###################################################
### FullFact8Block4
###################################################
Gen <- matrix(c(
  0,1,1,1,1,0,0,0,
  1,0,1,1,0,1,0,0), 
              nrow=2, 
              byrow=TRUE)

head(
  fac.design(nlevels=2, 
             nfactors=8, 
             blocks=4, 
             block.gen=Gen))

rm(Gen)


###################################################
### ResponseSurface
###################################################
library(rsm)

s  <- c(0.0075, 0.01, 0.0125, 0.015, 0.0175)

v0 <- c(0.0050, 0.00625, 0.0075, 0.00875, 0.0100)

k  <- c(1000, 2000, 3000, 4000, 5000)

t0 <- c(340, 345, 350, 355, 360)

Ccd <- ccd(basis=4, 
           n0=4, 
           alpha=2, 
           coding=list(x1 ~ -5 + s*400,
                       x2 ~ -6 + v0*800,
                       x3 ~ -3 + k*0.001,
                       x4 ~ -70 + t0*0.2),
           randomize=FALSE)

head(Ccd)

Levels <- as.data.frame(
  decode.data(Ccd))[, c("s", "v0", "k", "t0")]

Ps <- pistonSimulation(m=rep(60, nrow(Levels)),
                       s=Levels$s,
                       v0=Levels$v0,
                       k=Levels$k,
                       p0=rep(110000, nrow(Levels)),
                       t=rep(296, nrow(Levels)),
                       t0=Levels$t0,
                       each=30,
                       seed=123)

Ps <- simulationGroup(Ps, 30)

Ccd$meantime <- aggregate(Ps["seconds"], 
                          by=Ps["group"], 
                          FUN=mean)$seconds

Rsm <- rsm(meantime ~ SO(x1, x2, x3, x4), 
           data= Ccd)

summary(Rsm)

rm(k, s, t0, v0)


###################################################
### PlotMainEffectsResponseSurface
###################################################
layout(matrix(1:4, 2, byrow=TRUE))

with(Ccd,
     interaction.plot(
       x.factor=x1, 
       trace.factor=rep(0, length(x1)), 
       response=meantime, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(0.45, 0.75)))

with(Ccd,
     interaction.plot(
       x.factor=x2, 
       trace.factor=rep(0, length(x2)), 
       response=meantime, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(0.45, 0.75)))

with(Ccd,
     interaction.plot(
       x.factor=x3, 
       trace.factor=rep(0, length(x3)), 
       response=meantime, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(0.45, 0.75)))

with(Ccd,
     interaction.plot(
       x.factor=x4, 
       trace.factor=rep(0, length(x4)), 
       response=meantime, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(0.45, 0.75)))

layout(1)


###################################################
### PlotResponseSurface
###################################################
contour(Rsm, ~ SO(x1+x3), 
        image = FALSE)


###################################################
### SteepestDescent
###################################################
steepest(Rsm, 
         dist=seq(0, 2.5, by=0.5), 
         descent=TRUE)


###################################################
### ResponseSurfaceCanonical
###################################################
canonical.path(Rsm, 
               dist=seq(0, 2.5, by=0.5), 
               descent=TRUE)


###################################################
### Chap11End
###################################################
rm(HADPAS, KEYBOARDS, Levels, Ps, STRESS, Ccd, Rsm)
detach(package:rsm)
detach(package:FrF2)
detach(package:DoE.base)
detach(package:mistat)
