## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----loadlibrary---------------------------------------------------------
library(DLMtool)

## ----assignObjs----------------------------------------------------------
for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]])

## ----sfinit, eval=FALSE--------------------------------------------------
#  sfInit(parallel=TRUE, cpus=2)
#  

## ---- eval=FALSE---------------------------------------------------------
#  sfExportAll()

## ----seed----------------------------------------------------------------
set.seed(1) 

## ------------------------------------------------------------------------
OM <- new('OM', Blue_shark, Generic_fleet, Imprecise_Biased)

## ------------------------------------------------------------------------
slotNames(OM)

## ---- eval=FALSE---------------------------------------------------------
#  class?OM

## ------------------------------------------------------------------------
avail('DLM_output')

## ----MPs-----------------------------------------------------------------
MPs <- c("Fratio", "DCAC", "Fdem", "DD")    


## ---- eval=FALSE---------------------------------------------------------
#  ?Fratio
#  ?DBSRA

## ------------------------------------------------------------------------
Fratio

## ----SNAPMSE, eval=FALSE-------------------------------------------------
#  SnapMSE <- runMSE(OM, MPs, nsim=16, reps=1, proyears=30, interval=10)

## ---- fig.show='asis', fig.width=7, fig.height=7-------------------------
plot(SnapMSE)

## ------------------------------------------------------------------------
avail('DLM_data')

## ------------------------------------------------------------------------
slotNames(China_rockfish)

## ---- eval=FALSE---------------------------------------------------------
#  class?DLM_data

## ------------------------------------------------------------------------
Can(China_rockfish)
Cant(China_rockfish)
Needed(China_rockfish)

## ---- fig.width=7--------------------------------------------------------
RockReal<-TAC(China_rockfish)
plot(RockReal)

## ---- fig.width=7, fig.height=7------------------------------------------
RockReal <- Sense(RockReal,"DCAC")

## ------------------------------------------------------------------------
avail('Stock')

## ------------------------------------------------------------------------
ourstock <- Snapper

## ------------------------------------------------------------------------
ourstock@M <- c(0.2,0.25)
ourstock@maxage <- 18
ourstock@D <- c(0.05,0.3)
ourstock@Frac_area_1 <- c(0.05,0.15)
ourstock@Prob_staying <- c(0.4,0.99)

## ------------------------------------------------------------------------
ourfleet <- Generic_FlatE
ourfleet@Vmaxlen <- c(0.5, 1)

## ------------------------------------------------------------------------
ourOM <- new('OM',ourstock, ourfleet, Imprecise_Biased)

## ----availOut------------------------------------------------------------
avail("DLM_output")

## ----availIn-------------------------------------------------------------
avail("DLM_output")

## ----chooseMPs-----------------------------------------------------------
MPs <- c("BK", "CC1", "CompSRA", "DBSRA", "DBSRA4010", "DCAC", "DCAC4010", "DepF", "DynF",
         "EDCAC", "Fratio", "Itarget1", "Itarget4", "MCD", "MCD4010", "SBT1")

## ---- eval=FALSE---------------------------------------------------------
#  ourMSE <- runMSE(ourOM, MPs=MPs, proyears=20, interval=5, nsim=16,reps=1)

## ---- fig.width=7, fig.height=7------------------------------------------
Tplot(ourMSE)

## ------------------------------------------------------------------------
Results <- summary(ourMSE) 
head(Results)
Targetted <- subset(Results, Results$Yield>30 & Results$POF<50 & Results$P50<20)
Targetted

## ----ourMSE2, eval=FALSE-------------------------------------------------
#  TargMP <- Targetted$MP[grep("FMSYref",Targetted$MP,invert=T)]
#  ourMSE2 <- runMSE(ourOM, TargMP, proyears=20, interval=5, nsim=16, reps=1)

## ----Sub, eval=TRUE, echo=FALSE------------------------------------------
TargMP <- Targetted$MP[grep("FMSYref",Targetted$MP,invert=T)]
ourMSE2 <- Sub(ourMSE, MPs=TargMP)

## ---- fig.width=7, fig.height=5------------------------------------------
CheckConverg(ourMSE2)

## ---- fig.width=7, fig.height=7------------------------------------------
Pplot(ourMSE2)

## ---- fig.width=7, fig.height=5------------------------------------------
Kplot(ourMSE2)

## ---- fig.width=7, fig.height=3.5----------------------------------------
Tplot2(ourMSE2)

## ---- TradePlot, fig.width=7, fig.height=6-------------------------------
TradePlot(ourMSE2, XThresh=c(0,0), YThresh=c(0,0), ShowLabs = TRUE)

## ---- fig.width=7, fig.height=6------------------------------------------
VOI(ourMSE2)

## ---- fig.width=7, fig.height=3.5----------------------------------------
summary(ourReefFish)

## ----ourReefFish, eval=FALSE---------------------------------------------
#  ourReefFish <- TAC(ourReefFish)

## ----Load, echo=FALSE----------------------------------------------------
data(ourReefFish)

## ----plotOurReefFish, fig.width=7, fig.height=7--------------------------
plot(ourReefFish)

## ------------------------------------------------------------------------
ourOM@Btcv
ourOM@Btbias

## ---- fig.width=7, fig.height=7------------------------------------------
Tplot2(ourMSE2)

## ----sense, fig.width=7, fig.height=7------------------------------------
ourReefFish <- Sense(ourReefFish,'MCD')

## ------------------------------------------------------------------------
sapply(1,Fdem_CC,Red_snapper,reps=5)

## ------------------------------------------------------------------------
AvC <-function(x, DLM_data, reps)rlnorm(reps, log(mean(DLM_data@Cat[x,], na.rm=T)), 0.1) 

## ------------------------------------------------------------------------
class(AvC) <-"DLM_output"
environment(AvC) <-asNamespace('DLMtool')

## ---- eval=FALSE---------------------------------------------------------
#  sfExport("AvC")

## ------------------------------------------------------------------------
THC<-function(x,DLM_data, reps){
  rlnorm(reps,log(DLM_data@Cat[x,order(DLM_data@Cat[x,],decreasing=T)[3]]),0.1)
}
class(THC)<-"DLM_output"
environment(THC) <- asNamespace('DLMtool')


## ---- eval=FALSE---------------------------------------------------------
#  sfExport("THC")

## ------------------------------------------------------------------------
matlenlim <- function (x, DLM_data, ...) {
    dependencies = "DLM_data@LFC, DLM_data@LFS"
    Allocate <- 1
    Effort <- 1
    Spatial <- c(1, 1)
    newLFC <- DLM_data@L50[x] * 0.95
    newLFS <- DLM_data@L50[x]
    Vuln <- c(newLFC, newLFS)
    c(Allocate, Effort, Spatial, Vuln)
}
class(matlenlim) <- "DLM_input"
environment(matlenlim) <- asNamespace("DLMtool")


## ---- eval=FALSE---------------------------------------------------------
#  sfExport("matlenlim")

## ------------------------------------------------------------------------
area1_50<-function(x,DLM_data, ...){ 
  Allocate<-0 # Fraction of effort reallocated to open area
  Effort<-1  # Fraction of effort in last historical year
  Spatial<-c(0.5,1) # Fraction of effort found in each area
  Vuln<-rep(NA,2) # Length vulnerability is not specified   
  c(Allocate, Effort, Spatial, Vuln) # Input controls stitched togther
}
class(area1_50)<-"DLM_input"
environment(area1_50) <- asNamespace('DLMtool')

## ---- eval=FALSE---------------------------------------------------------
#  sfExport("area1_50")

## ---- eval=FALSE---------------------------------------------------------
#  new_MPs <- c("AvC","THC","matlenlim","area1_50")
#  OM <- new('OM',Porgy, Generic_IncE, Imprecise_Unbiased)
#  PorgMSE <- runMSE(OM,new_MPs,maxF=1,nsim=20,reps=1,proyears=20,interval=5)

## ---- fig.width=7, fig.height=7------------------------------------------
Tplot(PorgMSE)  

## ---- eval=FALSE---------------------------------------------------------
#  OM@D
#  OM@D <- c(0.05,0.3)
#  PorgMSE2 <- runMSE(OM, new_MPs, maxF=1, nsim=20, reps=1, proyears=20, interval=5)

## ----SubPorg, eval=TRUE, echo=FALSE--------------------------------------
ind <- which(PorgMSE@OM$D < 0.3)
PorgMSE2 <- Sub(PorgMSE, sims=ind)


## ---- fig.width=7,fig.height=7-------------------------------------------
Tplot(PorgMSE2)  

## ----]-------------------------------------------------------------------
slotNames('DLM_data')

## ------------------------------------------------------------------------
DLMDataDir()

## ------------------------------------------------------------------------
Madeup<-new('DLM_data')                             #  Create a blank DLM object
Madeup@Name<-'Test'                                 #  Name it
Madeup@Cat<-matrix(20:11*rlnorm(10,0,0.2),nrow=1)   #  Generate fake catch data
Madeup@Units<-"Million metric tonnes"               #  State units of catch
Madeup@AvC<-mean(Madeup@Cat)                        #  Average catches for time t (DCAC)
Madeup@t<-ncol(Madeup@Cat)                          #  No. yrs for Av. catch (DCAC)
Madeup@Dt<-0.5                                      #  Depletion over time t (DCAC)
Madeup@Dep<-0.5                                     #  Depletion relative to unfished 
Madeup@vbK<-0.2                                     #  VB maximum growth rate
Madeup@vbt0<-(-0.5)                                 #  VB theoretical age at zero length
Madeup@vbLinf<-200                                  #  VB maximum length
Madeup@Mort<-0.1                                    #  Natural mortality rate
Madeup@Abun<-200                                    #  Current abundance
Madeup@FMSY_M<-0.75                                 #  Ratio of FMSY/M
Madeup@L50<-100                                     #  Length at 50% maturity
Madeup@L95<-120                                     #  Length at 95% maturity
Madeup@BMSY_B0<-0.35                                #  BMSY relative to unfished

## ---- fig.width=7, fig.height=3.5----------------------------------------
summary(Atlantic_mackerel)

## ------------------------------------------------------------------------
Can(Atlantic_mackerel)
Cant(Atlantic_mackerel)
Needed(Atlantic_mackerel)

## ---- eval=TRUE----------------------------------------------------------
Atlantic_mackerel <- TAC(Atlantic_mackerel,reps=48)

## ---- fig.width=6, fig.height=8------------------------------------------
plot(Atlantic_mackerel)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  sfStop()

