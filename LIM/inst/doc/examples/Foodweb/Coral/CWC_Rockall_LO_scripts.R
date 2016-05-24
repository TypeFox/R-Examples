

#########################################################################
##                                                                     ##
## INVERSE FOODWEB MODELLING - COLDWATER CORAL WEB                     ##
##    Dick van Oevelen  -  Karline Soetaert                            ##
##                                                                     ##
##    This file contains the runs as used in the                       ##
##    Limnology and Oceanography paper                                 ##
##                                                                     ##
##                                                                     ##
## van Oevelen, Dick, Gerard Duineveld, Marc Lavaleye, Furu Mienis,    ##
## Karline Soetaert, and Carlo H. R. Heip, 2009.                       ##
## The cold-water coral community as hotspot of carbon cycling on      ##
## continental margins: A food web analysis from Rockall Bank          ##
## (northeast Atlantic). Limnology and Oceangraphy 54:1829-1844.       ##
##                                                                     ##
##                                                                     ##
##    do not run this script if you do not have plenty of time !       ##
##                                                                     ##
## Scripts written by Dick van Oevelen, with some tweaking by Karline  ##
##                                                                     ##
## Before trying to run these macros, make sure to set the working     ##
## directory to the one containing this file                           ##
#########################################################################

require(LIM)

## =============================================================================
## Solve the cold-water coral food web
## =============================================================================
# Name of the file with the inverse model
# File      <- "CWC_Rockall4 LOrev.input"
# LIM       <- Setup(file=File)

LIM       <- LIMCoralRockall

# The dimensionality of the problem:
LIM[c("NUnknowns", "NEquations", "NConstraints", "NComponents",
      "NExternal", "NVariables")]

SollseiDef<- Lsei(LIM,parsimonious = TRUE)

plotweb(Flowmatrix(LIM))
RangesDef <- Xranges(LIM)

# This takes about 3 seconds for ONE iteration!
# For the paper we did 10000... here 500 are saved in a file
# By notrun, it is ignored
notrun <- function ()
    SolBayDef <- Xsample(lim  = LIM,
                     exact= NULL,
                     iter = 500,
                     type = "mirror",
                     jmp  = c(RangesDef[,2]-RangesDef[,1])/4,
                     x0   = SollseiDef$X)
# 500 runs are saved in a file
load(file = "MCR.Rdata")

## =============================================================================
## Markov chain sampling
## =============================================================================

# This takes about 3 seconds for ONE iteration!
# For the paper we did 10000... here 500
# By function notrun, it is ignored
notrun <- function ()
    SolBayDef <- Xsample(lim  = LIM,
                     exact= NULL,
                     iter = 500,
                     type = "mirror",
                     jmp  = c(RangesDef[,2]-RangesDef[,1])/4,
                     x0   = SollseiDef$X)
# load 500 samples instead
load(file = "MCR.Rdata")
## =============================================================================
## Run the sensitivity analysis
## =============================================================================

## Again a very long simulation - skip it - but then fig 4 cannot be displayed..
# set counters
notrun <- function ()  {
 i <- 1 # counts the succesful initializations
 j <- 1 # counts the initialization attempts (not all are succesful)
 SolBaySens <- NULL # this is where the runs will be saved
 Pars       <- NULL # this is where the feasible parameter values will be stored
 PARS       <- c("FracCWC","FracEUN","FracHES","FracSPO","FracHYD","FracCRI",
                "FracPOL","FracBIV","FracLIM","FracASP","FracSUS",
                "FracOMN","FracCRA","FracURC","FracSTA","FracFIS")
                # list with parameters that are initialized
 while(i <= 200) #run until 200 succesful initializations are achieved
  {
  readLIM$pars$val[which(readLIM$pars$name %in% PARS)]    <-rnorm(length(PARS),3.4,0.98)  # init. fractionation for each organism
  readLIM$pars$val[which(readLIM$pars$name == "FracBIO")] <-runif(1, min=1, max=3.4)      # init. fractionation for the biofilm
  fracbio <- readLIM$pars$val[which(readLIM$pars$name == "FracBIO")]
  readLIM$pars$val[which(readLIM$pars$name == "IsoBIO")]  <-runif(1, min=4+fracbio, max=8.7+fracbio) # init. isotope value biofilm 
  LIM <- Setup(readLIM)
  j<-j+1
  Sollsei <- Lsei(LIM,parsimonious=TRUE)
  if(Sollsei$IsError==TRUE) next   # if initilization results in inconsistent model -> retry
  i <- i+1
  Ranges <- Xranges(LIM)
  Solxsam<- Xsample(lim  = LIM,
                    exact= NULL,
                    iter = 250, 
                    type = "mirror",
                    jmp  = c(Ranges[,2]-Ranges[,1])/4,
                    x0   = Sollsei$X)
  SolBaySens <- rbind(SolBaySens,Solxsam)
  }
} # end not run

## =============================================================================
##
## Fig 2 - four structures
##
## =============================================================================

# pdf(file = "C:/DataWerk/Hermes/wpColdwaterCorals/Final L&O/Fig2 + Fig5.pdf",
#  onefile=TRUE, family="serif", paper="a4", colormodel="cmyk")
 windows()
 par(mfrow=c(2,2),mgp=c(2,0.5,0),family="serif")

 #color assignment to flows
 BIO <- c("CWC","EUN","HES","URC","SPO","STA","CRA","HYD","CRI","LIM","ASP",
   "POL","BIV","SUS","OMN","FIS","INF")
 FM <- CM <- Flowmatrix(LIM,colMeans(SolBayDef))
 CM[which(CM>0.0)] <- colors()[c(170)] #darkgrey
 CM["PHY_w",]      <- colors()[c(497)] #green
 CM["DET_w",]      <- colors()[c(572)] #brown
 CM["ZOO_w",]      <- colors()[c(62)]  #purple grey
 CM[,"DIC"]        <- colors()[c(421)] #pinkish
 CM[BIO,BIO]       <- colors()[c(77)]  #darkorange
 # upper flow values
 maxflow_   <- c(50,5,0.5,0.05)
 for (i in 1:4){
   plotweb(flowmat=FM,nullflow=c(7.e-05,maxflow_[i]),
           minflow = 0.00001, maxflow = maxflow_[i],lab.size = 1.4,
           minarrow=0.5,maxarrow=10,fig.size = 1.3,leg.title="",bty="n",
           mar = c(0,0,0,0),
           arr.col = CM,lcol=CM,legend=TRUE)
   mtext(LETTERS[i],adj=0,cex=1.3)
   }
#   dev.off()
# dev.print(pdf,"C:/DataWerk/Hermes/wpColdwaterCorals/Final L&O/Fig2.pdf")


#Original function cannot plot 1x10-5 as log scale in legend. To achieve this
# replace the original legend function (in plotweb) with the following and
# load modified plotweb function to R.
# Moreover, make sure that package diagram 1.4 is loaded for full
# functionality of arr.col
#
# MODIFIED VERSION:
#  legend("bottomright", legend = c(format.pval(tmax, leg.digit),
#         expression(paste("1x") * group("","" * 10^{-05}, ""))),
#         cex = sizeleg, title = title,
#         lwd = c(maxarrow, minarrow), bty = bty)


## =============================================================================
##
## Fig 3 - flow ranges
##
## =============================================================================

# Original, simpler figure
 A4()
 par(mar=c(4, 8, 0, 10),family="serif")
 web    <- colMeans(SolBayDef)
 std    <- apply(SolBayDef, 2, sd)
 ord    <- order(web,decreasing = TRUE)
 xlab_  <- expression(paste("Flow value ") * group("(", "mmol C " * m^{-2} * d^{-1}, ")"))
 Plotranges(min=rev((web[ord]-std[ord])), max=rev((web[ord]+std[ord])), value=rev((web[ord])),
            lab.cex=0.3,log="x", ,
            xlab=xlab_,labels=rep("",length(web))  ,
            main="",pch=20,pch.col=c("black","grey"),
            xlim=c(1e-5,100))
 even   <- seq(2,140,2)
 uneven <- seq(1,139,2) 
 mtext(rev((LIM$Unknowns[ord][even]  )), side = 2, line=0.5,at=even, adj = 0, las = 2,cex=.5)
 mtext(rev((LIM$Unknowns[ord][uneven])), side = 4, line=0.5,at=uneven, adj = 0, las = 2,cex=.5) 
 for (i in 0:3)
  {
  ticks <- 1e-5*10^i*rep(1:9)
  segments(ticks,-4.5,ticks,-6)
  }

# Version for revision

 par(mar=c(4, 10, 0.1, 10),mgp=c(2.,0.6,0),family="serif",yaxs="i")
 web    <- colMeans(SolBayDef)
 std    <- apply(SolBayDef, 2, sd)
 ord     <- rev(order(web,decreasing = TRUE))
 even   <- seq(2,140,2)
 uneven <- seq(1,139,2) 
 plot(c(1:2,web[ord][even],1:2),c(-1:0,even,141:142),log="x",axes=FALSE,pch=19,
      xlab=xlab_,ylab="",las=2, cex.lab=0.85,
      xlim=c(1e-05,1e+02),type="n")
 labs <- c(expression(paste("1x") * group("","" * 10^{-05}, "")),
          expression(paste("1x") * group("","" * 10^{-03}, "")),
          expression(paste("1x") * group("","" * 10^{-01}, "")),
          expression(paste("1x") * group("","" * 10^{1}, "")))
 axis(1,labels=labs,at=c(1e-05,1e-03,1e-01,1e+01),cex.axis=0.8)
 axis(2,labels=rep("",70),at=even,las=2,cex.axis=0.5,tick=TRUE,lwd.ticks=0.5)
 axis(2,labels=LIM$Unknowns[ord][even],at=even,las=2,cex.axis=0.6,tick=FALSE) 
 axis(4,labels=rep("",70),at=uneven,las=2,cex.axis=0.5)
 axis(4,labels=LIM$Unknowns[ord][uneven],at=uneven,las=2,cex.axis=0.6,tick=FALSE) 
 segments(1e-06,even,1e+03,even,col="grey",lty="dotted")
 segments((web-std)[ord],1:140 ,(web+std)[ord],1:140,col="black") 
 points(web[ord],1:140,  pch=20,cex=0.75,col="black")
 box()

message("CoV <50% for:",length(which(std/web<0.5 ))/LIM$NUnknowns,"%")
message("CoV <75% for:",length(which(std/web<0.75))/LIM$NUnknowns,"%")


## =============================================================================
##
## Fig 4 - uncertainty analysis
##
##
## =============================================================================

## Can only be displayed if SolBaySens has been made; as this is toggled off,
## so are these figures

notrun <- function(){
 windows()
 par(mar=c(5, 5, 5, 5),family="serif")
 par(mgp=c(3, 0.4, 0))
 X  <- (colMeans(SolBaySens)-colMeans(SolBayDef))/colMeans(SolBayDef)
 Y  <- ( colMeans(SolBaySens)-colMeans(SolBayDef))
 XY <- cbind(X,Y)
 xlim_ <- c(1e-5,1e3)
 ylim_ <- c(1e-5,1e2)
 plot(subset(XY, X>1e-5 & X<10),
     log="xy",
     xlim=xlim_,
     ylim=ylim_,
     pch="+",
     ylab="Absolute flow change",
     axes=FALSE,
     las=2,
     cex.lab=1,
     xlab="Relative flow change")
 points(abs(subset(XY, X<0)),pch="-",cex=1.1,col="black")
 deviations <- subset(XY, X>10)
 ord <- order(deviations[,"X"])
 deviations <- deviations[ord,]
 arrows(deviations[1,1]+5,
       deviations[1,2]+0.05,
       deviations[1,1],
       deviations[1,2], length = 0.0, angle = 30)
 points(deviations[1,1]+7,
       deviations[1,2]+0.08,pch="a")
 arrows(deviations[2,1]+15,
       deviations[2,2]+0.03,
       deviations[2,1],
       deviations[2,2], length = 0.0, angle = 30)
 points(deviations[1,1]+22,
       deviations[1,2]+0.05,pch="b")
 arrows(deviations[3,1],
       deviations[3,2],
       deviations[3,1]+2,
       deviations[3,2]-0.015, length = 0.0, angle = 30)
 points(deviations[3,1]+2,
       deviations[3,2]-0.018,pch="c")
 points(deviations[4:nrow(deviations),],pch=letters[4:nrow(deviations)],cex=0.95)
 points(abs(subset(XY, X< -10)),pch=LETTERS[1:nrow(subset(XY, X< -10))],cex=0.95)
#axis(2,at=10^(seq(-5,3,1)),labels=TRUE,las=2)
 labs <- c(expression(paste("1x") * group("","" * 10^{-05}, "")),
          expression(paste("1x") * group("","" * 10^{-04}, "")),
          expression(paste("1x") * group("","" * 10^{-03}, "")),
          expression(paste("1x") * group("","" * 10^{-02}, "")),
          expression(paste("1x") * group("","" * 10^{-01}, "")),
          expression(paste("1x") * group("","" * 10^{0}, "")),
          expression(paste("1x") * group("","" * 10^{1}, "")),
          expression(paste("1x") * group("","" * 10^{2}, "")),
          expression(paste("1x") * group("","" * 10^{3}, "")))

 axis(1,at=10^(seq(-5,3,1)),labels=labs,tck=-0.01,cex.axis=0.9)
 axis(2,at=10^(seq(-5,3,1)),labels=labs,tck=-0.01,las=2,cex.axis=0.9)
 legend("bottomright",
       legend=rownames(deviations),
       pch=letters[1:nrow(deviations)],
       bty="n",
       cex=.75,
       pt.cex=0.85)
 for (i in 0:7)
  {
  ticks<-1e-5*10^i*rep(1:9)
  segments(ticks,5e-6,ticks,6.2e-6)
  segments(5e-6,ticks,5.7e-6,ticks)
  }
 box()

# data analysis
# difference for <100% and <1000%
  c(length(which(abs(X)<1)), length(which(abs(X)<1))/LIM$NUnknowns,
    length(which(abs(X)<10)), length(which(abs(X)<10))/LIM$NUnknowns)
  # max absolute and relative difference
  c(max(abs(Y)),max(abs(X)))
  #max absolute difference for 100<flowdifference<1000%
  sort(Y[which(abs(X)>10)])
  max(Y[which(abs(X)>10)])
}


###############################
#                             #
#  DATA OUTPUT                # 
#                             #
#                             #
###############################
#food sources of CWC
web<-colMeans(SolBayDef)
fm <- Flowmatrix(LIM,web)
SOURCE <- c("DET_w")
#food sources of CWC reef
  Idet_w<-sum(fm[which(rownames(fm)%in%SOURCE,arr.ind = TRUE),])
  Odet_w<-sum(fm[,which(colnames(fm)%in%SOURCE,arr.ind = TRUE)])
  Icwcreef<-c(sum(fm["ZOO_w",]),sum(fm["PHY_w",]),Idet_w-Odet_w)
  sum(Icwcreef)
  data.frame(rbind(c("ZOO","PHY","DET"),round(Icwcreef,2),c(round(Icwcreef/sum(Icwcreef),2)*100)))
#food sources of CWC
  Icwc <- c(fm["ZOO_w","CWC"],fm["PHY_w","CWC"],fm["DET_w","CWC"])
  data.frame(rbind(c("ZOO","PHY","DET"),
                   c(round(Icwc,2)),
                   c(round(Icwc/sum(Icwc),2))))
  sum(Icwc)
  sum(Icwc)/sum(Icwcreef)

#% input consumption by CWC versus CWC reef
sum(Icwc)/sum(Icwcreef)

#Total respiration
round(sum(fm[,"DIC"]),2)

# difference between Icwcreef and DIC
sum(Icwcreef)-round(sum(fm[,"DIC"]),2)
#Total export
sum(fm[,"EXP"])
sum(fm[,"EXP"])/(sum(Icwcreef)-sum(fm[,"DIC"]))

#%contributions of compartments to export
round(fm[,"EXP"]/sum(fm[,"EXP"])*100,1)

# comparison in food uptake by sponges between default and sensitivity analysis (sens analysis toggled off)
c(colMeans(SolBayDef)["ZOO_w->SPO"],colMeans(SolBayDef)["PHY_w->SPO"],colMeans(SolBayDef)["DET_w->SPO"])
#c(colMeans(SolBaySens)["ZOO_w->SPO"],colMeans(SolBaySens)["PHY_w->SPO"],colMeans(SolBaySens)["DET_w->SPO"])
sum(colMeans(SolBayDef)["ZOO_w->SPO"],colMeans(SolBayDef)["PHY_w->SPO"],colMeans(SolBayDef)["DET_w->SPO"])
#sum(colMeans(SolBaySens)["ZOO_w->SPO"],colMeans(SolBaySens)["PHY_w->SPO"],colMeans(SolBaySens)["DET_w->SPO"])



