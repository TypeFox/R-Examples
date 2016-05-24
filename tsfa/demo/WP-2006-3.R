### Generating postscript files for  Figures in paper can be changed with
PSFILES <- FALSE

######################################################################
######################################################################

# Generate figures and tables for paper  "Money and Credit Factors"
# Paul Gilbert and Erik Meijer, Bank of Canada Working Paper 2006-3

######################################################################

options(PlotTitles=FALSE, ModSeriesNames=FALSE)

data("CanadianMoneyData.asof.28Jan2005", package="CDNmoney")
data("CanadianCreditData.asof.28Jan2005", package="CDNmoney")

require("tsfa") 


##############################################################################

cat("Starting the construction of the data...\n")

######################################################################

cpi <- 100 * M1total / M1real
seriesNames(cpi) <- "CPI"
popm <- M1total / M1PerCapita
seriesNames(popm) <- "Population of Canada"

z <- tframed(tbind(
    MB2001,
    MB486 + MB452 + MB453 ,
    NonbankCheq,
    MB472 + MB473 + MB487p,
    MB475,
    NonbankNonCheq + MB454 + NonbankTerm + MB2046 + MB2047 + MB2048 +
    MB2057 + MB2058 + MB482),
    names=c("Currency", "Personal chequing", "Non-bank chequing",
    "N-P demand & notice", "N-P term", "Investment" )
  )


TotalMoney <- tframed(rowSums(z), tframe(z))

z <- tbind (z, ConsumerCredit, ResidentialMortgage,
    ShortTermBusinessCredit, OtherBusinessCredit)

# tfplot(tfwindow(z, start=c(1981,11), end=c(2004,11)), graphs.per.page=3)

# 28Jan2005 data has NAs after Nov 2004
# investment series goes back only to Nov 1981.
z <-tfwindow(z, start=c(1981,11), end=c(2004,11))

# 1e8 * gives real $ per person
#(Credit aggregates, B and MB numbers in millions,
#   CPI in fraction*100, popm in units.)

scale <- tfwindow(1e8 /(popm * cpi), tf=tframe(z))

MBandCredit <- sweep(z, 1, scale, "*")

tfplot(MBandCredit, graphs.per.page=3)
tfplot(diff(MBandCredit), graphs.per.page=3)

#############################################################################
cat("##########################\n")
cat(" Various sample statistics\n")
cat("##########################\n")

start(MBandCredit)
end(MBandCredit)
periods(MBandCredit)

DX <- diff(MBandCredit, lag=1)

cat("number of observations ", periods(MBandCredit), "\n")

cat("number of series ", nseries(MBandCredit), "\n")

cat("means of differenced series ", colMeans(DX), "\n")

cat("stds of differenced series ", sqrt(diag(cov(DX))), "\n")


###################################################################################
#  Factors from the money and credit data
###################################################################################

###########################################
# check for number of factors
###########################################

zz <- eigen(cor(diff(MBandCredit, lag=1)), symmetric=TRUE)$values
cat("\nEigenvalues for scree plot\n")
print(zz)

cat("##########################\n")
cat("    Figure 1: scree plot\n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/screeplot.eps", onefile=FALSE,
              horizontal=FALSE, width=4, height=4)}
par(omi=c(0.1,0.1,0.1,0.1),mar=c(4.1,4.1,0.6,0.1))
# mar: bottom,left,top,right
plot(zz, ylab="Value", xlab="Eigenvalue number",
       pch=20:20,cex=1,type="o")
if (PSFILES) dev.off()


z <- FAfitStats(MBandCredit)    
#print(z, digits=4)


cat("\n##########################\n")
cat("   Table 1 \n")
cat("##########################\n")
print(z$fitStats[c("chisq", "df", "pval", "RMSEA", "CFI"),], digits=4)
print(z$seqfitStats, digits=4)

c2withML <- estTSF.ML(MBandCredit, 2)

# Arrange so 1=currency, 2=investment, 3=per cheq, 4=consumer credit , 5=N-P term
# (which is the order they appear as factors are added.)
# BpermuteTarget=z in the next is just to reorder and fix signs
# These values are determined by an initial run (which was not in the order
# that works best for presentation).


z <- matrix(0,10,3)
z[matrix(c( 1,6,2,1:3),3,2)] <- c(10, 56, 41)
c3withML <- estTSF.ML(MBandCredit, 3, BpermuteTarget=z)

z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(13, 54, 37, 24)
c4withML  <- estTSF.ML(MBandCredit, 4, BpermuteTarget=z)

z <- matrix(0,10,5)
z[matrix(c( 1,6,2,7,5,1:5),5,2)] <- c(13, 67, 34, 30, 72)
c5withML <- estTSF.ML(MBandCredit, 5, BpermuteTarget=z)

cat("\n##########################\n")
cat("     Table 2 \n")
cat("standardized Loadings for 4 factor model\n")
cat("##########################\n")
print(DstandardizedLoadings(c4withML) )

cat("\n   Table 2 Factor correlations\n")
print(c4withML$model$Phi, digits=3)

cat("\n   Table 2 communalities for 4 factor model\n")
print(1 - c4withML$estimates$uniquenesses)

cat("\n##########################\n")
cat(" Table 3 \n")
cat("Loadings for 4 factor model\n")
cat("##########################\n")
print(loadings(c4withML) )

FactorNames4 <-  c("Transaction factor","Long-term factor",
                   "Potential spending factor","Consumer credit factor")
cat("##########################\n")
cat("       Figure 2 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/f4growth.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(factors(c4withML)),
       Title="Factors from 4 factor model (year-to-year growth rate)",
       lty=c("solid"),
       col=c("black"),
       xlab=c(""),ylab=FactorNames4,
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

cat("##########################\n")
cat("    Figure 3 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/f4.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(factors(c4withML),
       Title= "Factors from 4 factor model",
       lty=c("solid"),
       col=c("black"),
       xlab=c(""),ylab=FactorNames4,
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()


z <- explained(c4withML)

cat("##########################\n")
cat("     Figure 4  \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/e4growtha-corr.eps", onefile=FALSE,
        horizontal=FALSE, paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(MBandCredit), ytoypc(z), series=1:5, graphs.per.page=5,
       lty=c("solid", "dashed"), 
       col=c("black", "black"),
#       col=c("black", "red"),
       ylab=c("Currency", "Personal chequing", "Non-bank chequing",
              "N-P demand & notice", "N-P term"),
       ylim=list(NULL,NULL,c(-100,100),NULL,NULL),
       Title=
   "Explained money indicator 1-5 (year-to-year growth rate) using 4 factors",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

cat("##########################\n")
cat("      Figure 5 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/e4growthb.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(MBandCredit), ytoypc(explained(c4withML)), series=6:10,
       graphs.per.page=5,
       lty=c("solid", "dashed"), 
       col=c("black", "black"),
#       col=c("black", "red"),
       ylab=c("","","","","",
              "Investment","Consumer credit", "Residential mortgage",
              "Short-term business credit", "Other business credit"),
       Title=
   "Explained money indicator 6-10 (year-to-year growth rate)using 4 factors",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

cat("##########################\n")
cat("      Figure 6 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/e4a.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(  MBandCredit,  explained(c4withML),  series=1:5, graphs.per.page=5,
        lty=c("solid", "dashed"), 
       col=c("black", "black"),
#	col=c("black", "red"),
        Title=
  "Observered money indicators 1-5 and their explained portion, using 4 factors",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

cat("##########################\n")
cat("     Figure 7 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/e4b.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(  MBandCredit,  explained(c4withML),  series=6:10, graphs.per.page=5,
        lty=c("solid", "dashed"), 
       col=c("black", "black"),
#	col=c("black", "red"),
        Title=
 "Observered money indicator 6-10 and their explained portion, using 4 factors",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

#tfplot(  diff(MBandCredit), diff(explained(c4withML)),   graphs.per.page=3)


###################################################################################
#  Two, three and five factors models
###################################################################################

cat("\n##########################\n")
cat("\n      Table 4 \n")
cat("##########################\n")
print(DstandardizedLoadings(c2withML), digits=3)

cat("\n   Table 4 communalities for 2 factor model\n")
print(1 - c2withML$estimates$uniquenesses)

cat("\n   Table 4  Factor correlations\n")
print(c2withML$model$Phi, digits=3)

cat("##########################\n")
cat("\n    Table 5 \n")
DstandardizedLoadings(c3withML)

cat("\n   Table 5 Factor correlations\n")
print(c3withML$model$Phi, digits=3)
cat("\n   Table 5 communalities for 3 factor model\n")
print(1 - c3withML$estimates$uniquenesses)

cat("##########################\n")
cat("\n    Table 6 \n")
DstandardizedLoadings(c5withML)

cat("     Table 6 Factor correlations\n")
print(c5withML$model$Phi, digits=3)

cat("\n   Table 6 communalities for 5 factor model\n")
print(1 - c5withML$estimates$uniquenesses)


cat("##########################\n")
cat("    Figure 8  part 1 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/fAll2growth.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=4)}
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c2withML)),
       ytoypc(factors(c3withML)),
       ytoypc(factors(c5withML)), series=1:2,
       xlab=c(""),ylab=c("Factor 1","Factor 2"),
       lty=c("solid", "dotdash", "dashed", "dotted"),
#       col=c("black","green","red","blue"),
       col=c("black","black","black","black"),
       Title= 
 "Factors transaction and long (year-to-year growth rate) using 2-, 3-, 4-, and 5-factor models",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

cat("##########################\n")
cat("     Figure 8 part 2 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/f3Allgrowth.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=2)}
tfplot(ytoypc(factors(c4withML)),
       ytoypc(factors(c3withML)),
       ytoypc(factors(c5withML)), series=3,
       lty=c("solid", "dashed", "dotted"),
       xlab=c(""),ylab=c("","","Factor 3"),
#       col=c("black","red","blue"),
       col=c("black","black","black"),
       Title=
   "Factor near (year-to-year growth rate) using 3-, 4-, and 5-factor models",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

###################################################################################
########## Rotation method sensitivity
###################################################################################

# BpermuteTarget just helps put the factors in the same order, and with the
#  same signs as c4withML. It does not affect the estimation or rotation.

# this is a bit similar but not so interesting. 3 and 4 don't load on anything.
# Plotted factor growth rates are very similar
c4withMLg0.5 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
     rotation="oblimin", rotationArgs=list(gam=0.5))
loadings(c4withMLg0.5)
DstandardizedLoadings(c4withMLg0.5)
DstandardizedLoadings(c4withMLg0.5) - DstandardizedLoadings(c4withML)
summary(c4withMLg0.5)
tfplot(factors(c4withML), factors(c4withMLg0.5))
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLg0.5)))

# very similar
c4withMLgneg0.5 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
     rotation="oblimin", rotationArgs=list(gam=-0.5))
loadings(c4withMLgneg0.5)
DstandardizedLoadings(c4withMLgneg0.5)
DstandardizedLoadings(c4withMLgneg0.5) - DstandardizedLoadings(c4withML)
summary(c4withMLgneg0.5)
tfplot(factors(c4withML), factors(c4withMLgneg0.5))
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLgneg0.5)))

# very similar
c4withMLgneg1.0 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
     rotation="oblimin", rotationArgs=list(gam=-1.0))
loadings(c4withMLgneg1.0)
DstandardizedLoadings(c4withMLgneg1.0)
DstandardizedLoadings(c4withMLgneg1.0) - DstandardizedLoadings(c4withML)
summary(c4withMLgneg1.0)
tfplot(factors(c4withML), factors(c4withMLgneg1.0))
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLgneg1.0)))

# this does not converge
# c4withMLg1.0 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
#      rotation="oblimin", rotationArgs=list(gam=1.0))
# loadings(c4withMLg1.0)
# DstandardizedLoadings(c4withMLg1.0)
# DstandardizedLoadings(c4withMLg1.0) - DstandardizedLoadings(c4withML)
# summary(c4withMLg1.0)


# very similar
c4withMLbQ <- estTSF.ML(MBandCredit, 4, rotation="bentlerQ",
     BpermuteTarget=loadings(c4withML))
loadings(c4withMLbQ)
DstandardizedLoadings(c4withMLbQ)
DstandardizedLoadings(c4withMLbQ) - DstandardizedLoadings(c4withML)
summary(c4withMLbQ)
### tfplot(factors(c4withML), factors(c4withMLbQ))
### tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLbQ)))

cat("##########################\n")
cat("     Figure 9 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/f4Manygrowth.eps",
              onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLg0.5)),
       ytoypc(factors(c4withMLgneg0.5)), ytoypc(factors(c4withMLgneg1.0)),
       ytoypc(factors(c4withMLbQ)),
       xlab=c(""),ylab=FactorNames4,
       lty=c("solid", "dashed", "dotted", "dotdash", "longdash"),
#      col=c("black","red","blue","green","pink"),
       col=c("black","black","black","black","black"),
       Title=
 "Factors from various 4-factor models (year-to-year growth rate)\n and oblimin with gam=0 (solid)",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()


#geomin  2 and 3 each have one modestly different loading
# factor 2 has personal chequing mixed in with investment and credit
# factor 3 loads on currency, personal chequing, and investment
#  so the separation is not so interesting.
c4withMLgm <- estTSF.ML(MBandCredit, 4, rotation="geominQ",
     BpermuteTarget=loadings(c4withML))
loadings(c4withMLgm)


cat("\n##########################\n")
cat("     Table 7 \n")
cat("##########################\n")
print(DstandardizedLoadings(c4withMLgm), digits=3)


DstandardizedLoadings(c4withMLgm) - DstandardizedLoadings(c4withML)
summary(c4withMLgm)
#tfplot(factors(c4withML), factors(c4withMLgm))

cat("##########################\n")
cat("     Figure 10 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/fgmgrowth.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLgm)),
       xlab=c(""),ylab=FactorNames4,
       lty=c("solid", "dashed"),
       col=c("black", "black"),
#       col=c("black","red"),
       Title=
 "Factors from geomin (dashed) 4-factor model (year-to-year growth rate) \n and oblimin with gam=0 (solid)",
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()


c4withMLnotNorm  <- estTSF.ML(MBandCredit, 4, normalize=FALSE,
       BpermuteTarget=loadings(c4withML))


# There is only a qualitative statement about this in the paper
DstandardizedLoadings(c4withML)
DstandardizedLoadings(c4withMLnotNorm)
DstandardizedLoadings(c4withML) - DstandardizedLoadings(c4withMLnotNorm)

###################################################################################
########## Sensitivity to sample period
###################################################################################

# BpermuteTarget=loadings(c4withML) is not good enough in some cases

# there are difficulties interpretating factors 2 and 3 in tis case
z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(11, 104, 20, 13)
c4withMLbefore90 <- estTSF.ML(tfwindow(MBandCredit, end=c(1989,12)), 4,
         BpermuteTarget=z)

c4withMLafter95 <- estTSF.ML(tfwindow(MBandCredit, start=c(1995,1)), 4,
         BpermuteTarget=loadings(c4withML))

z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(11, 104, 20, 13)
c4withMLbefore95 <- estTSF.ML(tfwindow(MBandCredit, end=c(1994,12)), 4,
         BpermuteTarget=z)

c4withMLafter00 <- estTSF.ML(tfwindow(MBandCredit, start=c(2000,1)), 4,
         BpermuteTarget=loadings(c4withML))

c4withML90to00 <- estTSF.ML(tfwindow(MBandCredit, start=c(1990,1), end=c(2000,1)), 4,
         BpermuteTarget=loadings(c4withML))


cat("##########################\n")
cat("     Figure 11 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/fSubSamplesGrowth-corr.eps",
              onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(factors(c4withML)),         ytoypc(factors(c4withMLbefore90)),
       ytoypc(factors(c4withMLbefore95)), ytoypc(factors(c4withMLafter95)),
       ytoypc(factors(c4withMLafter00)),  ytoypc(factors(c4withML90to00)),
       xlab=c(""),
       ylab=FactorNames4,
#       ylim=list(NULL,c(-50,50),c(-50,50),NULL),
       ylim=list(NULL,c(-20,20),c(-25,40),NULL),
       graphs.per.page=4,
       lty=c("dashed", "dotted", "dotdash", "longdash", "dotted",
             "twodash"),
#       col=c("red","blue","green","pink","violet","brown"),
       col=c("black","black","black","black","black","black"),
       Title=
    paste("Factors (year to year growth) using full sample and sub-samples\n",
           "ML estimation with quartimin rotation objective", sep=""),
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

cat("##########################\n")
cat("     Figure 12 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/eSubSamplesGrowtha-corr.eps",
              onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(MBandCredit),                 ytoypc(explained(c4withML)),
       ytoypc(explained(c4withMLbefore90)), ytoypc(explained(c4withMLbefore95)),
       ytoypc(explained(c4withMLafter95)),  ytoypc(explained(c4withMLafter00)),
       ytoypc(explained(c4withML90to00)), series=1:5, graphs.per.page=5,
       ylab=c("Currency", "Personal chequing", "Non-bank chequing",
              "N-P demand & notice", "N-P term"),
#       ylim=list(NULL,NULL,c(-100,100),NULL,c(-100,100)),
       ylim=list(NULL,NULL,c(-70,70),NULL,c(-70,70)),
       lty=c("solid", "dashed", "dotted", "dotdash", "longdash", "dotted",
             "twodash"),
#       col=c("black", "red","blue","green","pink","violet","brown"),
       col=c("black","black","black","black","black","black","black"),
       Title= paste("Explained money indicators 1-5 (year to year growth)\n",
                     "using 4 factors, full sample and sub-samples", sep=""),
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

cat("##########################\n")
cat("     Figure 13 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/eSubSamplesGrowthb.eps",
              onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=8)}
tfplot(ytoypc(MBandCredit), ytoypc(explained(c4withML)),
       ytoypc(explained(c4withMLbefore90)), ytoypc(explained(c4withMLbefore95)),
       ytoypc(explained(c4withMLafter95)),  ytoypc(explained(c4withMLafter00)),
       ytoypc(explained(c4withML90to00)), series=6:10, graphs.per.page=5,
       ylab=c("","","","","",
              "Investment","Consumer credit", "Residential mortgage",
              "Short-term business credit", "Other business credit"),
       lty=c("solid", "dashed", "dotted", "dotdash", "longdash", "dotted",
             "twodash"),
#       col=c("black", "red","blue","green","pink","violet","brown"),
       col=c("black","black","black","black","black","black","black"),
       Title= paste("Explained money indicators 6-10 (year to year growth)\n",
               "using 4 factors, full sample and sub-samples", sep=""),
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,4.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()



###################################################################################
########## Comparison with Aggregates
###################################################################################

M1   <-  tfwindow(M1total, start=c(1981,11), end=c(2004,11)) * scale
seriesNames(M1) <- "Real per Capita M1"

z <-  tframed(MB2001 + MB486 + MB487p + MB452 + MB452adj + MB472 + NonbankCheq)

M1p   <-  tfwindow(z, start=c(1981,11), end=c(2004,11)) * scale
seriesNames(M1p) <- "Real per Capita M1+"

M2pp <-  tfwindow(M1total
                 + MB472 + MB473 + MB452 + MB453 + MB454
                 + NonbankCheq + NonbankNonCheq + NonbankTerm +
                 + MB2046 + MB2047 + MB2048
                 + MB2057 + MB2058, start=c(1981,11), end=c(2004,11))* scale

seriesNames(M2pp) <- "Real per Capita M2++"

#transaction and investment
f <- tframed(factors(c4withML)[,1:2], tf=tframe(factors(c4withML)))
mnF <- colMeans(f)
mnM <- colMeans(cbind(M1p, M2pp))
f <- sweep(f, 2, mnM/mnF, "*")

cat("##########################\n")
cat("    Figure 14 \n")
cat("##########################\n")
if (PSFILES) {postscript(file="figsDynamic/M1M2growth.eps", onefile=FALSE, horizontal=FALSE,
              paper="special",pagecentre=FALSE, width=6, height=4)}
tfplot(ytoypc(f), ytoypc(cbind(M1, M2pp)), graphs.per.page=2,
       lty=c("dashed", "solid"),
#       col=c("red","black"),
       col=c("black","black"),
       Title=
    paste("(year to year growth) M1+ and M2++ (solid) and scaled Bartlett Predictors\n",
           "computed using ordinary ML parameters (dashed)", sep=""),
       ylab=c("M1+ vs.\ntransaction factor", "M2++ vs.\nlong-term factor" ),
       par=list(omi=c(0.1,0.1,0.1,0.1),mar=c(3.1,5.1,0.6,0.1)),
# mar: bottom,left,top,right
       reset.screen=TRUE)
if (PSFILES) dev.off()

