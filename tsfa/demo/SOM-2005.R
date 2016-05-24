### Generating postscript files for  Figures in paper can be changed with
PSFILES <- FALSE

######################################################################

# Generate figures and tables for paper
# Time Series Factor Analysis with an Application to Measuring Money
# Erik Meijer & Paul Gilbert,  2005 
# University of Groningen, Research School SOM
# Research Report 05F10

######################################################################

options(PlotTitles=FALSE, ModSeriesNames=FALSE)

data("CanadianMoneyData.asof.6Feb2004", package="CDNmoney")

require("tsfa")

##############################################################################

cat("Starting the construction of the data...\n")

######################################################################

# This grouping of the data is slightly different than that use in
# Gilbert & Pichette (2003) because there is going to be a small change in the
# reporting forms and in the future we will not be able to do the grouping
# as in Gilbert & Pichette. You can see the old grouping with
# > help("CanadianMoneyData", package="CDNmoney")

#  MB453 is with personal chequing rather than investment, as this will 
#   be necessary if the proposal to combine demand and notice proceeds.
# Data extends back to Nov 1981, and could go back further except for MB482.

# Earlier work (prior to 2004) was done with a grouping that has 
#   demand/notice split for bank deposits:

#                      Previously                  Now

# currency           MB2001                    MB2001
# personal cheq.     MB486 + MB452 	       MB486 + MB452 + MB453
# NonbankCheq        NonbankCheq	       NonbankCheq
# N-P demand&notice  MB472 + MB473 + MB487p    MB472 + MB473 + MB487p
# N-P term           MB475		       MB475
# Investment         NonbankNonCheq + MB454    NonbankNonCheq + MB454
#                     + NonbankTerm + MB2046    + NonbankTerm + MB2046
#       	      + MB2047 + MB2048         + MB2047 + MB2048
#       	      + MB2057 + MB2058         + MB2057 + MB2058
#       	      + MB482 + MB453	        + MB482


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
    names=c("currency", "personal cheq.", "NonbankCheq",
    "N-P demand & notice", "N-P term", "Investment")
    )

# 6Feb2004 data has NAs after Nov 2003
z <-tfwindow(z, start=c(1986,1), end=c(2003,11))

# 1e8 * gives real $ per person
#(B and MB numbers in millions, CPI in fraction*100, popm in units.)
scale <- 1e8 /tfwindow(popm * cpi, start=c(1986,1), end=c(2003,11))

MBcomponents <- sweep(z, 1, scale, "*")

M1   <-  tfwindow(M1total, start=c(1986,1), end=c(2003,11)) * scale
seriesNames(M1) <- "Real per Capita M1"

M2pp <-  tfwindow(M1total
                 + MB472 + MB473 + MB452 + MB453 + MB454 
                 + NonbankCheq + NonbankNonCheq + NonbankTerm +
                 + MB2046 + MB2047 + MB2048
                 + MB2057 + MB2058, start=c(1986,1), end=c(2003,11))* scale

seriesNames(M2pp) <- "Real per Capita M2++"

#par (ask=T)
#tfplot(MBcomponents, graphs.per.page=3)
##############################################################################

cat("Starting the construction of true parameters for simulation...\n")

#write(t(MBcomponents),file="figsDynamic/money.dat",ncolumns=ncol(MBcomponents))

######################################################################
# Compute the differenced series.

DX <- diff(MBcomponents, lag=1)
#write(t(DX),file="figsDynamic/dmoney.dat",ncolumns=ncol(DX))


######################################################################
# The various population parameters

cat("number of observations ", periods(MBcomponents), "\n")

cat("number of series ", nseries(MBcomponents), "\n")

Nfact   <- 2
SDX     <- cov(DX)
stds    <- sqrt(diag(SDX))           # standard deviations

meanDX  <- matrix(apply(DX,2, mean), 1, nseries(MBcomponents))

T       <- nrow(DX)
M       <- ncol(DX)

cat("number of rows in DX: ", T, "\n")
cat("number of columns in DX: ", M, "\n")


# Chi-square of null model. Null model: Sigma = Omega = diagonal.
# ML solution: Omega = diag(diag(SDX)).
# So chi-square = 2 * log LR
#               = 2 * (T/2) * [log det Sigma + tr(Sigma^{-1} S)
#                      - log det S - tr(S^{-1}S)]
#               = T * [sum(i=1..M) log S_ii - log det S]
# because tr(Sigma^{-1} S) = tr(S^{-1}S) = M.

chisqnull <- T * ( sum(log(diag(SDX))) - log(det(SDX)) )
cat("chisqnull: ", chisqnull, "\n")

#test <- T * ( log(det(diag(diag(SDX)))) - log(det(SDX)) )
#print(test)
 
 
######################################################################
# What are we estimating?

poppar <- factanal(factors = Nfact, covmat = SDX, n.obs = periods(DX),
            scores = "none", rotation = "none")

# rotation

# row norms of B
rownorms  <- sqrt(apply(poppar$loadings^2, 1,sum))

# Quartimin (= oblimin with parameter 0) rotation with Kaiser normalization
KBorth    <- diag(1/rownorms) %*% poppar$loadings
rotpoppar <- GPFoblq(A = KBorth, Tmat = diag(1, Nfact), method="quartimin")
Boblq     <- diag(stds * rownorms) %*% rotpoppar$loadings
PhiOblq   <- rotpoppar$Phi
Omega     <- diag(stds * poppar$uniquenesses * stds)
Psi       <- 0.5 * Omega

# So the true parameters in the simulation are
cat("\ntrue value Boblq\n")
print(Boblq)
cat("\ntrue value PhiOblq\n")
print(PhiOblq)
cat("\ntrue value Omega\n")
print(Omega)
cat("\ntrue value Psi\n")
print(Psi)
SigmaOblq <- Boblq %*% PhiOblq %*% t(Boblq) + Omega
cat("\ntrue value SigmaOblq\n")
print(SigmaOblq)


cat("##########################\n")
cat("\nCommunalities (difference scale)\n")
cat("##########################\n")
print(1 - matrix(poppar$uniquenesses,ncol=1))


# ######################################################################
# # How well do the means fit?
 
kappaOblq <- solve(t(Boblq) %*% solve(SigmaOblq) %*% Boblq) %*%
                   t(Boblq) %*% solve(SigmaOblq) %*% t(meanDX)
meanhat <- t(Boblq %*% kappaOblq)
 
cat("\nmeanhat\n")
print(meanhat)
cat("\nmeanDX\n")
print(meanDX)
cat("\n100*(meanhat-meanDX)/meanDX\n")
print(100*(meanhat-meanDX)/meanDX) 

######################################################################
# Predict the factor scores
# (Bartlett predictor)

etaBart <- MBcomponents %*% solve(Omega) %*% Boblq %*% (
              solve( t(Boblq) %*% solve(Omega) %*% Boblq ) )

# Differenced factor scores and their covariance matrix
DetaBart <- diff(etaBart, lag=1)
SDE      <- cov(DetaBart)       

# Generally, SDE != PhiOblq. Therefore, transform factors so that
# they have the correct covariance matrix.
# The idea is that if (P SDE P' = PhiOblq), then the sample
# covariance matrix of (P DetaBart) is PhiOblq, just as desired.
#
# Compute Choleski roots of SDE and Phi

RR1 <- chol(SDE)      # upper triangular: SDE = RR1' RR1
RR2 <- chol(PhiOblq)  # ditto
PP  <- t(RR2) %*% solve(t(RR1))

# So now PP * SDE * PP' = RR2' inv(RR1') [RR1'RR1] inv(RR1)] RR2
#                       = RR2'RR2 = PhiOblq

# The true factors in the simulation:
etaTrue <- tframed(etaBart %*% t(PP), tf=tframe(MBcomponents))
#write(t(etaTrue),file="figsDynamic/truefacs.dat",ncolumns=ncol(etaTrue))
              
# Check correctness
DetaTrue <- diff(etaTrue, lag=1)
SigDeta  <- cov(DetaTrue)
CheckFax <- 100 * (SigDeta - PhiOblq)/PhiOblq   # relative difference in %
CheckFax


## Further analytical computations

SystX <- diag(Boblq %*% cov(etaTrue) %*% t(Boblq)) 
CommX <- SystX / (SystX + diag(Psi))

cat("##########################\n")
cat("\nCommunalities (original scale - undifferenced)\n")
cat("##########################\n")
print(CommX)

poppar
corDX = cor(DX)
corDX
cat("Eigen values of the sample correlation matrix of differenced indicators\n")
zz <- eigen(x=corDX,symmetric=TRUE,only.values=TRUE)
print(zz$values)

# Figure not used in paper
#postscript(file="figsDynamic/screeplot.eps", onefile=FALSE)
 plot(zz$values)
#dev.off()

#######################################################################
#
#        Now do simulations.
#
#######################################################################

cat("#############################################################\n")

cat("Starting the single simulation and estimation ...\n")
cat("#############################################################\n")

#  generate  simulated data

rngValue10 <- list(seed=10, kind="Mersenne-Twister", normal.kind="Inversion")

oldrng <- setRNG(rngValue10) # do this to be able to reproduce the result
simBoblq  <- etaTrue %*% t(Boblq) + matrix(rnorm(215*6),215,6) %*% Psi^0.5 
tframe(simBoblq)  <- tframe(etaTrue)

# The above can be done with
#simBoblq  <- simulate(TSFmodel(Boblq, f=etaTrue, 
#             positive.measures=FALSE), Cov=Psi, rng=rngValue10)
# but this produces a slightly different sequence because a few random numbers
# are reserved for initial conditions in dynamic models.


cov(simBoblq)

# par(ask=TRUE)
# Figure not used in paper
#postscript(file="figsDynamic/simplot03.eps", onefile=FALSE)
#tfplot(simBoblq, graphs.per.page=3)
#dev.off()

ML  <- estTSF.ML (simBoblq, 2, BpermuteTarget=Boblq)
QML <- estTSF.ML(simBoblq, 2, BpermuteTarget=Boblq, normalize=FALSE)

cat("mean difference between factors(ML) and factors(QML)")
print(apply(factors(ML) - factors(QML), 2, "mean"))


# Figure not used in paper
 tfplot(factors(ML),  ylab=c("Factor 1", "Factor 2"))

tfplot(diff(factors(ML)),  ylab=c("Factor 1", "Factor 2"))


# Figure 3
if(PSFILES) postscript(file="figsDynamic/haty01a.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
tfplot(explained(ML),
   series=1:3, graphs.per.page=3,
   Title= paste("Observed vs. predicted values\n",
            "(using sample estimates of the parameters)",
            "of indicators 1-3 in the first replication.", sep=""),
   ylab=c("Currency", "Personal cheq.", "NonbankCheq"),
   col=c("black", "black",  "black",   "black",   "black",  "black" ))
if (PSFILES) dev.off()

# Figure not used in paper
#postscript(file="figsDynamic/haty01b.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
tfplot(explained(ML),
   series=4:6, graphs.per.page=3,
   Title= paste("Observed vs. predicted values\n",
            "(using sample estimates of the parameters)",
            "of indicators 4-6 in the first replication.", sep=""),
   ylab=c("N-P demand & notice", "N-P term", "Investment"))
#dev.off()



LBtrue <- solve(t(Boblq) %*% solve(Omega) %*% Boblq) %*%  
                t(Boblq) %*% solve(Omega)

hatftrue <- tframed(simBoblq %*% t(LBtrue), tframe(etaTrue))

# Figure not used in paper
#postscript(file="figsDynamic/dfscore02.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
tfplot(hatftrue, etaTrue, ylab=c("Factor 1", "Factor 2"))
#dev.off()


###################################################################################
# Compute growth rates

cat("Starting plots of growth rates ...\n")

growth <- percentChange(etaTrue )

# Figure not used in paper
  tfplot(growth, percentChange(factors(ML)),
       Title= paste("Growth Rate of True Factors and Bartlett Predictors\n",
	            "computed using ordinary ML parameters", sep=""))

# Figure 2
if(PSFILES) postscript(file="figsDynamic/diff01sub.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
tfplot(diff(etaTrue ), factors(diff(ML)), 
       Title= paste( "Differenced True Factors and Bartlett Predictors\n",
                     "computed using ordinary ML parameters", sep=""),
       ylab=c("Factor 1", "Factor 2"),
       col=c("black", "black",  "black",   "black",   "black",  "black" ),
       start=c(1995,1))
if (PSFILES) dev.off()

# test Figure 2
##png(file="figsDynamic/variationdiff01.png",width = 480, height = 640, pointsize=12, bg = "white")
tf <-tframe(diff(etaTrue))
par(mfcol=c(4,1), mar=c(4, 6, 1, 2))
tfOnePlot(tframed(cbind(factors(diff(ML))[,1], diff(etaTrue)[,1]), tf=tf),
       Title= 
    "Predictors and True Factors\n computed using ordinary ML parameters",
       ylab=c("Factor 1 (differenced)\nPredictor and True"),
       lty= c(2:1),
       col=c("black", "black"))

tfOnePlot(tframed(factors(diff(ML))[,1]- diff(etaTrue)[,1], tf=tf),
       Title=
        "Predictors minus True Factors\ncomputed using ordinary ML parameters",
       ylab=c("Factor 1 (differenced)\nPredictor minus True"),
       col="black")

tfOnePlot(tframed(cbind(factors(diff(ML))[,2], diff(etaTrue)[,2]), tf=tf),
       Title=
        "Predictors and True Factors\ncomputed using ordinary ML parameters",
       ylab=c("Factor 2 (differenced)\nPredictor and True"),
       lty= c(2:1),
       col=c("black", "black"))

tfOnePlot(tframed(factors(diff(ML))[,2]- diff(etaTrue)[,2], tf=tf),
       Title=
	"Predictors minus True Factors\ncomputed using ordinary ML parameters",
       ylab=c("Factor 2 (differenced)\nPredictor minus True"),
       col="black")
#dev.off()

tfplot(growth, percentChange(hatftrue),
  Title= paste("Growth Rate of True Factors and Bartlett Predictors\n",
	       "computed using true parameters", sep=""))

###################################################################################

#                         Monte Carlo 
 
###################################################################################

cat("#############################################################\n")
cat("Starting the main monte carlo experiment...\n")
cat("#############################################################\n")

rngValue10 <- list(seed=10, kind="Mersenne-Twister", normal.kind="Inversion")

#print(summary(mc100))


#EE.ML <- EstEval(TSFmodel(Boblq, f=etaTrue, positive.measures=FALSE),
#    replications=100, rng=rngValue10, quiet=FALSE,
#    simulation.args=list(Cov=Psi, noIC=TRUE),
#    estimation="estTSF.ML", 
#    estimation.args=list(2, BpermuteTarget=Boblq),
#    criterion ="TSFmodel")


EE.ML.1000 <- EstEval(TSFmodel(Boblq, f=etaTrue, positive.measures=FALSE),
    replications=1000, rng=rngValue10, quiet=FALSE,
    simulation.args=list(Cov=Psi),
    estimation="estTSF.ML", 
    estimation.args=list(2, BpermuteTarget=Boblq),
    criterion ="TSFmodel")

# Figure not used in paper
#postscript(file="figsDynamic/sc01sim.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
tfplot(EE.ML.1000, ylab=c("Factor 1", "Factor 2"))
#dev.off()

# Figure 1
if(PSFILES) postscript(file="figsDynamic/Fig-sim01.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
 sm <- summaryStats(EE.ML.1000)  
 tfplot(etaTrue, factors(ML), hatftrue, sm$meanhatf, 
       sm$meanhatf + 1.96 * sm$sdhatf, sm$meanhatf - 1.96 * sm$sdhatf,
   ylab=c("Factor 1", "Factor 2"),
   lty=c("solid", "dashed", "dashed", "dotdash", "dotted","dotted"), 
   lwd=c(   .1,      2,        .1,       1,         1,       1 ), 
   col=c("black", "black",  "black",   "black",   "black",  "black" )) 
if (PSFILES) dev.off()

# Figure not used in paper
#postscript(file="figsDynamic/Dsc01sim.eps", onefile=FALSE)
 tfplot(EE.ML.1000, diff.=TRUE) 
#dev.off()

# Figure not used in paper
#postscript(file="figsDynamic/PCsc01sim.eps", onefile=FALSE)
 tfplot(EE.ML.1000, percentChange.=TRUE)
#dev.off()

# Figure not used in paper
#postscript(file="figsDynamic/Bsc01sim.eps", onefile=FALSE)
 tfplot(EE.ML.1000, PCcentered.=TRUE)
#dev.off()


summary(EE.ML.1000)
#  one rotation is not converging  !!!!!!!!!
i <-0
ii <- NULL
    for (m in EE.ML.1000$result) {
        i <- i+1
        if (!is.null(TSFmodel(m)$dots$rotationConverged) && !TSFmodel(m)$dots$rotationConverged)
            ii <- c(ii, i)
    }
ii
# [1] 626

cat("###############################\n")
cat("      Table 1\n")
cat("###############################\n")

if (is.null(ii)) {
   cat("all rotations converged\n")
   z <- summary(EE.ML.1000)
   print(cbind(Boblq, z$meanhatB.error, z$SDhatB), digits=2)
 }else {
   cat("rotation number which did not converge: ", ii, "\n")
   if (! any(626 == ii ))cat("rotation is not the same as original problem\n")
   zz <- EE.ML.1000
   zz$result[[626]] <- NULL

   # Table 1 
   z <- summary(zz)
   print(cbind(Boblq, z$meanhatB.error, z$SDhatB), digits=2)
   }
   

######################################################
cat("#############################################################\n")
cat("          Comparison with Canadian Monetary aggregates.\n")
cat("#############################################################\n")

# Comparison with Canadian Monetary aggregates.
# (Factors from real, per capita MB numbers.)

######################################################

MBwithML <- estTSF.ML(MBcomponents, 2)

cat("###############################\n")
cat("      Table 2\n")
cat("###############################\n")

z <- summary(MBwithML)
print(cbind(z$B.estimate, z$DstdB.estimate, 1-MBwithML$estimates$uniquenesses), digits=3)

# Figure 4
if(PSFILES) postscript(file="figsDynamic/Fig-app01a.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app01a.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(explained(MBwithML), 
   series=1:3, graphs.per.page=3,
   col=c("black", "black",  "black",   "black",   "black",  "black" ),
   Title= "Explained money indicator data using 2 factors"
   )
if (PSFILES) dev.off()

# Figure not used in paper
#postscript(file="figsDynamic/Fig-app01b.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app01b.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(explained(MBwithML), 
   series=4:6, graphs.per.page=3,
   Title= "Explained money indicator data using 2 factors"
   )
#dev.off()



# Figure 5
if(PSFILES) postscript(file="figsDynamic/Fig-app02a.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app02a.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(explained(diff(MBwithML)),  
   series=1:3, graphs.per.page=3,
   col=c("black", "black",  "black",   "black",   "black",  "black" ),
   Title= "Explained money indicator differenced data using 2 factors"
   )
if (PSFILES) dev.off()

#  this works too
#tfplot(diff(explained(MBwithML)),  
#   series=1:3, graphs.per.page=3,
#   col=c("black", "black",  "black",   "black",   "black",  "black" ),
#   Title= "Explained money indicator differenced data using 2 factors"
#   )

# Figure not used in paper
#postscript(file="figsDynamic/Fig-app02b.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app02b.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(explained(diff(MBwithML)),
   series=4:6, graphs.per.page=3,
   Title="Explained money indicator differenced data using 2 factors"
   )
#dev.off()


#mnF <- apply(factors(MBwithML), 2, mean)
#mnM <- apply(cbind(M1, M2pp),   2, mean)
mnF <- colMeans(factors(MBwithML))
mnM <- colMeans(cbind(M1, M2pp))
#f   <- factors(MBwithML) %*%  diag(mnM/mnF)   this looses tframe but
f <- sweep(factors(MBwithML), 2, mnM/mnF, "*") # this does not


# Figure 6
if(PSFILES) postscript(file="figsDynamic/Fig-app03.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app03.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(cbind(M1, M2pp), f, graphs.per.page=2,
   Title= paste("M1 and M2++ (solid) and scaled Bartlett Predictors\n",
                "computed using ordinary ML parameters (dashed)", sep=""),
   ylab=c("Real perCapita M1", "Real per Capita M2++" ),
   lty=c(     "solid",              "dashed"), 
   col=c(     "black",              "black")) 
if (PSFILES) dev.off()

# Figure 7
if(PSFILES) postscript(file="figsDynamic/Fig-app04.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app04.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(diff(cbind(M1, M2pp)), diff(f), graphs.per.page=2,
   Title=paste("Differenced M1 and M2++ (solid) and scaled Bartlett Predictors\n",
               "computed using ordinary ML parameters (dashed)", sep=""),
   ylab=c("Real perCapita M1", "Real per Capita M2++" ),
   lty=c(     "solid",              "dashed"), 
   col=c(     "black",              "black")) 
if (PSFILES) dev.off()

tfplot(tfwindow(f, start=c(1991,1), end=c(2000,12)))
tfplot(diff(tfwindow(f, start=c(1991,1), end=c(2000,12))))
dim(tfwindow(f, start=c(1991,1), end=c(2000,12)))

# average per period growth over this sub sample
(tfwindow(f, start=c(1991,2), end=c(2000,12))[119,] -
  tfwindow(f, start=c(1991,2), end=c(2000,12))[1,])/119
#  factor 1   factor 2
#  1.334144 430.661481

(tfwindow(f, start=c(1991,1), end=c(2000,12))[120,] -
  tfwindow(f, start=c(1991,1), end=c(2000,12))[1,])/120
#   factor 1    factor 2
#  0.8939283   428.0244655

# average of growth per period over this sub sample
colMeans(diff(tfwindow(f, start=c(1991,1), end=c(2000,12))))
#   factor 1    factor 2
#  0.9014403 431.6213097

#tfplot(diff(factors(MBwithML)), 
#   graphs.per.page=2, Title= "Estimated factors (differenced) using estTSF.ML")

#tfplot(diff(factors(MBwithML)),
#   graphs.per.page=2, Title= "Estimated factors using estTSF.ML")


print(summary(MBwithML))

cat("#############################################################\n")

cat("#######   Sample Sensitivity  - sample starting in 1995\n")

cat("#############################################################\n")

MBwithMLbefore90 <- estTSF.ML(tfwindow(MBcomponents, end=c(1989,12)), 2)
# there is possibly a sign reversal on the next
# MBwithMLafter95 <- estTSF.ML(tfwindow(MBcomponents, start=c(1995,1)), 2)
#  next fixes sign on savings put then seems to reverse sign on trans.
MBwithMLafter95 <- estTSF.ML(tfwindow(MBcomponents, start=c(1995,1)), 2,
                      BpermuteTarget=TSFmodel(MBwithML)$B)
MBwithMLafter95 <- estTSF.ML(tfwindow(MBcomponents, start=c(1995,1)), 2,
                      BpermuteTarget=loadings(MBwithML))
MBwithMLbefore95 <- estTSF.ML(tfwindow(MBcomponents, end=c(1994,12)), 2)
MBwithMLafter01 <- estTSF.ML(tfwindow(MBcomponents, start=c(2001,1)), 2)



# Figure 8
if(PSFILES) postscript(file="figsDynamic/Fig-app05a.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app05a.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(MBcomponents,                explained(MBwithML),   
       explained(MBwithMLbefore90), explained(MBwithMLbefore95), 
       explained(MBwithMLafter95),  explained(MBwithMLafter01), 
   lty=c("solid", "dashed", "dashed", "dashed", "dashed","dashed"), 
   col=c("black", "black",  "black",  "black",  "black", "black" ), 
   series=1:3, graphs.per.page=3,
   Title= paste("Explained money indicator data using 2 factors\n",
               "full sample and sub-samples", sep="")
   )
if (PSFILES) dev.off()

# Figure not used in paper
#postscript(file="figsDynamic/Fig-app05b.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app05b.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(MBcomponents,                explained(MBwithML),   
       explained(MBwithMLbefore90), explained(MBwithMLbefore95), 
       explained(MBwithMLafter95),  explained(MBwithMLafter01), 
   series=4:6, graphs.per.page=3,
   Title= paste("Explained money indicator data using 2 factors\n",
               "full sample and sub-samples", sep="")
   )
#dev.off()


# Figure 9
if(PSFILES) postscript(file="figsDynamic/Fig-app06a.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app06a.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(diff(MBcomponents),                diff(explained(MBwithML)),
       diff(explained(MBwithMLbefore90)), diff(explained(MBwithMLbefore95)), 
       diff(explained(MBwithMLafter95)),  diff(explained(MBwithMLafter01)), 
   lty=c("solid", "dashed", "dashed", "dashed", "dashed","dashed"), 
   col=c("black", "black",  "black",  "black",  "black", "black" ), 
   lwd=c(   .2,      1,        1,        1,        1,       1), 
   series=1:3, graphs.per.page=3,
   Title= paste("Explained money indicator differenced data using 2 factors\n",
              "full sample and sub-samples", sep="")
   )
if (PSFILES) dev.off()

# Figure not used in paper
#postscript(file="figsDynamic/Fig-app06b.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app06b.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(diff(MBcomponents), diff(explained(MBwithML)), 
       diff(explained(MBwithMLbefore90)), diff(explained(MBwithMLbefore95)), 
       diff(explained(MBwithMLafter95)),  diff(explained(MBwithMLafter01)), 
   series=4:6, graphs.per.page=3,
   Title= paste("Explained money indicator differenced data using 2 factors\n",
               "full sample and sub-samples", sep="")
   )
#dev.off()



# Figure 10
if(PSFILES) postscript(file="figsDynamic/Fig-app07.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app07b.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(factors(MBwithML),         factors(MBwithMLbefore90), 
       factors(MBwithMLbefore95), factors(MBwithMLafter95), 
       factors(MBwithMLafter01),
   lty=c("dotdash", "solid", "solid", "solid","solid"), 
   col=c("black",   "black", "black", "black","black" ), 
   lwd=c(   .2,       .2,      .2,      .2,     .2), 
   graphs.per.page=2,
   Title= paste("Bartlett Predictors based on full sample and sub-samples\n",
               "computed using ordinary ML parameters (dashed)", sep="")
   )
if (PSFILES) dev.off()
#summary (MBwithMLafter95)
#loadings(MBwithMLafter95)
#loadings(MBwithML)
tfplot(predict(MBwithMLafter95))
tfplot(predict(MBwithMLafter95, newdata=MBcomponents))
tfplot(factors(MBwithML), predict(MBwithMLafter95, newdata=MBcomponents))
tfplot(diff(factors(MBwithML)), diff(predict(MBwithMLafter95, newdata=MBcomponents)))


# colors make next much easier to see, but below is B&W
tfplot(diff(factors(MBwithML)),         diff(factors(MBwithMLbefore90)), 
       diff(factors(MBwithMLbefore95)), diff(factors(MBwithMLafter95)), 
       diff(factors(MBwithMLafter01)),
   graphs.per.page=2,
   Title= paste("Differenced Bartlett Predictors based on full sample and sub-samples\n",
               "computed using ordinary ML parameters (dashed)", sep="")
   )


# This graph suggests things are very good on all subsamples in differenced 
#  terms,but it is too difficult to make clear without windowing to 
#  small sections(change start), which would mean multiple figures in the paper.
#  It also works much better in colour.
# Figure not used in paper
#postscript(file="figsDynamic/Fig-app08.eps", onefile=FALSE, horizontal=FALSE, width=6, height=8)
#png(file="figsDynamic/Fig-app08b.png",width = 480, height = 480, pointsize=12, bg = "white")
tfplot(diff(factors(MBwithML)),         diff(factors(MBwithMLbefore90)), 
       diff(factors(MBwithMLbefore95)), diff(factors(MBwithMLafter95)), 
       diff(factors(MBwithMLafter01)),
   lty=c("dotdash", "solid", "solid", "solid","solid"), 
#   col=c("black",   "black", "black", "black","black" ), 
   lwd=c(   .2,       .2,      .2,      .2,     .2),
#   start=c(1997,1),
   graphs.per.page=2,
   Title= paste("Differenced Bartlett Predictors based on full sample and sub-samples\n",
               "computed using ordinary ML parameters (dashed)", sep="")
   )
## Figure not used in paper
dev.off()

###################################

cat("#############################################################\n")
cat("finished calculations.R \n")
cat("#############################################################\n")
