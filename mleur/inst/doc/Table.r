#Source: TableSnow.R
#CPU CUSTOMIZATION
#Using intel i7 with 6 cores and 12 compute nodes
nslaves<-12
#
#about 2.8 hrs for 25,000
NSIM <- 25*10^3
#
#move to data directory
require("rwm")
savews("mleur/snow")
#load required libraries
require("mleur")
require("snow")
require("rlecuyer")
require("xtable")
#
#System Requirements
#requires package snow for multicore computing
#multicore computer. See below CPU CUSTOMIZATION - number of cores.
#
#OUTPUT FILES. 
#User dependent, where to put the output?
#Set base name - must end in /
#This is system dependent - Windows or Mac
#Windows or MAC versions
if (.Platform$OS.type=="windows") 
  BASE <- "d:/r/2011/mleur/snow/" else 
  BASE <- "/Users/aim/R/2011/mleur/snow/"
#Names 
TextOutput <- "TableMLEUR.txt" 
TexOutput  <- "TableMLEUR.tex"
outName    <- "out"
Title <- "Unit Root Simulations"
ShortTitle <- "MLEUR-Aug13"
OutFileTxt  <- paste(BASE, "Table",  ShortTitle, ".txt", sep="")
OutFileTex  <- paste(BASE, "Table",  ShortTitle, ".tex", sep="")
ansFile  <- paste(BASE, "ans",  ShortTitle, ".rda", sep="")
outFile  <- paste(BASE, "out",  ShortTitle, ".rda", sep="")
GantChartFile <- paste(BASE, "GantChart", ShortTitle, ".pdf", sep="")
#
wnoises <- c("normal","t","stable","GARCH11")
phis <- c(0.65, 0.85, 0.9, 0.95, 1.0)
NS <- c(60, 70, 80, 90, 100)
NumPars <- length(wnoises)*length(phis)*length(NS)
x <- vector("list", NumPars)
it <- 0
for (i in 1:length(wnoises))
    for (j in 1:length(phis))
        for (k in 1:length(NS)){
            it <- it+1
            x[[it]] <- list(wnoise = wnoises[i], phi=phis[j], n=NS[k])
        }
nout <- length(x)
#
#one way to improve the efficiency is to 
irand <- sample(1:nout, size=nout, replace=FALSE)
irandPerm <- order(irand)
xcopy <- x
for (i in 1:nout)
    x[[i]] <- xcopy[[irand[i]]]
#
OneIt <- function(x){
  phi <- x$phi
  n <- x$n
  wn <- x$wnoise
  GetPower(phi=phi, n=n, NSIM=NSIM, tests=c("DF", "MLEp", "MLEn"), noiseDist=wn)
}
#
cl <- makeCluster(spec=nslaves, type="SOCK")
clusterSetupRNG(cl, seed=rep(775123,6))
#initialized seed to a fixed default
#Export variables
clusterExport(cl, list("x", "NSIM", "OneIt"))
#Export library
clusterEvalQ(cl, require(mleur))
date()
#parallel lapply
s <- snow.time(OUT<-parLapply(cl, x=x, fun=OneIt))
stopCluster(cl) #stop cluster
date()
#snow times
s
#Gantt chart
graphics.off()
plot(s)
title(sub=paste("Randomized Order - Total Elapsed Time: ", round(s[[1]],2), ", NSIM = ", NSIM))
#adjust raw output
OUT2 <- OUT
for (i in 1:nout)
    OUT[[i]] <- OUT2[[irandPerm[i]]]
#
nout <- length(OUT)
ans1 <- character(nout)
ans2 <- matrix(numeric(5*nout), nrow=nout, ncol=5)
ans2 <- as.data.frame.matrix(ans2)
names(ans2) <- c("n", "phi", "DF", "MLEp", "MLEn")
for (i in 1:nout){
    a <- OUT[[i]]
    ans1[i]       <- a$noiseDist 
    ans2[i,] <- c(a$n, a$phi, (a$power)*100)
}
ans <- data.frame(dist=ans1, n=ans2$n, phi=ans2$phi, DF=ans2$DF, MLEn=ans2$MLEn, MLEp=ans2$MLEp)
#output
save(NSIM, ans,  file=ansFile)
write.table(ans, file=OutFileTxt, row.names=FALSE)
#latex table
MOE <- 100*1.96*0.5/sqrt(NSIM)
#latex table
Title <- paste(Title, "0.95 level MOE for table percentage= ", round(MOE,2))
tb <- xtable(ans,  digits=c(0, 0, 0, 2, 1, 1, 1), caption=Title) 
print(tb, include.rownames=FALSE, include.colnames=TRUE, file=OutFileTex)
#Gantt chart
pdf(file=GantChartFile)
plot(s)
title(sub=paste("Total Elapsed Time: ", round(s[[1]],2), ", NSIM = ", NSIM))
dev.off()
