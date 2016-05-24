## ----knitr_options, echo=FALSE, results=FALSE----------------------------
library(knitr)
opts_chunk$set(fig.width = 12)

## ----loading, include=FALSE----------------------------------------------
library(aSPU)

## ----loading2------------------------------------------------------------
data(SAMD11)
attach(SAMD11)

## ----ZsPsM---------------------------------------------------------------
round(ZsM,3)
PsM

## ----corM----------------------------------------------------------------
round(corSNPM,2)
round(corPheM,2)

## ----ZsPsF---------------------------------------------------------------
round(ZsF,3)
PsF

## ----corF----------------------------------------------------------------
round(corSNPF,2)
round(corPheF,2)

## ----outFZ---------------------------------------------------------------
(outFZC <- MTaSPUsSetC(ZsF, corSNP=corSNPF, corPhe = corPheF,
           pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=100, Ps=FALSE))
(outFZ <- MTaSPUsSet(ZsF, corSNP=corSNPF, corPhe = corPheF,
           pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=100, Ps=FALSE))

## ----wrwr----------------------------------------------------------------
#write.table(ZsF, quote=FALSE, row.names=FALSE, col.names=FALSE, file="ZsF.txt")
#write.table(corPheF, quote=FALSE, row.names=FALSE, col.names=FALSE, file="corPheF.txt")
#write.table(corSNPF, quote=FALSE, row.names=FALSE, col.names=FALSE, file="corSNPF.txt")


## ----outFP---------------------------------------------------------------
(outFPC <- MTaSPUsSetC(PsF, corSNP=corSNPF, corPhe = corPheF,
           pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=100, Ps=TRUE))
(outFP <- MTaSPUsSet(PsF, corSNP=corSNPF, corPhe = corPheF,
           pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=100, Ps=TRUE))

## ----outMPZ--------------------------------------------------------------
(outMPC <- MTaSPUsSetC(PsM, corSNP=corSNPM, corPhe = corPheM,
           pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=100, Ps=TRUE))
(outMZC <- MTaSPUsSetC(ZsM, corSNP=corSNPM, corPhe = corPheM,
           pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=100, Ps=FALSE))

## ----Zsmcors-------------------------------------------------------------
round(ZsM,3)
round(corSNPM,2)
round(corPheM,2)

## ----plots, echo=FALSE---------------------------------------------------
plotG <- function(Ps, zlim = NULL, main = NULL, yt = NULL, title = "SNPs") {        
    log10P <- -log(Ps,10)  
    pos = 1:nrow(log10P)
    y = 1:ncol(log10P)
    log10P <- log10P
    val <- sqrt(seq(0, 1, len=251))
    col <- rgb(1, rev(val), rev(val))

    if(is.null(yt)) {
        yt = -length(pos)/15
    }

    if(is.null(zlim)) {
        maxP <- max(log10P, na.rm=TRUE)
        zlim <- c(0, maxP)
    }
    image.plot(pos, y, log10P, xaxt="n", yaxt="n", ylab="", xlab="",
                    zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n", main = main )
    title(xlab=title, mgp=c(1, 0, 0))
    text(yt,1,"BMI", xpd = TRUE)
    text(yt,2,"Height", xpd = TRUE)
    text(yt,3,"HIP", xpd = TRUE)
    text(yt,4,"WC", xpd = TRUE)
    text(yt,5,"Weight", xpd = TRUE)
    text(yt,6,"WHR", xpd = TRUE)
}

plotLD <- function(ldmatrix, zlim = NULL, main = NULL, yt = NULL, title = "SNPs") {
#    log10P <- -log(Ps,10)
    pos = 1:nrow(ldmatrix)
    y = 1:ncol(ldmatrix)
    val <- sqrt(seq(0, 1, len=251))
    col <- rgb(1, rev(val), rev(val))

    if(is.null(yt)) {
        yt = -length(pos)/15
    }

    if(is.null(zlim)) {
        maxP <- max(ldmatrix, na.rm=TRUE)
        zlim <- c(0, maxP)
    }
    image.plot(pos, y, ldmatrix, xaxt="n", yaxt="n", ylab="", xlab="",
        zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n", main = main )
    title(xlab=title, mgp=c(1, 0, 0))
    title(ylab=title, mgp=c(1, 0, 0))
}


## ----plot_MGAS, echo=FALSE, fig.width=7, fig.height=7--------------------
data(someGs)
par(mfrow = c(2,2))
plotG(someGs$LCORL[[1]], main = "LCORL (P-values)", zlim = c(0,18))
plotG(someGs$RASA2[[1]], main = "RASA2 (P-values)", zlim = c(0,12))
plotLD(abs(someGs$LCORL[[2]]), main = "LCORL (LDmatrix)")
plotLD(abs(someGs$RASA2[[2]]), main = "RASA2 (LDmatrix)")

## ----plot_MT, echo=FALSE, fig.width=7, fig.height=7----------------------
data(someGs)
par(mfrow = c(2,2))
plotG(someGs$STK33[[1]], main = "STK33 (P-values)", zlim = c(0,12))
plotG(someGs$RPGRIP1L[[1]], main = "RPGRIP1L (P-values)", zlim = c(0,12))
plotLD(abs(someGs$STK33[[2]]), main = "STK33 (LDmatrix)")
plotLD(abs(someGs$RPGRIP1L[[2]]), main = "RPGRIP1L (LDmatrix)")

