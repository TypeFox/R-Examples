## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcovNA))

alpha <- 0.55

data(phosphor); x <- y <- phosphor[,1:2]; x[10,2] <- NA; x[15,1] <- NA
mcdc <- CovMcd(y)                                           # mcdc - complete
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)    # mcd - robust sequentioal imputation + MCD
mcdna <- CovNAMcd(x)                                        # mcdna - norm + MCD
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)

data(heart); x <- y <- heart; x[10,2] <- NA; x[2,1] <- NA
mcdc <- CovMcd(y)
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)
mcdna <- CovNAMcd(x)
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)

data(starsCYG); x <- y <- starsCYG; x[10,2] <- NA; x[2,1] <- NA; x[33,1] <- NA; x[41,1] <- NA
mcdc <- CovMcd(y)
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)
mcdna <- CovNAMcd(x)
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)

data(stackloss); x <- y <- stack.x; x[10,2] <- NA; x[6,1] <- NA; x[13,3] <- NA
mcdc <- CovMcd(y)
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)
mcdna <- CovNAMcd(x)
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)

data(coleman); x <- y <- data.matrix(subset(coleman, select = -Y)); x[5,2] <- NA; x[8,4] <- NA; x[13,3] <- NA
mcdc <- CovMcd(y)
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)
mcdna <- CovNAMcd(x)
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)

data(salinity); x <- y <- data.matrix(subset(salinity, select = -Y)); x[1,2] <- NA; x[8,3] <- NA; x[13,3] <- NA
mcdc <- CovMcd(y)
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)
mcdna <- CovNAMcd(x)
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)

data(wood); x <- y <- data.matrix(subset(wood, select = -y)); x[1,2] <- NA; x[10,3] <- NA; x[13,4] <- NA
mcdc <- CovMcd(y)
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)
mcdna <- CovNAMcd(x)
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)

data(hbk); x <- y <- data.matrix(subset(hbk, select = -Y)); x[30,2] <- NA; x[40,3] <- NA; x[17,3] <- NA
mcdc <- CovMcd(y)
ximp <- impSeq(x); mcds <- CovMcd(ximp)                     # mcd - sequential imputation + MCD
ximp <- impSeqRob(x, alpha=alpha); mcd <- CovMcd(ximp$x)
mcdna <- CovNAMcd(x)
as.vector(which(mcdc@wt==0))
as.vector(which(mcds@wt==0))
as.vector(which(mcd@wt==0))
as.vector(which(mcdna@wt==0))
##cbind(mcdc@wt, mcdna@wt, mcd@wt)
