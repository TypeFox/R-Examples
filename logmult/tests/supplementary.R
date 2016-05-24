## Variations based on Becker & Clogg (1989), cf. Becker-Clogg1989.R

library(logmult)
data(color)

caithness <- color[,,"Caithness"]


## Adding the table again as supplementary rows
rsup <- caithness
csup <- caithness
rownames(rsup) <- paste(rownames(caithness), "2")
colnames(csup) <- paste(colnames(caithness), "2")

names(dimnames(rsup)) <- names(dimnames(csup)) <- names(dimnames(caithness))

mA1r <- rc(caithness, 1, rowsup=rsup, weighting="u", start=NA)
mA1c <- rc(caithness, 1, colsup=csup, weighting="u", start=NA)
mA1 <- rc(caithness, 1, rowsup=rsup, colsup=csup, weighting="u", start=NA)

rscoresA1 <- mA1$assoc$row[1:4, 1, 1] * sqrt(mA1$assoc$phi[1, 1])
cscoresA1 <- mA1$assoc$col[1:5, 1, 1] * sqrt(mA1$assoc$phi[1, 1])

rmA1 <- gnm(Freq ~ Eye + Hair + Eye:cscoresA1[col(rsup)],
               data=as.table(rsup), family=poisson)
cmA1 <- gnm(Freq ~ Eye + Hair + Hair:rscoresA1[row(csup)],
               data=as.table(csup), family=poisson)

rscA1 <- getContrasts(rmA1, pickCoef(rmA1, "cscores"), ref="mean")$qvframe[,1]
cscA1 <- getContrasts(cmA1, pickCoef(cmA1, "rscores"), ref="mean")$qvframe[,1]

stopifnot(all.equal(rscA1, rscoresA1, mA1$assoc$row[5:8,,] * sqrt(mA1$assoc$phi[1, 1]),
                    tolerance=1e-6, check.attributes=FALSE))
stopifnot(all.equal(cscA1, cscoresA1, mA1$assoc$col[6:10,,] * sqrt(mA1$assoc$phi[1, 1]),
                    tolerance=1e-6, check.attributes=FALSE))


## Merging some categories and re-adding existing ones
rsup <- rbind(Fairer=colSums(caithness[1:2,]),
              Darker=colSums(caithness[3:4,]),
              Dark2=caithness[4,])

csup <- cbind(FairRed=rowSums(caithness[,1:2]),
              Darker=rowSums(caithness[,3:5]),
              Black2=caithness[,5])

names(dimnames(rsup)) <- names(dimnames(csup)) <- names(dimnames(caithness))

# With one dimension
mB1r <- rc(caithness, 1, rowsup=rsup, weighting="u", start=NA)
mB1c <- rc(caithness, 1, colsup=csup, weighting="u", start=NA)
mB1 <- rc(caithness, 1, rowsup=rsup, colsup=csup, weighting="u", start=NA)

rscoresB1 <- mB1$assoc$row[1:4, 1, 1] * sqrt(mB1$assoc$phi[1, 1])
cscoresB1 <- mB1$assoc$col[1:5, 1, 1] * sqrt(mB1$assoc$phi[1, 1])

rmB1 <- gnm(Freq ~ Eye + Hair + Eye:cscoresB1[col(rsup)],
            data=as.table(rsup), family=poisson)
cmB1 <- gnm(Freq ~ Eye + Hair + Hair:rscoresB1[row(csup)],
            data=as.table(csup), family=poisson)

rscB1 <- getContrasts(rmB1, pickCoef(rmB1, "cscores"), ref="mean")$qvframe[,1]
cscB1 <- getContrasts(cmB1, pickCoef(cmB1, "rscores"), ref="mean")$qvframe[,1]

stopifnot(all.equal(rscB1, mB1$assoc$row[5:7,,] * sqrt(mB1$assoc$phi[1, 1]),
                    tolerance=1e-6, check.attributes=FALSE))
stopifnot(all.equal(cscB1, mB1$assoc$col[6:8,,] * sqrt(mB1$assoc$phi[1, 1]),
                    tolerance=1e-6, check.attributes=FALSE))

# With two dimensions
mB2r <- rc(caithness, 2, rowsup=rsup, weighting="u", start=NA)
mB2c <- rc(caithness, 2, colsup=csup, weighting="u", start=NA)
mB2 <- rc(caithness, 2, rowsup=rsup, colsup=csup, weighting="u", start=NA)

rscoresB2 <- sweep(mB2$assoc$row[1:4,, 1], 2, sqrt(mB2$assoc$phi[1,]), "*")
cscoresB2 <- sweep(mB2$assoc$col[1:5,, 1], 2, sqrt(mB2$assoc$phi[1,]), "*")

rmB2 <- gnm(Freq ~ Eye + Hair + Eye:cscoresB2[col(rsup), 1] + Eye:cscoresB2[col(rsup), 2],
            data=as.table(rsup), family=poisson)
cmB2 <- gnm(Freq ~ Eye + Hair + Hair:rscoresB2[row(csup), 1] + Hair:rscoresB2[row(csup), 2],
            data=as.table(csup), family=poisson)

rscB2 <- cbind(getContrasts(rmB2, pickCoef(rmB2, "cscores.*1\\]"), ref="mean")$qvframe[,1],
               getContrasts(rmB2, pickCoef(rmB2, "cscores.*2\\]"), ref="mean")$qvframe[,1])
cscB2 <- cbind(getContrasts(cmB2, pickCoef(cmB2, "rscores.*1\\]"), ref="mean")$qvframe[,1],
               getContrasts(cmB2, pickCoef(cmB2, "rscores.*2\\]"), ref="mean")$qvframe[,1])

stopifnot(all.equal(rscB2, sweep(mB2$assoc$row[5:7,,], 2, sqrt(mB2$assoc$phi[1,]), "*"),
                    tolerance=1e-6, check.attributes=FALSE))
stopifnot(all.equal(cscB2, sweep(mB2$assoc$col[6:8,,], 2, sqrt(mB2$assoc$phi[1,]), "*"),
                    tolerance=1e-6, check.attributes=FALSE))

# With two dimensions and marginal weighting
mB3r <- rc(caithness, 2, rowsup=rsup, weighting="m", start=NA)
mB3c <- rc(caithness, 2, colsup=csup, weighting="m", start=NA)
mB3 <- rc(caithness, 2, rowsup=rsup, colsup=csup, weighting="m", start=NA)

rscoresB3 <- sweep(mB3$assoc$row[1:4,, 1], 2, sqrt(mB3$assoc$phi[1,]), "*")
cscoresB3 <- sweep(mB3$assoc$col[1:5,, 1], 2, sqrt(mB3$assoc$phi[1,]), "*")

rmB3 <- gnm(Freq ~ Eye + Hair + Eye:cscoresB3[col(rsup), 1] + Eye:cscoresB3[col(rsup), 2],
            data=as.table(rsup), family=poisson)
cmB3 <- gnm(Freq ~ Eye + Hair + Hair:rscoresB3[row(csup), 1] + Hair:rscoresB3[row(csup), 2],
            data=as.table(csup), family=poisson)

rscB3 <- cbind(getContrasts(rmB3, pickCoef(rmB3, "cscores.*1\\]"),
                            ref=rowSums(rsup)/sum(rsup))$qvframe[,1],
               getContrasts(rmB3, pickCoef(rmB3, "cscores.*2\\]"),
                            ref=rowSums(rsup)/sum(rsup))$qvframe[,1])
cscB3 <- cbind(getContrasts(cmB3, pickCoef(cmB3, "rscores.*1\\]"),
                            ref=colSums(csup)/sum(csup))$qvframe[,1],
               getContrasts(cmB3, pickCoef(cmB3, "rscores.*2\\]"),
                            ref=colSums(csup)/sum(csup))$qvframe[,1])

stopifnot(all.equal(rscB3, sweep(mB3$assoc$row[5:7,,], 2, sqrt(mB3$assoc$phi[1,]), "*"),
                    tolerance=1e-6, check.attributes=FALSE))
stopifnot(all.equal(cscB3, sweep(mB3$assoc$col[6:8,,], 2, sqrt(mB3$assoc$phi[1,]), "*"),
                    tolerance=1e-6, check.attributes=FALSE))

