.proctime00 <- proc.time()
require(VLMC)

## Currently, there's also an extensive AIC example in
## ../man/OZrain.Rd, i.e.  example(OZrain)

data(bnrf1)
tit <- paste("VLMC for BNRF1 Epstein-Barr, N =", length(bnrf1EB))
N <- length(bnrf1EB)

nC <- length(cutoffs <- c(seq(2,6, by = .125),
                          seq(6.5, 10, by = 0.5),
                          seq(11, 15, by = 1)))

 ABIC.EB <- matrix(NA, 2, nC, dimnames = list(c("AIC", "BIC"), NULL))
sizes.EB <- matrix(NA, 4, nC)
for(ic in 1:nC) {
    cuto <- cutoffs[ic]
    v <- vlmc(bnrf1EB, cutoff = cuto)
    ABIC.EB [, ic] <- AIC(v, k = c(2, log(N)))
    sizes.EB[, ic] <- v$size
}
rownames(sizes.EB) <- names(v$size)

## Hmm... BIC seems completely off
## Table
cbind(cutoff = cutoffs, t(ABIC.EB), t(sizes.EB))

if(!dev.interactive(orNone=TRUE)) pdf("AIC-etc.pdf")
par(mfrow = c(2,1), mgp = c(1.5, .6, 0), mar = .1 + c(4,3,3,1))
plot(cutoffs, ABIC.EB[1,], type = "o", main = paste("AIC of", tit), log = 'xy',
     sub = paste("qchisq(0.95, 4 -1) / 2 = ",format(qchisq(0.95, 3) / 2)))
plot(cutoffs, ABIC.EB[2,], type = "o", main = paste("BIC of", tit), log = 'xy')


## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
