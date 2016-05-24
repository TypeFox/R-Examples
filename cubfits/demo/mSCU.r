suppressMessages(library(cubfits, quietly = TRUE))

b <- b.Init$roc
phi.Obs <- ex.train$phi.Obs
y <- ex.train$y
y.list <- convert.y.to.list(y)
mSCU <- calc_scu_values(b, y.list, phi.Obs)$mSCU
plot(mSCU, log10(phi.Obs), main = "Expression vs mSCU",
     xlab = "mSCU", ylab = "Expression (log10)")

### Compare with CAI with weights seqinr::cubtab$sc.
library(seqinr)
w <- caitab$sc
names(w) <- codon.low2up(rownames(caitab))
CAI <- calc_cai_values(y, y.list, w = w)$CAI

plot(mSCU, CAI, main = "CAI vs mSCU",
     xlab = "mSCU", ylab = "CAI")
