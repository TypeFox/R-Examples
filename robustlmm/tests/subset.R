## Test disabled
quit()

## testing subset
require(robustlmm)

fit1 <- rlmer(Yield ~ (1 | Batch), Dyestuff, subset=Batch != "F")
fit2 <- rlmer(Yield ~ (1 | Batch), subset(Dyestuff, Batch != "F"))

stopifnot(all.equal(getInfo(fit1)[-1], getInfo(fit2)[-1]))
