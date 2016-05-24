library("RUnit")
library("prc")

test.prc <- function() {


RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=5e-4 # i386 has a lower tolerance


# gnls fit
t1=Sys.time()
fit=prc (mtct.eg$V3_BioV3B_2500, 2500, mtct.eg$V3_BioV3B_500, 500, init.method="gnls", opt.method="gnls",verbose=TRUE); fit
checkEqualsNumeric(coef(fit), c(16.5845944, 28048.1001348,    -2.1822061,     0.4606717), tolerance=tolerance)
checkEqualsNumeric(diag(fit$Sigma.hat), c(2.345293e+00, 4.031360e+04, 1.237249e-02, 7.617220e-04), tolerance=tolerance)
t2=Sys.time()
print(t2-t1)

# optim fit
t3=Sys.time()
fit.optim=prc (mtct.eg$V3_BioV3B_2500, 2500, mtct.eg$V3_BioV3B_500, 500, init.method="gnls", opt.method="optim",verbose=TRUE); fit.optim
checkEqualsNumeric(coef(fit.optim), c(17.269514, 28150.094194,    -2.157746,     0.466794), tolerance=tolerance)
checkEqualsNumeric(diag(fit.optim$Sigma.hat), c(1.787694e+00, 4.392289e+04, 1.271869e-02, 8.203423e-04), tolerance=tolerance)
t4=Sys.time()
print(t4-t3)
# optim is only slightly faster

# naive fit
fit.naive=prc (mtct.eg$V3_BioV3B_2500, 2500, mtct.eg$V3_BioV3B_500, 500, init.method="optim", method="naive", verbose=TRUE); fit.naive
checkEqualsNumeric(coef(fit.naive), c(17.9585174, 29509.4983037,    -1.8251283,     0.5626997), tolerance=tolerance)

# prediction
logfi.1000 = predict(fit, new.dilution=1000)
checkEqualsNumeric(mean(logfi.1000), 9.212831, tolerance=tolerance)

# plot and slope
plot(fit,lcol=2)
theta=coef(fit)
x=log(500)
k=5
checkEqualsNumeric(four_pl_prc(theta["c"], theta["d"], theta["b"], theta["f"], x, k=k), 7.803386, tolerance=tolerance)

checkEquals(four_pl_prc(theta["c"], theta["d"], theta["b"], theta["f"], xx=log(15), k=k), NaN) # under log(c)
checkEqualsNumeric(four_pl_prc(theta["c"], theta["d"], theta["b"], theta["f"], xx=log(30000), k=k), 10.24355, tolerance=tolerance) # above log)(d)

slope=s.dot.f (theta["c"], theta["d"], theta["b"], theta["f"], x, k=k)
abline(reg=c(four_pl_prc(theta["c"], theta["d"], theta["b"], theta["f"], x, k=k) - slope*x,slope), col=4)
checkEqualsNumeric(slope, 1.022359, tolerance=tolerance)
lines(fit.optim)

# 3P prc fit
fit.3P=prc (mtct.eg$V3_BioV3B_2500, 2500, mtct.eg$V3_BioV3B_500, 500, verbose=TRUE, model="3P"); fit.3P
checkEqualsNumeric(coef(fit.3P), c(21.259201, 31946.684225,    -1.161139,     1.000000), tolerance=tolerance)



}
