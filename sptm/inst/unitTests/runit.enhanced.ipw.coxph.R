library("RUnit")
library("sptm")

test.enhanced.ipw.coxph <- function() {

tolerance=1e-3
# more stringent tolerance for one system to ensure algorithm accuracy
if (R.Version()$system %in% c("x86_64, mingw32")) {
    tolerance=1e-6
} 
RNGkind("Mersenne-Twister", "Inversion")

n=100
beta= c(log(.5), log(.7), log(1.2)) 
t0=2.9999
init = c(log(0.0373*t0),beta)   

dat = sim.fong(n, family="PH", beta, random.censoring="0", design="CC", auxiliary="weak", var.S=1, var.W=1, prevalence=0.1, non.adherence.ratio=0.15, seed=5) 
dat$indicators = !is.na(dat$s)
checkEqualsNumeric(mean(dat$ft), 49.44573, tol=tolerance)

fit = enhanced.ipw.coxph (Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, subset=dat$indicators, imputation.formulae=s~w, verbose=3)
checkEqualsNumeric(coef(fit), c(-0.5218921, -0.2329766, -3.4640868), tol=tolerance)


}
