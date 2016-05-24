library("RUnit")
library("sptm")

test.stm <- function() {

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



#################################################################################################
# PH

# seed=5 is chosen so that not all cases have phase 2 variable s due to non-adherence
dat = sim.fong(n, family="PH", beta, random.censoring="0", design="CC", auxiliary="weak", var.S=1, var.W=1, prevalence=0.1, non.adherence.ratio=0.15, seed=5) 
checkEqualsNumeric(mean(dat$ft), 49.44573, tol=tolerance)


# stm
est   = stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PH", t0=t0, init=init, var.est.type="1", verbose=3)
checkEqualsNumeric(mean(est), -0.1835283, tol=tolerance)

# null init
est.1 = stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PH", t0=t0, init=NULL, var.est.type="1", verbose=3)
checkEqualsNumeric(mean(est.1), -0.1835159, tol=tolerance)

# stm with var.est.type=2
est.2 = stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PH", t0=t0, init=init, var.est.type="2", verbose=3)
checkEqualsNumeric(mean(est.2), -0.1819624, tol=tolerance)

# the following call fails on my laptop: Matrix package is not available for R-devel
# stm.cal
est.3 = stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PH", t0=t0, init=init, var.est.type="1", verbose=3, imputation.formula=s~w)
checkEqualsNumeric(mean(est.3), -0.2612065, tol=tolerance)



#################################################################################################
# PO

dat = sim.fong(n, family="PO", beta, random.censoring="0", design="CC", auxiliary="weak", var.S=1, var.W=1, prevalence=0.1, non.adherence.ratio=0.15, seed=5) 
checkEqualsNumeric(mean(dat$ft), 36.63459, tol=tolerance)

# stm
est   = stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PO", t0=t0, init=init, var.est.type="1", verbose=3)
checkEqualsNumeric(mean(est), -0.1698263, tol=tolerance)



#################################################################################################
# P2

dat = sim.fong(n, family="P2", beta, random.censoring="0", design="CC", auxiliary="weak", var.S=1, var.W=1, prevalence=0.1, non.adherence.ratio=0.15, seed=5) 
checkEqualsNumeric(mean(dat$ft), 34.25107, tol=tolerance)

# stm
est   = stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="P2", t0=t0, init=init, var.est.type="1", verbose=3)
checkEqualsNumeric(mean(est), -0.12702, tol=tolerance)


}

##################################################################################################
## performance
#
#dat = sim.fong(n=200, family="PH", beta, random.censoring="0", design="CC", auxiliary="weak", var.S=1, var.W=1, prevalence=0.1, non.adherence.ratio=0.15, seed=1) 
#system.time({stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PH", t0=t0, init=init, var.est.type="1", verbose=FALSE, imputation.formula=s~w)})
##13
#
#dat = sim.fong(n=200, family="PO", beta, random.censoring="0", design="CC", auxiliary="weak", var.S=1, var.W=1, prevalence=0.1, non.adherence.ratio=0.15, seed=1) 
#system.time({stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PO", t0=t0, init=init, var.est.type="1", verbose=FALSE, imputation.formula=s~w)})
##12
#
#dat = sim.fong(n=200, family="P2", beta, random.censoring="0", design="CC", auxiliary="weak", var.S=1, var.W=1, prevalence=0.1, non.adherence.ratio=0.15, seed=1) 
#system.time({stm(Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="P2", t0=t0, init=init, var.est.type="1", verbose=FALSE, imputation.formula=s~w)})
##19
