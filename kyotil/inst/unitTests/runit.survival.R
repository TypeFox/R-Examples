library("RUnit")
library("survival")
library("kyotil")

test.timedep <- function() {

tolerance=1e-5
if(file.exists("D:/gDrive/3software/_checkReproducibility")) tolerance=1e-6
RNGkind("Mersenne-Twister", "Inversion")


## check actual incidence.density
#seed=2
#age.sim="tvaryinggroup" # default
#y=2; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.01, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=2; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.05, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.0015, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.01, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.05, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#
#seed=2
#age.sim="baselinegroup"
#y=2; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.01, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=2; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.05, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=.0015,age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.01, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.05, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#
#seed=2
#age.sim="continuous"
#y=2; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.01, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=2; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.05, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=.0015,age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.01, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)
#y=5; dat=sim.dat.tvarying(n=10000,followup.length=y, incidence.density=0.05, age.sim=age.sim, seed=seed); print(mean(subset(dat, for.non.tvarying.ana, d, drop=T))/y)


age.sim="tvaryinggroup"
dat=sim.dat.tvarying.three(n=200,followup.length=5, incidence.density=0.1, age.sim=age.sim, seed=1)
checkTrue(nrow(dat)==294)
checkEqualsNumeric(mean(dat$X), 3.30203, tolerance=tolerance)

age.sim="bt"
dat=sim.dat.tvarying.three(n=200,followup.length=5, incidence.density=0.1, age.sim=age.sim, seed=1)
checkTrue(nrow(dat)==294)
checkEqualsNumeric(mean(dat$X), 3.287375, tolerance=tolerance)

age.sim="baselinegroup"
dat=sim.dat.tvarying.three(n=200,followup.length=5, incidence.density=0.1, age.sim=age.sim, seed=1)
checkTrue(nrow(dat)==294)
checkEqualsNumeric(mean(dat$X), 3.14795, tolerance=tolerance)

age.sim="continuous"
dat=sim.dat.tvarying.three(n=200,followup.length=5, incidence.density=0.1, age.sim=age.sim, seed=1)
checkTrue(nrow(dat)==299)
checkEqualsNumeric(mean(dat$X), 3.23414, tolerance=tolerance)

dat=sim.dat.tvarying.three(n=6000,followup.length=3, incidence.density=0.05, age.sim=age.sim, seed=1)
f.tvarying = Surv(tstart,tstop,d) ~ trt*agegrp 
f =          Surv(X,d)            ~ trt*baseline.agegrp 
fits=list()
fits[["tvarying"]]=coxph(f.tvarying, dat)
fits[["baseline"]]=coxph(f, subset(dat, for.non.tvarying.ana))
checkEqualsNumeric(mean(coef(fits[[1]])), -0.5773688, tolerance=tolerance)
checkEqualsNumeric(mean(coef(fits[[2]])), -0.6146536, tolerance=tolerance)



#################################################################################
# check make.timedep.dataset

n=3000; followup.length=5; incidence.density=0.015; age.sim="continuous"

dat.0=sim.dat.tvarying.two(n, followup.length, incidence.density, age.sim, seed=1)
dat=subset(dat.0, for.non.tvarying.ana)
dat.timedep = make.timedep.dataset (dat, "X", "d", "baseline.age", 6)
checkTrue (nrow(dat.timedep) == nrow(dat.0))
fit.1=coxph(Surv(tstart,tstop,d) ~ trt*agegrp, dat.0); 
fit.2=coxph(Surv(tstart,tstop,d) ~ trt*.timedep.agegrp, dat.timedep); 
checkEqualsNumeric (coef(fit.1), coef(fit.2), tolerance=tolerance)


dat.0=sim.dat.tvarying.three(n, followup.length, incidence.density, age.sim, seed=1)
dat=subset(dat.0, for.non.tvarying.ana)
dat.timedep = make.timedep.dataset (dat, "X", "d", "baseline.age", 6, 12)
checkTrue (nrow(dat.timedep) == nrow(dat.0))
fit.1=coxph(Surv(tstart,tstop,d) ~ trt*agegrp, dat.0); 
fit.2=coxph(Surv(tstart,tstop,d) ~ trt*.timedep.agegrp, dat.timedep); 
checkEqualsNumeric (coef(fit.1), coef(fit.2), tolerance=tolerance)



# cox.zph.2
fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps, data=ovarian) 
temp.2 <- cox.zph.2(fit) 
checkEqualsNumeric (temp.2$table[,3], c(0.4035486, 0.1236237, 0.1625539), tolerance=tolerance)

}
