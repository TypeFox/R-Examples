library("RUnit")
library("krm")


test.krm.score.test <- function() {

tolerance=1e-3
# more stringent tolerance for one system to ensure algorithm accuracy
if (R.Version()$system %in% c("x86_64, mingw32")) {
    tolerance=1e-6
} 
RNGkind("Mersenne-Twister", "Inversion")


dat=sim.liu.2008(n=100, a=0, seed=1)
z=as.matrix(subset(dat, select=c(z.1,z.2,z.3,z.4,z.5)))
rho=1
K=kyotil::getK(z,kernel="rbf",para=rho^-2)

test = krm.score.test (y~x, dat, K, regression.type="logistic") 
checkEqualsNumeric(test, c(0.3708008, 0.2130736, 0.3552138, 0.2102588), tolerance = tolerance)

test = krm.score.test (y~x, dat, K, regression.type="linear") 
checkEqualsNumeric(test, c(0.3492788,        NA, 0.3327407,        NA), tolerance = tolerance)


# performance in time
#dat=sim.liu.2008(n=50, a=0, seed=1)
#z=as.matrix(subset(dat, select=c(z.1,z.2,z.3,z.4,z.5)))
#rho=1
#K=kyotil::getK(z,kernel="rbf",para=rho^-2)
#system.time({
#    krm.score.test (y~x, dat, K, regression.type="logistic", verbose=TRUE) 
#})
## 0.0-4 n=50: 0.008 second on gizmo; n=100: 0.028 second on gizmo; n=200: 0.164 second on gizmo
## 0.0-5 n=50: 0.008 second on gizmo; n=100: 0.012 second on gizmo; n=200: 0.076 second on gizmo
## a default krm.most call makes 10x2000 calls to krm.score.test 
#
#
}
