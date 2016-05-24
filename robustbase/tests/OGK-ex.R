library(robustbase)

## minimal testing only
data(ruspini, package = "cluster")

rub1 <- covOGK(ruspini, 1, scaleTau2, covGK, hard.rejection, consistency=FALSE)
rub2 <- covOGK(ruspini, 2, scaleTau2, covGK, hard.rejection, consistency=FALSE)

AE <- function(x,y) all.equal(x,y, tolerance = 2e-15)
## The following test is already fulfilled by Kjell Konis'  original code:
stopifnot(AE(c(rub1$wcov)[c(1,3:4)],
             c(917.99893333333, 94.9232, 2340.319288888888)),
          all.equal(rub1$wcov, rub2$wcov, tolerance=0)
          ,
          AE(c(rub1$cov)[c(1,3:4)],
             c(923.5774514441657, 91.5385216376565, 2342.4556232436971))
          ,
          AE(c(rub2$cov)[c(1,3:4)],
             c(927.2465953711782, 91.8009184487779, 2346.5790105548940))
          )

data(milk)
cM1 <- covOGK(milk, 1, sigmamu = scaleTau2, weight.fn = hard.rejection)
cM2 <- covOGK(milk, 2, sigmamu = scaleTau2, weight.fn = hard.rejection)

symnum(cov2cor(cM1 $cov))
symnum(cov2cor(cM2 $cov))
symnum(cov2cor(cM1 $wcov))
symnum(cov2cor(cM2 $wcov))

cMQn  <- covOGK(milk, sigmamu = s_Qn, weight.fn = hard.rejection)
cMSn  <- covOGK(milk, sigmamu = s_Sn, weight.fn = hard.rejection)
cMiqr <- covOGK(milk, sigmamu = s_IQR, weight.fn = hard.rejection)
cMmad <- covOGK(milk, sigmamu = s_mad, weight.fn = hard.rejection)

as.dist(round(cov2cor(cMQn$wcov), 3))
as.dist(round(cov2cor(cMSn$wcov), 3))
as.dist(round(cov2cor(cMiqr$wcov), 3))
as.dist(round(cov2cor(cMmad$wcov), 3))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
