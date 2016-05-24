library(robustX)

L1.obj <- function(m, X) sum(sqrt(rowSums((X - rep(m, each = nrow(X)))^2)))

## Really 'simple' 3-D example
X <- rbind(cbind(1:12, round(abs(-5:6)^1.5) %% 7, rep(1:4,3)),
           c(3,0, 20), c(6,1, 200))## 'outliers'

c0 <- system.time(for(i in 1:10) m0 <- L1median(X, method="nlminb",tol =3e-16))
c1 <- system.time(for(i in 1:10) m1 <- L1median(X, method="nlm",  tol = 3e-16))
c2 <- system.time(for(i in 1:10) m2 <- L1median(X, method="HoCr", tol = 3e-16))
c3 <- system.time(for(i in 1:10) m3 <- L1median(X, method="Vardi",tol = 3e-16))
(cM <- rbind(nlminb = c0, nlm = c1, HoCr = c2, Vardi = c3))
mNms <- rownames(cM) # method names

pX <- matrix(NA, ncol(X), length(mNms),
             dimnames = list(NULL, mNms))
## FIXME, once we have common return value
pX[,"nlminb"] <- m0$par
pX[,"nlm"   ] <- m1$estimate
pX[,"HoCr"  ] <- m2
pX[,"Vardi" ] <- m3
t(pX)

apply(pX, 2, function(m) L1.obj(m,X)) - 259.7393299943
## hmm, "nlminb" is slightly "worse" than the others ..

for(j in 1:ncol(pX))
    stopifnot(all.equal(pX[,j], pX[,1 + j %% ncol(pX)], tol = 1e-6))
## hmm:  1e-6  currently needed 'cause of "nlminb"

##--- even much simpler: Examples where  coord.wise.median == data point
(x3 <- cbind(c(0,1,8),0:2))
x5 <- rbind(x3,
            c(1,0),
            c(2,1))
(x5 <- x5[order(x5[,1],x5[,2]), ])
m3.0 <- L1median(x3, method="nlminb", trace=2, tol=3e-16)
m3.1 <- L1median(x3, method="nlm",    trace=2, tol=3e-16)
m3.2 <- L1median(x3, method="HoCrJo", trace=2)
m3.3 <- L1median(x3, method="Vardi",  trace=2)

pX <- matrix(NA, ncol(x3), length(mNms), dimnames = list(NULL, mNms))
pX[,"nlminb"] <- m3.0$par
pX[,"nlm"   ] <- m3.1$estimate
pX[,"HoCr"  ] <- m3.2
pX[,"Vardi" ] <- m3.3
t(pX)
for(j in 1:ncol(pX))
    stopifnot(all.equal(pX[,j], pX[,1 + j %% ncol(pX)], tol = 1e-6))


m5.0 <- L1median(x5, method="nlminb", trace=2, tol=3e-16)
m5.1 <- L1median(x5, method="nlm",    trace=2, tol=3e-16) # fails !!
m5.2 <- L1median(x5, method="HoCrJo", trace=2)
m5.3 <- L1median(x5, method="Vardi",  trace=2)

pX <- matrix(NA, ncol(x5), length(mNms), dimnames = list(NULL, mNms))
pX[,"nlminb"] <- m5.0$par
pX[,"nlm"   ] <- m5.1$estimate
pX[,"HoCr"  ] <- m5.2
pX[,"Vardi" ] <- m5.3
t(pX)

## FIXME; nlm() currently fails here:
pX <- pX[, - which("nlm" == mNms)]

for(j in 1:ncol(pX))
    stopifnot(all.equal(pX[,j], pX[,1 + j %% ncol(pX)], tol = 1e-6))



stopifnot(require(MASS))
##-> lazy loading  data sets

(sl.HC <- L1median (stackloss, trace = TRUE, method = "HoCrJo"))

system.time(rr0 <- L1median(stackloss, method="HoCrJo", tol = 1e-14))
system.time(rr1 <- L1median(stackloss, method="Nelder-M", tol = 1e-14))

system.time(rr2 <- L1median(stackloss, method="BFGS", tol = 1e-14))
## MM: Hmm, this ("CG") now takes MUCH longer (factor 5 - 10 !):
system.time(rr3 <- L1median(stackloss, method="CG", tol = 1e-14))

system.time(rr3.2 <- L1median(stackloss, method="CG", tol = 1e-14, type = 2))

system.time(rr3.3 <- L1median(stackloss, method="CG", tol = 1e-14, type = 3))

## nlm with gradient:
system.time(rr4 <- L1median(stackloss, method="nlm", tol = 1e-16))

##--> fastest! {faster than rr0 by almost factor of two!}
system.time(rr4 <- L1median(stackloss, method="nlm", tol = 1e-16, trace = 2))

mm <- rbind(HoCrJo= rr0, "NelderM"=rr1$par, BFGS=rr2$par,
            CG = rr3$par, CG.2 = rr3.2$par, CG.3= rr3.3$par,
           "nlm.grad" = rr4$estimate)
L1.obj <- function(m, Xmat) sum(sqrt(rowSums((Xmat - rep(m, each = nrow(Xmat)))^2)))
mm <- cbind(mm, Obj = apply(mm, 1, L1.obj, Xmat = stackloss))
print(mm[,"Obj"] - min(mm[,"Obj"]), digits = 4)
## HoCrJo is best (since it iterates too long!) *together*  with nlm;
##   but only significantly wrt "NelderM"
op <- options(digits=12)# more than usual to see (irrelevant) differences:
mm
options(op)

