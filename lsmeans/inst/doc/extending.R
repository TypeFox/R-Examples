### R code from vignette source 'extending.rnw'

###################################################
### code chunk number 1: extending.rnw:40-42
###################################################
options(show.signif.stars=FALSE, prompt="R> ", continue="   ")
set.seed(271828)


###################################################
### code chunk number 2: extending.rnw:53-54 (eval = FALSE)
###################################################
## help("extending-lsmeans", package="lsmeans")


###################################################
### code chunk number 3: extending.rnw:62-66
###################################################
fake = expand.grid(rep = 1:5, A = c("a1","a2"), B = c("b1","b2","b3"))
fake$y = c(11.46,12.93,11.87,11.01,11.92,17.80,13.41,13.96,14.27,15.82,
           23.14,23.75,-2.09,28.43,23.01,24.11,25.51,24.11,23.95,30.37,
           17.75,18.28,17.82,18.52,16.33,20.58,20.55,20.77,21.21,20.10)


###################################################
### code chunk number 4: extending.rnw:72-77
###################################################
library(MASS)
fake.rlm = rlm(y ~ A * B, data = fake)

library(lsmeans)
lsmeans(fake.rlm, ~B | A)


###################################################
### code chunk number 5: extending.rnw:83-84
###################################################
fake.lts = ltsreg(y ~ A * B, data = fake)


###################################################
### code chunk number 6: extending.rnw:89-90
###################################################
lsmeans:::recover.data.lm


###################################################
### code chunk number 7: extending.rnw:93-94
###################################################
recover.data.lqs = lsmeans:::recover.data.lm


###################################################
### code chunk number 8: extending.rnw:97-99
###################################################
rec.fake = recover.data(fake.lts)
head(rec.fake)


###################################################
### code chunk number 9: extending.rnw:113-114
###################################################
args(lsmeans:::lsm.basis.lm)


###################################################
### code chunk number 10: extending.rnw:121-122
###################################################
MASS:::predict.lqs


###################################################
### code chunk number 11: extending.rnw:126-137
###################################################
lsm.basis.lqs = function(object, trms, xlev, grid, ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = coef(object)
    Xmat = model.matrix(trms, data=object$model)
    V = rev(object$scale)[1]^2 * solve(t(Xmat) %*% Xmat)
    nbasis = matrix(NA)
    dfargs = list(df = nrow(Xmat) - ncol(Xmat))
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs)
}


###################################################
### code chunk number 12: extending.rnw:141-142
###################################################
lsmeans(fake.lts, ~ B | A)


###################################################
### code chunk number 13: extending.rnw:153-154 (eval = FALSE)
###################################################
## nbasis = estimability::nonest.basis(Xmat)


###################################################
### code chunk number 14: extending.rnw:194-197
###################################################
form = ~ data$x + data[[5]]
base::all.vars(form)
lsmeans::.all.vars(form)


###################################################
### code chunk number 15: extending.rnw:204-205
###################################################
.get.offset(terms(~ speed + offset(.03*breaks)), head(warpbreaks))


