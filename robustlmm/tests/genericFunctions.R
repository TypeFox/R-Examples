## test the availability of generic functions.
ae <- function(target, current, ...) {
    ret <- all.equal(target, current, ...)
    if (isTRUE(ret)) return(ret)
    print(ret)
    stop("Objects not equal")
}
chgClass <- function(object) {
    class(object) <- sub("(merMod|mer)", "rlmerMod", class(object))
    object
}

require(robustlmm)


set.seed(3)
sleepstudy2 <- within(sleepstudy, {
    Group <- letters[1:4]
    Covar <- rnorm(180)
})
rfm <- rlmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
             rho.e = cPsi, rho.b = cPsi, doFit=FALSE)
fm <- lmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2)    
## all three print the same:
print(rfm)
show(rfm)
summary(rfm)

## object information
## ae(df.residual(fm), df.residual(rfm))
ae(formula(fm), formula(rfm))
stopifnot(isLMM(rfm))
stopifnot(isREML(rfm))
ae(model.frame(fm), model.frame(rfm))
ae(model.matrix(fm), model.matrix(rfm), check.attributes=FALSE)
nobs(rfm)
## ae(getInitial(fm), getInitial(rfm))
ae(terms(fm), terms(rfm))
weights(rfm)

## basic accessors for the results
ae(chgClass(coef(fm)), coef(rfm))
## dummy.coef(rfm)
stopifnot(inherits(try(deviance(rfm), silent=TRUE), "try-error"))
stopifnot(inherits(try(extractAIC(rfm), silent=TRUE), "try-error"))
family(rfm)
## after version 1.1 fitted values are named
if (packageVersion("lme4") > "1.1") ae(fitted(fm), fitted(rfm))
ae(fixef(fm), fixef(rfm))
stopifnot(inherits(try(logLik(rfm), silent=TRUE), "try-error"))
ae(chgClass(ranef(fm)), ranef(rfm))
## after version 1.1 fitted values are named
if (packageVersion("lme4") > "1.1") ae(resid(fm), resid(rfm))
ae(sigma(fm), sigma(rfm))
## weighted.residuals(rfm)

## var-covar methods
## VarCorr(rfm)
ae(chgClass(VarCorr(fm)), VarCorr(rfm))
ae(vcov(fm), vcov(rfm), tolerance=1e-4)
## vcov(rfm)

## confidence intervals
## confint(rfm)

## other (deprecated?)
theta(rfm)

## other methods
getInfo(rfm)
compare(fm, rfm)
update(rfm, ~ . + Covar)

## predict method (various examples)
ae(predict(fm), predict(rfm))
if (packageVersion("lme4") > "1.1") {
    ae(predict(fm,re.form=NA), predict(rfm,re.form=NA))
    newdata <- with(sleepstudy, expand.grid(Subject=unique(Subject),
                                            Days=3:5, Group=letters[1:2]))
    ae(predict(fm,newdata), predict(rfm,newdata))
    ae(predict(fm,newdata,re.form=NA), predict(rfm,newdata,re.form=NA))
    ae(predict(fm,newdata,re.form=~(1|Subject)), predict(rfm,newdata,re.form=~(1|Subject)))
}
