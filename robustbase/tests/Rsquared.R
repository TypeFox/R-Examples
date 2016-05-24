require(robustbase)

set.seed(17)# reproducibility!
## to check:
## - for the empty model
summary(lmrob(Y ~ 0, coleman))
## - with and without an intercept in the  model
summary(lmrob(Y ~ 1, coleman))
writeLines(sfm <- capture.output(
                       summary(lmrob(Y ~ ., coleman)))) # and this must be "identical":
sfm2 <- capture.output(summary(lmrob(Y ~ ., coleman, model=FALSE, x=FALSE, y=FALSE)))
iCall <- grep("lmrob.*coleman", sfm)# the only line that differs
stopifnot(sfm[-iCall] == sfm2[-iCall])
## w/o intercept:
summary(lmrob(Y ~ . - 1, coleman, model=FALSE, x=FALSE, y=FALSE))

## - when prior-weights are included
wts <- c(rep(0.05, 10), rep(2, 10))
summary(lmrob(Y ~ . - 1, coleman, model=FALSE, x=FALSE, y=FALSE,
              weights = wts))
## - should work for object with NA in the coefficients, and
## - should work for object with NA in the observations --> both in ./NAcoef.R

## check equality with lm() for classical model
test <- function(formula, data,
                 items=c("coefficients", "residuals", "df", "scale",
                         "r.squared", "adj.r.squared", "weights"),
                 tol = 1e-4, ...)
{
    lmrCtrl <- lmrob.control(psi = "hampel", tuning.psi = c(1000, 2000, 3000),
                             method="SMDM", ...)
    sc <- summary(lm   (formula, data))
    sr <- summary(lmrob(formula, data, control= lmrCtrl))
    names(sc)[names(sc) == "sigma"] <- "scale"
    ret <- all.equal(sc[items], sr[items], tolerance=tol)
    if (!isTRUE(ret)) {
        print(sr)
        for (i in seq_along(items)) {
            print(sc[items[i]])
            print(sr[items[i]])
        }
        print(ret)
        stop(sprintf("all.equal(sc[items], sr[items], tol.. = %g) are not all TRUE",
                     tol))
    }
    ret
}

set.seed(101)

test(Y ~ 0, coleman, c("residuals", "df", "coefficients",
                       "r.squared", "adj.r.squared"), tol=1e-10)
test(Y ~ 1,     coleman, tol = 2e-4)
test(Y ~ .,     coleman, tol = 4e-4)
test(Y ~ . - 1, coleman, tol = 4e-4)

