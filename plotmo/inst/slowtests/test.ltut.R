# test.ltut.R: test modified version of linmod from
#              Friedrich Leisch "Creating R Packages: A Tutorial"
#
# Contains three version of linmod (grep for "###" or "source")
#   1. original code from Leisch tutorial
#   2. minimal changes for Milborrow tutorial (and for plotmo)
#   3. production version (includes error checks)

options(warn=2) # treat warnings as errors

if(!interactive())
    postscript(paper="letter")

# test that we got an error as expected from a try() call
expect.err <- function(object, expected.msg="")
{
    if(class(object)[1] == "try-error") {
        msg <- attr(object, "condition")$message[1]
        if(length(grep(expected.msg, msg, fixed=TRUE)))
            cat("Got error as expected from ",
                deparse(substitute(object)), "\n", sep="")
        else
            stop(sprintf("Expected: %s\n  Got:      %s",
                         expected.msg, substr(msg, 1, 1000)))
    } else
        stop("did not get expected error ", expected.msg)
}
almost.equal <- function(x, y)
{
    length(x) == length(y) && max(abs(x - y)) < 1e-10
}
check.lm <- function(fit, ref, check.coef.names=TRUE)
{
    cat("check ", deparse(substitute(fit)), " vs ",
        deparse(substitute(ref)), "\n", sep="")

    stopifnot(coef(fit) == coef(ref))
    if(check.coef.names)
        stopifnot(names(coef(fit)) == names(coef(ref)))

    stopifnot(identical(dim(fit$coefficients), dim(ref$coefficients)))
    stopifnot(length(fit$coefficients) == length(ref$coefficients))
    stopifnot(almost.equal(fit$coefficients, ref$coefficients))

    stopifnot(identical(dim(fit$residuals), dim(ref$residuals)))
    stopifnot(length(fit$residuals) == length(ref$residuals))
    stopifnot(almost.equal(fit$residuals, ref$residuals))

    stopifnot(identical(dim(fit$fitted.values), dim(ref$fitted.values)))
    stopifnot(length(fit$fitted.values) == length(ref$fitted.values))
    stopifnot(almost.equal(fit$fitted.values, ref$fitted.values))

    if(!is.null(fit$vcov) && !is.null(ref$vcov)) {
        stopifnot(identical(dim(fit$vcov), dim(ref$vcov)))
        stopifnot(length(fit$vcov) == length(ref$vcov))
        stopifnot(almost.equal(fit$vcov, ref$vcov))
    }
    ref.sigma <- ref$sigma
    if(is.null(ref.sigma)) # in lm models, sigma is only available from summary()
        ref.sigma <- summary(ref)$sigma
    stopifnot(almost.equal(fit$sigma, ref.sigma))

    stopifnot(almost.equal(fit$df, ref$df))

    stopifnot(almost.equal(fitted(fit), fitted(ref)))
    stopifnot(names(fitted(fit)) == names(fitted(ref)))

    stopifnot(almost.equal(residuals(fit), residuals(ref)))
    stopifnot(names(residuals(fit)) == names(residuals(ref)))

    # TODO this doesn't test predict with newdata
    stopifnot(almost.equal(predict(fit), predict(ref)))
    stopifnot(names(predict(fit)) == names(predict(ref)))
}
### Version1: original code from tutorial

source("linmod.leisch.tutorial.R")

cat("==example issues with predict with functions in the tutorial\n")
data(trees)
tr <- trees
rownames(tr) <- paste("tree", 1:nrow(trees), sep="")
fit1 <- linmod(Volume~., data=tr)
expect.err(try(predict(fit1, newdata=data.frame(Girth=10, Height=80))), "object 'Volume' not found")
expect.err(try(predict(fit1, newdata=as.matrix(tr[1:3,]))), "'data' must be a data.frame, not a matrix or an array")
library(plotmo)
expect.err(try(plotmo(fit1)), "object 'Volume' not found")
fit2 <- linmod(cbind(1, tr[,1:2]), tr[,3])
stopifnot(coef(fit1) == coef(fit2))
# following fail because newdata is a data.frame not a matrix
expect.err(try(predict(fit2, newdata=tr[,1:2])), "requires numeric/complex matrix/vector arguments")
expect.err(try(predict(fit2, newdata=data.frame(Girth=10, Height=80))), "requires numeric/complex matrix/vector arguments")
expect.err(try(predict(fit2, newdata=as.matrix(data.frame(Girth=10, Height=80)))), "non-conformable arguments")
expect.err(try(plotmo(fit2)), "requires numeric/complex matrix/vector arguments")

cat("==a plotmo method function can deal with the issues\n")
plotmo.predict.linmod <- function(object, newdata, ...)
{
    if(is.null(object$formula))                                # x,y interface?
        plotmo:::plotmo.predict.defaultm(object, newdata, ...) # pass matrix not data.frame
    else {
        # add dummy response column to newdata
        newdata[[as.character(as.list(object$formula)[[2]])]] <- 1
        plotmo:::plotmo.predict.default(object, newdata, ...)
    }
}
plotmo(fit1, pt.col=2, caption="fit1 with original tutorial code and plotmo.predict.linmod")
plotmo(fit2, pt.col=2, caption="fit2 with original tutorial code and plotmo.predict.linmod")
remove(plotmo.predict.linmod)

### Version2: minimal changes version for vignette "Guidelines for S3 Regression Models"

source("linmod.milbo.tutorial.R")

cat("==check that example issues with functions in the tutorial have gone\n")
fit1.form <- linmod(Volume~., data=tr)
cat("==print(summary(fit1.form))\n")
print(summary(fit1.form))
stopifnot(abs(predict(fit1.form, newdata=data.frame(Girth=10, Height=80)) - 16.234045) < 1e-5)
stopifnot(sum(abs(predict(fit1.form, newdata=as.matrix(tr[1:3,])) - c(4.8376597, 4.5538516, 4.8169813))) < 1e-5)

lm.tr <- lm(Volume~., data=tr)
check.lm(fit1.form, lm.tr)

fit1.mat <- linmod(tr[,1:2], tr[,3]) # note no need for intercept term
cat("==print(summary(fit1.mat))\n")
print(summary(fit1.mat))
stopifnot(abs(predict(fit1.mat, newdata=data.frame(Girth=10, Height=80)) - 16.234045) < 1e-5)
stopifnot(sum(abs(predict(fit1.mat, newdata=tr[1:3,1:2]) - c(4.8376597, 4.5538516, 4.8169813))) < 1e-5)
stopifnot(abs(predict(fit1.mat, newdata=as.matrix(data.frame(Girth=10, Height=80))) - 16.234045) < 1e-5)

check.lm(fit1.mat, lm.tr)

cat("==example plots\n")

library(plotmo)
data(trees)

fit1.form <- linmod(Volume~., data=trees)
print(fit1.form)
print(summary(fit1.form))

fit1.mat <- linmod(trees[,1:2], trees[,3])
print(fit1.mat)
print(summary(fit1.mat))

plotmo(fit1.form)
plotmo(fit1.mat)

plotres(fit1.form)
plotres(fit1.mat)

cat("==test model building with different numeric args\n")

x <- tr[,1:2]
y <- tr[,3]
fit2.mat <- linmod(x, y)
check.lm(fit2.mat, lm.tr)

# check consistency with lm
expect.err(try(linmod(y~x)), "invalid type (list) for variable 'x'")
expect.err(try(lm(y~x)), "invalid type (list) for variable 'x'")

fit3.mat <- linmod(as.matrix(x), as.matrix(y))
check.lm(fit3.mat, lm.tr)

fit4.form <- linmod(y ~ as.matrix(x))
lm4 <- linmod(y ~ as.matrix(x))
check.lm(fit4.form, lm4)
stopifnot(coef(fit4.form)  == coef(lm.tr),
          gsub("as.matrix(x)", "", names(coef(fit4.form)), fixed=TRUE)  == names(coef(lm.tr)))

xm <- as.matrix(x)
fit5.form <- linmod(y ~ xm)
lm5 <- linmod(y ~ xm)
check.lm(fit5.form, lm5)
stopifnot(coef(fit5.form)  == coef(lm.tr),
          gsub("xm", "", names(coef(fit5.form)), fixed=TRUE)  == names(coef(lm.tr)))

cat("==test correct use of global x1 and y1\n")
x1 <- tr[,1]
y1 <- tr[,3]
lm1 <- linmod(y1~x1)

fit6.mat <- linmod(x1, y1)
check.lm(fit6.mat, lm1, check.coef.names=FALSE)
# production version only:
# stopifnot(coef(fit6.mat) == coef(lm1),
#           names(coef(fit6.mat)) == c("(Intercept)", "V1")) # names(coef(lm1) are "(Intercept)" "x1"

fit6.form <- linmod(y1~x1)
check.lm(fit6.form, lm1)

cat("==check integer input (sibsp is an integer) \n")

library(earth) # for etitanic data
data(etitanic)
tit <- etitanic[seq(1, nrow(etitanic), by=60), ] # small set of data for tests (18 cases)
tit$survived <- tit$survived != 0 # convert to logical
rownames(tit) <- paste("pas", 1:nrow(tit), sep="")
cat(paste(colnames(tit), "=", sapply(tit, class), sep="", collapse=", "), "\n")

fit7.mat <- linmod(tit$age, tit$sibsp)
lm7 <- lm.fit(cbind(1, tit$age), tit$sibsp)
stopifnot(coef(fit7.mat) == coef(lm7)) # coef names will differ

fit7.form <- linmod(sibsp~age, data=tit)
lm7.form  <- lm(sibsp~age, data=tit)
check.lm(fit7.form, lm7.form)

fit8.mat <- linmod(tit$sibsp, tit$age)
lm8 <- lm.fit(cbind(1, tit$sibsp), tit$age)
stopifnot(coef(fit8.mat) == coef(lm8)) # coef names will differ

fit8.form <- linmod(age~sibsp, data=tit)
lm8.form  <- lm(age~sibsp, data=tit)
check.lm(fit8.form, lm8.form)

# drop=FALSE so response is a data frame
fit1a.mat <- linmod(trees[,1:2], trees[, 3, drop=FALSE])
print(fit1a.mat)
print(summary(fit1.mat))
plotres(fit1a.mat) # plot caption shows response name "Volume"

cat("==test model building with different non numeric args\n")

library(earth) # for etitanic data
data(etitanic)
tit <- etitanic[seq(1, nrow(etitanic), by=60), ] # small set of data for tests (18 cases)
tit$survived <- tit$survived != 0 # convert to logical
rownames(tit) <- paste("pas", 1:nrow(tit), sep="")
cat(paste(colnames(tit), "=", sapply(tit, class), sep="", collapse=", "), "\n")

lm9 <- lm(survived~., data=tit)
fit9.form <- linmod(survived~., data=tit)
check.lm(fit9.form, lm9)

options(warn=2) # treat warnings as errors
expect.err(try(linmod(tit[,c(1,3,4,5,6)], tit[,"survived"])), "NAs introduced by coercion")
options(warn=1)
expect.err(try(linmod(tit[,c(1,3,4,5,6)], tit[,"survived"])), "NA/NaN/Inf in foreign function call (arg 1)")

options(warn=2) # treat warnings as errors
expect.err(try(lm(pclass~., data=tit)), "using type = \"numeric\" with a factor response will be ignored")
# minimal version
expect.err(try(linmod(pclass~., data=tit)), "(converted from warning) NAs introduced by coercion")
expect.err(try(linmod(tit$pclass, tit$survived)), "(converted from warning) NAs introduced by coercion")
# # production version
# expect.err(try(linmod(pclass~., data=tit)), "'y' is not numeric or logical")
options(warn=1)

lm10 <- lm(pclass~., data=tit) # will give warnings
fit10.form <- linmod(as.numeric(pclass)~., data=tit)
stopifnot(coef(fit10.form) == coef(lm10))
stopifnot(names(coef(fit10.form)) == names(coef(lm10)))
# check.lm(fit10.form, lm10) # fails because lm10 fitted is all NA

# production version: (minimal version just gives warnings and builds lousy model)
# expect.err(try(linmod(pclass~., data=tit)), "'y' is not numeric or logical")
# expect.err(try(linmod(tit[,-1], tit[,1])), "'y' is not numeric or logical")
# expect.err(try(linmod(1:10, paste(1:10))), "'y' is not numeric or logical")

fit10a.form <- linmod(survived~pclass, data=tit)
lm10a <- lm(survived~pclass, data=tit)
check.lm(fit10a.form, lm10a)

expect.err(try(linmod(paste(1:10), 1:10)), "requires numeric/complex matrix/vector arguments")

lm11 <- lm(as.numeric(pclass)~., data=tit)
fit11.form <- linmod(as.numeric(pclass)~., data=tit)
check.lm(fit11.form, lm11)

cat("==data.frame with strings\n")

df.with.string <-
    data.frame(1:5,
               c(1,2,-1,4,5),
               c("a", "b", "a", "a", "b"),
               stringsAsFactors=FALSE)
colnames(df.with.string) <- c("num1", "num2", "string")

fit30.form <- linmod(num1~num2, df.with.string)
lm30       <- lm(num1~num2, df.with.string)
check.lm(fit30.form, lm30)

fit31.form <- linmod(num1~., df.with.string)
lm31       <- lm(num1~., df.with.string)
check.lm(fit31.form, lm31)

expect.err(try(linmod(string~., df.with.string)), "non-numeric argument to binary operator")
# production version
# expect.err(try(linmod(string~., df.with.string)), "'y' is not numeric or logical")

vec <- c(1,2,3,4,3)
options(warn=2) # treat warnings as errors
expect.err(try(linmod(df.with.string, vec)), "NAs introduced by coercion")
options(warn=1)
# minimal version
expect.err(try(linmod(df.with.string, vec)), "NA/NaN/Inf in foreign function call (arg 1)")
# production version
# expect.err(try(linmod(df.with.string, vec)), "NA in 'x'")

options(warn=2) # treat warnings as errors
expect.err(try(linmod(df.with.string, vec)), "NAs introduced by coercion")
options(warn=1)
# minimal version
expect.err(try(linmod(df.with.string, vec)), "NA/NaN/Inf in foreign function call (arg 1)")
# production version
# expect.err(try(linmod(df.with.string, vec)), "NA in 'x'")

cat("==more variables  than cases\n")

set.seed(1)
x2 <- matrix(rnorm(6), nrow=2)
y2 <- c(1,2)
# production version
# expect.err(try(linmod(y2~x2)), "more variables than cases")
# minimal version
expect.err(try(linmod(y2~x2)), "'size' cannot exceed nrow(x) = 2")

x3 <- matrix(1:10, ncol=2)
y3 <- c(1,2,9,4,5)
# production version will give a better error message
expect.err(try(linmod(y3~x3)), "singular matrix 'a' in 'solve'")

cat("==nrow(x) does not match length(y)\n")
# note that the production version gives better error messages

x4 <- matrix(1:10, ncol=2)
y4 <- c(1,2,9,4)
expect.err(try(linmod(x4, y4)), "singular matrix 'a' in 'solve'")

x5 <- matrix(1:10, ncol=2)
y5 <- c(1,2,9,4,5,9)
expect.err(try(linmod(x5, y5)), "singular matrix 'a' in 'solve'")

cat("==y has multiple columns\n")

vec <- c(1,2,3,4,3)
y2 <- cbind(c(1,2,3,4,9), vec^2)
expect.err(try(linmod(vec, y2)), "'qr' and 'y' must have the same number of rows")
# following does not issue any error message, it should
# expect.err(try(linmod(y2~vec)), "error message")

### Version 3: production version of linmod

source("linmod.R")

cat("==check that example issues with functions in the tutorial have gone\n")
fit1.form <- linmod(Volume~., data=tr)
cat("==print(summary(fit1.form))\n")
print(summary(fit1.form))
lm.tr <- lm(Volume~., data=tr)
check.lm(fit1.form, lm.tr)
stopifnot(abs(predict(fit1.form, newdata=data.frame(Girth=10, Height=80)) - 16.234045) < 1e-5)
stopifnot(sum(abs(predict(fit1.form, newdata=as.matrix(tr[1:3,])) - c(4.8376597, 4.5538516, 4.8169813))) < 1e-5)

fit1.mat <- linmod(tr[,1:2], tr[,3]) # note no need for intercept term
cat("==print(summary(fit1.mat))\n")
print(summary(fit1.mat))
check.lm(fit1.mat, lm.tr)
stopifnot(abs(predict(fit1.mat, newdata=data.frame(Girth=10, Height=80)) - 16.234045) < 1e-5)
stopifnot(sum(abs(predict(fit1.mat, newdata=tr[1:3,1:2]) - c(4.8376597, 4.5538516, 4.8169813))) < 1e-5)
stopifnot(abs(predict(fit1.mat, newdata=as.matrix(data.frame(Girth=10, Height=80))) - 16.234045) < 1e-5)

cat("==print.default(fit1.form)\n")
print.default(fit1.form)

cat("==check single x variable\n")
fit1a.form <- linmod(Volume~Height, data=tr)
cat("==print(summary(fit1a.form))\n")
print(summary(fit1a.form))
lma.tr <- lm(Volume~Height, data=tr)
check.lm(fit1a.form, lma.tr)

stopifnot(abs(predict(fit1a.form, newdata=data.frame(Height=80)) - 36.34437) < 1e-5)
stopifnot(abs(predict(fit1a.form, newdata=data.frame(Girth=99, Height=80)) - 36.34437) < 1e-5)
stopifnot(sum(abs(predict(fit1a.form, newdata=as.matrix(tr[1:3,])) - c(20.91087, 13.19412, 10.10742))) < 1e-5)

fit1a.mat <- linmod(tr[,2,drop=FALSE], tr[,3])
cat("==print(summary(fit1a.mat))\n")
print(summary(fit1a.mat))
check.lm(fit1a.mat, lma.tr)
stopifnot(abs(predict(fit1a.mat, newdata=data.frame(Height=80)) - 36.34437) < 1e-5)
stopifnot(sum(abs(predict(fit1a.mat, newdata=tr[1:3,2]) - c(20.91087, 13.19412, 10.10742))) < 1e-5)
stopifnot(abs(predict(fit1a.mat, newdata=as.matrix(data.frame(Height=80))) - 36.34437) < 1e-5)

# check that rownames got propagated
stopifnot(names(fit1.form$residuals)[1] == "tree1")
stopifnot(names(fit1.form$fitted.values)[3] == "tree3")
stopifnot(names(fit1.mat$residuals)[1] == "tree1")
stopifnot(names(fit1.mat$fitted.values)[3] == "tree3")
stopifnot(!is.null(names(fit1.mat$residuals)))
stopifnot(!is.null(names(fit1.mat$fitted.values)))
cat("==print.default(fit1.mat)\n")
print.default(fit1.mat)

# check that we don't artificially add rownames when no original rownames
fit1a.mat <- linmod(trees[,1:2], trees[,3])
stopifnot(is.null(names(fit1a.mat$residuals)))
stopifnot(is.null(names(fit1a.mat$fitted.values)))

cat("==example plots\n")

library(plotmo)
data(trees)

fit1.form <- linmod(Volume~., data=trees)
print(fit1.form)
print(summary(fit1.form))

fit1.mat <- linmod(trees[,1:2], trees[,3])
print(fit1.mat)
print(summary(fit1.mat))

plotmo(fit1.form)
plotmo(fit1.mat)

plotres(fit1.form)
plotres(fit1.mat)

# test that keep arg works correctly for plotmo and plotres
fit1.mat.keep <- linmod(trees[,1:2], trees[,3], keep=TRUE)
fit1.mat.keep$call <- NULL # trick to force use of x and y in plotmo
plotmo(fit1.mat.keep, pt.col=3)
plotres(fit1.mat.keep)

cat("==test model building with different numeric args\n")

x <- tr[,1:2]
y <- tr[,3]
fit2.mat <- linmod(x, y)
check.lm(fit2.mat, lm.tr)

# check consistency with lm
expect.err(try(linmod(y~x)), "invalid type (list) for variable 'x'")
expect.err(try(lm(y~x)), "invalid type (list) for variable 'x'")

fit3.mat <- linmod(as.matrix(x), as.matrix(y))
check.lm(fit3.mat, lm.tr)

fit4.form <- linmod(y ~ as.matrix(x))
lm4 <- linmod(y ~ as.matrix(x))
check.lm(fit4.form, lm4)
stopifnot(coef(fit4.form)  == coef(lm.tr),
          gsub("as.matrix(x)", "", names(coef(fit4.form)), fixed=TRUE)  == names(coef(lm.tr)))

xm <- as.matrix(x)
fit5.form <- linmod(y ~ xm)
lm5 <- linmod(y ~ xm)
check.lm(fit5.form, lm5)
stopifnot(coef(fit5.form)  == coef(lm.tr),
          gsub("xm", "", names(coef(fit5.form)), fixed=TRUE)  == names(coef(lm.tr)))

# test that non-intercept models give an error message (for now,
# until we test them, and remove -1 in predict.linmod)

cat("==test no intercept model gives an error message\n")
expect.err(try(linmod(Volume~.-1, data=tr)), "the first column of 'x' is not an intercept column (all 1s)")

cat("==test correct use of global x1 and y1\n")
x1 <- tr[,1]
y1 <- tr[,3]
lm1 <- linmod(y1~x1)

fit6.mat <- linmod(x1, y1)
check.lm(fit6.mat, lm1, check.coef.names=FALSE)
# production version only:
stopifnot(coef(fit6.mat) == coef(lm1),
          names(coef(fit6.mat)) == c("(Intercept)", "V1")) # names(coef(lm1) are "(Intercept)" "x1"

fit6.form <- linmod(y1~x1)
check.lm(fit6.form, lm1)

cat("==check integer input (sibsp is an integer) \n")

library(earth) # for etitanic data
data(etitanic)
tit <- etitanic[seq(1, nrow(etitanic), by=60), ] # small set of data for tests (18 cases)
tit$survived <- tit$survived != 0 # convert to logical
rownames(tit) <- paste("pas", 1:nrow(tit), sep="")
cat(paste(colnames(tit), "=", sapply(tit, class), sep="", collapse=", "), "\n")

fit7.mat <- linmod(tit$age, tit$sibsp)
lm7 <- lm.fit(cbind(1, tit$age), tit$sibsp)
stopifnot(coef(fit7.mat) == coef(lm7)) # coef names will differ

fit7.form <- linmod(sibsp~age, data=tit)
lm7.form  <- lm(sibsp~age, data=tit)
check.lm(fit7.form, lm7.form)

fit8.mat <- linmod(tit$sibsp, tit$age)
lm8 <- lm.fit(cbind(1, tit$sibsp), tit$age)
stopifnot(coef(fit8.mat) == coef(lm8)) # coef names will differ

fit8.form <- linmod(age~sibsp, data=tit)
lm8.form  <- lm(age~sibsp, data=tit)
check.lm(fit8.form, lm8.form)

# drop=FALSE so response is a data frame
fit1a.mat <- linmod(trees[,1:2], trees[, 3, drop=FALSE])
print(fit1a.mat)
print(summary(fit1.mat))
plotres(fit1a.mat) # plot caption shows response name "Volume"

cat("==test model building with different non numeric args\n")

library(earth) # for etitanic data
data(etitanic)
tit <- etitanic[seq(1, nrow(etitanic), by=60), ] # small set of data for tests (18 cases)
tit$survived <- tit$survived != 0 # convert to logical
rownames(tit) <- paste("pas", 1:nrow(tit), sep="")
cat(paste(colnames(tit), "=", sapply(tit, class), sep="", collapse=", "), "\n")

lm9 <- lm(survived~., data=tit)
fit9.form <- linmod(survived~., data=tit)
check.lm(fit9.form, lm9)

expect.err(try(linmod(tit[,c(1,3,4,5,6)], tit[,"survived"])),
           "non-numeric column in 'x'")
fit9a.form <- lm.fit(data.matrix(cbind("(Intercept)"=1, tit[,c(1,3,4,5,6)])), tit[,"survived"])
lm9 <- lm.fit(data.matrix(cbind("(Intercept)"=1, tit[,c(1,3,4,5,6)])), tit[,"survived"])
stopifnot(coef(fit9a.form) == coef(lm9))
stopifnot(names(coef(fit9a.form)) == names(coef(lm9)))

options(warn=2) # treat warnings as errors
expect.err(try(lm(pclass~., data=tit)), "using type = \"numeric\" with a factor response will be ignored")
expect.err(try(linmod(pclass~., data=tit)), "'y' is not numeric or logical")
options(warn=1)

lm10 <- lm(pclass~., data=tit) # will give warnings
fit10.form <- linmod(as.numeric(pclass)~., data=tit)
stopifnot(coef(fit10.form) == coef(lm10))
stopifnot(names(coef(fit10.form)) == names(coef(lm10)))
# check.lm(fit10.form, lm10) # fails because lm10 fitted is all NA

expect.err(try(linmod(pclass~., data=tit)), "'y' is not numeric or logical")
expect.err(try(linmod(tit[,-1], tit[,1])), "non-numeric column in 'x'")
expect.err(try(linmod(1:10, paste(1:10))), "'y' is not numeric or logical")

fit10a.form <- linmod(survived~pclass, data=tit)
lm10a <- lm(survived~pclass, data=tit)
check.lm(fit10a.form, lm10a)

expect.err(try(linmod(tit[,"pclass"], tit[,"age"])), "non-numeric column in 'x'")

expect.err(try(linmod(paste(1:10), 1:10)), "non-numeric column in 'x'")

lm11 <- lm(as.numeric(pclass)~., data=tit)
fit11.form <- linmod(as.numeric(pclass)~., data=tit)
check.lm(fit11.form, lm11)

# logical data (not numeric)
bool.data <- data.frame(x=rep(c(TRUE, FALSE, TRUE), length.out=10),
                        y=rep(c(TRUE, FALSE, FALSE), length.out=10))
lm12 <- lm(y~x, data=bool.data)
fit12.form <- linmod(y~x, data=bool.data)
check.lm(fit12.form, lm12)
fit12.xy <- linmod(bool.data$x, bool.data$y)
# hack: delete mismatching names so check.lm() doesn't fail
names(lm12$coefficients) <- NULL     # were "(Intercept)" "xTRUE"
names(fit12.xy$coefficients) <- NULL # were "(Intercept)" "V1"
check.lm(fit12.xy, lm12)

cat("==data.frame with strings\n")

df.with.string <-
    data.frame(1:5,
               c(1,2,-1,4,5),
               c("a", "b", "a", "a", "b"),
               stringsAsFactors=FALSE)
colnames(df.with.string) <- c("num1", "num2", "string")

fit30.form <- linmod(num1~num2, df.with.string)
lm30       <- lm(num1~num2, df.with.string)
check.lm(fit30.form, lm30)

fit31.form <- linmod(num1~., df.with.string)
lm31       <- lm(num1~., df.with.string)
check.lm(fit31.form, lm31)

expect.err(try(linmod(string~., df.with.string)), "'y' is not numeric or logical")

vec <- c(1,2,3,4,3)
expect.err(try(linmod(df.with.string, vec)), "non-numeric column in 'x'")
expect.err(try(linmod(tit$pclass, tit$survived)), "non-numeric column in 'x'")

cat("==x is singular\n")

set.seed(1)
x2 <- matrix(rnorm(6), nrow=2)
y2 <- c(1,2)
expect.err(try(linmod(y2~x2)), "'x' is singular (it has 4 columns but its rank is 2)")

x3 <- matrix(1:10, ncol=2)
y3 <- c(1,2,9,4,5)
expect.err(try(linmod(y3~x3)), "'x' is singular (it has 3 columns but its rank is 2)")

expect.err(try(linmod(trees[1,1:2], trees[1,3])), "'x' is singular (it has 3 columns but its rank is 1)")

cat("==nrow(x) does not match length(y)\n")

x4 <- matrix(1:10, ncol=2)
y4 <- c(1,2,9,4)
expect.err(try(linmod(x4, y4)), "nrow(x) is 5 but length(y) is 4")

x5 <- matrix(1:10, ncol=2)
y5 <- c(1,2,9,4,5,9)
expect.err(try(linmod(x5, y5)), "nrow(x) is 5 but length(y) is 6")

cat("==y has multiple columns\n")

vec <- c(1,2,3,4,3)
y2 <- cbind(c(1,2,3,4,9), vec^2)
expect.err(try(linmod(vec, y2)), "nrow(x) is 5 but length(y) is 10")
expect.err(try(linmod(y2~vec)), "nrow(x) is 5 but length(y) is 10")

cat("==NA in x\n")

x <- tr[,1:2]
y <- tr[,3]
x[2,2] <- NA
expect.err(try(linmod(x, y)), "NA in 'x'")

x <- tr[,1:2]
y <- tr[,3]
y[9] <- NA
expect.err(try(linmod(x, y)), "NA in 'y'")

cat("==misc tests with different kinds of data\n")

data3 <- data.frame(s=c("a", "b", "c", "a", "a"), num=c(1,9,4,2,6), y=c(1,2,3,5,3), stringsAsFactors=F)
stopifnot(sapply(data3, class) == c("character", "numeric", "numeric"))
a40 <- linmod(y~., data=data3)
stopifnot(sum(abs(a40$coef - c(2.571, -1.857, -0.143, 0.143))) < 0.001)
stopifnot(sum(abs(predict(a40, newdata=data3[1:3,]) - c(2.71429, 2.00000, 3.00000))) < 0.001)

data4 <- data.frame(s=c("a", "b", "c", "a", "a"), num=c(1,9,4,2,6), y=c(1,2,3,5,3), stringsAsFactors=T)
stopifnot(sapply(data4, class) == c("factor", "numeric", "numeric"))
expect.err(try(linmod(data4[,1:2], data4[,3])), "non-numeric column in 'x'")

# following gives no error (and matches lm)
a41 <- linmod(y~., data=data4)
stopifnot(sum(abs(predict(a41, newdata=data3[1:3,]) - c(2.71429, 2.00000, 3.00000))) < 0.001)

data5 <- data.frame(s=c("a", "b", "c", "a", "a"), num=c(1,9,4,2,6), y=c(1,2,3,5,3), stringsAsFactors=F)
stopifnot(sum(abs(predict(a41, newdata=data5[1:3,1:2]) - c(2.71429, 2.00000, 3.00000))) < 0.001)

data6 <- data.frame(s=c("a", "b", "c", "a9", "a"), num=c(1,9,4,2,6), y=c(1,2,3,5,3), stringsAsFactors=T)
expect.err(try(predict(a41, newdata=data6[1:3,1:2])), "ncol(newdata) is 4 but should be 3")

expect.err(try(predict(a41, newdata=data.frame(s=1, num=2, y=3))), "ncol(newdata) is 2 but should be 3")

expect.err(try(predict(a41, newdata=1:9)), "object 's' not found") # issued by model.matrix in predict.linmod

expect.err(try(predict(a41, newdata=data.frame())), "'newdata' is empty")

tr.na <- trees
tr.na[9,3] <- NA
expect.err(try(linmod(Volume~.,data=tr.na)), "NA in 'y'")
expect.err(try(linmod(tr.na[,1:2], tr.na[,3])), "NA in 'y'")

tr.na <- trees
tr.na[10,1] <- NA
expect.err(try(linmod(Volume~.,data=tr.na)), "NA in 'x'")
expect.err(try(linmod(tr.na[,1:2], tr.na[,3])), "NA in 'x'")

a42 <- linmod(trees[,1:2], trees[, 3])
newdata1 <- data.frame(Girth=20)
expect.err(try(predict(a42, newdata=newdata1)), "ncol(newdata) is 1 but should be 2")
expect.err(try(predict(a42, newdata=data.frame())), "'newdata' is empty")
newdata.with.NA <- data.frame(Girth=20, Height=NA)
expect.err(try(predict(a42, newdata=newdata.with.NA)), "NA in 'newdata'")

a43 <- linmod(Volume~.,data=trees)
expect.err(try(predict(a43, newdata=newdata.with.NA)), "NA in 'newdata'")

y6 <- 1:5
x6 <- data.frame()
expect.err(try(linmod(x6, y6)), "'x' is empty")

y7 <- data.frame()
x7 <- 1:5
expect.err(try(linmod(x7, y7)), "'y' is empty")

# duplicated column names
data7 <- matrix(1:25, ncol=5)
colnames(data7) <- c("y", "x1", "x1", "x3", "x4")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"x1\" in 'x' is duplicated")

colnames(data7) <- c("y", "x1", "x2", "x2", "x4")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"x2\" in 'x' is duplicated")

colnames(data7) <- c("y", "x1", "x2", "x2", "x2")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"x2\" in 'x' is duplicated")

# column name V2 will be created but it clashes with the existing column name
colnames(data7) <- c("y", "V2", "", "V3", "V4")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"V2\" in 'x' is duplicated")

# missing column names
trees1 <- trees
colnames(trees1) <- NULL
cat("a52\n")
a52 <- linmod(trees1[,1:2], trees1[,3])
print(summary(a52))

trees1 <- trees
colnames(trees1) <- c("", "Height", "Volume") # was Girth Height Volume
cat("a53\n")
a53 <- linmod(trees1[,1:2], trees1[,3])
print(summary(a53))
cat("a53.formula\n")
expect.err(try(linmod(Volume~., data=trees1)), "attempt to use zero-length variable name")

# very long names to test formatting in summary.linmod
trees1 <- trees
colnames(trees1) <- c("Girth.a.very.long.name.in.fact.an.exceptionally.long.name",
                      "Height.a.very.long.name.in.fact.an.exceptionally.long.name",
                      "Volume.a.very.long.name.in.fact.an.exceptionally.long.name")
cat("a55\n")
a55 <- linmod(Volume.a.very.long.name.in.fact.an.exceptionally.long.name~
              Girth.a.very.long.name.in.fact.an.exceptionally.long.name+
              Height.a.very.long.name.in.fact.an.exceptionally.long.name,
              data=trees1)
print(summary(a55))
options(show.signif.stars=FALSE) # god i hate significance stars
print(summary(a55))

# intercept-only model
a56.form <- linmod(Volume~1, data=trees)
print(summary(a56.form))
stopifnot(length(coef(a56.form)) == 1)
plotres(a56.form)
expect.err(try(plotmo(a56.form)), "x is empty")
expect.err(try(linmod(rep(1, length.out=nrow(trees)), trees$Volume)), "'x' is singular (it has 2 columns but its rank is 1)")

# various tests for bad args
expect.err(try(linmod(trees[,1:2])), "no 'y' argument")

expect.err(try(linmod(Volume~., data=trees, nonesuch=99)), "unused argument (nonesuch = 99)")
expect.err(try(linmod(trees[,1:2], trees[,3], nonesuch=linmod)), "unused argument (nonesuch = function (...)")
expect.err(try(summary(linmod(trees[,1:2], trees[,3]), nonesuch=linmod)), "unused argument (nonesuch = function (...)")
expect.err(try(print(linmod(trees[,1:2], trees[,3]), nonesuch=linmod)), "unused argument (nonesuch = function (...)")

if(!interactive()) {
    dev.off()        # finish postscript plot
    q(runLast=FALSE) # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
