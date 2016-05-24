# test.plotmo.x.R: test plotmo_x and related functions

options(warn=1) # print warnings as they occur

library(plotmo)
library(earth)
data(ozone1)
data(etitanic)

get.tit <- function()
{
    tit <- etitanic
    pclass <- as.character(tit$pclass)
    # change the order of the factors so not alphabetical
    pclass[pclass == "1st"] <- "first"
    pclass[pclass == "2nd"] <- "class2"
    pclass[pclass == "3rd"] <- "classthird"
    tit$pclass <- factor(pclass, levels=c("class2", "classthird", "first"))
    # log age is so we have a continuous predictor even when model is age~.
    set.seed(2015)
    tit$logage <- log(tit$age) + rnorm(nrow(tit))
    tit$parch <- NULL
    # by=12 gives us a small fast model with an additive and a interaction term
    tit <- tit[seq(1, nrow(etitanic), by=12), ]
}

if(!interactive())
    postscript(paper="letter")

printf <- function(format, ...) cat(sprintf(format, ...), sep="") # like c printf

strip.space <- function(s) gsub("[ \t\n]", "", s)

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
X <- X1 <- X2 <- Y <- DF <- NULL
get.data <- function()
{
    X <<- matrix(c(1,2,3,4,5,6,7,8,9,
                   2,3,3,5,6,7,8,9,9), ncol=2)
    colnames(X) <- c("xx1", "xx2")
    X1 <<- X[,1]
    X2 <<- X[,2]
    Y  <<- c(1,2,7,4,5,6,6,6,6)
    DF <<- data.frame(Y=Y, X1=X1, X2=X2)
}
stopifnot1 <- function(x, y){
    xname <- deparse(substitute(x))
    yname <- deparse(substitute(y))
    if(!all(x == y))
        stop(sprintf("%s == %s failed\n", xname, yname, call.=FALSE))
    printf("%s == %s passed\n", xname, yname)
}
printf("====== standard earth.formula model with a data frame\n")

get.data()
earth.form.df.dot <- earth(Y~., data=DF)
plotmo(earth.form.df.dot, caption="test basic use of DF")
printf("-- test basic use of DF\n")
rv <- plotmo(earth.form.df.dot, trace=100)
stopifnot1(rv, X)

printf("-- test use same DF even when other variables change\n")
get.data()
earth.form.df.dot <- earth(Y~., data=DF)
X1 <- "rubbish"
rv <- plotmo(earth.form.df.dot, trace=100)
stopifnot1(rv, X)

printf("-- test detect that DF is now trashed\n")
get.data()
earth.form.df.dot <- earth(Y~., data=DF)
DF <- "rubbish"
X1 <- "rubbish" # DF is corrupt and will treated as NULL by plotmo, so make sure plotmo doesn't find the global X1
# invalid 'envir' argument of type 'character'
expect.err(try(plotmo(earth.form.df.dot, trace=100)), "cannot get the original model predictors")

# Removed this test because this no longer fails, because we get the formula using formula(object)
# printf("-- DF is NULL so will get '.' in formula and no 'data' argument\n")
# get.data()
# earth.form.df.dot <- earth(Y~., data=DF)
# DF <- NULL
# # '.' in formula and no 'data' argument
# expect.err(try(plotmo(earth.form.df.dot, trace=100)), "cannot get the original model predictors")

printf("-- DF is NULL so will pick up X1 with same values from global environment\n")
get.data()
earth.form.df <- earth(Y~X1+X2, data=DF)
DF <- NULL
rv <- plotmo(earth.form.df, trace=100)
stopifnot1(rv, X)

printf("-- DF is NULL so will will pick up trashed X1 from global environment\n")
earth.form.df <- earth(Y~X1+X2, data=DF)
DF <- NULL
X1 <- "rubbish"
# variable lengths differ (found for 'X1')
expect.err(try(plotmo(earth.form.df, trace=100)), "cannot get the original model predictors")

printf("-- DF has only one column, so will pick up X1 from it and X2 from global environment\n")
get.data()
earth.form.df <- earth(Y~X1+X2, data=DF)
DF <- data.frame(Y=Y, X1=X1)
DF[1,2] <- 99
X2[1] <- 98
rv <- plotmo(earth.form.df, trace=100)
stopifnot1(rv[1,1], 99)
stopifnot1(rv[1,2], 98)

printf("-- sanity check, make sure we are back to normal\n")
get.data()
earth.form.df <- earth(Y~X1+X2, data=DF)
rv <- plotmo(earth.form.df, trace=100)
stopifnot1(rv, X)

printf("-- change the data frame, make sure we pick up the changed value\n")
get.data()
earth.form.df <- earth(Y~X1+X2, data=DF)
DF[1,2] <- 99
rv <- plotmo(earth.form.df, trace=100)
stopifnot1(rv[1,1], 99)

printf("-- change order of columns in the data frame, should be ok\n")
get.data()
earth.form.df <- earth(Y~X1+X2, data=DF)
DF <- data.frame(X2=X2, X1=X1)
rv <- plotmo(earth.form.df, trace=100)
stopifnot1(rv, X)

printf("======= standard earth.formula model with a data frame and keepxy\n")

get.data()
earth.form.df.keepxy <- earth(Y~., data=DF, keepxy=TRUE)
printf("-- test basic use of DF\n")
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("-- test use same DF even when other variables change\n")
earth.form.df.keepxy <- earth(Y~., data=DF, keepxy=TRUE)
X1 <- "rubbish"
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("-- DF is now trashed but it doesn't matter because keepxy=T\n")
DF <- "rubbish"
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("-- DF is NULL but it doesn't matter because keepxy=T\n")
get.data()
earth.form.df.keepxy <- earth(Y~., data=DF, keepxy=TRUE)
DF <- NULL
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("-- DF and X1 are NULL but it doesn't matter because keepxy=T\n")
DF <- NULL
X1 <- "rubbish"
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("-- sanity check, make sure we are back to normal\n")
get.data()
earth.form.df.keepxy <- earth(Y~., data=DF, keepxy=TRUE)
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("-- change the data frame, but it doesn't matter because keepxy=T\n")
get.data()
earth.form.df.keepxy <- earth(Y~., data=DF, keepxy=TRUE)
DF[1,2] <- 99
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("-- change order of columns in the data frame, should be ok\n")
get.data()
earth.form.df.keepxy <- earth(Y~., data=DF, keepxy=TRUE)
DF <- data.frame(X2=X2, X1=X1)
rv <- plotmo(earth.form.df.keepxy, trace=100)
stopifnot1(rv, X)

printf("======= standard lm model with a data frame but with model=FALSE\n")

get.data()
lm.form.df.model.false.with.dot <- lm(Y~., data=DF, model=FALSE)
printf("-- test basic use of DF\n")
rv <- plotmo(lm.form.df.model.false.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- test use same DF even when other variables change\n")
get.data()
lm.form.df.model.false.with.dot <- lm(Y~., data=DF, model=FALSE)
X1 <- "rubbish"
rv <- plotmo(lm.form.df.model.false.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- test detect that DF is now trashed\n")
DF <- "rubbish"
# invalid 'envir' argument of type 'character'
expect.err(try(plotmo(lm.form.df.model.false.with.dot, trace=100)), "cannot get the original model predictors")

printf("-- DF is NULL so will pick up X1 with same values from global environment\n")
get.data()
lm.form.df.model.false <- lm(Y~X1+X2, data=DF, model=FALSE)
DF <- NULL
rv <- plotmo(earth.form.df, trace=100)
stopifnot1(rv, X)

printf("-- DF is NULL so will will pick up trashed X1 from global environment\n")
get.data()
lm.form.df.model.false <- lm(Y~X1+X2, data=DF, model=FALSE)
DF <- NULL
X1 <- "rubbish"
# variable lengths differ (found for 'X1')
expect.err(try(plotmo(lm.form.df.model.false, trace=100)), "cannot get the original model predictors")

printf("-- sanity check, make sure we are back to normal\n")
get.data()
lm.form.df.model.false <- lm(Y~X1+X2, data=DF, model=FALSE)
rv <- plotmo(lm.form.df.model.false, trace=100)
stopifnot1(rv, X)

printf("-- change the data frame, make sure we pick up the changed value\n")
get.data()
lm.form.df.model.false <- lm(Y~X1+X2, data=DF, model=FALSE)
DF[1,2] <- 99
rv <- plotmo(lm.form.df.model.false, trace=100)
stopifnot1(rv[1,1], 99)

printf("-- change order of columns in the data frame, should be ok\n")
get.data()
lm.form.df.model.false <- lm(Y~X1+X2, data=DF, model=FALSE)
DF <- data.frame(X2=X2, X1=X1)
rv <- plotmo(lm.form.df.model.false, trace=100)
stopifnot1(rv, X)

printf("======= standard lm with a data frame and model=TRUE (the default)\n")

get.data()
lm.form.df.with.dot <- lm(Y~., data=DF)
printf("-- test basic use of DF\n")
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- test use same DF even when other variables change\n")
get.data()
lm.form.df.with.dot <- lm(Y~., data=DF)
X1 <- "rubbish"
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- DF is now trashed but it doesn't matter because keepxy=T\n")
DF <- "rubbish"
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- DF is NULL but it doesn't matter because keepxy=T\n")
get.data()
lm.form.df.with.dot <- lm(Y~., data=DF)
DF <- NULL
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- DF and X1 are NULL but it doesn't matter because keepxy=T\n")
DF <- NULL
X1 <- "rubbish"
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- sanity check, make sure we are back to normal\n")
get.data()
lm.form.df.with.dot <- lm(Y~., data=DF)
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- change the data frame, but it doesn't matter because keepxy=T\n")
get.data()
lm.form.df.with.dot <- lm(Y~., data=DF)
DF[1,2] <- 99
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("-- change order of columns in the data frame, should be ok\n")
get.data()
lm.form.df.with.dot <- lm(Y~., data=DF)
DF <- data.frame(X2=X2, X1=X1)
rv <- plotmo(lm.form.df.with.dot, trace=100)
stopifnot1(rv, X)

printf("======= standard lm with a data frame and model=FALSE but x=TRUE\n")

get.data()
lm.form.df.model.false.x.true <- lm(Y~., data=DF, model=FALSE, x=TRUE)
printf("-- test basic use of DF\n")
rv <- plotmo(lm.form.df.model.false.x.true, trace=100)
stopifnot1(rv, X)

printf("-- test DF not available (shouldn't matter)\n")
DF <- "rubbish"
rv <- plotmo(lm.form.df.model.false.x.true, trace=100)
stopifnot1(rv, X)

printf("-- test $x trashed causes failure\n")
get.data()
lm.form.df.model.false.x.true <- lm(Y~., data=DF, model=FALSE, x=TRUE)
DF <- "rubbish"
X2 <- "rubbish1"
lm.form.df.model.false.x.true[["x"]] <- "nonesuch"
expect.err(try(plotmo(lm.form.df.model.false.x.true, trace=100)), "cannot get the original model predictors")

printf("-- test ok with $x trashed but DF ok\n") # although with trace!=100 will get downstream failures in predict.lm, that's ok
get.data()
lm.form.df.model.false.x.true[["x"]] <- "nonesuch"
# Warning: object$x may be corrupt
rv <- plotmo(lm.form.df.model.false.x.true, trace=100)
stopifnot1(rv, X)

printf("-- test \"warning: object$x may be corrupt\", same as above but set options(warn=2)\n")
old.warn <- options(warn=2)
get.data()
lm.form.df.model.false.x.true[["x"]] <- "nonesuch"
# Warning: object$x may be corrupt
expect.err(try(plotmo(lm.form.df.model.false.x.true, trace=100)), "x may be corrupt")
options(warn=old.warn[[1]])
stopifnot1(rv, X)

printf("====== strings in the data.frame\n")

tit1 <- get.tit()

tit1$char.pclass <- as.character(tit1$pclass)

earth.survived.vs.pclass <- earth(survived~pclass, data=tit1, linpreds=TRUE)
x.earth.survived.vs.pclass <- plotmo(earth.survived.vs.pclass, trace=100, linpreds=TRUE)
stopifnot(is.factor(x.earth.survived.vs.pclass[[1]]))

earth.survived.vs.char.pclass <- earth(survived~char.pclass, data=tit1)
x.earth.survived.vs.char.pclass <- plotmo(earth.survived.vs.char.pclass, trace=100)
stopifnot(is.factor(x.earth.survived.vs.char.pclass[[1]]))

stopifnot(x.earth.survived.vs.pclass == x.earth.survived.vs.char.pclass)

lm.survived.vs.pclass <- earth(survived~pclass, data=tit1, linpreds=TRUE)
x.lm.survived.vs.pclass <- plotmo(lm.survived.vs.pclass, trace=100, linpreds=TRUE)
stopifnot(is.factor(x.lm.survived.vs.pclass[[1]]))

lm.survived.vs.char.pclass <- earth(survived~char.pclass, data=tit1)
x.lm.survived.vs.char.pclass <- plotmo(lm.survived.vs.char.pclass, trace=100)
stopifnot(is.factor(x.lm.survived.vs.char.pclass[[1]]))

stopifnot(x.lm.survived.vs.pclass == x.lm.survived.vs.char.pclass)

stopifnot(x.lm.survived.vs.pclass == x.earth.survived.vs.pclass)

printf("-- test.plotmo.x done\n")

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
