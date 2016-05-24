library(compare)

# Utility function
testCompare <- function(FUN, ..., result=TRUE) {
    if (interactive()) {
        FUN(...)
    } else {
        stopifnot(FUN(...)$result == result)
    }
}

###############
# compareIdentical
###############

# identical integers
testCompare(compareIdentical, 1:10, 1:10)
# non-identical integers
testCompare(compareIdentical, 1:10, 10:1, result=FALSE)

# identical floats
testCompare(compareIdentical, 1:10/10, 1:10/10)
# non-identical floats
testCompare(compareIdentical, 1:10/10, 3:12/10 - 2/10, result=FALSE)

# identical strings
testCompare(compareIdentical, letters, letters)
# non-identical strings
testCompare(compareIdentical, letters, LETTERS, result=FALSE)

# identical data frames
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=1:26, y=letters, z=factor(letters),
                         stringsAsFactors=FALSE)
testCompare(compareIdentical, model, comparison)
# non-identical data frames
comparison <- data.frame(x=26:1, y=letters, z=factor(letters),
                         stringsAsFactors=FALSE)
testCompare(compareIdentical, model, comparison, result=FALSE)

# identical lists
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- model
testCompare(compareIdentical, model, comparison)
# non-identical lists
comparison <- list(a=1:26, b=letters,
              c=list(x=factor(LETTERS)))
testCompare(compareIdentical, model, comparison, result=FALSE)

###############
# compareEqual
###############

# equal integers
testCompare(compareEqual, 1:10, as.numeric(1:10))
# non-equal integers
testCompare(compareEqual, 1:10, as.numeric(10:1), result=FALSE)

# equal floats
testCompare(compareEqual, 1:10/10, 3:12/10 - 2/10)
# non-equal floats
testCompare(compareEqual, 1:10/10, 10:1/10, result=FALSE)

# equal floats (rounding)
testCompare(compareEqual, 1:10, 1:10 + .1, round=TRUE)
# non-equal floats (rounding)
testCompare(compareEqual, 1:10, 1:10 + .6, round=TRUE, result=FALSE)

# equal floats (round=signif)
testCompare(compareEqual, c(.1, 1, 10), c(.1, 1, 10),
            round=function(x) { signif(x, 1) })
# non-equal floats (round=signif)
testCompare(compareEqual, c(.1, 1, 10), c(.13, 1.3, 13),
            round=function(x) { signif(x, 1) })
# non-equal floats (rounding)
testCompare(compareEqual, c(.1, 1, 10), c(.13, 1.3, 13),
            round=TRUE, result=FALSE)

# equal strings
testCompare(compareEqual, letters, paste(" ", letters, " "),
            trim=TRUE)
# non-equal strings
testCompare(compareEqual, letters, paste(" ", LETTERS, " "),
            trim=TRUE, result=FALSE)

# equal strings
testCompare(compareEqual, letters, LETTERS,
            ignoreCase=TRUE)
# non-equal strings
testCompare(compareEqual, letters, rev(LETTERS),
            ignoreCase=TRUE, result=FALSE)

# equal factors
testCompare(compareEqual, factor(letters), factor(letters))
# non-equal factors
testCompare(compareEqual, factor(letters), factor(LETTERS), result=FALSE)

# equal factors
testCompare(compareEqual, factor(letters[1:10]),
            factor(letters[1:10], levels=letters),
            dropLevels=TRUE)
# non-equal factors
testCompare(compareEqual, factor(letters[1:10]),
            factor(LETTERS[1:10], levels=LETTERS),
            dropLevels=TRUE, result=FALSE)

# equal factors
testCompare(compareEqual, factor(letters),
            factor(letters, levels=rev(letters)),
            ignoreLevelOrder=TRUE)
# non-equal factors
testCompare(compareEqual, factor(letters),
            factor(rev(letters), levels=rev(letters)),
            ignoreLevelOrder=TRUE, result=FALSE)

# equal matrices/arrays/tables
testCompare(compareEqual, cbind(1:10, 10:1), cbind(1:10, 10:1))
# non-equal matrices/arrays/tables
testCompare(compareEqual, cbind(1:10, 10:1), cbind(10:1, 1:10), result=FALSE)
# equal matrices/arrays/tables
testCompare(compareEqual, array(1:8, dim=rep(2, 3)),
            array(1:8, dim=rep(2, 3)))
# non-equal matrices/arrays/tables
testCompare(compareEqual, array(1:8, dim=rep(2, 3)),
            array(8:1, dim=rep(2, 3)), result=FALSE)
# equal matrices/arrays/tables
testCompare(compareEqual, table(a=rep(1:2, 5), b=rep(1:5, 2)),
            table(a=rep(1:2, 5), b=rep(1:5, 2)))
# non-equal matrices/arrays/tables
testCompare(compareEqual, table(a=rep(1:2, 5), b=rep(1:5, 2)),
            table(a=rep(3:4, 5), b=rep(1:5, 2)), result=FALSE)
# non-equal matrices/arrays/tables (because of dimnames)
testCompare(compareEqual, table(a=rep(1:2, 5), b=rep(1:5, 2)),
            table(b=rep(1:2, 5), a=rep(1:5, 2)), result=FALSE)
# equal matrices/arrays/tables
testCompare(compareEqual, cbind(1:10, 10:1), cbind(1:10, 10:1) + .1,
            round=TRUE)
# non-equal matrices/arrays/tables
testCompare(compareEqual, cbind(1:10, 10:1), cbind(10:1, 1:10) + .1,
            round=TRUE, result=FALSE)
# equal matrices/arrays/tables (ignoring dim order)
testCompare(compareEqual, table(a=1:2, b=3:4),
            table(b=3:4, a=1:2),
            ignoreDimOrder=TRUE)
# non-equal matrices/arrays/tables (ignoring dim order)
testCompare(compareEqual, table(a=1:2, b=3:4),
            table(b=3:4, a=1:2), result=FALSE)
# equal matrices/arrays/tables
testCompare(compareEqual, cbind(letters, rev(letters)),
            cbind(letters, rev(letters)))
# non-equal matrices/arrays/tables
testCompare(compareEqual, cbind(letters, rev(letters)),
            cbind(rev(letters), letters), result=FALSE)
# equal matrices/arrays/tables
# NOTE cbind() creates col names so I use matrix() instead
testCompare(compareEqual, matrix(c(letters, rev(letters)), ncol=2),
            matrix(c(LETTERS, rev(LETTERS)), ncol=2),
            ignoreCase=TRUE)
# non-equal matrices/arrays/tables
testCompare(compareEqual, matrix(c(letters, rev(letters)), ncol=2),
            matrix(c(rev(LETTERS), LETTERS), ncol=2),
            ignoreCase=TRUE, result=FALSE)
                                      
# equal data frames
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=1:26, y=letters, z=factor(letters),
                         stringsAsFactors=FALSE)
testCompare(compareEqual, model, comparison)
# non-equal data frames
comparison <- data.frame(x=1:26, y=LETTERS, z=factor(letters),
                         stringsAsFactors=FALSE)
testCompare(compareEqual, model, comparison, result=FALSE)

# equal data frames (columns NOT identical)
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=as.numeric(1:26), y=LETTERS,
                         z=factor(letters, levels=rev(letters)),
                         stringsAsFactors=FALSE)
testCompare(compareEqual, model, comparison,
            ignoreCase=TRUE, ignoreLevelOrder=TRUE)
# non-equal data frames
testCompare(compareEqual, model, comparison, result=FALSE)

# equal data frames (allow different col order)
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(y=letters, z=factor(letters), x=1:26, 
                         stringsAsFactors=FALSE)
testCompare(compareEqual, model, comparison,
            ignoreColOrder=TRUE)
# non-equal data frames
testCompare(compareEqual, model, comparison, result=FALSE)

# non-equal data frames with zero columns
model <- data.frame()
comparison <- model
attr(comparison, "difference") <- 1
testCompare(compareEqual, model, comparison, result=FALSE)

# equal lists
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- model
testCompare(compareEqual, model, comparison)
# equal lists (allow different component order)
comparison <- list(a=1:26, 
                   c=list(x=factor(letters)),
                   b=letters)
testCompare(compareEqual, model, comparison,
            ignoreComponentOrder=TRUE)
# equal lists (allow different component order and different name case)
comparison <- list(A=1:26, 
                   C=list(x=factor(letters)),
                   B=letters)
testCompare(compareEqual, model, comparison,
            ignoreComponentOrder=TRUE,
            ignoreNameCase=TRUE)
# non-equal lists
comparison <- list(a=1:26, b=letters,
                   c=list(x=1:26))
testCompare(compareEqual, model, comparison, result=FALSE)

# equal "lm" objects
x <- 1:10
y <- rnorm(10)
lm1 <- lm(y ~ x)
testCompare(compareEqual, lm1, lm1)

# non-equal "lm" objects
y2 <- rnorm(10)
lm2 <- lm(y2 ~ x)
testCompare(compareEqual, lm1, lm2, result=FALSE)

# equal expressions (also tests "call")
testCompare(compareEqual, expression(x + y), expression(x + y))

# non-equal expressions
testCompare(compareEqual, expression(x + y), expression(y + x), result=FALSE)

###############
# compareCoerce
###############

# identical float from integer
testCompare(compareCoerce, as.numeric(1:10), 1:10)
# non-identical float from integer
testCompare(compareCoerce, as.numeric(1:10), 10:1, result=FALSE)

# identical float from string
testCompare(compareCoerce, 1:10, as.character(1:10))
# non-identical float from string
testCompare(compareCoerce, 1:10, as.character(10:1), result=FALSE)

# equal float from string
testCompare(compareCoerce, 1:10/10,
            c(".1", ".2", ".3", ".4", ".5",
              ".6", ".7", ".8", ".9", "1.0"))
# non-equal float from string
testCompare(compareCoerce, 10:1/10,
            c(".1", ".2", ".3", ".4", ".5",
              ".6", ".7", ".8", ".9", "1.0"), result=FALSE)

# identical string from factor
testCompare(compareCoerce, letters, factor(letters))
# non-identical string from factor
testCompare(compareCoerce, letters, factor(LETTERS), result=FALSE)

# identical factor from string
testCompare(compareCoerce, factor(letters), letters)
# non-identical factor from string 
testCompare(compareCoerce, factor(letters), LETTERS, result=FALSE)

# equal data frames
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=1:26, y=factor(letters), z=letters,
                         stringsAsFactors=FALSE)
testCompare(compareCoerce, model, comparison)
# non-equal data frames
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=26:1, y=factor(LETTERS), z=LETTERS,
                         stringsAsFactors=FALSE)
testCompare(compareCoerce, model, comparison, result=FALSE)

# equal lists
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- model
testCompare(compareCoerce, model, comparison)
# equal lists (after coercion)
model <- list(a=1:26, b=letters)
comparison <- list(a=1:26, 
                   b=letters)
class(comparison) <- "dummy"
testCompare(compareCoerce, model, comparison)
# equal lists (names in different order)
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- list(a=1:26, 
                   c=list(x=factor(letters)),
                   b=letters)
testCompare(compareCoerce, model, comparison,
            ignoreComponentOrder=TRUE)
# equal lists (names in different order and different case)
comparison <- list(A=1:26, 
                   C=list(x=factor(letters)),
                   B=letters)
testCompare(compareCoerce, model, comparison,
            ignoreComponentOrder=TRUE,
            ignoreNameCase=TRUE)
# non-equal lists
comparison <- list(a=1:26, b=letters,
                   c=list(x=26:1))
testCompare(compareCoerce, model, comparison, result=FALSE)

# Some tests to make sure that simple comparisons pass through
# equal "lm" objects
x <- 1:10
y <- rnorm(10)
lm1 <- lm(y ~ x)
testCompare(compareCoerce, lm1, lm1)

# non-equal "lm" objects
y2 <- rnorm(10)
lm2 <- lm(y2 ~ x)
testCompare(compareCoerce, lm1, lm2, result=FALSE)

# equal expressions (also tests "call")
testCompare(compareCoerce, expression(x + y), expression(x + y))

# non-equal expressions
testCompare(compareCoerce, expression(x + y), expression(y + x), result=FALSE)

###############
# compareShorten
###############

# identical integer
testCompare(compareShorten, 1:5, 1:5)
# non-identical integer
testCompare(compareShorten, 1:5, 5:1, result=FALSE)
# identical integer, model longer
testCompare(compareShorten, 1:10, 1:5)
# identical integer, comparison longer
testCompare(compareShorten, 1:5, 1:10)
# non-identical integer, comparison longer
testCompare(compareShorten, 1:5, 10:1, result=FALSE)

# identical matrix
testCompare(compareShorten,
            matrix(1:10, ncol=2),
            matrix(1:10, ncol=2))
# non-identical matrix
testCompare(compareShorten,
            matrix(1:10, ncol=2),
            matrix(10:1, ncol=2), result=FALSE)
# identical matrix, different number of dims
testCompare(compareShorten,
            matrix(1:10, ncol=2),
            array(1:100, dim=c(5, 2, 10)))
# identical array, different dims
testCompare(compareShorten,
            array(1:100, dim=c(5, 2, 10)),
            array(1:100, dim=c(5, 2))[, , drop=FALSE])
# identical table, different dims
testCompare(compareShorten, Titanic,
            apply(Titanic, 1:2, sum))

# equal lists
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- model
testCompare(compareShorten, model, comparison)
# equal lists (model longer)
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- list(a=1:26, b=letters)
testCompare(compareShorten, model, comparison)
# equal lists (comparison longer)
model <- list(a=1:26, b=letters)
comparison <- list(a=1:26, b=letters,
                   c=list(x=factor(letters)))
testCompare(compareShorten, model, comparison)

# Some tests to make sure that simple comparisons pass through
# equal "lm" objects
x <- 1:10
y <- rnorm(10)
lm1 <- lm(y ~ x)
testCompare(compareShorten, lm1, lm1)

# non-equal "lm" objects
y2 <- rnorm(10)
lm2 <- lm(y2 ~ x)
testCompare(compareShorten, lm1, lm2, result=FALSE)

# equal expressions (also tests "call")
testCompare(compareShorten, expression(x + y), expression(x + y))

# non-equal expressions
testCompare(compareShorten, expression(x + y), expression(y + x), result=FALSE)

###############
# compareIgnoreOrder
###############

# identical integer, different order
testCompare(compareIgnoreOrder, 1:10, 10:1)
# non-identical integer (regardless of order)
testCompare(compareIgnoreOrder, 1:10, 11:2, result=FALSE)

# identical float, different order
testCompare(compareIgnoreOrder, 10:1/10, 3:12/10 - 2/10)
# non-identical float, different order
testCompare(compareIgnoreOrder, 11:2/10, 3:12/10 - 2/10, result=FALSE)

# identical table, different order WITHIN dimensions
testCompare(compareIgnoreOrder,
            Titanic,
            do.call("[", c(list(Titanic),
                           lapply(dimnames(Titanic),
                                  order))))

# identical data frame, different order
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=1:26, y=letters, z=factor(letters),
                         stringsAsFactors=FALSE)[26:1, ]
testCompare(compareIgnoreOrder, model, comparison)
# non-identical data frame, different order
comparison <- data.frame(x=1:26, y=LETTERS, z=factor(LETTERS),
                         stringsAsFactors=FALSE)[26:1, ]
testCompare(compareIgnoreOrder, model, comparison, result=FALSE)

# equal lists
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- model
testCompare(compareIgnoreOrder, model, comparison)
# equal lists (different component order)
comparison <- list(a=1:26, 
                   c=list(x=factor(letters)),
                   b=letters)
testCompare(compareIgnoreOrder, model, comparison)
# equal lists (different component order and case)
comparison <- list(A=1:26, 
                   C=list(x=factor(letters)),
                   B=letters)
testCompare(compareIgnoreOrder, model, comparison,
            ignoreNameCase=TRUE)

# Some tests to make sure that simple comparisons pass through
# equal "lm" objects
x <- 1:10
y <- rnorm(10)
lm1 <- lm(y ~ x)
testCompare(compareIgnoreOrder, lm1, lm1)

# non-equal "lm" objects
y2 <- rnorm(10)
lm2 <- lm(y2 ~ x)
testCompare(compareIgnoreOrder, lm1, lm2, result=FALSE)

# equal expressions (also tests "call")
testCompare(compareIgnoreOrder, expression(x + y), expression(x + y))

# non-equal expressions
testCompare(compareIgnoreOrder,
            expression(x + y), expression(y + x), result=FALSE)


###############
# compareIgnoreNameCase
###############

# identical float, identical names
model <- 1:10
names(model) <- letters[1:10]
comparison <- 1:10
names(comparison) <- letters[1:10]
testCompare(compareIgnoreNameCase, model, comparison)
# non-identical float, identical names
comparison <- 10:1
names(comparison) <- letters[1:10]
testCompare(compareIgnoreNameCase, model, comparison, result=FALSE)

# identical float, names differ by case
comparison <- 1:10
names(comparison) <- LETTERS[1:10]
testCompare(compareIgnoreNameCase, model, comparison)
# non-identical float, names differ by case
comparison <- 10:1
names(comparison) <- LETTERS[1:10]
testCompare(compareIgnoreNameCase, model, comparison, result=FALSE)

# identical float, names differ
comparison <- 1:10
names(comparison) <- letters[10:1]
testCompare(compareIgnoreNameCase, model, comparison, result=FALSE)
# non-identical float, names differ
comparison <- 10:1
names(comparison) <- letters[10:1]
testCompare(compareIgnoreNameCase, model, comparison, result=FALSE)

# identical data frame, colnames differ by case
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(X=1:26, Y=letters, Z=factor(letters),
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNameCase, model, comparison)
# non-identical data frame, colnames differ by case
comparison <- data.frame(X=1:26, Y=LETTERS, Z=factor(LETTERS),
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNameCase, model, comparison, result=FALSE)

# identical data frame, rownames differ by case
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    row.names=letters,
                    stringsAsFactors=FALSE)
comparison <- data.frame(X=1:26, Y=letters, Z=factor(letters),
                         row.names=LETTERS,
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNameCase, model, comparison,
            colsOnly=FALSE)
# non-identical data frame, rownames differ by case
comparison <- data.frame(X=1:26, Y=LETTERS, Z=factor(LETTERS),
                         row.names=LETTERS,
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNameCase, model, comparison,
            colsOnly=FALSE, result=FALSE)

# equal lists
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- model
testCompare(compareIgnoreNameCase, model, comparison)
# equal lists (different component order)
comparison <- list(a=1:26, 
                   c=list(x=factor(letters)),
                   b=letters)
testCompare(compareIgnoreNameCase, model, comparison,
            ignoreComponentOrder=TRUE)
# equal lists (different component order and case)
comparison <- list(A=1:26, 
                   C=list(x=factor(letters)),
                   B=letters)
testCompare(compareIgnoreNameCase, model, comparison,
            ignoreComponentOrder=TRUE)

###############
# compareIgnoreNames
###############

model <- 1:10
names(model) <- letters[1:10]

# identical float, identical names
comparison <- 1:10
names(comparison) <- letters[1:10]
testCompare(compareIgnoreNames, model, comparison)
# non-identical float, identical names
comparison <- 10:1
names(comparison) <- letters[1:10]
testCompare(compareIgnoreNames, model, comparison, result=FALSE)

# identical float, names differ
comparison <- 1:10
names(comparison) <- letters[10:1]
testCompare(compareIgnoreNames, model, comparison)
# identical float, names differ
comparison <- 10:1
names(comparison) <- letters[10:1]
testCompare(compareIgnoreNames, model, comparison, result=FALSE)

# identical data frame, names differ
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(a=1:26, b=letters, c=factor(letters),
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNames, model, comparison)
# non-identical data frame, names differ
comparison <- data.frame(a=1:26, b=LETTERS, c=factor(LETTERS),
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNames, model, comparison, result=FALSE)

# identical data frame, rownames differ
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    row.names=letters,
                    stringsAsFactors=FALSE)
comparison <- data.frame(a=1:26, b=letters, c=factor(letters),
                         row.names=LETTERS,
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNames, model, comparison,
            colsOnly=FALSE)
# non-identical data frame, rownames differ
comparison <- data.frame(a=1:26, b=LETTERS, c=factor(LETTERS),
                         row.names=LETTERS,
                         stringsAsFactors=FALSE)
testCompare(compareIgnoreNames, model, comparison,
            colsOnly=FALSE, result=FALSE)

# equal lists 
model <- list(a=1:26, b=letters,
              c=list(x=factor(letters)))
comparison <- model
testCompare(compareIgnoreNames, model, comparison)
# equal lists (different component order)
comparison <- list(a=1:26, 
                   c=list(x=factor(letters)),
                   b=letters)
testCompare(compareIgnoreNames, model, comparison,
            ignoreComponentOrder=TRUE)
# equal lists (different component order and case)
comparison <- list(A=1:26, 
                   C=list(x=factor(letters)),
                   B=letters)
testCompare(compareIgnoreNames, model, comparison,
            ignoreComponentOrder=TRUE,
            ignoreNameCase=TRUE)
# equal lists (same component order BUT different names)
comparison <- list(x=1:26, 
                   y=letters,
                   z=list(x=factor(letters)))
testCompare(compareIgnoreNames, model, comparison,
            ignoreComponentOrder=TRUE,
            ignoreNameCase=TRUE)

###############
# compareIgnoreAttrs
###############

model <- matrix(1:10, ncol=2)

# identical float, identical attrs
comparison <- matrix(1:10, ncol=2)
testCompare(compareIgnoreAttrs, model, comparison)
# non-identical float, identical attrs
comparison <- matrix(10:1, ncol=2)
testCompare(compareIgnoreAttrs, model, comparison, result=FALSE)

# identical float, attrs differ
comparison <- matrix(1:10, ncol=5)
testCompare(compareIgnoreAttrs, model, comparison)
# non-identical float, identical attrs
comparison <- matrix(10:1, ncol=5)
testCompare(compareIgnoreAttrs, model, comparison, result=FALSE)

# identical data frame, attrs differ
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(a=1:26, b=letters, c=factor(letters),
                         stringsAsFactors=FALSE)
attr(comparison, "test") <- "test"
testCompare(compareIgnoreAttrs, model, comparison)
# non-identical data frame, attrs differ
comparison <- data.frame(a=1:26, b=LETTERS, c=factor(LETTERS),
                         stringsAsFactors=FALSE)
attr(comparison, "test") <- "test"
testCompare(compareIgnoreAttrs, model, comparison, result=FALSE)

# equal lists
# (after dropping attributes;  coercion does not take care of rownames!)
model <- list(a=1:26, b=letters)
comparison <- model
attr(comparison, "test") <- "test"
testCompare(compareIgnoreAttrs, model, comparison)

###############
# Combining comparisons
###############

# coerce AND sort
model <- 10:1/10
comparison <- c(".1", ".2", ".3", ".4", ".5",
                ".6", ".7", ".8", ".9", "1.0")
# fails then succeeds
comp <- compareCoerce(model, comparison)
if (!comp$result) {
    compareIgnoreOrder(comp$tM, comp$tC, comp$transform)
}
# Other way around fails BOTH TIMES
# (sorting character produces different order than sorting numeric)
# NOTE though that the sort order is system-dependent
# (on Windows, this fails then succeeds!)
comp <- compareIgnoreOrder(model, comparison)
if (!comp$result) {
    compareCoerce(comp$tM, comp$tC, comp$transform)
}

###############
# compare()
###############

# Obviously identical
compare(1:10, 1:10)
# Sensibly identical
compare(1:10/10, 3:12/10 - 2/10)
# Being picky
compare(1:10/10, 3:12/10 - 2/10, equal=FALSE)
# Sensibly different
compare(1:10/10,
        c(".1", ".2", ".3", ".4", ".5",
          ".6", ".7", ".8", ".9", "1.0"))
# Being generous
compare(1:10/10,
        c(".1", ".2", ".3", ".4", ".5",
          ".6", ".7", ".8", ".9", "1.0"),
        coerce=TRUE)
# Bending over backwards
compare(10:1/10,
        c(".1", ".2", ".3", ".4", ".5",
          ".6", ".7", ".8", ".9", "1.0"),
        coerce=TRUE, ignoreOrder=TRUE)
compare(10:1/10, 1:10/10,
        coerce=TRUE, ignoreOrder=TRUE)
# You're crazy !!
inTheBallPark <- matrix(c(".1", ".2", ".3", ".4", ".5",
                          ".6", ".7", ".8", ".9", "1.0"),
                        ncol=2)
compare(1:10/10,
        inTheBallPark,
        coerce=TRUE, ignoreOrder=TRUE, ignoreAttrs=TRUE)
# Testing passing special argument, 'ignoreCase', through to compareEqual()
model <- letters
comparison <- matrix(LETTERS, ncol=2)
# Should fail
compare(model, comparison)
# Should still fail
compare(model, comparison, ignoreAttrs=TRUE)
# Should now succeed
compare(model, comparison, ignoreCase=TRUE, ignoreAttrs=TRUE)

# Doing some data frame compares ...
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=1:26, y=letters, z=factor(letters),
                         stringsAsFactors=FALSE)
compare(model, comparison)
# Allowing any transformation at all
comparison <- data.frame(X=1:26, Y=letters, Z=factor(letters),
                         row.names=letters,
                         stringsAsFactors=FALSE)
compare(model, comparison, allowAll=TRUE)
# Allowing any transformation of columns
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    row.names=letters,
                    stringsAsFactors=FALSE)
comparison <- data.frame(X=1:26, Y=factor(LETTERS), Z=letters,
                         row.names=LETTERS,
                         stringsAsFactors=FALSE)
compare(model, comparison, allowAll=TRUE)
# Allowing any transformation of columns
# BUT the data frames ARE DIFFERENT!
# NOTE that tests for dropping names and dropping
# attributes are considered, but not performed
# because having done previous checks ignoring
# case of names and case of row names, these
# comparisons are unnecessary(!)
model <- data.frame(x=1:26, y=letters, z=factor(letters),
                    row.names=letters,
                    stringsAsFactors=FALSE)
comparison <- data.frame(X=1:26, Y=factor(LETTERS), Z=rev(letters),
                         row.names=LETTERS,
                         stringsAsFactors=FALSE)
compare(model, comparison, allowAll=TRUE)

# Testing whether compareEqual transforms are retried if they
# fail and are only recorded if they succeed
# equal floats (rounding)
compare(1:10, 1:10 + .1, round=TRUE)
# non-equal floats (rounding)
compare(1:10, 1:10 + .6, round=TRUE)
# equal after coercion
# Should ONLY report rounding AFTER coercion
compare(as.numeric(1:10), as.character(1:10 + .1),
        round=TRUE, coerce=TRUE)
# non-equal after coercion
compare(as.numeric(1:10), as.character(1:10 + .6),
        round=TRUE, coerce=TRUE)
# data frame with column that needs sorting THEN rounding
model <- data.frame(x=as.numeric(1:10))
comparison <- data.frame(x=10:1 + .1)
compare(model, comparison, round=TRUE, allowAll=TRUE)

# Tests of recurseFun for data frame comparisons
# Standalone compare should report all transforms tried on columns
# whether succeeds or fails
model <- data.frame(x=1:26, 
                    y=letters, 
                    z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=1:26, 
                         y=factor(letters), 
                         z=letters,
                         stringsAsFactors=FALSE)
compareCoerce(model, comparison)
comparison <- data.frame(X=1:26, 
                         Y=factor(letters), 
                         Z=letters,
                         stringsAsFactors=FALSE)
compareCoerce(model, comparison)
# Overall compare should only report each transform that was tried ONCE
model <- data.frame(x=1:26, 
                    y=letters, 
                    z=factor(letters),
                    stringsAsFactors=FALSE)
comparison <- data.frame(x=1:26, 
                         y=factor(letters), 
                         z=letters,
                         stringsAsFactors=FALSE)
compare(model, comparison, allowAll=TRUE)
comparison <- data.frame(X=1:26, 
                         Y=factor(letters), 
                         Z=letters,
                         stringsAsFactors=FALSE)
compare(model, comparison, allowAll=TRUE)

# Example where columns need to be coerced before being sorted
# (and the coercion needs to be persistent)
model <- data.frame(x=1:10)
comparison <- data.frame(x=as.character(10:1))
testCompare(compare, model, comparison, allowAll=TRUE)

# Example where columns need to be put in the right order
# BEFORE coercing the columns.
model <- data.frame(x=1:26, 
                    y=letters, 
                    z=factor(letters),
                    row.names=letters,
                    stringsAsFactors=FALSE)
comparison <- data.frame(X=1:26, 
                         Z=letters,
                         Y=factor(LETTERS), 
                         row.names=LETTERS,
                         stringsAsFactors=FALSE)
testCompare(compare, model, comparison, allowAll=TRUE)

# List versus data frame (need to drop attributes after coercion
# to get rid of row names)
model <- list(a=1:26, b=letters)
comparison <- data.frame(a=1:26, 
                         b=letters)
testCompare(compare, model, comparison, allowAll=TRUE)

# Proper "lm" comparison
lm1CO2 <- lm(uptake ~ Treatment + conc, data=CO2)
lm2CO2 <- lm(uptake ~ conc + Treatment, data=CO2)
# NOTE: overall result FALSE,
# BUT 'residuals' and 'fitted.values' components are the same !
compare(lm1CO2, lm2CO2)
# With allowAll=TRUE, still overall FALSE, 
# BUT now 'coefficients' and 'model' are ALSO the same !
compare(lm1CO2, lm2CO2, allowAll=TRUE)

###############
# compareName()
###############

# identical 
model <- 1:10
x <- 1:10
testCompare(compareName, model, "x")
testCompare(compareName, model, "x", ignore.case=FALSE)
# non-identical 
x <- 10:1
testCompare(compareName, model, "x", result=FALSE)
testCompare(compareName, model, "x", ignore.case=FALSE, result=FALSE)

# identical, ignore case
model <- 1:10
X <- 1:10
rm("x") # !!!!!
testCompare(compareName, model, "x")
testCompare(compareName, model, "x", ignore.case=FALSE, result=FALSE)
# non-identical, ignore case
X <- 10:1
testCompare(compareName, model, "x", result=FALSE)
testCompare(compareName, model, "x", ignore.case=FALSE, result=FALSE)

# identical, ignore case, local environment
model <- 1:10
tempenv <- new.env()
with(tempenv, X <- 1:10)
testCompare(compareName, model, "x", compEnv=tempenv)
testCompare(compareName, model, "x", ignore.case=FALSE, compEnv=tempenv, result=FALSE)

###############
# compareFile
###############

# Allow rounding of one answer out of several
# No rounding; both FALSE
compList <- compareFile(textConnection("{ x <- y <- 1:10 + .1 }"),
                        modelNames=c("x", "y"),
                        modelAnswers=list(x=1:10, y=1:10))
stopifnot(all(sapply(compList, function(x) { !isTRUE(x) })))
# All rounding; both TRUE
compList <- compareFile(textConnection("{ x <- y <- 1:10 + .1 }"),
                        modelNames=c("x", "y"),
                        modelAnswers=list(x=1:10, y=1:10),
                        round=TRUE)
stopifnot(all(sapply(compList, isTRUE)))
# Only round 'y'; only 'y' TRUE
compList <- compareFile(textConnection("{ x <- y <- 1:10 + .1 }"),
                        modelNames=c("x", "y"),
                        modelAnswers=list(x=1:10, y=1:10),
                        round=list(y=TRUE))
stopifnot(!isTRUE(compList[[1]]) && isTRUE(compList[[2]]))

# All rounding; both TRUE
compList <- compareFile(textConnection("{ x <- y <- 1:10 + .1 }"),
                        modelNames=c("x", "y"),
                        modelAnswers=list(x=1:10, y=1:10),
                        round=floor)
stopifnot(all(sapply(compList, isTRUE)))
# Only round 'y'; only 'y' TRUE
compList <- compareFile(textConnection("{ x <- y <- 1:10 + .1 }"),
                        modelNames=c("x", "y"),
                        modelAnswers=list(x=1:10, y=1:10),
                        round=list(y=function(x) signif(x, 1)))
stopifnot(!isTRUE(compList[[1]]) && isTRUE(compList[[2]]))


