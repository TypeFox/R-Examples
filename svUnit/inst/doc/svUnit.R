### R code from vignette source 'svUnit.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: svUnit.Rnw:229-230
###################################################
options(useFancyQuotes = TRUE)


###################################################
### code chunk number 2: svUnit.Rnw:233-234 (eval = FALSE)
###################################################
## install.packages("svUnit")


###################################################
### code chunk number 3: svUnit.Rnw:245-256
###################################################
library(svUnit)
Square <- function (x) return(x^2)
test(Square) <- function () {
    checkEqualsNumeric(9, Square(3))
    checkEqualsNumeric(10, Square(3))   # This intentionally fails
    checkEqualsNumeric(9, SSSquare(3))  # This raises error
    checkEqualsNumeric(c(1, 4, 9), Square(1:3))
    checkException(Square("xx"))
}
clearLog()
(runTest(Square))


###################################################
### code chunk number 4: svUnit.Rnw:279-302
###################################################
library(svUnit)
## Create two R functions that include their own test cases
Square <- function (x) return(x^2)
test(Square) <- function () {
    checkEqualsNumeric(9, Square(3))
    checkEqualsNumeric(c(4, 9), Square(2:3))
    checkException(Square("xx"))
}

Cube <- function (x) return(x^3)
test(Cube) <- function () {
    checkEqualsNumeric(27, Cube(3))
    checkEqualsNumeric(c(8, 28), Cube(2:3))
    checkException(Cube("xx"))
}

## Add a separate test case
test_Integrate <- svTest(function () {
    checkTrue(1 < 2, "check1")
    v <- c(1, 2, 3)  # The reference
    w <- 1:3         # The value to compare to the reference
    checkEquals(v, w)
})


###################################################
### code chunk number 5: svUnit.Rnw:315-319
###################################################
clearLog()
runTest(Square)
runTest(test_Integrate)
Log()


###################################################
### code chunk number 6: svUnit.Rnw:342-344
###################################################
runTest(Cube)
Log()


###################################################
### code chunk number 7: svUnit.Rnw:407-412
###################################################
clearLog()
checkEqualsNumeric(1, log(exp(1)))
checkException(log("a"))
checkTrue(1 == 2)
Log()


###################################################
### code chunk number 8: svUnit.Rnw:471-485
###################################################
## Clear test exclusion list for running all test suites
options(svUnit.excludeList = NULL)
## Clear the logger
clearLog()
## Run all currently defined tests
runTest(svSuiteList(), name = "AllTests")
## Get some statistics
stats(Log())[, 1:3]
## A slightly different presentation than with print
summary(Log())
## Metadata collected on the machine where tests are run
metadata(Log())
## List content of the log
ls(Log())


###################################################
### code chunk number 9: svUnit.Rnw:495-500
###################################################
myTest <- Log()$testCube
class(myTest)
myTest
summary(myTest)
stats(myTest)


###################################################
### code chunk number 10: svUnit.Rnw:508-511
###################################################
ls(Log())
rm(test_R, envir = Log())
ls(Log())


###################################################
### code chunk number 11: svUnit.Rnw:532-541
###################################################
test_function <- function () {
    checkTrue(1 < 2, "check1")
    v <- c(1, 2, 3)  # The reference
    w <- 1:3         # The object to compare to the reference
    checkEqualsNumeric(v, w)
}
## Turn this function into a test
test_function <- as.svTest(test_function)
is.svTest(test_function)


###################################################
### code chunk number 12: svUnit.Rnw:554-557
###################################################
clearLog()
runTest(test_function)
Log()


###################################################
### code chunk number 13: svUnit.Rnw:567-577
###################################################
## A very simple function
Square <- function (x) return(x^2)

## A test case to associate with the Square() function
test(Square) <- function () {
    checkEqualsNumeric(9, Square(3))
    checkEqualsNumeric(c(1, 4, 9), Square(1:3))
    checkException(Square("xx"))
}
is.test(Square)  # Does this object contain tests?


###################################################
### code chunk number 14: svUnit.Rnw:582-583
###################################################
test(Square)


###################################################
### code chunk number 15: svUnit.Rnw:589-591
###################################################
runTest(Square)
Log()  # Remember we didn't clear the log!


###################################################
### code chunk number 16: svUnit.Rnw:625-628 (eval = FALSE)
###################################################
## # Create a test unit on disk and view its content
## unit <- makeUnit(Square)
## file.show(unit, delete.file = TRUE)


###################################################
### code chunk number 17: svUnit.Rnw:847-848
###################################################
example(unitTests.svUnit)


###################################################
### code chunk number 18: svUnit.Rnw:903-907
###################################################
## Reset default exclusion list
options(svUnit.excludeList = c("package:sv", "package:RUnit"))
## List all currently available tests 
svSuiteList()


###################################################
### code chunk number 19: svUnit.Rnw:916-919
###################################################
## Clear exclusion list
options(svUnit.excludeList = NULL)
svSuiteList()


###################################################
### code chunk number 20: svUnit.Rnw:928-929
###################################################
(mySuite <- svSuiteList())


###################################################
### code chunk number 21: svUnit.Rnw:934-936 (eval = FALSE)
###################################################
## myUnit <- makeUnit(mySuite, name = "ExampleTests")
## file.show(myUnit, delete.file = TRUE)


###################################################
### code chunk number 22: svUnit.Rnw:947-950
###################################################
clearLog()
runTest(mySuite)
summary(Log())


