#' # pyExec
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

expect_that(pyExec("import os"), equals(0))

#' ## get/set current working directory
## by default
expect_that(pyCall("os.chdir", args=list(getwd())), equals(NULL))
expect_that(gsub("\\", "/", pyCall("os.getcwd"), fixed=TRUE), equals(getwd()))

#' ## test builtins functions (__builtins__)
pyImport("sys")
expect_that(pyCall("abs", args=list(-5)), equals(5))
expect_that(pyCall("sum", args=list(1:5)), equals(sum(1:5)))
#' NOTE: since all integer variables are translated to long by default
#' the following code produces an error
expect_that(
    expect_that(
        pyCall("hex", args=list(255), namespace=builtinNsp),
        prints_text("TypeError")
    ),
    throws_error()
)

#' ## Test pyExec
#' pyExec can execute multiple lines printed values are shown in the R stdout
#' evaluated values not.
cmd <- "
x = 3
y = 4
z = x * y
print(z)
"
expect_that(pyExec(cmd), prints_text("12"))

#' ### Test error handling
expect_that(expect_that(pyExec("x=4/0"), prints_text("ZeroDivisionError")), throws_error())

#' ## Test pyExecp
#' pyExecp executes a single line, every evaluated statement is printed in the R stdout
expect_that(pyExecp("z"), prints_text("12"))
#' ### Test error handling
expect_that(expect_that(pyExecp("x=4/0"), prints_text("ZeroDivisionError")), throws_error())

#' ## Differences between pyExec and pyExecp
#' pyExecp executes only a the first line of the code the others are obmitted
#' pyExecp is intended to behave more like the Python terminal
pyExec('
def fun():
    return("Hello R!")
')
expect_that(pyExec("fun()"), prints_text("^([A-Z]+)?$", perl = TRUE))
expect_that(pyExecp("fun()"), prints_text("Hello R"))

#' ## Test pyExecg
#' pyExecg executes the provided code and returns the during the execution assigned values
expect_that(pyExecg('x = 5*5')[['x']], equals(25))
#' ### Test error handling
expect_that(expect_that(pyExecg("x=4/0"), prints_text("ZeroDivisionError")), throws_error())  

#' ### Test different options of pyExecg
expect_that(pyExecg("x=fun()")[['x']], equals("Hello R!"))
#' #### returnToR
expect_that(pyExecg("x=4", returnToR=FALSE)[['x']], equals(NULL))
#' #### mergeNamespaces 
expect_that(pyExecg("some_new_variable=4", mergeNamespaces=TRUE, override = TRUE)[[1]], equals(4))
expect_that(pyPrint("some_new_variable"), prints_text("4"))
#' #### override
expect_that(pyExecg("some_new_variable=1", mergeNamespaces=TRUE, override=FALSE)[[1]], equals(1))
# should be still 4 since override is FALSE
expect_that(pyPrint("some_new_variable"), prints_text("4"))
expect_that(pyExecg("some_new_variable2=5", mergeNamespaces=TRUE, override=FALSE)[[1]], equals(5))
# show that the variable get's assigned when it doesn't already exits
expect_that(pyPrint("some_new_variable2"), prints_text("5"))
expect_that(pyExecg("some_new_variable=1", mergeNamespaces=TRUE, override=TRUE)[[1]], equals(1))
# should be 1 since override is TRUE
expect_that(pyPrint("some_new_variable"), prints_text("1"))

#' ### Test error handling
cmd <- '
a = "this is a multi line test"
b = 3
b = 4/0
b = 5
'
expect_that(
    expect_that(
        pyExecg(cmd, mergeNamespaces=TRUE),
        prints_text("ZeroDivisionError")
        ),
    throws_error()
    )

#' ***NOTE:*** since all the commands are executed in a new namespace
#' non of the script will have any effect since the function will exit before
#' the namespaces are merged
expect_that(
    expect_that(
        pyPrint(b),
        prints_text("NameError")
        ),
    throws_error())

## -----------------------------------------------
##
#' ## pyExecfile
##
## -----------------------------------------------
pyExecfile(file.path(path.package("PythonInR"), "testing/Test_cases.py"))

myInt = 6
myDouble = 3.14
myString = "Test String!"
# <<< NOTE: This is necessary since my test cases are written in Linux (utf-8)
#           when I run them on Windows, Windows will assume that the file has
#           the local encoding and produce an error (since the encoding of the
#           reference variable is messed up) even when the encoding in Python
#           is correct. (This worked on Windows with latin1 as default encoding) >>>
myUnicode = iconv('Äöüß 945 hdfji', from="UTF-8")
    
myList = list(2, 3, "Hallo")
myTuple = list(1, 2, "Hallo")
mySet = list(myTuple)

expect_that(pyGet("myInt"), equals(myInt))
expect_that(pyGet("myDouble"), equals(myDouble))
expect_that(pyGet("myString"), equals(myString))
expect_that(pyGet("myUnicode"), equals(myUnicode))
expect_that(pyGet("myList"), equals(myList))
expect_that(pyGet("myTuple"), equals(myTuple))

