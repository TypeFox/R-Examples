#' # pyObject
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

## Object
expect_that(pyExec("import os"), equals(0))
expect_that(class(pyObject("os", regFinalizer = FALSE))[1], equals("module.os"))
expect_that(class(pyObject("os", regFinalizer = TRUE))[1], equals("module.os"))

## Tuple
expect_that(pyExec('myPyTuple = (1, 2, 5, "Hello R!")'), equals(0))
# create a virtual Python tuple for an existing tuple
myTuple <- pyTuple("myPyTuple")
expect_that(class(myTuple)[1], equals("PythonInR_Tuple"))
expect_that(myTuple[1], equals(2))
expect_that(myTuple[1] <- "", throws_error())
# create a new Python tuple and virtual tuple
newTuple <- pyTuple('myNewTuple', list(1:3, 'Hello Python'))
expect_that(newTuple[0], equals(1:3))

## List
expect_that(pyExec('myPyList = [1, 2, 5, "Hello R!"]'), equals(0))
# create a virtual Python list for an existing list
myList <- pyList("myPyList")
expect_that(myList[1], equals(2))
# create a new Python list and virtual list
myNewList <- pyList('myNewList', list(1:3, 'Hello Python'))
expect_that(myNewList[0], equals(1:3))

## Dict
expect_that(pyExec('myPyDict = {"a":1, "b":2, "c":3}'), equals(0))
# create a virtual Python dictionary for an existing dictionary
myDict <- pyDict("myPyDict")
expect_that(myDict["a"], equals(1))
# create a new Python dict and virtual dict
myNewDict <- pyDict('myNewDict', list(p=2, y=9, r=1))
expect_that(myNewDict["r"], equals(1))
