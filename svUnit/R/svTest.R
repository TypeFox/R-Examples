svTest <- function (testFun)
{
	## Create a 'svTest' object, using testFun: a function without arguments
	## that contains one or more checkXX() assertions
	if (!is.function(testFun))
		stop("'testFun' must be a function or a 'svTest' object")
	## Check that there are no arguments
	if (length(formals(testFun)) > 0)
		stop("'testFun' must be a function without any arguments")
	## This is a S3 object of class 'svTest', subclassing 'function'
	class(testFun) <- c("svTest", "function")
	return(testFun)
}

print.svTest <- function (x, ...)
{
	cat("svUnit test function:\n")
	print(body(x))
	return(invisible(x))
}

as.svTest <- function (x)
{
	## Coercion to a 'svTest' object
	return(svTest(x))
}

is.svTest <- function (x)
{
	## It this a svTest object
	return(inherits(x, "svTest"))
}

is.test <- function (x)
{
	## Is this a 'svTest'object
	## or do this object contain a non NULL 'test' attribute?
	return(is.svTest(x) || !is.null(attr(x, "test")))
}

test <- function (x)
{
	## If x is a 'svTest' object, return it, otherwise,
	## get the 'test' attribute from the object, if it exists
	if (is.svTest(x)) {
		return(x)
	} else {
		res <- attr(x, "test")
		if (is.null(res)) {
			## Create a dummy test with only a DEACTIVATED entry
			res <- svTest(function() DEACTIVATED("Object has no tests!"))
		}
		return(res)
	}
}

`test<-` <- function (x, value)
{
	## Add 'value' as a 'test' attribute to 'x' after coercing it to 'svTest'
	attr(x, "test") <- as.svTest(value)
    return(x)
}

makeUnit <- function (x, ...)
	UseMethod("makeUnit")

makeUnit.default <- function (x, name = make.names(deparse(substitute(x))),
dir = tempdir(), objfile = "", codeSetUp = NULL, codeTearDown = NULL, ...)
{
	## Take an object and make a unit from the tests it contains
	## It is saved in a file runit<name>.R in 'dir'
	name <- as.character(name)[1]
	name <- sub("^test\\.(.+)\\.$", "\\1", name)
	## Under Windows, we transform \\ into /
	dir <- gsub("\\\\", "/", as.character(dir)[1])
	Unit <- .prepareUnit(name, dir)
	## Just get the test from the object
	Test <- test(x)
	## Make required initialisation to allow locating objects
	.writeSetUp(unit = Unit, file = objfile, code = codeSetUp)
	.writeTearDown(unit = Unit, code = codeTearDown)
	## Write the test function in the file
	.writeTest(unit = Unit, objname = name, obj = x)
	## Return the name of the test function
	return(Unit)
}

makeUnit.svTest <- function (x, name = make.names(deparse(substitute(x))),
dir = tempdir(), objfile = "", codeSetUp = NULL, codeTearDown = NULL, ...)
{
	## I know: this is not needed, but it is there in case additional work
	## would be needed in the future, and also to show that makeUnit is
	## designed to work on 'svTest' objects
	return(makeUnit.default(x, name = name, dir = dir, objfile = objfile,
		codeSetUp = codeSetUp, codeTearDown = codeTearDown, ...))
}

runTest <- function (x, ...)
	UseMethod("runTest")

runTest.default <- function (x, name = deparse(substitute(x)), objfile = "",
tag = "", msg = "", ...)
{
	## Run the test for the 'test' attribute of this object
	name <- paste("test(", name, ")", sep = "")
	return(runTest(test(x), name = name, objfile = objfile, tag = tag,
		msg = msg, ...))
}

runTest.list <- function(x, ...) {
  ## Run each test in x, giving each test the name it has in x
  lapply(names(x), function(name, item=x[[name]]) {
    unit <- ifelse(is.null(attr(item, "unit")), "**root**", attr(item, "unit"))
    runTest(item, name=name, unit=unit, ...)
  })
}

runTest.svTest <- function (x, name = deparse(substitute(x)), objfile = "",
tag = "", msg = "", ...)
{
	if (!is.svTest(x))
		stop("'x' must be a 'svTest' object")
	## Names of object and test
	test <- as.character(name)[1]
	test <- .runTest(x, test = test, objfile = objfile, tag = tag, msg = msg, ...)
	.Log <- Log()
	return(invisible(.Log[[test]]))
}
