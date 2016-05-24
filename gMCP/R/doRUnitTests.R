#' Run the R unit (and optional the JUnit) test suite for gMCP
#' 
#' Runs the R unit (and optional the JUnit) test suite for gMCP and prints the
#' results.
#' 
#' The environment variable GMCP_UNIT_TESTS may be used to specify which
#' unit tests should run: "extended", "interactive", "java" or a combination
#' of these separated by comma (without blanks). A short cut for all three is
#' "all".
#' 
#' @param extended If \code{TRUE} (or if the environment variable
#' GMCP_UNIT_TESTS equals "extended" or "all") an extended version of the R
#' unit test suite for gMCP will be used.  The run will take significantly
#' longer time.
#' @param java If \code{TRUE} (or if the environment variable
#' GMCP_UNIT_TESTS equals "java" or "all") the GUI and its logic is tested with JUnit tests.
#' You need JUnit 4 classes in the classpath or specify the path to a JUnit 4
#' jar file via the parameter junitLibrary.
#' @param interactive If \code{TRUE} (or if the environment variable
#' GMCP_UNIT_TESTS equals "interactive" or "all") the interactive part of the RUnit tests is
#' run.  The user have to look at results and answer questions.
#' @param junitLibrary A character String specifying the path to a JUnit 4 jar
#' file to run the JUnit tests.  You can download it from
#' \url{http://www.junit.org/}.  Alternatively you can use the environment
#' variable GMCP_JUNIT_LIBRARY to specify the path.
#' @param outputPath During the RUnit tests files maybe produced at this
#' location.  If missing the current working directory is used if nothing else
#' is specified in the environment variable GMCP_UNIT_TEST_OPATH.
#' @return None of interest so far - the function prints the results to the
#' standard output.  (Perhaps in future versions a value will be returned that
#' can be processed by the GUI.)
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @examples
#' 
#' 
#' \dontrun{
#' unitTestsGMCP()
#' unitTestsGMCP(extended=TRUE, java=TRUE, interactive=TRUE, outputPath="~/RUnitTests")
#' 
#' }
#' 
#' 
#' @export unitTestsGMCP
unitTestsGMCP <- function(extended=FALSE, java=FALSE, interactive=FALSE, junitLibrary, outputPath) {
  # Adapted and extended from the code from http://rwiki.sciviews.org/doku.php?id=developers:runit
	if(!requireNamespace("RUnit", quietly=TRUE)) {
		stop("Please install package RUnit to run the unit tests.")
	}
	if (extended) Sys.setenv(GMCP_UNIT_TESTS=paste(Sys.getenv("GMCP_UNIT_TESTS"),"extended"), sep=",")
	if (interactive) Sys.setenv(GMCP_UNIT_TESTS=paste(Sys.getenv("GMCP_UNIT_TESTS"),"interactive"), sep=",")
	if (java) Sys.setenv(GMCP_UNIT_TESTS=paste(Sys.getenv("GMCP_UNIT_TESTS"),"java"), sep=",")
	if (missing(outputPath)) {
		if (Sys.getenv("GMCP_UNIT_TEST_OPATH")=="") {
			Sys.setenv(GMCP_UNIT_TEST_OPATH=getwd())
		}
	} else {
		Sys.setenv(GMCP_UNIT_TEST_OPATH=outputPath)
	}
	pkg <- "gMCP" 
	path <- system.file("unitTests", package=pkg)
	cat("\nRunning unit tests\n")
	print(list(pkg=pkg, getwd=getwd(), pathToUnitTests=path))
	
	library(package=pkg, character.only=TRUE)
	
	## If desired, load the name space to allow testing of private functions
	## if (is.element(pkg, loadedNamespaces()))
	##     attach(loadNamespace(pkg), name=paste("namespace", pkg, sep=":"), pos=3)
	##
	## or simply call PKG:::myPrivateFunction() in tests
	
	## --- Testing ---
  
	# Yes, these functions always exist, since we stopped if RUnit could not be required:
	defineTestSuite <- get("defineTestSuite", envir = loadNamespace("RUnit"))
	runTestSuite <- get("runTestSuite", envir = loadNamespace("RUnit"))
	printTextProtocol <- get("printTextProtocol", envir = loadNamespace("RUnit"))
	checkEquals <-  get("checkEquals", envir = loadNamespace("RUnit"))
	checkTrue <-  get("checkTrue", envir = loadNamespace("RUnit"))
	
	## Define tests
	testSuite <- defineTestSuite(name=paste(pkg, "unit testing"), dirs=path)
	
	## Run
	testRuns <- runTestSuite(testSuite)
	
	## Default report name
	pathReport <- file.path(path, "report")
	
	## Report to stdout and text files
	cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
	printTextProtocol(testRuns, showDetails=FALSE)
	printTextProtocol(testRuns, showDetails=FALSE,
			fileName=paste(pathReport, "Summary.txt", sep=""))
	printTextProtocol(testRuns, showDetails=TRUE,
			fileName=paste(pathReport, ".txt", sep=""))
	
	if (java || tests("java")) {
		# Test whether junit*.jar is in classpath
		if (!missing(junitLibrary)) {
			.jaddClassPath(junitLibrary)
		}
		if (Sys.getenv("GMCP_JUNIT_LIBRARY")!="") {
			.jaddClassPath(Sys.getenv("GMCP_JUNIT_LIBRARY"))
		}
		#testClass <- .jcall(.jnew("tests/RControlTest"), "Ljava/lang/Class;", method="getClass")
		testClasses <- .jcall(.jnew("tests/TestSuite"), "[Ljava/lang/Class;", method="getClasses", evalArray=FALSE)
		result <- try(.jcall("org.junit.runner.JUnitCore", "Lorg/junit/runner/Result;", method="runClasses", testClasses))
		if (("try-error" %in% class(result))) {
			cat("JUnit 4 is needed for JUnit tests (See http://www.junit.org/).")
			stop("Please specify the path to junit 4 jar file via junitLibrary.")
		}
		if (.jcall(result, "I", "getFailureCount")>0) {
			cat("------------------- JUNIT TEST SUMMARY --------------------\n\n")
			cat(.jcall(.jnew("tests/TestSuite"), "S", method="getResultString", result))
			stop(paste(.jcall(result, "I", "getFailureCount"),"failures in JUnit tests!"))
		} else {
			cat(.jcall(result, "I", "getRunCount"), " Java Unit Tests successful! (Runtime: ",.jcall(result, "J", "getRunTime")/1000," sec)")
		}
	}
}

# tests("java") gives back a logical specifying whether java tests should be run.
tests <- function(type) {
  return (any(type, "all") %in% strsplit(Sys.getenv("GMCP_UNIT_TESTS"),",")[[1]])  	
}

equals <- function(graph1, graph2, checkAttributes=FALSE, verbose=FALSE) {
	if (length(getNodes(graph1))!=length(getNodes(graph2))) {
		if (verbose) cat("Wrong number of hypotheses.\n")
		return(FALSE);
	}
	if (any(getNodes(graph1)!=getNodes(graph2))) {
		if (verbose) cat("Names of nodes differ.\n")
		return(FALSE);
	}
	if ("entangledMCP" %in% class(graph1) != "entangledMCP" %in% class(graph2)) {
		if (verbose) cat("Only one graph is of class entangledMCP.\n")
		return(FALSE);
	}
	# Call this function recursivly for entangled graphs.	
	if ("entangledMCP" %in% class(graph1)) {
		equal <- TRUE
		for(i in 1:length(graph1@subgraphs)) {
			if (verbose) cat("Checking subgraphs at position ",i,".\n")
			if (!equals(graph1@subgraphs[[i]], graph2@subgraphs[[i]])) {
				equal <- FALSE
			}
		}
		return(equal)
	}
	# Real function:
	if (any(graph1@m!=graph2@m)) {
		if (verbose) cat("Transition matrices differ.\n")
		return(FALSE);
	}
	if (any(graph1@weights!=graph2@weights)) {
		if (verbose) cat("Node weights differ.\n")
		return(FALSE);
	}
	return(TRUE)
	#TODO Implement checkAttributes=TRUE	
}