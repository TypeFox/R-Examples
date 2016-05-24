## Adapted and extended from the code from http://rwiki.sciviews.org/doku.php?id=developers:runit
unitTestsCrossover <- function(extended=FALSE, java=FALSE, interactive=FALSE, junitLibrary, outputPath) {
	if(!requireNamespace("RUnit", quietly=TRUE)) {
		stop("Please install package RUnit to run the unit tests.")
	}
	if (extended) Sys.setenv(CROSSOVER_UNIT_TESTS=paste(Sys.getenv("CROSSOVER_UNIT_TESTS"),"extended"))
	if (interactive) Sys.setenv(CROSSOVER_UNIT_TESTS=paste(Sys.getenv("CROSSOVER_UNIT_TESTS"),"interactive"))
	if (missing(outputPath)) {
		if (Sys.getenv("CROSSOVER_UNIT_TEST_OPATH")=="") {
			Sys.setenv(CROSSOVER_UNIT_TEST_OPATH=getwd())
		}
	} else {
		Sys.setenv(CROSSOVER_UNIT_TEST_OPATH=outputPath)
	}
	pkg <- "Crossover" 
	path <- system.file("unitTests", package=pkg)
	cat("\nRunning unit tests\n")
	print(list(pkg=pkg, getwd=getwd(), pathToUnitTests=path))
	
	library(package=pkg, character.only=TRUE)
	
	## --- Testing ---
	
	# Yes, these functions always exist, since we stopped if RUnit could not be required:
	defineTestSuite <- get("defineTestSuite")
	runTestSuite <- get("runTestSuite")
	printTextProtocol <- get("printTextProtocol")
  
	## Define tests
	testSuite <- defineTestSuite(name=paste(pkg, "unit testing"), dirs=path)
	
	## Run
	tests <- with(loadNamespace(pkg), runTestSuite(testSuite))
	
	## Default report name
	pathReport <- file.path(path, "report")
	
	## Report to stdout and text files
	cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
	printTextProtocol(tests, showDetails=FALSE)
	printTextProtocol(tests, showDetails=FALSE,
			fileName=paste(pathReport, "Summary.txt", sep=""))
	printTextProtocol(tests, showDetails=TRUE,
			fileName=paste(pathReport, ".txt", sep=""))
	
	if (java || "java" %in% strsplit(Sys.getenv("CROSSOVER_UNIT_TESTS"),",")[[1]]) {
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
