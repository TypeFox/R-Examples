# Unit tests utilities
# 
# Author: Renaud Gaujoux
# Creation: 25 Apr 2012
###############################################################################

#' @include utils.R
#' @include logging.R 
NULL

#' Load RUnit Compatible Package
#' 
#' Loads the package responsible for the implementation of the RUnit framework,
#' choosing amongst \sQuote{RUnitX}, \sQuote{svUnit} and \sQuote{RUnit}.
#' 
#' @param ... arguments passed to \code{\link{requirePackage}}.
#' 
#' @return nothing
#' @export
#' 
requireRUnit <- local({
			
	.cache <- NULL
	function(...){
		
		if( !is.null(.cache) ) return(.cache)
		
		has_pkg <- function(x) length(find.package(x, quiet=TRUE)) > 0L
		
		ruf <- c('RUnit', 'svUnit')
		runit <- NULL
		for( pkg in ruf){
			if( require.quiet(pkg, character.only=TRUE) ){
				runit <- pkg
				break
			}
		}
		
		if( is.null(runit) )
			stop("Cannot find any package providing RUnit framework.")
		message("Using RUnit framework provider: ", runit)
		
		.cache <<- runit
		# return name of the loaded framework 
		invisible(runit)
	}
	
})


# Borrowed from RUnit::.existsTestLogger
.existsTestLogger <- function(envir = .GlobalEnv){
    exists(".testLogger", envir = envir) && inherits(.testLogger, "TestLogger")
}

#' Enhancing RUnit Logger
#' 
#' Adds a function or a local variable to RUnit global logger.
#' 
#' @param name name of the function or variable to add 
#' @param value object to append to the logger.
#' If \code{value} is a function it is added to the list and is accessible via 
#' \code{.testLogger$name}.
#' If \code{value} is a variable it is added to the local environment and is 
#' therefore accessible in all logging functions.
#' @param logger an optional RUnit logger object. 
#' If missing or \code{NULL}, the object \code{.testLogger} is searched in  
#' \code{.GlobalEnv} -- and an error is thrown if it does not exist. 
#' 
#' @return the modified logger object. Note that the global object is also 
#' modified if \code{logger} is \code{NULL}.
#' 
addToLogger <- function(name, value, logger=NULL){
	
	
	logobj <- 
		if( !is.null(logger) ) logger
		else{
			if( !.existsTestLogger() )
				stop("No global logger exists")
			
			get('.testLogger', envir=.GlobalEnv)
		}
	
	# get local logger environment
	logenv <- environment(logobj$incrementCheckNum)

	if( is.function(value) ){# add function to logger
		if( is.null(logobj[[name]]) ){
			environment(value) <- logenv 
			logobj[[name]] <- value
			
			# update global logger if necessary
			if( is.null(logger) ){
				ge <- .GlobalEnv
				assign('.testLogger', logobj, envir=ge)
			}
		}
	}else{ # assign object in logger's local environment if not already there
		if( !exists(name, envir=logenv) )
			assign(name, value, envir=logenv)
	}
	
	# return modified logger object
	logobj
}

#' Plot in Unit Tests
#' 
#' Saves a plot in a PNG file that will be included in unit test HTML reports.
#' 
#' @param expr expression that generate th eplot
#' @param width plot width
#' @param height plot height (not used)  
#' @param msg plot msg explaining the plot . It will be used as the caption
#' 
#' @export
#' @keywords internal
checkPlot <- function(expr, msg=NULL, width=1000, height=NULL){
	
	# get stuff from RUnit
	uf <- requireRUnit()
	if( is.null(uf) || uf != 'RUnit' ) return(TRUE)
	#.existsTestLogger <- RUnit:::.existsTestLogger	
	.testLogger <- if( .existsTestLogger() ) .GlobalEnv$.testLogger
	
	if (missing(expr)) {
		stop("'expr' is missing.")
	}
	
	
	plotfile <- 
	if (.existsTestLogger()) {
		
		.testLogger$incrementCheckNum()
		
		if( is.null(.testLogger$setPlot) ){
			# add .plot list to logger environment
			addToLogger('.plots', NULL)
			
			
			# add function setPlot to logger
			.plots <- NULL # to trick R CMD check 
			addToLogger('setPlot', 
				function(name, msg=''){
					##@bdescr
					## add a plot to a test function.
					##@edescr
					##@in testFuncName : [character] name of test function
					##@in name : [character] filename
					##@in msg : [character] message string
					##@edescr	
					
					if( is.null(.plots) ) .plots <<- list()
					.plots[[name]] <<- msg
				}
			)
			
			# add function setPlot to logger
			.getTestData <- 
				.currentTestSuiteName <- 
				.currentSourceFileName <- 
				.getCheckNum <- NULL # not to get NOTES is R CMD check
			addToLogger('getPlotfile', 
				function(name, msg=''){
					
					td <- .getTestData()
					# TODO from test function name
					#fname <- tail(names(td[[.currentTestSuiteName]]$sourceFileResults[[.currentSourceFileName]]), 1L)
					fname <- basename(tempfile(paste(.currentTestSuiteName, '_', .currentSourceFileName, '_', sep='')))
					paste(fname, .getCheckNum(), sep='_')
					
				}
			)
			
			# update local object with modified global logger
			.testLogger <- .GlobalEnv$.testLogger
		}
		
		.testLogger$getPlotfile()
	}
	else tempfile(tmpdir='.')
	
	# add extension to plot file
	plotfile <- paste(plotfile, 'png', sep='.')
	
	# reset the msg if none was provided
	if( is.null(msg) ) msg <- plotfile

	#plot in the PNG file
	png(filename=plotfile, width=width)
	
	# evaluate the expression that generates the plot
	res <- try( eval(expr, envir = parent.frame()) )
	# close the graphic device
	dev.off()
	
	# test if everything went alright
	fileinfo <- file.info(plotfile)	
	if( inherits(res, "try-error") || is.na(fileinfo$size[1]) || fileinfo$size[1] == 0 ){
		#make sure the plot file is removed
		unlink(plotfile)
		
		if (.existsTestLogger()) {
			.testLogger$setFailure()
		}
		stop("Problem when generating plot:", res, msg)
	}
	
	if (.existsTestLogger()) {
		.testLogger$setPlot(plotfile, msg)
	}
	return(TRUE)
	
}

if( FALSE ){
	
	library(NMF, lib='build/lib')
	utest('pkg/inst/tests/runit.algorithms.r', fun='test.brunet', framework='RUnit')
	
}

#' Extra Check Functions for RUnit
#' 
#' \code{checkWarning} checks if a warning is generated by an expression, and 
#' optionally follows an expected regular expression pattern.
#' 
#' @param expr an R expression
#' @param expected expected value as regular expression pattern.
#' If a logical, then it specifies if a warning is expected or not.
#' 
#' For backward compatibility, a \code{NULL} value is equivalent to \code{TRUE}.
#' @param msg informative message to add to the error in case of failure
#' 
#' @export
#' @rdname uchecks
#' 
#' @examples 
#' 
#' # check warnings
#' checkWarning({ warning('ah ah'); 3})
#' checkWarning({ warning('ah oh ah'); 3}, 'oh')
#' try( checkWarning(3) )
#' try( checkWarning({ warning('ah ah'); 3}, 'warn you') )
#' 
checkWarning <- function(expr, expected=TRUE, msg=NULL){
	
	# get stuff from RUnit
	uf <- requireRUnit()
	#.existsTestLogger <- RUnit:::.existsTestLogger	
	.testLogger <- if( .existsTestLogger() ) .GlobalEnv$.testLogger
	
	if (missing(expr)) {
		stop("'expr' is missing")
	}
#	if (is.null(silent)) {
#		silent <- FALSE
#		warning("'silent' has to be of type 'logical'. Was NULL. Set to FALSE.")
#	}
#	
	if (.existsTestLogger()) {
		.testLogger$incrementCheckNum()
	}
	
	pf <- parent.frame()
	warns <- NULL
	withCallingHandlers(eval(expr, envir = pf)
		, warning = function(w){
			warns <<- c(warns, w$message)
		}
	)
	
	# check that some warning was thrown
	if( length(warns) == 0L ){
        if( isFALSE(expected) ) return( TRUE )
		if (.existsTestLogger()) {
			.testLogger$setFailure()
		}
		stop("Warning not generated as expected\n", msg)
	}
	if( isFALSE(expected) ){
        if (.existsTestLogger()) {
			.testLogger$setFailure()
		}
		stop("Warning generated while none was expected:\n"
            , "  - Warning(s): ", if(length(warns)>1)"\n    * ",  str_out(warns, Inf, sep="\n    * ") ,"\n"
            , msg)
    }
	# check warnings
	if( is.null(expected) || isTRUE(expected) ) return(TRUE)
	if( any(grepl(expected, warns)) ) return(TRUE)
	
	# throw error
	if (.existsTestLogger()) {
		.testLogger$setFailure()
	}
	stop("Warning does not match expected pattern:\n"
		, "  - Warning(s): ", if(length(warns)>1)"\n    * ",  str_out(warns, Inf, sep="\n    * ") ,"\n"
		, "  - Pattern: '", expected,"'\n"
		, msg)
	
	TRUE
}

#' Make Vignette for Unit Tests
#' 
#' Builds a vignette for unit tests in a package using the \code{\link{utest}} 
#' and a template vignette file. 
#' 
#' @param pkg Package name
#' @param file Output file (.Rnw, .tex, or .pdf)
#' @param ... extra arguments passed to \code{\link{utest}}.
#' @param check logical that indactes the cal was made from R CMD check, in which case the vignette
#' is updated only if results of unit tests can be found in the unit test output directory, 
#' where they would have been generated by \code{\link{utest}}.
#' 
#' @return Result of running unit test suite
#'  
#' @export
#' 
makeUnitVignette <- function(pkg, file=paste(pkg, '-unitTests.pdf', sep=''), ..., check=FALSE){
	
	package <- pkg
	pkg <- sub("^package:", "", pkg)
	# generate the vignette for unit test on exit
	if( !is.null(file) )
		on.exit( writeUnitVignette(pkg, file, check=check) )
	# load this package
	if( !require(pkg, character.only = TRUE ) ){
		stop("Could not load package '", pkg, "' for testing [libPath= ", str_out(.libPaths(), Inf), "]")
	}
	
	# run unit tests if not check or if the test results are not there (e.g., R CMD build)
#	if( userIs('renaud') ){
#		env <- str_trim(capture.output(system('env', intern=TRUE)))
#		if( check )	write(env, file="~/check_env.txt")
#		else write(env, file="~/make_env.txt")
#	}
	
#	if( !check || !is.dir(utestPath(package=package)) ){
	if( !check ){
		
		# force running all tests 
		utestCheckMode(FALSE)
		
		# run unit tests
		tests <- utest(package, ...)
		
		# check for errors
		err <- getErrors(tests)
		errMsg <- NULL
		if( err$nFail > 0) {
			errMsg <- c(errMsg, sprintf( "unit test problems: %d failures\n", err$nFail))
		}
		if( err$nErr > 0) {
			errMsg <- c(errMsg, sprintf( "unit test problems: %d errors\n", err$nErr))
		}
		# stop if any failure or error occured
		if( length(errMsg) > 0L )
			stop(errMsg)
		
		# return result of unit test suite
		err
	}else{
		# do nothing: tests should have been already run by R CMD check
	}
	
}

#' Writes Unit Tests Vignette 
#' 
#' Writes a vignette that contains the results from running unit test suites.
#' 
#' @param pkg Package name
#' @param file Output Sweave (.Rnw) file
#' @param results result file or output character vector
#' @param check logical that indactes the cal was made from R CMD check, 
#' in which case the vignette is updated only if results of unit tests can 
#' be found in the unit test output directory, where they would have been 
#' generated by \code{\link{utest}}.
#' 
#' @export
#' 
writeUnitVignette <- function(pkg, file, results=NULL, check=FALSE){
	Rnw.template <- 
"
\\documentclass[10pt]{article}
%\\VignetteDepends{knitr}
%\\VignetteIndexEntry{@pkg@-unitTests}
%\\VignetteCompiler{knitr}
%\\VignetteEngine{knitr::knitr}
\\usepackage{vmargin}
\\setmargrb{0.75in}{0.75in}{0.75in}{0.75in}

<<setup, include=FALSE>>=
pkg <- '@pkg@'
require( pkg, character.only=TRUE )
prettyVersion <- packageDescription(pkg)$Version
prettyDate <- format(Sys.Date(), '%B %e, %Y')
authors <- packageDescription(pkg)$Author
@

\\usepackage[colorlinks]{hyperref}
\\author{\\Sexpr{authors}}
\\title{\\texttt{\\Sexpr{pkg}}: Unit testing results@resNote@}
\\date{\\texttt{\\Sexpr{pkg}} version \\Sexpr{prettyVersion} as of \\Sexpr{prettyDate}}
\\begin{document}
\\maketitle

@results@

\\section*{Session Information}
@sessionInfo@

\\end{document}
"
	verbatim_wrap <- function(...){
		c("\\\\begin{verbatim}\n", ..., "\n\\\\end{verbatim}")
	}
	# default is to load the unit test results from the global output directory
	if( is.null(results) ){
		upath <- utestPath(package=pkg)
		results <- list.files(upath, pattern="\\.txt$", full.names=TRUE)
		if( !length(results) ){
			results <- verbatim_wrap('Could not find any unit test result in "', upath, '"')
		}
	}
	
	if( is.file(results[1L]) ){
		resFile <- results[1L]
		name <- str_match(resFile, "([^.]+)\\.[^.]+$")[,2L]
		results <- c(str_c("\\\\section{", name, "}"), verbatim_wrap(readLines(resFile)))
	}else{
		resFile <- NULL
	}
	results <- paste(results, collapse="\n")
	
	# substitute template variables
	contents <- Rnw.template
	# package name
	contents <-	gsub("@pkg@", pkg, contents)
	# unit test results
	contents <-	gsub("@results@", results, contents)
	# session info (as when calling this function)
	contents <-	gsub("@sessionInfo@", gsub("\\", "\\\\", paste(toLatex(sessionInfo()), collapse="\n"), fixed=TRUE), contents)
	# note on how tests were performed
	resnote <- str_c("\\footnote{Vignette computed ", if( check ) ' via R CMD check/build ', ' on ', date(),"}")
	if( check ){ 
		# add path to included file if compiled from R CMD check (for debug purposes)
		lfile <- gsub("([_$])", "\\\\\\1", paste(resFile, collapse="\\\\"))
		resnote <- str_c(resnote, " \\footnote{File: '", lfile, "'}")
	}
	contents <-	gsub("@resNote@", gsub("\\", "\\\\", resnote, fixed=TRUE), contents)
	
	fileext <- toupper(file_extension(file))
	fileext <- charmatch(fileext, c('RNW', 'TEX', 'PDF'))
	if( is_NA(fileext) )
		stop("Invalid output file extension [",fileext,"] from file '", file, "'")
	
	fileRNW <- if( fileext == 1L ) file else str_c(pkg, '-unitTests.Rnw')
	fileTEX <- if( fileext == 2L ) file else str_c(pkg, '-unitTests.tex')
	filePDF <- if( fileext == 3L ) file else str_c(pkg, '-unitTests.pdf')
	
	# write into Rnw file
	writeLines(contents, fileRNW)	
	if( fileext == 1L ) return()
	
	# compile vignette
	rnw(fileRNW, fileTEX)
	if( fileext == 2L ) return()
	
	# Run texi2dvi tex file
	res <- tools::texi2dvi(fileTEX, pdf = TRUE, clean = TRUE )
	
	# copy file in main check directory
	if( check )	file.copy(filePDF, '../../..')
	res
}

# Unit test frameworks data
.UFdata <- list(
	RUnit = list(
		file_pattern="^runit.*\\.[rR]$"
		, fun_pattern="^test\\."
		, check_pattern = "^check.+"
		, check_functions = c(
				'checkTrue'
				, 'checkIdentical'
				, 'checkEquals'
				, 'checkEqualsNumeric'
				, 'checkException'
		)
	)
	, testthat = list(
		file_pattern="^test.*\\.[rR]$"
		, check_pattern = "^(expect_.+)|(test_that$)" 
		, check_functions = c(
				"test_that"
				, "expect_equal"
				, "expect_equivalent"
				, "expect_error"
				, "expect_false"
				, "expect_identical"
				, "expect_is"
				, "expect_match"
				, "expect_message"
				, "expect_output"
				, "expect_that"
				, "expect_true"
				, "expect_warning"   
		)
	)
)

#' Inferring Unit Test Framework
#' 
#' @param x an filename, a function or the body of a function
#' @param eval a logical that indicates if the value of \code{x} should be used.
#' 
#' @return the name of the framework as a character string or NULL if
#' it could not be detected.
#' 
#' @import codetools
#' @export
utestFramework <- function(x, eval=FALSE){
	
	# check if one should detect within an expression
	expr <- if( missing(eval) || !eval ) substitute(x) 
			else if( is.function(x) ) body(x)
	
	# walk code using codetools looking up for known test functions
	if( !is.null(expr) ){
		cw <- makeCodeWalker(leaf= function(e, w) if( is.symbol(e) ) cat(e, "\n"))
		s <- str_trim(capture.output(walkCode(expr, cw)))
		if( length(s) > 1L ){
			for( f in names(.UFdata) ){
				if( any(s %in% .UFdata[[f]]$check_functions) ){
					return(f)
				}
			}
		}
		# not found without evaluating
		if( !missing(eval) && !eval ) return()
		if( missing(eval) ){ # try evaluating
			return(utestFramework(x, eval=TRUE))
		}
	}
	
	if( !is.character(x) )
		stop("Invalid argument `x`: expecting a character string")
	path <- x
	framework <- NULL
	tf <- if( is.dir(path) ) list.files(path, "\\.[rR]$") else path
	for( f in names(.UFdata) ){
		if( any(grepl(.UFdata[[f]]$file_pattern, tf)) ){
			return(f)
		}
	}
	
	if( is.null(framework) )
		stop("Could not determine unit test framework used in directory: '", path, "'")
	framework
}

#' Embedded Unit Tests
#' 
#' The function \code{unit.test} provides a way to write unit tests embedded within
#' package source files.
#' These tests are stored and organised in the package namespace, and can be run using 
#' the unified interface provided by the function \code{link{utest}}.
#' Both Runit and testthat tests are supported -- and automatically detected.
#' 
#' 
#' @param x single character string used as test identifier/label
#' @param expr expression containing the actual test commands.
#' It is not evaluated, but only stored in the package namespace.
#' @param framework Unit test framework
#' @param envir the definition environment of object \code{x}.
#' 
#' @return a test function with no arguments that wrapping around \code{expr} 
#' 
#' @import digest
#' @export
#' 
unit.test <- function(x, expr, framework=NULL, envir=parent.frame()){
	
	sid <- as.character(deparse(substitute(x)))	
	hash <- suppressWarnings(digest(x))
	# get test environment
	eTest <- packageTestEnv()
	# wrap test into a function
	f <- function(){}
	environment(f) <- eTest
	body(f) <- substitute({expr})
	
	if( !grepl('"', sid) )
	{
		lmessage('Creating unit test for object: `', sid, '`')
		eval(substitute(attr(x, 'test') <- f, list(x=substitute(x), f=f)), envir)
	}else
		lmessage('Creating unit test: ', sid)
	
	# add the test to the package test environment
	eTest[[str_c(sid, ':', hash)]] <- list(test=f, name=sid, object=is.name(x))
	# return the test function
	f
}

#' Returns the package internal environment where unit tests are stored.
#' 
#' @param pkg package name.
#' If missing the caller's package is assumed. 
#' 
#' @export
packageTestEnv <- function(pkg){
	
	if( !missing(pkg) && !is.null(pkg) ){
		e <- packageEnv(pkg)
		return( e$.packageTest )
	}
	
	e <- packageEnv()
	# create test environment if necessary
	if( is.null(e$.packageTest) )
		e$.packageTest <- new.env(parent=e)
	e$.packageTest
}


list.tests <- function(x, pattern=NULL){
	
}

#unit.test(packageEnv, {print('test for packageEnv')})
#unit.test('lmlm', {print('test for something else')})

#utest <- function(x, ..., framework="RUnit", PACKAGE=NULL){
#		
#	if( missing(x) )
#		x <- packagePath('unitTests', PACKAGE=PACKAGE)
#	else if( class(x)[1] != 'character')
#		return( UseMethod('utest', x) )
#	
#	if( is.null(framework) ){
#		stop("Not implemented")
#	}else{
#		# change directory to run tests
#		owd <- setwd(x)
#		on.exit(setwd(owd))
#		# run tests under selected framework
#		class(x) <- framework
#		utest(x, ..., PACKAGE=PACKAGE)
#		# output test result
#	}
#}

#' Running Unit Tests
#' 
#' Run unit tests in a variety of settings.
#' This is still \strong{very} experimental.
#' 
#' @param x object to which a unit test is attached
#' @param ... extra arguments to allow extensions and are passed to 
#' the unit framework running funcitons. 
#'
#' @inline
#' @export
setGeneric('utest', function(x, ...) standardGeneric('utest'))
#' Run the unit test assoicated to a function. 
#' 
#' @param run a logical that indicates if the unit test should be run
setMethod('utest', 'function',
	function(x, run = TRUE){
		# get actual name of the function
		sid <- as.character(deparse(substitute(x, parent.frame())))
		# remove leading namespace specifications
		sid <- sub("^[^:]+:::?", "", sid)
		# get the package's  
		pkg <- attr(x, 'package')
		eTest <- packageTestEnv(pkg)
		if( is.null(eTest) ) return()
		tfun <- ls(eTest, pattern=str_c("^", sid, ":"))		
	}
)
#' Run a package test suite
#' 
#' @param filter pattern to match files that contain the definition of 
#' the unit tests functions to run.
#' @param fun patter to match the test functions to run.
#' @param testdir directory where to look for the test files
#' @param framework unit test framework
#' @param quiet a logical that indicates if the tests should be run silently
#' @param lib.loc path to a library where installed packages are searched for.
#' Used is of the form \code{x='package:*'}.
#'  
setMethod('utest', 'character', 
		function(x, filter="^runit.+\\.[rR]$", fun="^test\\.", ...
				, testdir='tests', framework=c('RUnit', 'testthat')
				, quiet = Sys.getenv("RCMDCHECK") != "FALSE"
				, lib.loc = NULL){
			
			cat("#########################\n")
			
			# define environment variable that identifies a test run (if not already defined) 
			if( is.na(utestCheckMode(raw = TRUE)) ){
				oldUM <- utestCheckMode(TRUE)
				on.exit( utestCheckMode(oldUM), add=TRUE)
			}
			
			#print(system('env'))
			# detect type of input string
			path <- 
					if( grepl("^package:", x) ){# installed package
						pkg <- sub("^package:", "", x)
						if( is.null(path <- path.package(pkg, quiet=TRUE)) ){
							library(pkg, character.only=TRUE, lib.loc=lib.loc)
							path <- path.package(pkg)
						}
						file.path(path, testdir)
					}else{
						# try to find a corresponding development package
						if( require.quiet(devtools) 
								&& is.package(pkg <- as.package(x, quiet=TRUE)) ){
							load_all(pkg, TRUE)
							file.path(pkg$path, 'inst', testdir)
						}else{ # assume x is a path  
							x
						}
					}
			
			# check that the path exists
			if( !file.exists(path) ){
				if( !hasArg(testdir) ){ # try another default
					opath <- path
					path <- file.path(dirname(path), 'unitTests')
					if( !file.exists(path) )
						stop("Could not find any default unit test directory ['", opath, "' nor '", path, "'].")
				} else {
					stop("Unit test directory '", path, "' does not exist")
				}
			}
			
			message("Running unit tests in: '", path, "'")
			# detect unit test framework: RUnit or testthat?
			framework <- 
					if( missing(framework) ) utestFramework(path)
					else match.arg(framework)
			message("Using unit test framework: ", framework)
			
			# load default patterns
			up <- .UFdata[[framework]]
			if( missing(filter) ) filter <- up$file_pattern
			if( missing(fun) ) fun <- up$fun_pattern
			
			# run tests
			path <- normalizePath(path)
			# remove/create output directory
			opath <- utestPath(package=x)
			if( file.exists( opath ) ){
				unlink(opath, recursive=TRUE) 
			}
			dir.create(opath, recursive=TRUE)
			# copy results in working directory on exit
			on.exit(
				{ if( file.exists(opath) )
					file.copy(opath, '.', recursive=TRUE)
				}
			, add=TRUE)
			#
			
			if( is.dir(path) ){ # all tests in a directory
				if( framework == 'RUnit' ){ # RUnit
					
					requireRUnit("Running RUnit test suites")
					s <- defineTestSuite(x, path
							, testFileRegexp=filter
							, testFuncRegexp=fun, ...)
					str(s)
					utest(s, quiet=quiet, outdir=opath)
					
				}else if( framework == 'testthat' ){ # testthat
					
					requirePackage('testthat', "Running testthat unit test suites")
					test_dir(path, filter=filter, ...)
					
				}
			}else{ # single test file
				if( framework == 'RUnit' ){ # RUnit
					
					requireRUnit("Running RUnit unit test file")
					runTestFile(path, testFuncRegexp=fun, ...)
					
				}else if( framework == 'testthat' ){ # testthat
					
					requirePackage('testthat', "Running testthat unit test file")
					test_file(path, ...)
					
				}
			}
			
		}
)

setOldClass('RUnitTestSuite')
#' Runs a RUnit test suite
#' 
#' @param outdir output directory
setMethod('utest', 'RUnitTestSuite',
	function(x, ..., quiet=FALSE, outdir=NULL){
		requireRUnit("Running RUnit test suites")
		
		pathReport <- file.path(outdir, str_c("utest.", sub("[:]", "_", x$name)))
		
		t <- system.time({
			if( quiet ){
				suppressWarnings(suppressMessages(out <- capture.output(
					tests <- runTestSuite(x, ...)
				)))
			}else 
				tests <- runTestSuite(x, ...)
		})
		
		## Report to stdout and text files
		cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
		summary_file <- paste(pathReport, ".Summary.txt", sep="")
		printTextProtocol(tests, showDetails=FALSE,	fileName=summary_file)
		# append timing
		st <- c("\nTotal execution time\n***************************"
				, paste(capture.output(print(t)), collapse="\n"))
		write(st, file=summary_file, append=TRUE)
		# detailed report
		details_file <- paste(pathReport, ".Details.txt", sep="")
		printTextProtocol(tests, showDetails=TRUE, fileName=details_file)
		write(st, file=details_file, append=TRUE)
		#
		
		## Report to HTML file
		printHTMLProtocol(tests, fileName=paste(pathReport, ".html", sep=""))
		
		## Return stop() to cause R CMD check stop in case of
		##  - failures i.e. FALSE to unit tests or
		##  - errors i.e. R errors
        tmp <- getErrors(tests)
        if(tmp$nFail > 0 | tmp$nErr > 0) {
            # filter failures and errors
            test_summary <- capture.output(.summaryRUnit(tests, c('error', 'failure'), print = TRUE))
            stop(paste0(test_summary, collapse = "\n"))
#			stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
#							", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
		}
		
		tests
	}
)

.summaryRUnit <- function(testData, type = c('error', 'failure', 'deactivated'), print = FALSE){
    
    # taken from printTextProtocol
    res <- testData
    for (tsName in names(testData)) {
        
        if( print ){
            cat("Functions:", testData[[tsName]]$nTestFunc, "|"
                , "Errors:", testData[[tsName]]$nErr, "|"
                , "Failures:", testData[[tsName]]$nFail, "\n")
        }
        srcFileRes <- testData[[tsName]][["sourceFileResults"]]
        
        for (i in seq_along(srcFileRes)) {
            fname <- basename(names(srcFileRes)[i])
            testFuncNames <- names(srcFileRes[[i]])
            keep <- integer()
            for (j in seq_along(testFuncNames)) {
                  funcList <- srcFileRes[[i]][[testFuncNames[j]]]
                  if (funcList$kind %in% type){
                      keep <- c(keep, j)
                      if( print ){
                          if (funcList$kind == "error") {
                                cat("ERROR in ", fname, "::", testFuncNames[j], ": ", funcList$msg, sep = "")
                          }
                          else if (funcList$kind == "failure") {
                                cat("FAILURE in ", fname, "::", testFuncNames[j], ": ", funcList$msg, sep = "")
                          }
                          else if (funcList$kind == "deactivated") {
                                cat("DEACTIVATED ", fname, "::", testFuncNames[j], ": ", funcList$msg, sep = "")
                          }
                      }
                }
            }
            if( length(keep) ) res[[tsName]][["sourceFileResults"]][[i]] <- srcFileRes[[i]][keep]
            else res[[tsName]][["sourceFileResults"]] <- res[[tsName]][["sourceFileResults"]][-i]
        }
    }
    invisible(res)
}

#' Unit Tests Result Directory
#' 
#' Returns the path to the directory where the results of unit tests are stored.
#' This path is used by \code{\link{utest}} to save unit test results, which are 
#' read by \code{\link{makeUnitVignette}} to update the unit test vignette when 
#' runnning R CMD check.  
#' 
#' @param ... extra arguments passed to \code{\link{packagePath}}, e.g., \code{package}.
#' 
#' @export
utestPath <- function(...){
	packagePath('tests-results', ...)
}

#uTest <- function(file, fun, ...){
#	
#	library(RUnit)
#	tdir <- packagePath('unitTests')
#	ufiles <- list.files(tdir)
#	
#	get.tfile <- function(file){
#		i <- grep(paste(file,"(\\.[rR])?",sep=''), ufiles)
#		if( length(i) > 0L ) ufiles[i[1L]]
#		else NULL
#	}
#	
#	tfile <- file
#	if( is.null(tfile <- get.tfile(tfile)) ){
#		tfile <- paste('runit.', file, sep='')
#		if( is.null(tfile <- get.tfile(tfile)) ){
#			tfile <- paste('testthat.', file, sep='')
#			if( is.null(tfile <- get.tfile(tfile)) )
#				stop("Could not find test file '", file, "' (nor runit.% or testthat.% versions) in '", tdir, "'")
#		}
#	}
#	tfile <- file.path(tdir, tfile)
#	
#	if( !missing(fun) ){
#		e <- new.env()
#		source(tfile, local=e)
#		tfun <- fun
#		if( !exists(tfun, e, inherits=FALSE) ){
#			tfun <- paste('test.', fun, sep='')
#			if( !exists(tfun, e, , inherits=FALSE) )
#				stop("Could not find test function '", fun, "' (not test.% version) in '", tfile, "'")
#		}
#		tfun <- gsub(".", "\\.", tfun, fixed=TRUE)
#		runTestFile(tfile, testFuncRegexp=str_c("^", tfun, "$"), ...)
#	}else 
#		runTestFile(tfile, ...)
#}

