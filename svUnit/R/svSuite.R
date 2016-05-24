svSuite <- function (tests)
{
	## Check provided tests and build a 'svSuite' object
	tests <- as.character(tests)
	## Remove NAs and empty strings ("") from tests
	tests <- tests[!is.na(tests) & !(tests == "")]
	if (length(tests) > 0) {
		## Tests must be character strings like:
		## * package:PKG
		## * package:PKG (TESTSUITE)
		## * dir:MYDIR
		## * test(OBJ) where OBJ is any object with a 'test' attribute
		## * OBJ being a 'svTest' object (with non "exotic" name!),
		## Syntax is checked, but not existence/validity of corresponding tests!
		check1 <- (regexpr("^package:[a-zA-Z][a-zA-Z._0-9]*$", tests) > -1)
		check2 <- (regexpr("^package:[a-zA-Z][a-zA-Z._0-9]* *\\(.+\\)$", tests) > -1)
		check3 <- (regexpr("^dir:.+", tests) > -1)
		check4 <- (regexpr("^test\\(.+\\)$", tests) > -1)
		check5 <- (regexpr("^[a-zA-Z0-9_.]+$", tests) > -1)
		wrong <- ((check1 + check2 + check3 + check4 + check5) == 0)
		if (any(wrong))
			stop("Wrong 'tests' data: must be 'package:PKG', 'package:PKG (SUITE)',\n'dir:MYDIR', 'test(OBJ)' or 'OBJ'")
	}
	## This is a 'svSuite' object subclassing 'character'
	class(tests) <- c("svSuite", "character")
	return(tests)
}

as.svSuite <- function (x)
	return(svSuite(x))

is.svSuite <- function (x)
	return(inherits(x, "svSuite"))

print.svSuite <- function (x, ...)
{
	if (!is.svSuite(x))
		stop("'x' must be a 'svSuite' object")
	if (length(x) < 1) {
		cat("An empty svUnit test suite\n")
	} else {
		cat("A svUnit test suite definition with:\n")
		## Separate unit tests from tests embedded in objects
		isSuite <- regexpr("^[package:|dir:]", x) > -1
		if (any(isSuite)) {
			Suites <- x[isSuite]
			msg <- ifelse (length(Suites) == 1, "\n- Test suite:\n",
				"\n- Test suites:\n")
			cat(msg)
			print(Suites)
		}

		if (any(!isSuite)) {
			Objs <- x[!isSuite]
			msg <- ifelse (length(Objs) == 1, "\n- Test function:\n",
				"\n- Test functions:\n")
			cat(msg)
			print(Objs)
		}
	}
	return(invisible(x))
}

svSuiteList <- function (packages = TRUE, objects = TRUE, dirs = getOption("svUnit.dirs"),
	excludeList = getOption("svUnit.excludeList"), pos = .GlobalEnv,
	loadPackages = FALSE)
{
	## List unit test (1) in loaded packages (2) in objects in pos and (3) in
	## directories, possibly filtering them using an exclusion list
	## Note: Komodo should list test unit files in loaded projects too!
	if (length(packages) < 1)
		stop("'package' cannot have zero length")
	if (length(objects) < 1)
		stop("'objects' cannot have zero length")

	items <- character()

	## 1) Unit test files in loaded packages
	if (packages[1] != FALSE) {
		if (is.character(packages)) {	# We assume it is a list of packages
			Pkgs <- packages
		} else {	# We use the list of all loaded packages
			Pkgs <- .packages()
		}
		for (Pkg in Pkgs) {
			## Look for test units in the package
			path <- system.file(package = Pkg, "unitTests")
			if (path != "" && file.info(path)$isdir) {
				pkgname <- paste("package", Pkg, sep = ":")
				items <- c(items, pkgname)
				Files <- list.files(path = path, full.names = TRUE)
				for (File in Files) { # Add all subdirectories too
					if (file.info(File)$isdir)
						items <- c(items, paste(pkgname, " (", basename(File),
							")", sep = ""))
				}
			}
		}
	}

	## 2) Tests embedded in objects located in 'pos' environment
	if (objects[1] != FALSE) {
		envir = as.environment(pos)
		if (is.character(objects)) {
			tests <- character()
			for (Oname in objects) {
				if (exists(Oname, envir = envir, inherits = FALSE)) {
					Obj <- get(Oname, envir = envir, inherits = FALSE)
					if (is.svTest(Obj)) {
						tests <- c(tests, Oname)
					} else if (is.test(Obj)) {
						tests <- c(tests, paste("test(", Oname, ")", sep = ""))
					}
				}
			}
		} else {	# We list all objects in pos
			Objs <- mget(ls(envir = envir), envir = envir)
			Onames <- names(Objs)
			tests <- character()
			if (length(Objs) > 0) {
				for (i in 1:length(Objs)) {
					if (is.svTest(Objs[[i]])) {
						tests <- c(tests, Onames[i])
					} else if (is.test(Objs[[i]])) {
						tests <- c(tests, paste("test(", Onames[i], ")", sep = ""))
					}
				}
			}
		}
		items <- c(items, sort(tests))
	}

	## 3) Additional directories (check that they are valid and existing dirs)
	if (!is.null(dirs)) {
		## Check if each entry exists as a directory, exclude it if not
		## Prepend "dir:" to tag them as additional directories
		Dirs <- character()
		for (Dir in dirs)
			if (file.exists(Dir) && file.info(Dir)$isdir)
				Dirs <- c(Dirs, paste("dir", Dir, sep = ":"))
		items <- c(items, sort(Dirs))
	}

	## Filter the resulting list with 'excludeList'
	if (!is.null(excludeList)) {
		for (pattern in excludeList)
			items <- items[regexpr(pattern, items) == -1]
	}

	## Do we load the package?
	if (loadPackages) {
		## Get a list of packages we need for the suite
		Pkgs <- items[regexpr("^package:", items)]
		PkgsSrch <- unique(sub(" +\\(.+$", "", Pkgs))
		l <- length(PkgsSrch)
		if (l > 0) {
			PkgsName <- sub("^package:", "", PkgsSrch)
			Search <- search()
			for (i in 1:l) {
				if (!PkgsSrch[i] %in% Search) {
					res <- try(library(PkgsName[i], character.only = TRUE),
						silent = TRUE)
					if (inherits(res, "try-error"))
						warning("Cannot load package '", PkgsName[i], "'")
				}
			}
		}
	}

	## Make it a 'svSuite' object subclassing 'character'
	class(items) <- c("svSuite", "character")
	return(items)
}

makeUnit.svSuite <- function (x, name = make.names(deparse(substitute(x))),
dir = tempdir(), objfile = "", codeSetUp = NULL, codeTearDown = NULL,
pos = .GlobalEnv, ...)
{
	## Take an 'svSuite' object and make a unit from its function tests
	## that are not written yet in a test unit in a file
	## They are saved in a file named runit<name>.R in 'dir'
	if (!is.svSuite(x))
		stop("'x' must be a 'svSuite' object")
	name <- as.character(name)[1]
	## Under Windows, we transform \\ into /
	dir <- gsub("\\\\", "/", as.character(dir)[1])
	## Collect all items that are not 'package:...' or 'dir:...'
	isObj <- regexpr("^[package:|dir:]", x) == -1
	Objs <- sub("^test[(](.+)[)]$", "\\1", x[isObj])
	if (length(Objs) == 0) { # No objects, return NULL
		return(NULL)
	} else {	# Make a sourceable test unit file with tests collected in Objs
		Unit <- .prepareUnit(name, dir)
		.writeSetUp(unit = Unit, file = objfile, code = codeSetUp)
		.writeTearDown(unit = Unit, code = codeTearDown)
		for (objname in Objs)
			.writeTest(unit = Unit, objname = objname, pos = pos)
	}
	return(Unit)
}

runTest.svSuite <- function (x, name = make.names(deparse(substitute(x))),
unitname = NULL, ...)
{
	## Compile and run the test for this 'svSuite' object
	if (!is.svSuite(x))
		stop("'x' must be a 'svSuite' object")
	name <- as.character(name[1])

	## Decode tests contained in x
	tests <- as.character(x)
	dirs <- character()
	## Package suites...
	isPkg <- regexpr("^package:", tests) > -1
	if (any(isPkg)) {
		Pkgs <- tests[isPkg]
		Subdirs <- sub("^.+[(](.+)[)] *$", "\\1", Pkgs)
		Subdirs[Subdirs == Pkgs] <- ""
		Pkgs <- sub("^package:([^ ]+).*$", "\\1", Pkgs)
		for (i in 1:length(Pkgs)) {
			if (Subdirs[i] == "") {
			dir <- system.file(package = Pkgs[i], "unitTests")
		} else {
			dir <- system.file(package = Pkgs[i], "unitTests", Subdirs[i])
		}
			if (dir != "") dirs <- c(dirs, dir)
		}
	}

	## Add directories, and possibly make a temporary unit for test objects
	if (any(!isPkg)) {
		tests <- tests[!isPkg]
		## Directories
		isDir <- regexpr("^dir:", tests) > -1
		if (any(isDir))
			dirs <- c(sub("^dir:", "", tests[isDir]), dirs)
		## Objects
		if (any(!isDir)) {
			## Make a temporary unit for the tests of these objects
			if (!is.null(Unit <- makeUnit(x, name = name))) {
				## Add this path to dirs
				dirs <- c(dirname(Unit), dirs)
			}
		}
	}

	## Now, list all files in these dirs with name being runit*.R
	files <- character()
	for (dir in dirs)
		files <- c(files, list.files(dir, pattern = "^runit.+\\.[rR]$",
			full.names = TRUE))
	if (length(files) == 0) return(NULL)	# Nothing to run!	
	## Under Windows, transform all \\ into / in the file names
	files <- gsub("\\\\", "/", files)
	## Added by Thomas Wurtzler to control which unit test to run
	if (!is.null(unitname)) {
		unitname <- deparse(substitute(unitname))
		testNames <- gsub("^.*runit(.+)\\.[rR]$", "\\1", files)
		keep <- which(testNames == unitname)
		files <- files[keep]
		if (length(files) == 0) {
			warning("Test unit ", unitname, " not found")
			return(NULL)	
		}
	}
	## Run this test suite now, that is, source each file in .TestSuiteEnv
	## and run each testxxx function in it, using .setUp and .tearDown too
	## Record the list of tests
	.lastSuite <- list()
	for (file in files)
		.lastSuite[[basename(file)]] <- list(file = file)
	.Log <- Log()
	.Log$.lastSuite <- .lastSuite

	## Source each runit*.R file in turn
	for (unit in names(.lastSuite)) {
		## Create a new environment for this suite.
		.ThisTestSuiteEnv <- new.env(parent = .GlobalEnv)
        ## Store it in SciViews:TempEnv so that we can inspect it in case of
		## stop on error. But please do not remove the local alias.  #1327
        .assignTemp(".TestSuiteEnv", .ThisTestSuiteEnv)
		## Source the corresponding file
		Unit <- .lastSuite[[unit]]$file
		sys.source(Unit, envir = .ThisTestSuiteEnv)
		## Make sure there are .setUp() and .tearDown() functions
		if (!exists(".setUp", envir = .ThisTestSuiteEnv, mode = "function",
			inherits = FALSE))
			.ThisTestSuiteEnv$.setUp <- function() {}
		if (!exists(".tearDown", envir = .ThisTestSuiteEnv, mode = "function",
			inherits = FALSE))
			.ThisTestSuiteEnv$.tearDown <- function() {}
		## List all test files in the unit
		tests <- ls(.ThisTestSuiteEnv, pattern = "^test.+$")
		## Keep only 'test*' objects that are function
		keep <- unlist(lapply(tests, function(n) exists(n,
			envir = .ThisTestSuiteEnv, mode = "function", inherits = FALSE)))
		tests <- tests[keep]
		.Log$.lastSuite[[unit]]$tests <- tests
		## Run each test in turn
		for (test in tests) {
			.runTest(envir = .ThisTestSuiteEnv, test = test, unit = Unit)
		}
	}
	return(invisible(files))
}
