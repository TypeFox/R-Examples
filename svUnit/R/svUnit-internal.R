.onLoad <- function (lib, pkg)
{
	## The default exclusion list, if it is not defined yet
	## Although there are unit tests defined in these packages (as examples),
	## we don't want to include them, by default, in our test suite!
	if (is.null(getOption("svUnit.excludeList")))
		options(svUnit.excludeList = c("package:sv", "package:RUnit"))
	## Look if the SciViews-K Unit Komodo extension is installed
	## TODO: this causes more problems than solutions => temporarily deactivated!
	#.installUpgradeKomodoExtension()
	## Install a callback to update the list of units automatically in the GUI
	## Use the mechanism introduced in svSocket 0.9-48 to allow execution of
	## this task callback from a socket client too
	h <- .getTemp(".svTaskCallbackManager", default = NULL)
	if (!is.null(h))
		h$add(guiSuiteAutoList, name = "guiSuiteAutoList")
}

.onUnload <- function (libpath)
{
	## Delete the task callback
	h <- .getTemp(".svTaskCallbackManager", default = NULL)
	if (!is.null(h))
		h$remove("guiSuiteAutoList")
	## Clear the list of units in the GUI client
	if (exists("koCmd", mode = "function"))
		get("koCmd")('sv.r.unit.getRUnitList_Callback("");')
}

.packageName <- "svUnit"

.komodoExtensionMinVersion <- "0.7.3"

.installUpgradeKomodoExtension <- function ()
{
	if (!exists("koCmd", mode = "function")) return()
	## Look if the SciViews-K Unit Komodo extension is installed and is of the
	## right version. Otherwise, propose to install, or update it
	xpiFile <- system.file("komodo", "sciviewskunit-ko.xpi", package = "svUnit")
	## Bug: sometimes this fails, preventing svUnit to load, despite it would
	## work well past this point. So, I put this in a try() to silently catch
	## the error and continue loading svUnit anyway (thanks Claudia Beleites)
	koVersion <- try(
		get("koCmd")('sv.socket.serverWrite(sv.r.unit.version + "." + sv.r.unit.release);'),
		silent = TRUE)
	if (inherits(koVersion, "try-error")) return()
	if (koVersion == "undefined.undefined") {
		## We need to install the extension
		cmd <- 'var res = ko.dialogs.okCancel("The SciViews-K Unit extension is required by svUnit",'
		cmd <- paste(cmd, '"OK", "Do you want to install the SciViews-K Unit extension now in Komodo?')
		cmd <- paste(cmd, 'You will be prompted for confirmation (click \'Install Now\')')
		cmd <- paste(cmd, 'and will have to restart Komodo at the end of the installation')
		cmd <- paste(cmd, '(click \'Restart Komodo\').", "svUnit");')
		cmd <- paste(cmd, ' if (res == "OK") { ko.open.URI("<<<data>>>"); }', sep = "")
		get("koCmd")(cmd, data = xpiFile)
	} else if (.compareVersion(koVersion, .komodoExtensionMinVersion) == -1) {
		## We need to upgrade the extension
		cmd <- 'var res = ko.dialogs.okCancel("A newer SciViews-K Unit extension is required by svUnit",'
		cmd <- paste(cmd, '"OK", "Your SciViews-K Unit extension is too old for this version of svUnit.')
		cmd <- paste(cmd, 'Do you want to upgrade it now?')
		cmd <- paste(cmd, 'You will be prompted for confirmation (click \'Install Now\')')
		cmd <- paste(cmd, 'and will have to restart Komodo at the end of the installation')
		cmd <- paste(cmd, '(click \'Restart Komodo\').", "svUnit");')
		cmd <- paste(cmd, ' if (res == "OK") { ko.open.URI("<<<data>>>"); }', sep = "")
		get("koCmd")(cmd, data = xpiFile)
	}
}

.compareVersion <- function (a, b)
{
    ## This is the same as utils::compareVersion(), but we don't want
	## dependencies on utils
	if (is.na(a))
        return(-1)
    if (is.na(b))
        return(1)
    a <- as.integer(strsplit(a, "[\\.-]")[[1]])
    b <- as.integer(strsplit(b, "[\\.-]")[[1]])
    for (k in 1:length(a)) {
        if (k <= length(b)) {
            if (a[k] > b[k])
                return(1)
            else if (a[k] < b[k])
                return(-1)
        }
        else {
            return(1)
        }
    }
    if (length(b) > length(a))
        return(-1)
    else return(0)
}

.kindLevels <- c("OK", "**FAILS**", "**ERROR**", "DEACTIVATED")

.kind <- function (val = TRUE)
{
    ## TRUE or 1 -> 1 = "OK"
    ## FALSE or 0 -> 2 = "**FAILS**"
    ## -1 -> 3 = "**ERROR**"
    ## -2 -> 4 = "DEACTIVATED"
    factor(.kindLevels[-(as.integer(val) - 2)], levels = .kindLevels)
}

.kindMax <- function (kinds)
{
    ## If there are no record, must be because all tests succeed!
    if (length(kinds) == 0) return(.kind(TRUE))
    Kinds <- as.numeric(kinds)
    if (sum(Kinds, na.rm = TRUE) == 0) return(.kind(NA))
    factor(.kindLevels[max(as.numeric(kinds), na.rm = TRUE)],
        levels = .kindLevels)
}

.formatTime <- function (x, secDigits = 0, minSec = 10^-secDigits,
prepend = " run in")
{
	## x is given in seconds, and it returns a pretty formatted string with time
	if (is.null(x) || is.na(x)[1]) return("")
	x <- as.numeric(x)
	Sec <- round(x %% 60, digits = secDigits)
	Min <- x %/% 60
	Hour <- Min %/% 60
	Min <- Min %% 60
	Time <- prepend
	Time <- ifelse (Hour > 0, paste(Time, Hour, "h"), Time)
	Time <- ifelse (Min > 0 | Time != prepend, paste(Time, Min, "min"), Time)
	Time <- ifelse (Sec > minSec | Time != prepend, paste(Time, Sec, "sec"), Time)
	Time <- ifelse (Time == prepend, paste(prepend, "less than", minSec, "sec"), Time)
	Time[is.na(Time)] <- ""
	return(Time)
}
## Test: .formatTime((0:10) * 400 + 0.56)

.formatResult <- function (result, level = getOption("svUnit.strLevel"))
{
	if (is.null(level)) level <- 1 else level <- as.integer(level[1])
	if (level < 1) return("")	# Return an empty string
	## Capture the report returned by the str() function
	capture.str <- function(data, level) {
		rval <- NULL
		file <- textConnection("rval", "w", local = TRUE)
		sink(file, type = "output")
		sink(file, type = "message")
		on.exit({
			sink(type = "output")
			sink(type = "message")
			close(file)
		})
		str(data, max.level = level)
        cat("\n")
		return(rval)
	}
	Str <- capture.str(result, level)
	return(paste(Str, collapse = "\n"))
}

.logTest <- function (timing, test, msg = "", description = NULL)
{
    .Log <- Log(description = description)
    ## Determine the name of the test
    if (missing(test)) {    # Is it defined globally?
        if (exists("..Test", envir = .Log, inherits = FALSE)) {
            test <- .Log$..Test
        } else {            # Try to guess it from the call
            ret <- try(test <- as.character(sys.call(-2))[1], silent = TRUE)
            if (inherits(ret, "try-error") || is.na(test)) {
                ## check...() probably called directly at the command line
                test <- "eval"
				## Convenient for collecting these tests together with tests run
				## inside runit*.R files (not embedded in test functions)
            } else if (test == "runTest") {
                ## Special case for runTest(myTest) or runTest(test(foo))
                test <- as.character(sys.call(-2))[2]
            } else if (test == "eval.with.vis") {
				test <- "eval"
			}
        }
    }
    ## Do we need to create 'test'?
    if (!exists(test, envir = .Log, inherits = FALSE)) {
        if (msg == "") msg <- .Log$..Msg
        .Log[[test]] <- structure(
            data.frame(
                msg = character(),
                call = character(),
                timing = numeric(),
                kind = .kind(logical()),
                res = character(),
                obj = character(),
                file = character(),
                tag = character(),
                stringsAsFactors = FALSE),
            time = Sys.time(),
            stats = c(tests = 1, timing = timing),
            context = c(unit = as.character(.Log$..Unit), test = test,
                msg = paste(msg, collapse = "\n")),
            class = c("svTestData", "data.frame"))
    } else {
        ## Just update stats
        attr(.Log[[test]], "stats") <-
			attr(.Log[[test]], "stats") + c(1, timing)
    }
    return(test)
}

.logTestData <- function (test, msg, call, timing, val, kind = .kind(val), res,
obj = .Log$..Obj, file = .Log$..File, tag = .Log$..Tag,
printTest = getOption("svUnit.printTest"))
{
    ## Add these data to .lastTest
    .Log$.lastTest <- structure(data.frame(
        msg = msg, call = call, timing = timing, kind = kind, res = res,
        obj = obj, file = file, tag = tag, stringsAsFactors = FALSE),
        class = c("svTestData", "data.frame"),
        row.names = as.character(attr(.Log[[test]], "stats")["tests"]))
    ## Add them also to the log of my test
    .Log[[test]][nrow(.Log[[test]]) + 1, ] <- .Log$.lastTest
    ## Do we print detailed results for this test?
	if (is.null(printTest)) printTest <- !interactive() # Guess it from context
	if (isTRUE(printTest)) print(.Log$.lastTest)
}

.prepareUnit <- function (name, dir)
{
	## Prepare for writing a test unit file
	dir <- gsub("\\\\", "/", as.character(dir)[1])
	## Check that dir exists (do not create it!)
	if (!file.exists(dir) || !file.info(dir)$isdir)
		stop("'dir' must be an existing directory")
	## If dir is tempdir(), then, make sure there are no other runit*.R files
	## left (should not!) - One can store only one unit at a time in tempdir()!
	if (dir == gsub("\\\\", "/", tempdir())) {
		runitfiles <- list.files(dir, pattern = "^runit.*\\.[r|R]$",
			full.names = TRUE)
		unlink(runitfiles)
	}
	Unit <- file.path(dir, paste("runit", name, ".R", sep = ""))
	cat("## Test unit '", name, "'\n", sep = "", file = Unit)
	return(Unit)
}

.writeSetUp <- function (unit, file = "", msg = "", tag = "", code = NULL)
{
	## Write the .setUp() function in the test unit file
	## Here, we write a context to localize tested objects and test unit files
	catUnit <- function(...) cat(..., sep = "", file = unit, append = TRUE)
	catUnit('\n.setUp <-\n')
	catUnit('function () {\n')
	catUnit('\t## Specific actions for svUnit: prepare context\n')
	catUnit('\tif ("package:svUnit" %in% search()) {\n')
	catUnit('\t\t.Log <- Log() ## Make sure .Log is created\n')
	catUnit('\t\t.Log$..Unit <- "', unit, '"\n')
	catUnit('\t\t.Log$..File <- "', as.character(file)[1], '"\n')
	catUnit('\t\t.Log$..Obj <- ""\n')
	catUnit('\t\t.Log$..Tag <- "', tag, '"\n')
	catUnit('\t\t.Log$..Msg <- "', paste(msg, collapse = "\n"), '"\n')
	catUnit('\t\trm(..Test, envir = .Log)\n')
	catUnit('\t}\n')
	if (!is.null(code)) catUnit(paste("\t", code, collapse = "\n"))
	catUnit('}\n')
}

.writeTearDown <- function (unit, code = NULL, rm.unit = TRUE, rm.file = TRUE)
{
	## Write the .tearDown() function in the test unit file
	## Here, we undo what was done in .setUp()
	catUnit <- function(...) cat(..., sep = "", file = unit, append = TRUE)
	catUnit('\n.tearDown <-\n')
	catUnit('function () {\n')
	if (!is.null(code)) catUnit(paste("\t", code, collapse = "\n"))
	catUnit('\t## Specific actions for svUnit: clean up context\n')
	catUnit('\tif ("package:svUnit" %in% search()) {\n')
	catUnit('\t\t.Log$..Unit <- ""\n')
	catUnit('\t\t.Log$..File <- ""\n')
	catUnit('\t\t.Log$..Obj <- ""\n')
	catUnit('\t\t.Log$..Tag <- ""\n')
	catUnit('\t\t.Log$..Msg <- ""\n')
	catUnit('\t\trm(..Test, envir = .Log)\n')
	catUnit('\t}\n')
	catUnit('}\n')
}

.writeTest <- function (unit, objname, pos = .GlobalEnv, obj = NULL)
{
	## Make sure that the name of a test function is syntactically correct
	## and starts with 'test'
	if (regexpr("^test", objname) > -1) {
		testname <- objname
	} else {
		testname <- paste("test", objname, sep = "")
	}
	## Write the first line in the file
	catUnit <- function(...) cat(..., file = unit, append = TRUE)
	catUnit('\n"', testname, '" <-\n', sep = "")
	## Get the object
	if (missing(obj)) {
		## Look for 'objname' in 'pos'
		if (!exists(objname, where = pos)) {
			Test <- ""
		} else {
			Test <- test(get(objname, pos = pos))
		}
	} else {
		Test <- test(obj)
	}
	if (is.character(Test)) {
		## Create a dummy test with DEACTIVATED entry indicating missing object
		body <- c(
			'function() {',
			paste('\tDEACTIVATED("Object', objname, 'not found!")'),
			'}\n')
	} else if (is.null(Test)) {
		## Create a dummy test with DEACTIVATED entry indicating missing test
		body <- c(
			'function() {',
			paste('\tDEACTIVATED("Object', objname, 'has no tests!")'),
			'}\n')
	} else {
		## Get the body of the test function
		capture.body <- function(Data) {
			rval <- NULL
			File <- textConnection("rval", "w", local = TRUE)
			sink(File)
			on.exit({ sink(); close(File) })
			dput(Data, file = File, control = "useSource")
			on.exit()
			sink()
			close(File)
			return(rval)
		}
		body <- capture.body(Test)
	}
	## Write the body of the test function in the file
	catUnit(body, sep = "\n")
}

.runTest <- function (x, envir, test, objfile = "", unit = "", tag = "",
msg = "")
{
	## Run one test in a protected environment catching errors and warnings
	## and preparing a suitable context
	name <- sub("^test\\.(.+)\\.$", "\\1", test)

	## The context is written in the .Log, but previous context is saved
	## and restored on return
	.Log <- Log()	# Make sure that .Log exists, or create a new one
	oContext <- c(Unit = .Log$..Unit, Obj = .Log$..Obj, File = .Log$..File,
		Msg = .Log$..Msg, Tag = .Log$..Tag)
	on.exit({
		.Log$..Unit <- as.character(oContext[1])
		.Log$..Obj <- as.character(oContext[2])
		.Log$..File <- as.character(oContext[3])
		.Log$..Msg <- as.character(oContext[4])
		.Log$..Tag <- as.character(oContext[5])
	})
	.Log$..Unit <- unit				# The unit file
	.Log$..Obj <- name 				# Name of the tested object
	.Log$..File <- objfile 			# Where is the code source of 'name'?
	.Log$..Msg <- paste(msg, collapse = "\n") # Message for this test
	.Log$..Tag <- tag 				# A tag in objfile to locate code
	## Define the test and save possible existing definition to restore it
	if (exists("..Test", envir = .Log, inherits = FALSE)) {
		oTest <- .Log$..Test
		on.exit(.Log$..Test <- oTest, add = TRUE)
	} else on.exit(rm("..Test", envir = .Log), add = TRUE)
	.Log$..Test <- test 			# Define the name of the test

	if (missing(envir)) {
		## The environment where the test is run
		envir <- new.env(parent = .GlobalEnv)
		envir[[test]] <- x 					# A copy of the test code
		envir$.setUp <- function() {}		# Fake .setUp
		envir$.tearDown <- function() {}	# Fake .tearDown
	}
	## We need this installed in our sandbox .TestEnv to run the test
	envir$.LogWarnings <- list() # A list to collect warnings

	## Clear the corresponding log, if it exists
	if (exists(test, envir = .Log, inherits = FALSE))
		rm(list = test, envir = .Log)

	## Evaluate the test function in .testEnv, catching errors
	owarn <- getOption("warn")
	on.exit(options(warn = owarn), add = TRUE)
	if (isTRUE(getOption("svUnit.StopOnWarning"))) nwarn <- 2 else nwarn <- -1
	options(warn = nwarn)

	## Evaluate the test function in the .TestEnv environment
	cmd <- paste("evalq(.LogRes <- try( { .setUp(); `", test,
		"`(); .tearDown() }, silent = TRUE), envir = envir)", sep = "")
	eval(parse(text = cmd))

	## Analyze error => is it a deactivation or error in the code?
	if (inherits(Res <- envir$.LogRes, "try-error")) {
		## We record this as a test returning **ERROR** or DEACTIVATED
		.logTest(0, test)
		## Did we encountered a DEACTIVATED() command or something else?
		if (regexpr("DEACTIVATED\\(", Res) > -1) {
			Msg <- sub("^[^:]+: *", "", as.character(Res))
			Msg <- sub("\n$", "", Msg)
			.logTestData(test, msg = Msg, call = "", timing = NA,
				val = -2, res = "\n")
		} else {
			## This is an error (something wrong with the code in the test fun)
			.logTestData(test, msg = "", call = deparse(sys.call()),
				timing = NA, val = -1, res = paste(Res, collapse = "\n"))
		}
	}
	return(test)
}

.assignTemp <- function (x, value)
    assign(x, value, envir = .TempEnv())

.getTemp <- function (x, default = character(0))
{
    if  (exists(x, envir = .TempEnv(), inherits = FALSE)) {
        return(get(x, envir = .TempEnv(), inherits = FALSE))
    } else { # Variable not found, return the default value
        return(default)
    }
}

.TempEnv <- function ()
{
    pos <-  match("SciViews:TempEnv", search())
    if (is.na(pos)) { # Must create it
        `SciViews:TempEnv` <- list()
        Attach <- function (...) get("attach", mode = "function")(...)
        Attach(`SciViews:TempEnv`, pos = length(search()) - 1)
        rm(`SciViews:TempEnv`)
        pos <- match("SciViews:TempEnv", search())
    }
    return(pos.to.env(pos))
}
