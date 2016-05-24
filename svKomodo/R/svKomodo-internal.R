.loadSvOptions <- function ()
{
	## Try loading SciViews options from ~/.SciViewsConfig.RData
	argOption <- function(arg)
	{
		args <- commandArgs(trailingOnly = TRUE)
		regexp <- paste("^", arg, "=", sep = "")
		where <- grepl(regexp, args)
		if (any(where)) {
			arg <- rev(args[where])[1]
			sub(regexp, "", arg)
		} else NULL
	}
	
	## For those items not defined yet
	.SciViewsConfig <- list() # Default value
	if (file.exists("~/.SciViewsConfig.RData"))
		load("~/.SciViewsConfig.RData")
	## Change options...
	ko.type <- argOption("ko.type")
	if (!is.null(ko.type)) options(ko.type = ko.type)
	if (is.null(getOption("ko.type")))
		options(ko.type = .SciViewsConfig$ko.type)
	if (is.null(getOption("ko.type")))
		options(ko.type = "http")
		
	ko.host <- argOption("ko.host")
	if (!is.null(ko.host)) options(ko.host = ko.host)
	if (is.null(getOption("ko.host")))
		options(ko.host = .SciViewsConfig$ko.host)
	if (is.null(getOption("ko.host")))
		options(ko.host = "localhost")
	
	ko.port <- argOption("ko.port")
	if (!is.null(ko.port)) options(ko.port = ko.port)
	if (is.null(getOption("ko.port")))
		options(ko.port = .SciViewsConfig$ko.port)
	if (is.null(getOption("ko.port")))
		options(ko.port = 7052)	
	
	ko.kotype <- argOption("ko.kotype")
	if (!is.null(ko.kotype)) options(ko.kotype = ko.kotype)
	if (is.null(getOption("ko.kotype")))
		options(ko.kotype = .SciViewsConfig$ko.kotype)
	if (is.null(getOption("ko.kotype")))
		options(ko.kotype = "socket")
	
	ko.serve <- argOption("ko.serve")
	if (!is.null(ko.serve)) options(ko.serve = ko.serve)
	if (is.null(getOption("ko.serve")))
		options(ko.serve = .SciViewsConfig$ko.serve)
	if (is.null(getOption("ko.serve")))
		options(ko.serve = 8888)	
	
	ko.activate <- argOption("ko.activate")
	if (!is.null(ko.activate)) options(ko.activate = ko.activate)
	if (is.null(getOption("ko.activate")))
		options(ko.activate = .SciViewsConfig$ko.activate)
	if (is.null(getOption("ko.activate")))
		options(ko.activate = 8888)
	
	ko.id <- argOption("ko.id")
	if (!is.null(ko.id)) options(ko.id = ko.id)
	if (is.null(getOption("ko.id")))
		options(ko.id = .SciViewsConfig$ko.id)
	if (is.null(getOption("ko.id")))
		options(ko.id = "SciViewsK")
	
	R.id <- argOption("R.id")
	if (!is.null(R.id)) options(R.id = R.id)
	if (is.null(getOption("R.id")))
		options(R.id = .SciViewsConfig$R.id)
	if (is.null(getOption("R.id")))
		options(R.id = "R")
	
	R.initdir <- argOption("R.initdir")
	if (!is.null(R.initdir)) options(R.initdir = R.initdir)
	if (is.null(getOption("R.initdir")))
		options(R.initdir = .SciViewsConfig$R.initdir)
	R.initdir <- getOption("R.initdir")
	if (is.null(R.initdir)) {
		options(R.initdir = NULL)
	} else {
		## Check that the dir exists...
		isfile <- try(file.exists(R.initdir))
		if (inherits(isfile, "try-error") || !isfile || !file.info(R.initdir)$isdir)
			options(R.initdir = NULL)
	}
		
	width <- argOption("width")
	if (!is.null(width)) options(width = width)
	if (is.null(getOption("width")))
		options(width = .SciViewsConfig$width)
	if (is.null(getOption("width")))
		options(width = 80)
	
	OutDec <- argOption("OutDec")
	if (!is.null(OutDec)) options(OutDec = OutDec)
	if (is.null(getOption("OutDec")))
		options(OutDec = .SciViewsConfig$OutDec)
	if (is.null(getOption("OutDec")))
		options(OutDec = ".")
		
	OutSep <- argOption("OutSep")
	if (!is.null(OutSep)) options(OutSep = OutSep)
	if (is.null(getOption("OutSep")))
		options(OutSep = .SciViewsConfig$OutSep)
	if (is.null(getOption("OutSep")))
		options(OutSep = ",")
}

.onAttach <- function (lib, pkg)
{
	## Make sure config is OK
	.loadSvOptions()
	
	## Check that usually inactivated rc.settings are set
	assignTemp(".old.rc.settings", rc.settings())
	rc.settings(ipck = TRUE)
	
	## Create our SciViews task callback manager
	## Should eliminate this???
	## PhG: inactivated for now, because it makes problems in R!!!
	#assignTemp(".svTaskCallbackManager", svTaskCallbackManager())
	koInstall()
	
	## Temporary change of a few R options that make problem with SciViews
	## Do not use fancy quotes
	assignTemp(".useFancyQuotes", getOption("useFancyQuotes"))
	options(useFancyQuotes = FALSE)
	## Limit output to 999 items
	assignTemp(".max.print", getOption("max.print"))
	options(max.print = 999)
	
	## If svStart function exists in.GlobalEnv, delegate the rest to it!
	if (exists("svStart", envir = .GlobalEnv, mode = "function",
		inherits = FALSE)) return()
	
	## Determine which server to start
	type <- getOption("ko.type")
	req <- require
	if (type == "http") {
		req("svHttp", character.only = TRUE, lib.loc = lib)
		if (interactive()) try(svHttp::startHttpServer())
	} else {
		req("svSocket", character.only = TRUE, lib.loc = lib)
		if (interactive()) try(svSocket::startSocketServer())
	}
		
	## This comes from svStart... and should really be placed here indeed!
	## Look if and where komodo is installed
	if (Sys.getenv("koAppFile") != "") {
		Komodo <- Sys.getenv("koAppFile")
	} else Komodo <- ""
	
	if (.Platform$OS.type == "unix") {
		if (Komodo == "")
			Komodo <- "/usr/local/bin/komodo"  # Default location
		if (!file.exists(Komodo))
			Komodo <- Sys.which("komodo")[1]
		if (length(Komodo) == 0 || Komodo == "")
			Komodo <- NULL
			
		## Just to avoid warnings while compiling outside of Windows...
		readRegistry <- function() return()
	} else {  # Windows
	    ## If komodo path was not passed in environment
		if (Komodo == "") {
			Komodo <- NULL
			
			err.null <- function (e) return(NULL)
			## On Windows, 'komodo' should be enough
			## But for reasons that escape me, Komodo seems to strip off its
			## own directory from the path variable. So, I have to restore
			## it from the Windows registry :-(

			## Try several ways to get Komodo path from registry.
			key <- "SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\App Paths\\komodo.exe"
			Komodo <-
				tryCatch(readRegistry(key, hive = "HLM")[["(Default)"]],
					error = err.null)

			if (is.null(Komodo) || !file.exists(Komodo)) {
				key <- "Applications\\komodo.exe\\shell\\open\\command"
				Komodo <-
					tryCatch(readRegistry(key, hive = "HCR")[["(Default)"]],
						error = err.null)
				if (!is.null(Komodo))
					Komodo <- sub(" *\\\"%[1-9\\*].*$", "", Komodo)
			}

			if (is.null(Komodo) || !file.exists(Komodo)) {
				key <- "SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment"
				Path <-
					tryCatch(readRegistry(key, hive = "HLM")$Path,
						error = err.null)
				if (!is.null(Path) && !is.na(Path) && Path != "") {
					Path <- strsplit(Path, ";")[[1]]
					Path <- Path[sapply(Path,
						function (x)
							file.exists(file.path(x, "komodo.exe")))][1]
					Komodo <- gsub("\\\\+", "\\\\", file.path(Path,
						"komodo.exe", fsep = "\\"))
				}
			}
		}
	}

	if (length(Komodo) && Komodo != "" && file.exists(Komodo)) {
		## Change the editor and the pager to Komodo
		options(pager2 = getOption("pager"))
		## A custom pager consists in displaying the file in Komodo
		svPager <- function (files, header, title, delete.file) {
			files <- gsub("\\", "\\\\", files[1], fixed = TRUE)
			res <- tryCatch(koCmd(sprintf('sv.r.pager("%s", "%s")', files, title)),
				error = function (e) return(FALSE))
			if (res == FALSE) {
				## Try using pager2 instead
				pager2 <- getOption("pager2")
				if (is.null(pager2)) {
					stop("You must start Komodo Edit for displaying this file!")
				} else if (is.function(pager2)) {
					pager2(files = files, header = header, title = title,
						delete.file = delete.file)
				} else {
					file.show(files, header = header, title = title,
						delete.file = delete.file, pager = pager2)
				}
			} else {
				if (delete.file)
					koCmd(sprintf('window.setTimeout("try { sv.tools.file.getfile(\\"%s\\").remove(false); } catch(e) {}", 10000);', files));
			}
		}

		options(browser2 = getOption("browser"))
		svBrowser <- function (url) {
			url <- gsub("\\", "\\\\", url, fixed = TRUE)
			## If the URL starts with '/', I could safely assume a file path
			## on Unix or Mac and prepend 'file://'
			url <- sub("^/", "file:///", url)
			res <- tryCatch(koCmd(sprintf("sv.command.openHelp(\"%s\")", url)),
				error = function (e) return(FALSE))
			if (res == FALSE) {
				## Try using browser2 instead
				browser2 <- getOption("browser2")
				if (is.null(browser2)) {
					stop("You must start Komodo Edit for browsing files")
				} else browseURL(url, browser = browser2)
			}
		}

		options(editor2 = getOption("editor"))
		options(editor = Komodo, browser = svBrowser, pager = svPager)
	} else {
		Komodo <- NULL
		
		## TODO: use packageStartupMessage instead!
		#cat("R cannot find Komodo.", file = stderr())
		#if (.Platform$OS.type == "unix") {
		#	cat("Please, follow instructions at",
		#		"http://www.sciviews.org/SciViews-K to install it correctly.",
		#		"In particular, you must create a symbolic link in /user/local/bin:",
		#		"sudo ln -s <KomodoBinLocation>/komodo /usr/local/bin/komodo",
		#		"otherwise, R cannot find it!", sep = "\n", file = stderr())
		#} else {
		#	cat("Please, make sure you install it correctly\n",
		#		"You can find it at http://www.activestate.com/Products/komodo_edit.\n",
		#		file = stderr())
		#}
	}

	## Make sure we use HTML help (required for Alt-F1 and Alt-Shift-F1)
	## to display R help in Komodo Edit
	options(help_type = "html")
	## Make sure the help server is started
	port <- try(tools::startDynamicHelp(), silent = TRUE)
	if (inherits(port, "try-error")) {
		## Dynamic help is already started
		if (R.Version()$`svn rev` >= 67550) {
			port <- tools::startDynamicHelp(NA)
		} else {
			port <- getNamespace("tools")$httpdPort
		}
	}
	## Record the home page for the help server in an option
	options(helphome = paste("http://127.0.0.1:", port,
		"/doc/html/index.html", sep = ""))

	## I need to get the help file URL, but help() does not provide it any
	## more! This is a temporary workaround for this problem
	assignTemp("getHelpURL", function (x, ...) {
		file <- as.character(x)
		if (length(file) == 0) return("")
		## Extension ".html" may be missing
		htmlfile <- basename(file)
		## Get the HTML help server port
		if (R.Version()$`svn rev` >= 67550) {
			port <- tools::startDynamicHelp(NA)
		} else {
			port <- getNamespace("tools")$httpdPort
		}
		if (length(file) > 1) {
			## If more then one topic is found
			paste("http://127.0.0.1:", port,
				"/library/NULL/help/", attr(x,"topic"), sep = "")
		} else {
			if(substring(htmlfile, nchar(htmlfile) -4) != ".html")
				htmlfile <- paste(htmlfile, ".html", sep="")
			paste("http://127.0.0.1:", port,
			"/library/", basename(dirname(dirname(file))),
			"/html/", htmlfile, sep = "")
		}
	})

## print() method of object returned by help() is very unflexible for R.app and
## does not allow in any way to use anything else than the R.app internal
## browser for help!!!
## That makes me very unhappy! Hey guys, I would like to use SciViews help
## browser here! So, no other solution than to be even harsher, and to force
## rewriting of the print function in base environment!!!
## (problem emailed to Simon Urbanek on 03/11/2009... I hope he will propose
## a work-around for this in R 2.12!!!)
#if (compareVersion(rVersion, "2.11.0") < 0) {
#	source("print.help_files_with_topic210.R")
#} else if (compareVersion(rVersion, "2.14.0") < 0) {
#	source("print.help_files_with_topic211.R")
#} else {
#	source("print.help_files_with_topic214.R")
#}
	
	R.initdir <- getOption("R.initdir")
	if (!is.null(R.initdir) && !is.na(R.initdir)) {
		## Change the working directory to the provided directory
		try(setwd(R.initdir), silent = TRUE)

		## Create a .Last.sys function that clears some variables in .GlobalEnv
		## and then, switch to R.initdir before closing R. The function is
		## stored in SciViews:TempEnv
		assignTemp(".Last.sys", function () {
			## Eliminate some known hidden variables from .GlobalEnv to prevent
			## saving them in the .RData file
			if (exists(".required", envir = .GlobalEnv, inherits = FALSE))
				rm(list = ".required", envir = .GlobalEnv, inherits = FALSE)
			## Note: .SciViewsReady is now recorded in SciViews:TempEnv instead
			## of .GlobalEnv, but we leave this code for old workspaces...
			if (exists(".SciViewsReady", envir = .GlobalEnv, inherits = FALSE))
				rm(list = ".SciViewsReady", envir = .GlobalEnv, inherits = FALSE)
			## If a R.initdir is defined, make sure to switch to it, so that
			## the session's workspace and command history are written at the
			## right place (in case of error, no change is made!)
			try(setwd(getOption("R.initdir")), silent = TRUE)
			## Clean up everything in Komodo
			tryCatch(
				svKomodo::koCmd("window.setTimeout(\"sv.r.closed();\", 1000);"),
				error = function (e) invisible(NULL))
		})

		msg <- paste("Session directory is", dQuote(R.initdir))
		msg2 <- NULL

		## Do we load .RData and .Rhistory now?
		args <- commandArgs()
		if (!"--vanilla" %in% args && !"--no-restore" %in% args &&
			!"--no.restore-data" %in% args) {
			if (file.exists(".RData")) {
				load(".RData", envir = .GlobalEnv)
				msg2 <- append(msg2, "data loaded")
			} else {
				msg2 <- append(msg2, "no data")
			}
	
			if (file.exists(".Rhistory")) {
				## On R Tk gui:
				## "Error in loadhistory(file) : no history mechanism available"
				## So, do it inside a try()
				history.loaded <- try(loadhistory(), silent = TRUE)
				if (inherits(history.loaded, "try-error"))  {
					msg2 <- append(msg2, "history cannot be loaded")
				} else {
					msg2 <- append(msg2, "history loaded")
				}
			} else {
				msg2 <- append(msg2, "no history")
			}
		} else {
			msg2 <- append(msg2, "data and history not loaded")
		}
	
		packageStartupMessage(paste(msg, " (", paste(msg2, collapse = ", "),
			")", sep = ""))
	
		## Do we reactivate Komodo now?
	#	koact <- getOption("ko.activate")
	#	debugMsg("Reactivate Komodo:", koact)
	#	if (getTemp(".SciViewsReady", FALSE) && koact) {
	#		if ((.Platform$pkgType == "mac.binary")) {
	#			system("osascript -e 'tell application \"Komodo\" to activate'",
	#				wait = FALSE)
	#		} else if (!is.null(Komodo)) {
	#			## TODO: The following starts komodo if not started yet,
	#			## but does not activate it!
	#			system(shQuote(Komodo), wait = FALSE)
	#		}
			## Indicate to Komodo that R is ready
			## and test also communication from R to Komodo!
		#	koCmd('sv.cmdout.message("<<<data>>>", 10000, true);',
		#		data = paste("'", getOption("R.id"), "' (R ",
		#		R.Version()$major, ".", R.Version()$minor,
		#		") connected. Session dir: ",
		#		path.expand(getOption("R.initdir")), sep = ""))
			## ... and refresh the object explorer
			## TODO!
				
			## Differ synching R <-> Komodo to avoid deadlock situation
			# That does not work!
			#koCmd('window.setTimeout("sv.r.objects.getPackageList(true, true, true);", 500)')
			#koCmd('window.setTimeout("sv.r.test(true, true);", 500)')
	#	}
		## Update info in Komodo
		#debugMsg("Contacting Komodo with koCmd")
		#invisible(koCmd(paste(
		#	"sv.socket.rUpdate()",
		#	"sv.cmdout.append('R is started')",
		#	"sv.command.updateRStatus(true)",
		#	sep = ";"))
		#)
		## Refreshing Komodo's GUI elements
#		try(koCmd(paste('sv.r.running = true; sv.socket.charset = "',  localeToCharset()[1],
#			'"; sv.cmdout.message("' , R.version.string, ' is ready!");',
#			' window.setTimeout("sv.r.objects.getPackageList(true, true, true);", 1000);',
#			sep = "")), silent = TRUE)
		#try(koRefresh(force = TRUE), silent = TRUE)
		#}
	
		## TODO: eliminate this first from svStart.R... otherwise, it is run twice!
		## Do we have a .Rprofile file to source?
	#	rprofile <- file.path(c(getwd(), Sys.getenv("R_USER")), ".Rprofile")
	#	rprofile <- rprofile[file.exists(rprofile)][1]
	#	if (!is.na(rprofile)) {
	#		source(rprofile)
	#		debugMsg("Loaded file:", rprofile)
	#	}
	
		## If there is a function called ".First" in .GlobalEnv, run it now
		.First <- NULL
		if (exists(".First", envir = .GlobalEnv, mode = "function",
			inherits = FALSE))
			eval(.First(), envir = .GlobalEnv)
	}
}

.onUnload <- function (libpath)
{
	## Restore R options that preexisted
	try(options(useFancyQuotes = getTemp(".useFancyQuotes", TRUE)),
		silent = TRUE)
	try(options(max.print = getTemp(".max.print", 99999)),
		silent = TRUE)
	
	## Restore editor, browser & pager
	editor2 <- getOption("editor2")
	if (!is.null(editor2)) {
		options(editor = editor2)
		options(editor2 = NULL)
	}
	browser2 <- getOption("browser2")
	if (!is.null(browser2)) {
		options(browser = browser2)
		options(browser2 = NULL)
	}
	pager2 <- getOption("pager2")
	if (!is.null(pager2)) {
		options(pager = pager2)
		options(pager2 = NULL)
	}
	#serve <- getOption("ko.serve")
	#if (!is.null(serve) && "package:svSocket" %in% search())
	#	stopSocketServer(port = as.integer(serve)[1])
	koUninstall()
	## Remove the SciViews task callback manager
	try(removeTaskCallback("SV-taskCallbackManager"), silent = TRUE)
	try(rmTemp(".svTaskCallbackManager"), silent = TRUE)
	
	## Restore rc.settings
	settings <- getTemp(".old.rc.settings", rc.settings())
	do.call("rc.settings", as.list(settings))
	rmTemp(".old.rc.settings")
}

.packageName <- "svKomodo"
