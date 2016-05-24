
#-----------------------------------------------------------------------------------------
#  miscellaneous routines including utilities used throughout the package.
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#  demoSets()			return names of text files with ID sets
#  buildDemoSets()		create the .rda of "biom" objects included as package data
#
#  set-1.tsv		7 metagenomes (in three groups)
#  set-2.tsv		24 metagenomes (in two groups)
#  set-3.tsv		32 metagenomes (in three groups)
#  set-4.tsv		16 metagenomes (in three groups)
#
#  do not use these for now:
#  set-5.tsv		1606 metagenomes (HMP)
#  set-6.tsv		4 projects of various sizes
#  set-7.tsv		10 metagenomes and 3 projects
#-----------------------------------------------------------------------------------------

demoSets <- function () {
	dir (file.path (path.package ("matR"), "extdata"), pattern="set-*", full.names=TRUE)
	}

buildDemoSets <- function (n=TRUE, file="demoObjects.rda") {
	setParams <- list(
		list (request="function", group_level="level2"),							# set-1
		list (request="organism", group_level="phylum"),							# set-2
		list (request="function", group_level="level1", evalue=1),					# set-3
		list (request="organism", group_level="domain", source="Greengenes"))		# set-4
	buildArgs <- mapply (
		c, 
		file = demoSets() [which], 
		setParams [which], 
		USE.NAMES=FALSE)
	objects <- lapply (
		buildArgs, 
		do.call, 
		what="biomRequest")
	names (objects) <- paste0 ("xx", 1:length (objects))

####  fix up non-ASCII characters in metadata

	for (j in 1:length (objects))
		objects [[j]] $ columns <- rapply (
			objects [[j]] $ columns, 
			iconv, 
			"character", 
			how='replace', 
			to='ASCII', 
			sub='?')

	save (list = names(objects), file = file, envir = list2env (objects))
	message(
		"Created file \"", file, "\" in ", getwd(),
		".\nFor package build move to \"matR/data\"")
	}


#-----------------------------------------------------------------------------------------
#  readSet()		read ids/metadata from a tsv file		(character or data.frame)
#  scrubSet()		clean up a specification of "ids"		(character)
#  scrapeSet()		return "resources" per specification	(character)
#  expandSet()		expand projects to metagenomes			(character or data.frame)
#-----------------------------------------------------------------------------------------

readSet <- function (file) {
	df <- read.table(file, header=F, sep="\t", colClasses="character", stringsAsFactors=TRUE)
	if(ncol(df) > 1) {
		colnames(df) <- df[1, , drop=TRUE]
		rownames(df) <- df[, 1, drop=TRUE]
		df[-1, -1, drop=FALSE]
	} else df[, , drop=TRUE]
	}

scrubSet <- function (x, resources = "metagenome") {
	if (is.data.frame (x)) {
		scrubSet (rownames (x), resources)
	} else {
		y <- strsplit (collapse (as.character (x)), "[[:space:]+]") [[1]]
		y <- y [y != ""]
		pfx <- match.arg(
			rep(resources, len=length(y)),
			c("metagenome", "project"),
			several.ok = TRUE)
		paste0(
			ifelse (substr(y,1,3) %in% c("mgm", "mgp"), 
				"",
				c (metagenome="mgm", project="mgp") [pfx]), 
			y)
		}
	}

scrapeSet <- function (x) {
	if (is.data.frame (x)) {
		scrapeSet (rownames (x))
	} else
		c("metagenome", "project") [match (substr(x, 1, 3), c("mgm", "mgp"), nomatch=1)]
	}

expandSet <- function (x) {
	if (is.data.frame (x)) {
		rownames (x) <- scrubSet (rownames (x))
		y <- scrapeSet (rownames (x))
		if (!any (y == "project")) return (x)
		z <- as.list (rownames (x))
		z [y == "project"] <- metadata (rownames (x) [y == "project"])
		j <- rep (1:length(z), sapply(z, length))
		x <- x [j, , drop=FALSE]
 		rownames (x) <- unlist (z)
 		x
	} else {
		x <- scrubSet (x)
		y <- scrapeSet (x)
		if (!any (y == "project")) return (x)
		z <- as.list (x)
		z [y == "project"] <- metadata (x [y == "project"])
		unlist (z)
		}
	}


#-----------------------------------------------------------------------------------------
#  handling of "suggested dependencies".
#  official dependencies are minimized to remove obstacles to installation.
#-----------------------------------------------------------------------------------------

hazPackages <- function() {
	me <- packageDescription ("matR")
	need <- unlist (strsplit (c (me$Imports, me$Suggests), "[^[:alnum:]\\.]+"))
	sapply (need, function (x) length (find.package (x, quiet = TRUE)) > 0)
	}

dependencies <- function (prompt = TRUE) {
	need <- !hazPackages()
	if (any (need)) {
		message("matR uses other software.. blah blah blah.. here is what we\'ll do")
		message("Suggested package(s) missing: ", collapse (names(need) [need]), "\n")
		if (prompt) {
			chooseBioCmirror (graphics = FALSE)
			cat ("\n")
			chooseCRANmirror (graphics = FALSE)
			cat ("\n")
			setRepositories (graphics = FALSE)
			cat ("\n")
			}
		install.packages (names (need) [need])
		haz <- hazPackages()
		if (all (haz)) {
			message("\nAll suggested packages have been installed.\nNow quit and restart R.")
		} else message("\nPackage(s) could not be installed: ", collapse (names(haz) [!haz]))
	} else 
		message ("All suggested packages appear to be installed.")
	}


#-----------------------------------------------------------------------------------------
#  'step.through()' to help with demos.
#
#  stepper():
#  read the file as text and echoes each line exactly.
#  each command must fit on a line, and blank lines are simply echoed.
#
#  stepper2():
#  parses the whole file first.  commands may span lines.
#  comments are not displayed.  commands are reformatted to standard appearance.
#-----------------------------------------------------------------------------------------

step.through <- function (demo, file) {
	stepper (if (!missing (file)) {
			file
		} else
			file.path (path.package ("matR"), "demo", paste0 (demo, ".R")))
	}

stepper <- function (file) {
	lines <- readLines (file)
	for (j in 1:length (lines))
		if (lines [j] != "") {
			cat (getOption ("prompt"), lines [j], sep = "")
			readLines (n = 1, warn = FALSE)
			R <- withVisible (eval.parent (parse (text = lines [j])))
			if (R$visible && !is.null (R$value)) print (R$value)
			}
		else cat ("\n")
	}

stepper2 <- function (file) {
	exprs <- parse (file = file)
	for (j in 1:length (exprs)) {
		cat (getOption ("prompt"), as.character (exprs [j]), sep = "")
		readLines (n = 1, warn = FALSE)
		eval (exprs [j])
		}
	}


#-----------------------------------------------------------------------------------------
#  utilities for docs.
#-----------------------------------------------------------------------------------------

seeDoc <- function (f) {
	outfile <- tempfile (fileext = ".html")
	browseURL (tools::Rd2HTML (f, outfile))
	}

buildHTMLDocs <- function (docs) {
	for (d in docs) tools::Rd2HTML (d, paste ("./html/", unlist (strsplit (d, ".", fixed = TRUE)) [1], ".html", sep=""), Links = tools::findHTMLlinks ())
	}


list2df <- function (li) {
#-----------------------------------------------------------------------------------------
#  support for bringing nested list structures into a data.frame with NAs where needed
#
#  unlist to depth one, while keeping names of top-level list entries
#  the df variables will be the union of all names
#  data.frame has one row per element of the original list; one column per variable
#  don't apply t() in case of a single variable
#  final result:  NA placed wherever a variable is missing from a list entry
#-----------------------------------------------------------------------------------------
	li <- sapply (li, unlist, simplify=FALSE)				
	vars <- sort (unique (unlist (lapply (li, names))))		
	y <- sapply (li, "[", vars)
	if (is.matrix (y)) y <- t(y)
	as.data.frame (y, stringsAsFactors=FALSE)
	}

tagline <- function () {
#-----------------------------------------------------------------------------------------
#  package name.
#  preprocess source with:
#	 sed s/XXXBUILDXXX/$commit/g matR/R/init.R > init.Rtemp
#	 mv init.Rtemp matR/R/init.R
#-----------------------------------------------------------------------------------------
	ss <- " c71ee9"
	if (substr (ss, 2, 9) == "XXXBUILD") ss <- ""
	paste0 ("matR: metagenomics analysis tools for R (", packageVersion("matR"), ss, ")")
	}

warning <- function (...) {
#-----------------------------------------------------------------------------------------
#  stylish warnings
#-----------------------------------------------------------------------------------------
	base::warning ("matR: ", ...,  call.=FALSE)
	}

collapse <- function (x, ..., sep = " ") {
#-----------------------------------------------------------------------------------------
#  make length-one string from a character vector
#-----------------------------------------------------------------------------------------
	paste(x, ..., sep = sep, collapse = sep)
	}

abbrev <- function (s, n, where = "right") {
#-----------------------------------------------------------------------------------------
#  abbreviates each element of a character vector to fit a given width
#  adds "..." where text is omitted (left, middle, or right)
#  -->still need abbrev from middle and left ...
#-----------------------------------------------------------------------------------------
	toolong <- function (s, n) sapply (s, function (x) { nchar (x) > n - 3 }, USE.NAMES = FALSE)
	paste (strtrim (s, width = n - 3), 
		ifelse (toolong (s, n), "...", ""), 
		sep = "")
	}
