### Make a .pot file with translatable strings found in a Komodo project/package
### Translatable strings are:
### 1) All names
### 2) In snippets:
###    [[%ask:R-desc:XXX]]
###    [[%ask:R-tip:XXX]]
###    [[%ask|pref:URL-help:XXX]] and [[%ask|pref:RWiki-help:XXX]]
###    [[%tr:XXX]
### 3) In macros:
###    Strings inside _("XXX"), with _() being a function returning its argument
kpf2pot <- function (kpfFile, potFile)
{
    if (missing(kpfFile)  || is.null(kpfFile) || is.na(kpfFile))
		stop("'kpfFile' must be provided")
	if (length(kpfFile) != 1)
		stop("You must provide only a single file name/path for 'kpfFile'")
	if (!file.exists(kpfFile))
		stop("The kpfFile is not found")
	if (missing(potFile)) 
        potFile <- sub("\\.kpf", ".pot", kpfFile)
	if (kpfFile == potFile)
		potFile <- paste(kpfFile, "pot", sep = ".")
	## Extract translatable strings from this file
	doc <- xmlRoot(xmlTreeParse(kpfFile))
	imax <- xmlSize(doc)
	if (imax < 1) stop("No node found in the file!")
	## Collect all strings to be translated
	s <- character(0)
	for (i in 1:imax) {
		n <- xmlGetAttr(doc[[i]], "name")
		if (!is.null(n)) {
			s <- c(s, n) 
			type <- xmlName(doc[[i]])
			## If this is a snippet, look for other translatable strings
			if (type == "snippet") {
				snip <- xmlValue(doc[[i]])
				chunks <- strsplit(snip, "\\[\\[|\\]\\]")[[1]]
				## Keep only chunks starting with R-desc:, R-tip:, URL-help:
				## RWiki-help: or %tr:
				chunks <- chunks[grep("^%ask:R-desc:|^%ask:R-tip:|^%ask:URL-help:|^%ask:RWiki-help:|^%pref:URL-help|^%pref:RWiki-help|%tr:", chunks)]
				## Are there any remaining chunks?
				l <- length(chunks)
				if (l > 0) {
					## Eliminate leading stuff
					chunks <- sub("^%ask:[^:]+:|%tr:", "", chunks)
					## and add to the list of items to translate
					s <- c(s, chunks)
				}
			} else if (type == "macro") {
				mac <- xmlValue(doc[[i]])
				## Collect tagged strings (i.e., strings inside _(...))
				repeat {
					str <- sub("^.*_\\(\"(.*[^\\\\])\"\\).*$", "\\1", mac)
					if (str == mac) break
					s <- c(s, gsub('\\\\"', '"', str))
					mac <- sub("^(.*)(_\\(\".*[^\\\\]\"\\))(.*)$",
						"\\1(...)\\3", mac)
				}
				repeat {
					str <- sub("^.*_\\('(.*[^\\\\])'\\).*$", "\\1", mac)
					if (str == mac) break
					s <- c(s, gsub("\\\\'", "'", str))
					mac <- sub("^(.*)(_\\('.*[^\\\\]'\\))(.*)$",
						"\\1(...)\\3", mac)
				}
			}
		}
	}
	## Keep only unique strings
	s <- unique(s)
	s <- s[nzchar(s)]
    tmp <- shQuote(encodeString(s), type = "cmd")
    con <- file(potFile, "wt")
    on.exit(close(con))
    writeLines(con = con, c("msgid \"\"", "msgstr \"\"",
		sprintf("\"Project-Id-Version: R %s.%s\\n\"", 
        R.version$major, R.version$minor),
		"\"Report-Msgid-Bugs-To: support@sciviews.org\\n\"", 
        paste("\"POT-Creation-Date: ", format(Sys.time(), "%Y-%m-%d %H:%M"), 
            "\\n\"", sep = ""),
		"\"PO-Revision-Date: YEAR-MO-DA HO:MI+ZONE\\n\"", 
        "\"Last-Translator: FULL NAME <EMAIL@ADDRESS>\\n\"", 
        "\"Language-Team: LANGUAGE <LL@li.org>\\n\"",
		"\"MIME-Version: 1.0\\n\"", 
        "\"Content-Type: text/plain; charset=utf-8\\n\"",
		"\"Content-Transfer-Encoding: 8bit\\n\"", 
        ""))
    for (e in tmp)
		writeLines(con = con, c("", paste("msgid", e), "msgstr \"\""))
	## Check that the .pot file is created
	return(invisible(file.exists(potFile)))	
}

kpz2pot <- function (kpzFile, potFile)
{
	if (missing(kpzFile)  || is.null(kpzFile) || is.na(kpzFile))
		stop("'kpzFile' must be provided")
	if (length(kpzFile) != 1)
		stop("You must provide only a single file name/path for 'kpzFile'")
	if (!file.exists(kpzFile))
		stop("The kpzFile is not found")
    if (missing(potFile)) 
        potFile <- sub("\\.kpz", ".pot", kpzFile)
	if (kpzFile == potFile)
		potFile <- paste(kpzFile, "pot", sep = ".")
    ## The kpz file is a zipped file containing package.kpf in a subdirectory
	f <- file.path(tempdir(), "package.kpf")
	unlink(f)  # Make sure the file does not exist yet
	unzip(kpzFile, junkpaths = TRUE, exdir = tempdir())
	if (!file.exists(file.path(tempdir(), "package.kpf")))
		stop("Impossible to extract the content of the .kpz file.")
	## Run kpf2pot() on this file
	kpf2pot(f, potFile)	
	## Delete extracted file
	unlink(f)
	## Check that the .pot file is created
	return(invisible(file.exists(potFile)))
}

kpfTranslate <- function (kpfFile, langs, poFiles, kpf2Files)
{
    if (missing(kpfFile)  || is.null(kpfFile) || is.na(kpfFile))
		stop("'kpfFile' must be provided")
	if (length(kpfFile) != 1)
		stop("You must provide only a single file name/path for 'kpfFile'")
	if (!file.exists(kpfFile))
		stop("The kpfFile is not found")
	proj <- sub("\\.kpf$", "", kpfFile)
	if (missing(poFiles) || is.na(poFiles)) poFiles <- NULL
	if (missing(langs)) {
		if (is.null(poFiles)) {
			## Try to get the list of suitable .po files in same dir as kpfFile
			pattern <- paste(basename(proj), ".+\\.po$", sep = "-")
			poFiles <- dir(dirname(kpfFile), pattern, full.names = TRUE)
			if (length(poFiles) < 1)
				stop("You must provide 'langs' (ex.: 'fr', or 'de'), or 'poFiles'")
		} else langs <- NULL
	}
	if (is.null(poFiles)) {
		if (is.null(langs))
			stop("You must provide 'langs' (ex.: 'fr', or 'de'), or 'poFiles'")
		## Try to guess poFiles from langs
        poFiles <- paste(proj, "-", langs, ".po", sep = "")
	}
	if (any(kpfFile == poFiles))
		stop("'poFiles' cannot be the same as 'kpfFile'")
	if (any(!file.exists(poFiles)))
		stop("One or more 'poFiles' not found!")
	if (missing(kpf2Files)) {
		## Guess kpf2Files from poFiles
		kpf2Files <- sub("\\.po$", ".kpf", poFiles)
	}
	if (any(kpfFile == kpf2Files))
		stop("'kpfFile' and 'kpf2Files' cannot be the same")
	if (any(poFiles == kpf2Files))
		stop("'poFiles' and 'kpf2Files' cannot be the same")
	if (length(poFiles) != length(kpf2Files))
		stop("Number of items must be the same in 'poFiles' and in 'kpf2Files'")
	## Make sure we create new resulting files
	unlink(kpf2Files)

	## Process each file in turn
	for (h in 1:length(poFiles)) {
		poFile <- poFiles[h]
		kpf2File <- kpf2Files[h]
		## Read the content of the .po file
		tr <- readLines(poFile, encoding = "UTF-8")
		## Keep only lines starting with msgid or msgstr
		trid <- tr[regexpr("^msgid ", tr) == 1]
		trid <- sub("^msgid ", "", trid)
		trmsg <- tr[regexpr("^msgstr ", tr) == 1]
		trmsg <- sub("^msgstr ", "", trmsg)
		## Check that both trid and trmsg have same length
		if (length(trid) != length(trmsg))
			stop("Unequal number of id and translated strings in the .po file '", poFile, "'")
		keep <- trid != "\"\""
		trid <- trid[keep]
		trmsg <- trmsg[keep]
	
		## We need to "unquote" the strings
		unquote <- function (s) {
			## Replace any \\\" by \"
			s <- gsub("\\\\\"", "\"", s)
			## Eliminate leading and trailing quotes
			s <- sub("^\"", "", s)
			s <- sub("\"$", "", s)
			return(s)
		}
		trid <- unquote(trid)
		trmsg <- unquote(trmsg)
		names(trmsg) <- trid
		
		## Extract translatable strings from the .kpf file
		doc <- xmlRoot(xmlTreeParse(kpfFile))
		imax <- xmlSize(doc)
		if (imax < 1) stop("No node found in the kpfFile!")
		## Collect all strings to be translated
		trans <- function (s) {
			tr <- as.character(trmsg[s])
			if (is.na(tr) || tr == "") return(s) else return(tr)
		}
		
		s <- character(0)
		for (i in 1:imax) {
			n <- xmlGetAttr(doc[[i]], "name")
			if (!is.null(n)) {
				## Replace name in attributes of this node
				node <- addAttributes(doc[[i]], name = trans(n), append = TRUE)
				type <- xmlName(node)
				## If this is a snippet, look for other translatable strings
				if (type == "snippet") {
					snip <- xmlValue(node)
					chunks <- strsplit(snip, "\\[\\[|\\]\\]")[[1]]
					## Translate chunks starting with R-desc:, R-tip:, URL-help:,
					## RWiki-help: or %tr:
					toTrans <- grep("^%ask:R-desc:|^%ask:R-tip:|^%ask:URL-help:|^%ask:RWiki-help:|^%pref:URL-help|^%pref:RWiki-help|%tr:", chunks)
					if (length(toTrans) > 0) {
						for (j in toTrans) {
							msg <- sub("^%ask:[^:]+:|%tr:", "", chunks[j])
							header <- sub("^(%ask:[^:]+:|%tr:).*$", "\\1", chunks[j])
							chunks[j] <- paste(header, trans(msg), sep = "")
						}
						## Reconstitute the snippet content using translated messages
						snip <- paste(chunks, c("[[", "]]"),
							sep = "", collapse = "")
						## We need to eliminate latest '[['
						snip <- sub("\\[\\[$", "", snip)
						xmlValue(node) <- snip
					}
				} else if (type == "macro") {
					mac <- xmlValue(node)
					## Translate tagged strings (i.e., strings inside _(...))
					repeat {
						str <- sub("^.*_\\(\"(.*[^\\\\])\"\\).*$", "\\1", mac)
						if (str == mac) break
						s <- trans(gsub('\\\\"', '"', str))
						s <- gsub('"', '\\"', s)
						mac <- sub("^(.*)(_\\(\".*[^\\\\]\"\\))(.*)$",
							paste("\\1%%%%%%(\"", s, "\")\\3", sep = ""), mac)
					}
					repeat {
						str <- sub("^.*_\\('(.*[^\\\\])'\\).*$", "\\1", mac)
						if (str == mac) break
						s <- trans(gsub("\\\\'", "'", str))
						s <- gsub("'", "\\'", s)
						mac <- sub("^(.*)(_\\('.*[^\\\\]'\\))(.*)$",
							paste("\\1%%%%%%('", s, "')\\3", sep = ""), mac)
					}
					mac <- gsub("%%%%%%", "_", mac)
					xmlValue(node) <- mac
				}
				## Replace the node with its translated version
				doc[[i]] <- node
			}
		}
		## In case the project has a name, we change it now
		projname <- xmlGetAttr(doc, "name")
		if (!is.null(projname))
			doc <- addAttributes(doc, name = basename(kpf2File), append = TRUE)
		
		## Make sure the directory where to place the kpf2File exists
		dir.create(dirname(kpf2File), showWarnings = FALSE, recursive = TRUE)
		
		## Save the translated XML content into the second .kpf file
		saveXML(doc, file = kpf2File, prefix = '<?xml version="1.0" encoding="UTF-8"?>\n<!-- Komodo Project File - DO NOT EDIT -->\n')
	}
	## We don't need .mo files => delete them
	moFiles <- sub("\\.po$", ".mo", poFiles)
	if (any(moFiles != poFiles))
		unlink(moFiles)
	## Check that all the translated .kpf files are produced
	res <- file.exists(kpf2Files)
	names(res) <- basename(kpf2Files)
	return(invisible(res))	
}

kpzTranslate <- function (kpzFile, langs, poFiles, kpz2Files)
{
    if (missing(kpzFile)  || is.null(kpzFile) || is.na(kpzFile))
		stop("'kpzFile' must be provided")
	if (length(kpzFile) != 1)
		stop("You must provide only a single file name/path for 'kpzFile'")
	if (!file.exists(kpzFile))
		stop("The kpzFile is not found")
	proj <- sub("\\.kpz$", "", kpzFile)
	if (missing(poFiles) || is.na(poFiles)) poFiles <- NULL
	if (missing(langs)) {
		if (is.null(poFiles)) {
			## Try to get the list of suitable .po files in same dir as kpfFile
			pattern <- paste(basename(proj), ".+\\.po$", sep = "-")
			poFiles <- dir(dirname(kpzFile), pattern, full.names = TRUE)
			if (length(poFiles) < 1)
				stop("You must provide 'langs' (ex.: 'fr', or 'de'), or 'poFiles'")
		} else langs <- NULL
	}
	if (is.null(poFiles)) {
		if (is.null(langs))
			stop("You must provide 'langs' (ex.: 'fr', or 'de'), or 'poFiles'")
		## Try to guess poFiles from langs
        poFiles <- paste(proj, "-", langs, ".po", sep = "")
	}
	if (any(kpzFile == poFiles))
		stop("'poFiles' cannot be the same as 'kpzFile'")
	if (any(!file.exists(poFiles)))
		stop("One or more 'poFiles' not found!")
	if (missing(kpz2Files)) {
		## Guess kpz2Files from poFiles
		kpz2Files <- sub("\\.po$", ".kpz", poFiles)
	}
	if (any(kpzFile == kpz2Files))
		stop("'kpzFile' and 'kpz2Files' cannot be the same")
	if (any(poFiles == kpz2Files))
		stop("'poFiles' and 'kpz2Files' cannot be the same")
	if (length(poFiles) != length(kpz2Files))
		stop("Number of items must be the same in 'poFiles' and in 'kpz2Files'")
	## Make sure we create new resulting files
	unlink(kpz2Files)
	
    ## The kpz file is a zipped file containing package.kpf in a subdirectory
	pack <- file.path(tempdir(), "package.kpf")
	unlink(pack)	# Make sure the file does not exist yet
	unzip(kpzFile, junkpaths = TRUE, exdir = tempdir())
	if (!file.exists(pack))
		stop("Impossible to extract the content of the .kpz file.")
	## kpf2Files are "package.kpf" files in respective subdirectories with names
	## of the packages, like mypack-fr, mypack-it, etc.
	kpf2Dirs <- file.path(tempdir(), sub("\\.kpz$", "", basename(kpz2Files)))
	kpf2Files <- file.path(kpf2Dirs, "package.kpf")
	
	## Call kpfTranslate on the created package.kpf file
	kpfTranslate(pack, poFiles = poFiles, kpf2Files = kpf2Files)
	
	## Eliminate the temporary extracted package.kpf file
	unlink(pack)
	odir <- getwd()
	on.exit(setwd(odir))
	## Compress the created files in kpz zipped archives named kpz2Files	
	for (h in 1:length(kpz2Files)) {
		kpz2File <- kpz2Files[h]
		kpf2Dir <- kpf2Dirs[h]
		if (file.exists(kpf2Dir)) {
			setwd(dirname(kpz2File))
			## Note: the 'zip' program must be accessible!
			cmd <- paste('zip -rqm9 "', basename(kpz2File), '" "', kpf2Dir, '"',
				sep = "")
			try(system(cmd, intern = TRUE, wait = TRUE), silent = TRUE)
		}
	}
	
	## Check that all the translated .kpz files are produced
	res <- file.exists(kpz2Files)
	names(res) <- basename(kpz2Files)
	return(invisible(res))	
}
