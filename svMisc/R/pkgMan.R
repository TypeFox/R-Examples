pkgManDescribe <- function (pkgName, print.it = TRUE)
{
	owarn <- getOption("warn")
	on.exit(options(warn = owarn))
	options(warn = -1)
	desc <- packageDescription(pkgName)
	options(warn = owarn)
	if (is.na(desc)) {
		## Package is apparently not installed... Try getting data from CRAN
		con <- url(file.path(getOption("repos")['CRAN'], "web", "packages",
			pkgName, 'DESCRIPTION', fsep = '/'))
        m <- try(open(con, "r"), silent = TRUE)
        if (!inherits(m, "try-error")) {
			on.exit(close(con))
			dcf <- try(read.dcf(con))
			## Build a 'packageDescription' object
			desc <- as.list(dcf[1, ])
			class(desc) <- "packageDescription"
		} else {
			return(invisible(NULL))
		}
	}
	if (isTRUE(print.it)) {
		write.dcf(as.data.frame.list(desc[!sapply(desc, is.na)],
			optional = TRUE), width = Inf)
		invisible(desc)
	} else desc
}

pkgManGetMirrors <- function ()
{
	## Cache the list of CRAN mirrors in SciViews:TempEnv
	tmpVar <- "pkgMan.CRANmirrors"
	if (existsTemp(tmpVar)) {
		mirrors <- getTemp(tmpVar)
	} else {
		mirrors <- getCRANmirrors()
		assignTemp(tmpVar, mirrors)
	}
	write.table(mirrors[, c("Name", "URL", "CountryCode")],
		row.names = FALSE, col.names = FALSE, sep = ';', quote = FALSE, na = "")
}

pkgManGetAvailable <- function (page = "next", pattern = "", n = 50,
keep = c("Package", "Version", "InstalledVersion", "Status"), reload = FALSE,
sep = ";", eol = "\t\n") {
	
	availablePkgs <- function (avpkg = available.packages(), installed = TRUE) {
		avpkg <- avpkg[order(toupper(avpkg[, "Package"])), , drop = FALSE]
		if (isTRUE(installed)) {
			inspkg <- installed.packages()
			ipkgnames <- unique(inspkg[, 'Package'])
	
			ipkgnames <- ipkgnames[ipkgnames %in% avpkg[, 'Package']]
			avpkg <- cbind(avpkg, InstalledVersion = NA, Status = NA)
			if (length(ipkgnames)) {
				pkgstatus <- sapply(ipkgnames, function (pkg) {
					compareVersion(avpkg[pkg, 'Version'], inspkg[pkg, 'Version'])
				})
				avpkg[ipkgnames, 'Status'] <- pkgstatus
				avpkg[ipkgnames, 'InstalledVersion'] <- inspkg[ipkgnames, 'Version']
			}
		}
		return(avpkg)
	}

	if (!existsTemp('avpkg.list') || isTRUE(reload)) {
		avpkg.list <- availablePkgs(available.packages(filters = c("R_version",
			"OS_type", "duplicates")), installed = FALSE)
		assignTemp('avpkg.list', avpkg.list)
	} else {
		avpkg.list <- getTemp('avpkg.list')
	}
	if (page == "first") {
		newSearch <- TRUE
		i0 <- 1
	} else {
		newSearch <- getTemp('avpkg.pattern', "") != pattern
		i0 <- getTemp('avpkg.idx', default = 1)
	}

	if (is.character(pattern) && pattern != "") {
		if (newSearch) {
			page <- "current"
			i0 <- 1
			idx <- grep(pattern, avpkg.list[,'Package'], ignore.case = TRUE)
			assignTemp('avpkg.pattern.idx', idx)
		} else {
			idx <- getTemp('avpkg.pattern.idx')
		}
		imax <- length(idx)
	} else {
		imax <- nrow(avpkg.list)
		idx <- seq(imax)
	}
	assignTemp('avpkg.pattern', pattern)

	if (page == "next") {
		i0 <- i0 + n
	} else if (page == "prev") {
		i0 <- i0 - n
	}
	outside <- i0 > imax || i0 < 1
	if (outside) return(NULL)
	assignTemp('avpkg.idx', i0)
	i1 <- min(i0 + n - 1, imax)
	i <- seq(i0, i1)
	cat(i0, i1, imax, "\t\n")
	write.table(availablePkgs(avpkg.list[idx[i], , drop = FALSE])[ ,
		keep, drop = FALSE], row.names = FALSE, col.names = FALSE, sep = sep,
		quote = FALSE, eol = eol, na = "")
}

pkgManGetInstalled <- function (sep = ";", eol = "\t\n")
{
	inspkg <- installed.packages(fields = "Description")
	inspkg <- inspkg[order(toupper(inspkg[ , "Package"])),
		c("Package", "Version", "Description")]

	inspkg[, 3] <- gsub("\n", " ", inspkg[, 3])
	inspkg <- cbind(inspkg, Installed = inspkg[, 'Package'] %in% .packages())
	write.table(inspkg, row.names = FALSE, col.names = FALSE, sep = sep,
		quote = FALSE, eol = eol, na = "")
}

pkgManSetCRANMirror <- function (url)
{
	repos <- getOption("repos")
	repos['CRAN'] <- url
	options(repos = repos)
}

pkgManInstall <- function (pkgs, install.deps = FALSE, ask = TRUE)
{
	dep <- suppressMessages(getNamespace("utils")$getDependencies(pkgs,
		available = getTemp('avpkg.list')))
	msg <- status <- ""
	if (!isTRUE(ask) && (isTRUE(install.deps) || all(dep %in% pkgs))) {
		msg <- captureAll(install.packages(dep))
		status <- "done"
	} else {
		l <- length(dep)
		msg <- sprintf(ngettext(l,
			"This will install package %2$s.",
			"This will install packages: %s and %s.",
		), paste(sQuote(dep[-l]), collapse = ", "), sQuote(dep[l]))
		status <- "question"
	}
	list(packages = dep, message = msg, status = status)
}

pkgManRemove <- function (pkgName)
{
	sapply(pkgName, function (pkgName) {
		packSearchName <- paste("package", pkgName, sep = ":")
		if (packSearchName %in% search()) detach(packSearchName,
			character.only = TRUE, unload = TRUE)
		if (pkgName %in% loadedNamespaces()) unloadNamespace(pkgName)
	
		dlli <- getLoadedDLLs()[[pkgName]]
		if (!is.null(dlli)) dyn.unload(dlli[['path']])

		pkgpath <- find.package(pkgName, quiet = TRUE)
		if (length(pkgpath) == 0L) return(FALSE)

		pkglib <- normalizePath(file.path(pkgpath, ".."))
		if (file.access(pkglib, 2) == 0) {
			remove.packages(pkgName, lib = pkglib)
			TRUE
		} else {
			#warning("No sufficient access rights to library", sQuote(pkglib))
			FALSE
		}
	}, simplify = FALSE)
}

pkgManLoad  <- function (pkgName)
{
	sapply(pkgName, library, character.only = TRUE, logical.return = TRUE,
		simplify = FALSE)
}

pkgManDetach <- function (pkgName)
{
	sapply(pkgName, function (pkgName) {
		tryCatch({
			packSearchName <- paste("package", pkgName, sep = ":")
			if (packSearchName %in% search())
				detach(packSearchName, character.only = TRUE, unload = TRUE)
			if(pkgName %in% loadedNamespaces()) unloadNamespace(pkgName)
			TRUE
		}, error = function(e) { conditionMessage(e) })
	}, simplify = FALSE)
}
