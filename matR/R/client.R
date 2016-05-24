
#-----------------------------------------------------------------------------------------
#  "Client" routines:  for data transfer using MG-RAST API.
#
#  General goals for API wrapping:
#    restricting access to API functionality is ok...
#    ...but it must be for the sake of creating clear concepts
#    give useful defaults to users
#    an R-intuitive interface
#    clear relationships between R terminology (concepts) and API terminology (concepts)
#    leave the API room to grow:
#		build in as little dependency on its current version, as possible
#-----------------------------------------------------------------------------------------


biomRequest <- function (x, request=c("function", "organism", "feature"), ..., 
	block, wait=TRUE, quiet=FALSE, file, outfile) {
#-----------------------------------------------------------------------------------------
#  Post and fulfill data requests.  Important capabilities:
#    file containing IDs only
#    file of IDs and metadata
#    specification of IDs as string or character vector
#    specification of IDs as data.frame with metadata
#    block, file output
#    varying API parameters
#
#  for "..." arguments see:
#    doc.MGRAST(2, head=c('matrix','organism','parameters','options'))
#    doc.MGRAST(2, head=c('matrix','function','parameters','options'))
#    doc.MGRAST(2, head=c('matrix','feature','parameters','options'))
#
#    asynchronous 	source 		result_type 	filter 			group_level 	grep 	length 
#    evalue 		identity	 filter_source 	hide_metadata 	id 				filter_level
#-----------------------------------------------------------------------------------------
	if (missing (request)) {
		request <- match.arg(request)
		warning ("\'request\' is defaulting to \'", request, "\'")
	} else
		request <- match.arg(request)

	if (!missing (file)) x <- readSet (file)
	x <- expandSet (x)
	add.metadata <- NULL
	if (is.data.frame (x)) {
		add.metadata <- x
		x <- rownames(x)									# look out for reordering!
	}

	if (missing (block)) block <- length(x)

	req <- new.env (parent = globalenv())
	param <- append (
		list (
			resource='matrix',
			request=request,
			asynchronous=1,
			verify=FALSE), 									# avoid misleading warings
		list(...))

	ledger <- data.frame (start = seq(1, length(x), block), stringsAsFactors=FALSE)
	ledger$stop 		<- c (ledger$start[-1] - 1, length(x))
	ledger$requested 	<- FALSE
	ledger$ticket		<- ""
	ledger$file			<- ""

	yy <- do.call (call.MGRAST, c (param, list (id = x [ledger[1,"start"]:ledger[1,"stop"]])))
	if ("id" %in% names (yy)) {
		ledger[1,"ticket"] <- yy ["id"]
		ledger[1,"requested"] <- TRUE
	} else {
		print (yy)
		stop ("can\'t interpret server response")
		}

	assign("IDs", x, req)
	assign("n", 1, req)
	assign("ledger", ledger, req)
	assign("param", param, req)
 	assign("add.metadata", add.metadata, req)
 	assign("outfile", if (missing (outfile)) NULL else outfile, req)

	if (wait) {
		biom(req, wait=TRUE, quiet=quiet)
	} else {
		message("returning ticket for queued request; apply biom() to fulfill")
		invisible(req)
		}
	}


biom.environment <- function (x, wait=TRUE, ..., quiet=FALSE) {
#-----------------------------------------------------------------------------------------
# data retrieval: function to fulfill requests
#
# "..." in prototype is just good practice for generics
# !quiet=TRUE is used to report on the download in complete detail
#-----------------------------------------------------------------------------------------
	assign("quiet", quiet, x)
	assign("wait", wait, x)
	with (x, {
		repeat {
			zz <- call.MGRAST("status", "instance", id=ledger[n, "ticket"], verify=FALSE)
			if("data" %in% names(zz)) {								# request may have arrived
				yy <- try (biom(zz$data))
				if (inherits (yy, "try-error")) {					# some other message from API
					print (zz$data)
					stop ("can\'t interpret server response as BIOM data")
					}
				tt <- tempfile()
				save (yy, file=tt)
				ledger[n, "file"] <- tt
				if (!quiet) print (ledger [n,])
				if (n == nrow(ledger)) break

				n <- n+1
				yy <- do.call (call.MGRAST, c (param, list (id = IDs [ledger[n,"start"]:ledger[n,"stop"]])))
				if ("id" %in% names (yy)) {
					ledger[n,"ticket"] <- yy ["id"]
					ledger[n,"requested"] <- TRUE
				} else {
					print (yy)
					stop ("can\'t interpret API response")
					}
				}
			if (!wait) {
				warning ("retrieval in blocks with wait=FALSE only checks one block per call")
				return()
				}
			Sys.sleep(5)
			}

		if (!quiet && nrow (ledger) > 1) message("assembling...")
		ll <- lapply(
				ledger[,"file"], 
				function (ff) { 
					load(ff); unlink(ff); yy
					})
		yy <- suppressWarnings (Reduce (merge.biom, ll))			# suppress dup row warning (not ideal)
		if (!quiet && nrow (ledger) > 1) message("done")

#		if(!is.null("add.metadata")) columns(yy) <- add.metadata
		if(!is.null(outfile)) writeLines(as.character(yy), file=outfile)
		invisible(yy)
		})
	}


metadata.character <- function (x, detail=NULL, ..., quiet=TRUE, file) {
#-----------------------------------------------------------------------------------------
# get metadata without data.
#
# detail = NULL			for projects, returns metagnomes		(list, because one-to-many)
#						for metagenomes, returns projects		(named vector, because one-to-one)
#
# detail = TRUE			equivalent to "verbosity=minimal"		(data.frame, for both samples and projects)
#
# detail = c("minimal","verbose","full")				for projects		(data.frame)
# 		 = c("minimal","metadata","stats","full")		for metagenomes		(data.frame)
# relayed directly as "verbosity" to call.MGRAST()
#
# suppressWarnings() is necessary here since API doesn't return everything it says it will
#-----------------------------------------------------------------------------------------
	if (!missing (file)) x <- readSet(file)
	x <- scrubSet(x)
	y <- scrapeSet(x) [1]

	if (is.null(detail) && y=="project") {
		f <- function (x) simplify2array (suppressWarnings (call.MGRAST (					# keep IDs, drop URLs
			"project", "instance", id=x, verbosity='full', quiet=quiet)) $ metagenomes,
			higher=FALSE) [1,]
		sapply (x, f, simplify=FALSE)									# sapply retains names

	} else if (is.null(detail) && y=="metagenome") {
		f <- function (x) suppressWarnings (call.MGRAST (
			"metagenome", "instance", id=x, verbosity='min', quiet=quiet)) $ project [1]
		sapply (x, f)

	} else {
		if (isTRUE (detail)) detail <- "minimal"
		f <- function (x) suppressWarnings (call.MGRAST (
			y, "instance", id=x, verbosity=detail, ..., quiet=quiet))
		z <- list2df (sapply (x, f, simplify=FALSE))
		z$id <- z$url <- NULL											# get rid of some annoying things
		z$library2 <- z$project2 <- z$sample2 <- NULL					# same, for metagenome metadata only
		z
		}
	}


dir.MGRAST <- function (from, to, length.out=0, ..., quiet=TRUE) {
#-----------------------------------------------------------------------------------------
#  here we translate just a bit:
#  arg names "from", "to", "length.out" are familiar to R people,
#  as are indices starting at 1.
#
#  for "..." arguments see:
#    doc.MGRAST(3, head=c('project','query','parameters','options'))}
#
#  verbosity = c("minimal", "verbose", "full")
#  order = c("id", "name")
#  limit = ["integer"]
#  offset = ["integer"]
#
#  --> allow "from" and "to" to contain names
#-----------------------------------------------------------------------------------------
	if (missing(from) && missing(to)) {								# length.out given or ==0
		from <- 1
	} else if (missing(from) && length.out) {						# to and length.out given
		from <- to - length.out + 1
	} else if (missing(from))	{									# to given, only
		from <- 1
		length.out <- to
		}
	else if (!missing(to))
		length.out <- to - from + 1
	args <- resolve (list (...), list(
		resource = "project",
		request = "query",
		verbosity = "minimal",
		order = "name",
		limit = length.out,
		offset = from - 1))
	y <- list2df (do.call (call.MGRAST, args) $ data)

####  make a data.frame.
####  content will vary according to "verbosity".
####  remove certain junk fields, and turn "status" (public/private) into a factor.
####  rownames of the data.frame will always be "mgpXX"

	rownames(y) <- y$id
	y$id <- y$created <- y$url <- y$version <- NULL
	y$status <- as.factor (y$status)
	y
	}


search.MGRAST <- function (public=NULL, detail=NULL, match.all=TRUE, ..., quiet=TRUE) {
#-----------------------------------------------------------------------------------------
#  for "..." arguments see:
#    doc.MGRAST(3, head=c('metagenome','query','parameters','options'))}
#
#  verbosity = c("minimal","mixs","metadata","stats","full")
#  status = c("both","public","private")
#  match = c("all","any")
#  offset = [integer]
#  limit = [integer]
#  order = [string]
#  direction = c("asc","desc")
#  ----
#  function		"search parameter: query string for function"
#  metadata		"search parameter: query string for any metadata field"
#  md5				"search parameter: md5 checksum of feature sequence"
#  organism		"search parameter: query string for organism"
#
#  --> would want no limit on search results, right?
#  --> maybe retrieve in chunks & provide assembly)
#-----------------------------------------------------------------------------------------
	args <- resolve (list (...), list(
		resource = "metagenome",
		request = "query",
		verbosity = if (is.null (detail)) "minimal" else if (isTRUE (detail)) "metadata" else detail,
		status = if (is.null (public)) "both" else if (isTRUE (public)) "public" else "private",
		match = if (match.all) "all" else "any",
		offset = 0,
		limit = 50))
	y <- list2df (do.call (call.MGRAST, args) $ data)
	if (args$limit && nrow(y) == args$limit)
		message ("limit reached; more results may be available")
	y
 	}
