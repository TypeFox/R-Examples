#
#  BIOM format comprises a simple standard for annotation of a two-dimensional
#  matrix:  http://biom-format.org/documentation/format_versions/biom-1.0.html.
#
#  Here, a corresponding "biom" class is implemented to emphasize integration of
#  data and metadata with R functionality, plus convenient import and export.
#  The objects have class(x) == c("biom", "list").  
#
#  The representation is:
#
#		BIOM element				"biom" class element
#       ------------				--------------------
#		id							id
#		format						format
#		format_url					---stored in biom_format_url exported variable---
#		type						type
#		generated_by				generated_by
#		date						date
#  
#		rows						 --
#			id						rownames (data)		OR		sparse $ dimnames [[1]]
#			metadata				rows
#		columns						 --
#			id						colnames (data) 	OR		sparse $ dimnames [[2]]
#			metadata				columns
#  
#		matrix_type					sparse
#		matrix_element_type			---implied by storage.mode() of data---
#		shape						sparse $ dim
#		data						data
#
#		comment						comment
#
#  Notes:
#
#		1) x$data is always "matrix" class
#
#		2) if the BIOM matrix_type is "dense" then 
#				dim(x$data)			is		BIOM shape
#				rownames(x$data)	has		BIOM row ids  (but dimnames is unnamed list)
#				colnames(x$data)	has		BIOM column ids  (but dimnames is unnamed list)
#
#		3) whereas if matrix_type is "sparse" then
#				dim(x$data)"		is 		(# nonzero entries), 3
#				dimnames(x$data)	is		NULL
#				x$sparse$dim		is		BIOM shape
#				x$sparse$dimnames	has		BIOM row and column ids (list of "rows=", "columns=")
#

##############################################################################
##############################################################################
##
##  BASIC CLASS UTILITIES
##
##############################################################################
##############################################################################

#-----------------------------------------------------------------------------
#  str() method to pretty-print object structure.
#  multiple "str" calls are necessary to handle the "list.len" argument properly.
#  important to unclass because ... something.
#-----------------------------------------------------------------------------

str.biom <- function (object, ...) {
	object <- unclass (object)
	str (object [c ("rows", "columns")], vec.len=1, nchar.max=35, list.len=3, 
		max.lev=3, no.list=T, give.attr=F)
	jj <- c("data", 
		if (!is.null (object$sparse)) "sparse", 
		"type", 
		"format", 
		"date", 
		"id", 
		"generated_by", 
		if (!is.null (object$comment)) "comment")
	str (object [jj], vec.len=1, nchar.max=35, no.list=T, give.attr=F)
	}

#-----------------------------------------------------------------------------
#  summary() method to nicely output meta-content (omitting the 'primary' data).
#-----------------------------------------------------------------------------

summary.biom <- function (object, ...) { 
	dd <- dim (object)
	ss <- paste (dd, collapse="x")
	if (!is.null (attr (dd, "nnz")))
		ss <- paste0 (ss, " of which ", attr (dd, "nnz"), " nonzero")
	xx <- with (object, { paste0(
			"[id:] ", id, "\n",

			"[generated_by:] ", generated_by, " on [date:] ", date, "\n",

			"[type:] ", if (exists ("sparse", inherits=F)) "sparse " else "dense ", 
			type, " (", ss, ")\n", 

			"[format:] ", format, "\n",
			
			if (exists ("comment", inherits=F)) paste("[comment:] ", comment, "\n")) })

	class (xx) <- c("biomsummary", "character")

#  returning a classed summary list, rather than a string, would be more proper.

	xx
	}

#-----------------------------------------------------------------------------
#  method to print the summary() object (also, return it invisibly).
#-----------------------------------------------------------------------------

print.biomsummary <- function (x, ...) {
	cat (x)
	invisible (x)
	}

#-----------------------------------------------------------------------------
#  print() method to output full contents (canonically returns object invisibly).
#-----------------------------------------------------------------------------

print.biom <- function (x, ...) {
	print (as.matrix (x, expand=TRUE))
	cat("\n")
	print (summary (x))
	invisible (x)
	}

#-----------------------------------------------------------------------------
#  head() method to show initial segment of data.
#  we avoid writing 1:nrow(x) to allow for 0 rows and/or columns.
#
#  dense: 		existing row and column names will be shown
#  sparse:		rows numbered and columns named by meaning
#-----------------------------------------------------------------------------

head.biom <- function (x, n=5, p=n, ...) {

# allow specifying each dim separately, here and in tail()

	if (is.null (x$sparse)) {
		j <- pmin (dim (x), c(n,p))
		x$data [seq (1, len=j[1]), seq (1, len=j[2]), drop=F]
	} else {
		j <- min (nrow (x$data), n)
		y <- x$data [seq (1, len=j), , drop=F]
		colnames (y) <- c('row', 'col', 'value')
		rownames (y) <- seq (1, len=j)
		y
		}
	}

#-----------------------------------------------------------------------------
#  same
#-----------------------------------------------------------------------------

tail.biom <- function (x, n=5, p=n, ...) {
	if (is.null (x$sparse)) {
		s <- pmax (dim(x) - c(n,p) + 1, c (1,1))
		k <- pmin (dim(x), c(n,p))
		x$data [seq (s[1], len=k[1]), seq (s[2], len=k[2]), drop=F]
	} else {
		s <- max (nrow (x$data) - n + 1, 1)
		k <- min (nrow (x$data), n)
		y <- x$data [seq (s, len=k), , drop=F]
		colnames (y) <- c('row', 'col', 'value')
		rownames (y) <- seq (s, len=k)
		y
		}
	}

#-----------------------------------------------------------------------------
#  dim(x) method returns the BIOM "shape",
#  with additional attribute "nnz" (number not zero) if matrix_type is sparse.
#  "nnz" equals the number of values in the sparse representation.
#  also, the two components of the returns value are named "rows" and "columns".
#-----------------------------------------------------------------------------

dim.biom <- function (x) {
	dd <- with (x, {
		if (exists ("sparse", inherits=FALSE)) {
			attr (sparse$dim, "nnz") <- nrow (data)
			sparse$dim
		} else
			dim (data)
		} )
	names (dd) <- c ('rows', 'columns')
	dd
	}

#-----------------------------------------------------------------------------
#  dimnames() method returns BIOM "ids" in a list of two components,
#  named "rows" and "columns".
#
#  those names need to be added to dimnames(), in the dense case.
#  otherwise, they are already present in sparse$dimnames.
#
#  because rownames() and colnames() are not generic, they are not implemented here.
#-----------------------------------------------------------------------------

dimnames.biom <- function (x) {
	with (x, {
		if (exists ("sparse", inherits=FALSE)) {
			sparse$dimnames
		} else {
			names (dimnames (data)) <- c ("rows", "columns")
			dimnames (data)
		} } )
	}

#-----------------------------------------------------------------------------
#  "metadata() returns BIOM "metadata" in a list of two components,
#  named "rows" and "columns".
#
#  metadata() is a new S3 generic function.
#
#  here as elsewhere, the return value is invisible, to avoid a massive screen dump.
#-----------------------------------------------------------------------------

metadata <- function (x, ...) UseMethod("metadata")

metadata.biom <- function (x, ...) {
	invisible (with (x, 
		list(
			rows = rows,
			columns = columns)))
	}

##############################################################################
##############################################################################
##
##  CONVERSION FROM
##
##############################################################################
##############################################################################

#-----------------------------------------------------------------------------
#  as.matrix() method returns:
#
#		for "sparse", the 3-column sparse representation "matrix" with,
#		also, attributes "rownames" and "colnames" attached, that
#		contain the corresponding BIOM ids.  But, on request, the full
#		"matrix" is returned with zeros filled in.
#
#		for "dense", a "matrix" of the object's data
#
#  the default value of "dense" is not "FALSE" to avoid the implication that
#  a dense matrix might be returned sparse.
#
#  in the biom context, we use the names "rows" and "columns", but here we
#  convert back to an R type, so names(dimnames(.)) of the returned matrix 
#  is NULL.  this is a piddly point.
#
#  return value is invisible to avoid dumping.
#-----------------------------------------------------------------------------

as.matrix.biom <- function (x, expand=NULL, ...) {
	invisible (with (x, {
		if (exists ("sparse", inherits=F)) {
			if (isTRUE (expand)) {
				if (is.character (data)) {
					jj <- data [,1:2]
					storage.mode (jj) <- "integer"
					data [,1:2] <- 1 + jj
				} else
					data [,1:2] <- 1 + data [,1:2]
				data <- sparse2dense (data, sparse$dim)
				dimnames (data) <- sparse$dimnames
				names (dimnames (data)) <- NULL

			} else {
				attr (data, "rownames") <- sparse$dimnames [[1]]
				attr (data, "colnames") <- sparse$dimnames [[2]]
				}
			}

		data } ))
	}

#-----------------------------------------------------------------------------
#  as.list() method reorganizes internally to resemble the BIOM specification.
#  and removes the class attribute.
#
#  in particular, the storage.mode() of "data" is used to assign matrix_element_type,
#  and matrix_type is assigned based on the presence of "sparse"
#
#  the exceptions are:
#		row and column "id"s are made into list elements, not part of "rows" and "columns"
#		"data" is not made into a list, but remains a matrix
#		"format_url" is missing
#
#  return value is invisible.
#-----------------------------------------------------------------------------

as.list.biom <- function (x, ...) {
	y <- list()
	y [c(
		'id',
		'format',
		'type',
		'generated_by',
		'date',
		'rows',
		'columns',
		'row.ids',
		'column.ids',
		'matrix_type',
		'matrix_element_type',
		'shape',
		'data',
		'comment')] <- list(
			x$id,
			x$format,
			x$type,
			x$generated_by,
			x$date,
			x$rows, 
			x$columns,
			dimnames (x) [[1]],
			dimnames (x) [[2]],
			if (is.null (x$sparse)) "dense" else "sparse",
			if (is.integer (x$data)) {
				"int"
			} else if (is.numeric (x$data)) {
				"float"
			} else "unicode",
			as.integer (dim (x)),			# clear "nnz" attribute if present
			x$data,
			x$comment)
	invisible (y)
	}

#-----------------------------------------------------------------------------
#  as.character() converts to JSON text, in a "character" mode object, or file.
#  if, to file, then the text is prettily-formatted, otherwise not.
#
#  try to respect 'unicode' character encoding when writing to a file.
#
#  the filename is returned (visibly) and, if none, then the JSON text (invisibly).
#-----------------------------------------------------------------------------

as.character.biom <- function (x, ..., file=NULL) {

	x <- within (as.list (x), {
		data <- matrix2list (data)
		row.ids <- row.ids
		column.ids <- column.ids
		rows <- mapply (list,
			id = as.list(row.ids),
			metadata = rows, 
			SIMPLIFY = F)
		columns <- mapply (list,
			id = as.list(column.ids), 
			metadata = columns,
			SIMPLIFY=F) 
		format_url <- biom_format_url
		} )
	x$row.ids <- x$column.ids <- NULL

	if (is.null (file)) {
		invisible (RJSONIO::toJSON (x, ...))
	} else {
		if (x$matrix_element_type == "unicode") {
			writeLines (enc2utf8 (RJSONIO::toJSON (x, pretty=TRUE, ...)), file, useBytes=T)
		} else
			writeLines (RJSONIO::toJSON (x, pretty=TRUE, ...), file)
		file
		}
	}


##############################################################################
##############################################################################
##
##  CONSTRUCTORS
##
##############################################################################
##############################################################################

biom <- function (x, ...) UseMethod("biom")

#-----------------------------------------------------------------------------
#  we assume that fromJSON() will return an appropriate list.
#  so biom.list() is called in turn, here.
#-----------------------------------------------------------------------------

biom.character <- function (x, ..., file=NULL, quiet=FALSE) {

#  there is some work to be done here, to use RJSONIO more intelligently
	if (is.null (file)) {
		biom (RJSONIO::fromJSON (x, asText=TRUE, simplify=TRUE, ...), quiet=quiet)
	} else
		biom (RJSONIO::fromJSON (file, simplify=TRUE, ...), quiet=quiet)
	}

#-----------------------------------------------------------------------------
#  from "matrix", just invent something appropriate for all fields.
#
#  argument "sparse" may be given in two formats:
#		c (integer, integer)		indicating matrix dimensions
#		list (character, character)	indicating matrix rownames and colnames
#
#  we insist on assigning unique dimnames, one way or another.
#  empty metadata is created.
#
#  note that constructing from a matrix does not allow specifying:
#		metadata
#		matrix_element_type (explicitly)
#		comment
#
#  storage.mode(x$data) implicitly specifies "matrix_element_type", and is used
#  if/when the object is converted back to JSON.
#-----------------------------------------------------------------------------

biom.matrix <- function (x, type=biom_table_types, sparse=NULL, ..., quiet=FALSE) {
	if (quiet) warning <- function (...) { }

	if (missing (type))
		warning ("unspecified type defaulting to \"", match.arg (type), "\"")

	if (is.null (sparse) && ncol(x) == 3)
		warning ("not interpreting three-column data as a sparse representation")

	y <- list()

	if (is.null (sparse)) {

		if (is.null (rownames (x))) {
			rownames (x) <- rownames (x, do.NULL=F, prefix="row")
		} else if (anyDuplicated (rownames (x)))
			stop ("non-unique rownames")

		if (is.null (colnames (x))) {
			colnames (x) <- colnames (x, do.NULL=F, prefix="column")
		} else if (anyDuplicated (colnames (x)))
			stop ("non-unique colnames")

		names (dimnames (x)) <- NULL

		y$rows <- replicate (nrow(x), list())
		y$columns <- replicate (ncol(x), list())

	} else {
		dimnames (x) <- NULL
		jj <- x [,1:2]
		storage.mode (jj) <- "integer"
		dd <- 1 + apply (jj, 2, max)

		y$sparse <- if (is.list (sparse)) {
				names(sparse) <- c("rows", "columns")
				list(
					dimnames = sparse,
					dim = c (
						length (sparse[[1]]), 
						length (sparse[[2]])))
			} else if (is.numeric (sparse)) {
				list(
					dim = sparse,
					dimnames = list(
						rows = paste0 ("row", 1:sparse[1]),
						columns = paste0 ("column", 1:sparse[2])))
			} else if (isTRUE (sparse)) {
				list(
					dim = dd,
					dimnames = list(
						rows = paste0 ("row", 1:dd[1]),
						columns = paste0 ("column", 1:dd[2])))
			} else
				stop ("bad specification of sparse data")

		if (any (dd > y$sparse$dim)) {
			stop ("sparse data does not fit in provided dimensions")
		} else if (any (y$sparse$dim > dd))
			warning ("sparse data could fit in smaller dimensions than provided")

		y$rows <- replicate (y$sparse$dim[1], list())
		y$columns <- replicate (y$sparse$dim[2], list())
		}		

	y$data			<-		x
	y$type			<-		match.arg (type)

	y$id			<-		""
	y$generated_by	<-		paste("BIOM.utils (", packageVersion("BIOM.utils"), ")", sep="")
	y$date			<-		strftime(Sys.time())
	y$format		<- 		biom_format

	class (y) <- c ("biom", "list")
	invisible (y)
	}

#-----------------------------------------------------------------------------
#  "list" method constructs using "matrix" method, then augments with any 
#  other provided data.
#
#  matrix_element_type, if present, is applied, superceding storage.mode(x$data).
#  so, note the storage.modef() may be integer, double, or character.
#  when the data is sparse and character, this is slightly strange, but still ok.
#-----------------------------------------------------------------------------

biom.list <- function (x, ..., quiet=FALSE) { 
	if (quiet) warning <- function (...) { }

	if (is.null (x$data)) {
		mm <- matrix (0, 0, 0)		# providing first "0" gives mode "numeric" rather than "logical"
	} else if (is.matrix (x$data)) {
		mm <- x$data
	} else {
		mm <- simplify2array (x$data)
		if (is.vector (mm)) {			#  whoops, one-column matrix oversimplified to vector
			mm <- as.matrix (mm)
		} else
			mm <- t(mm)
		}

	if (!is.null (x$matrix_element_type))
		storage.mode (mm) <- switch (x$matrix_element_type,
			int = "integer",
			float = "double",
			unicode = "character")

	rownames <- colnames <- NULL
	val.f <- function (x) all (c ('id', 'metadata') %in% names (x))	

	if (!is.null (x$rows)) {
		if (all (sapply (x$rows, val.f))) {	
			rownames <- sapply (x$rows, "[[", "id")
			x$rows <- lapply (x$rows, "[[", "metadata")
		} else
			warning ("ignoring malformed \"rows\"")
		}

	if (!is.null (x$columns)) {
		if (all (sapply (x$columns, val.f))) {
			colnames <- sapply (x$columns, "[[", "id")
			x$columns <- lapply (x$columns, "[[", "metadata")
		} else
			warning ("ignoring malformed \"columns\"")
		}

	if (isTRUE (x$matrix_type == "sparse")) {

		if (is.null (x$shape))
			stop ("\"shape\" required for sparse data")
		if (is.null (rownames))
			rownames <- paste0 ("row", 1:x$shape[1])
		if (is.null (colnames))
			colnames <- paste0 ("column", 1:x$shape[2])
		sparse <- list (rownames, colnames)

	} else {

		rownames (mm) <- 
			if (is.null (rownames)) {
				rownames (mm, do.NULL=F, prefix="row")
			} else rownames
		colnames (mm) <-
			if (is.null (colnames)) {
				colnames (mm, do.NULL=F, prefix="column")
			} else colnames
		sparse <- NULL
		}

	ll <- list (x=mm, sparse=sparse, quiet=quiet)
	ll$type <- x$type
	y <- do.call (biom, ll)					# allow "type" to be "missing" in the call

	if (!is.null (x$rows)) 				y$rows 			<- x$rows
	if (!is.null (x$columns)) 			y$columns 		<- x$columns
	if (!is.null (x$id)) 				y$id 			<- x$id
	if (!is.null (x$generated_by)) 		y$generated_by 	<- x$generated_by
	if (!is.null (x$date)) 				y$date 			<- x$date
										y$comment		<- x$comment
	invisible (y)
	}


##############################################################################
##############################################################################
##
##  PACKAGE UTILITIES
##
##############################################################################
##############################################################################

#-----------------------------------------------------------------------------
#  This routine builds the package example data:
#
#  BIOM.utils/data/examples.rda:
#		jtxt -- complete JSON text for BIOM object
#		smat -- sparse matrix of data, represented in three columns
#		dmat -- dense matrix, with row and column names
#		li1 -- minimal list, with data given as dense matrix with dimnames
#		li2 -- short list, data as above, but also with explicit row and column records
#		li3 -- short list as above, but with data formatted as list of rows
#		li4 -- complete list of BIOM components
#
#  BIOM.utils/inst/extdata/example-file.txt:
#		JSON text for a biom object
#
#-----------------------------------------------------------------------------

buildBiomExamples <- function (rdafile="examples.rda", jsonfile="example-json.txt") {
	require("MGRASTer")

	triple <- function (x) paste(x, x, x, sep="")

	jtxt <- iconv(
		readLines(
			MGRASTer::call.MGRAST (issue=FALSE, 'ma', 'or', id=c(4447943.3, 4447192.3, 4447102.3, 4447103.3), 
				gro='family', so='Ref', resu='ab', ev=15, quiet=TRUE),
			warn=FALSE),
		to="ASCII",
		sub="?")
	writeLines(jtxt, jsonfile)
	message ("Built ", jsonfile, " in: ", getwd(), 
		".  For package build, move to BIOM.utils/inst/extdata")

	dmat <- matrix(101:200, nrow=20, dimnames=list(letters[1:20], LETTERS[1:5]))
	li1 <- list(
		data=dmat,
		type="OTU")
	li2 <- list(
		data = dmat,
		type = "OTU table",
		rows = lapply(triple(rownames(dmat)),
			function (x) c(id=x, metadata=paste("metadata of", x))),
		columns = lapply(triple(colnames(dmat)),
			function (x) c(id=x, metadata=paste("metadata of", x))))
	li3 <- list(
		data=lapply(apply(unname(dmat), 1, list), `[[`, 1),
		type="OTU table",
		rows=li2$rows,
		columns=li2$columns)
	li4 <- RJSONIO::fromJSON (jtxt)
	li4 [c("matrix_element_value", "url")] <- NULL
	smat <- t(simplify2array(li4$data))

	save(smat, dmat, li1, li2, li3, li4, jtxt, file=rdafile)
	message ("Built ", rdafile, " in: ", getwd(), 
		".  For package build, move to BIOM.utils/data")

	biom_format <- "Biological Observation Matrix 1.0"
	biom_format_url <- "http://biom-format.org/documentation/format_versions/biom-1.0.html"
	names(biom_format_url) <- biom_format
	biom_fields <- c(
		"rows",
		"columns",
		"data",
		"shape",
		"matrix_type",
		"matrix_element_type",
		"type",
		"format",
		"format_url",
		"date",
		"id",
		"generated_by")
	biom_table_types <- c(
		"OTU table",
		"Pathway table",
		"Function table",
		"Ortholog table",
		"Gene table",
		"Metabolite table",
		"Taxon table")
	biom_matrix_element_types <- c(
		"int",
		"float",
		"unicode")
	biom_matrix_types <- c(
		"dense",
		"sparse")
	save(biom_format, biom_format_url, biom_fields, biom_table_types, 
		biom_matrix_element_types, biom_matrix_types, file="sysdata.rda")
	message ("Built sysdata.rda in: ", getwd(), 
		".  For package build, move to BIOM.utils/R")
	}

#-----------------------------------------------------------------------------
#  routine simply to apply all biom functions to an object.
#  for testing.
#-----------------------------------------------------------------------------

applyBiomMethods <- function (x) {
	cat ("======= str ============================\n")
	str(x)
	cat ("======= summary ========================\n")
	print (summary(x))
	cat ("======= print ==========================\n")
	print(x)
	cat ("======= head ===========================\n")
	print (head(x))
	cat ("======= tail ===========================\n")
	print (tail(x))
	cat ("======= dim ============================\n")
	print (dim(x))
	cat ("======= dimnames =======================\n")
	print (dimnames(x))
	cat ("======= metadata =======================\n")
	print (metadata(x))
	cat ("======= as.matrix ======================\n")
	print (as.matrix(x))
	cat ("======= as.matrix ======================\n")
	print (as.matrix(x, expand=T))
	cat ("======= as.list ========================\n")
	print (as.list(x))
	cat ("======= as.character ===================\n")
	print (as.character(x))
	}

#-----------------------------------------------------------------------------
#  return the example file of JSON text.  biom.character() may be applied to it.
#-----------------------------------------------------------------------------

exampleBiomFile <- function () {
	file.path (path.package ("BIOM.utils"), "extdata", "example-json.txt")
	}


##############################################################################
##############################################################################
##
##  LITTLE STUFF  (unexported)
##
##############################################################################
##############################################################################


#-----------------------------------------------------------------------------  
#  dense matrix to sparse
#-----------------------------------------------------------------------------  

dense2sparse <- function (x) {
	if (is.character (x)) {
		ss <- which (nchar(x) & !is.na(x), arr.ind=TRUE)
	} else
		ss <- which (x != 0, arr.ind=TRUE)
	cbind (unname (ss), unname (x) [ss])
	}

#-----------------------------------------------------------------------------  
#  sparse matrix to dense.  allowed to specify dimensions for result.
#
#  there is some subtlety here regarding the mode/type of the matrix
#  we need the routine to accept & preserve the matrix as "integer", "double" or 
#  "character", whence as.integer() below in assigning dimensions.
#  also the result value begins as NA, is promoted to whatever type by assignment
#  of the given data, and then is not change wrt type by replacing na's with zero.
#-----------------------------------------------------------------------------  

sparse2dense <- function (x, dim=NULL) {
	x <- unname (x)
	fill <- if (is.character (x)) "" else 0

	jj <- x [,1:2]
	storage.mode (jj) <- "integer"
	if (is.null (dim)) {
		if (nrow (x)) {
			dim <- apply (jj, 2, max)
		} else
			stop ("cannot omit dimensions for an empty sparse representation")
		}
	y <- matrix (fill, dim[1], dim[2])

	y [jj] <- x [,3]
	y
	}

#-----------------------------------------------------------------------------  
#  matrix to list of its rows (or columns).   a limited inverse of simplify2array().
#-----------------------------------------------------------------------------  

matrix2list <- function (x) {
	lapply(
		as.list (seq (1, len=nrow(x))),
		function (i, x) x[i,],
		unname (x))
	}

#-----------------------------------------------------------------------------
#  we give warnings in a slightly idiosyncratic way
#-----------------------------------------------------------------------------

warning <- function (...) base::warning ("BIOM.utils: ", ...,  call.=FALSE)

#-----------------------------------------------------------------------------
#  concatenate things with spaces between
#-----------------------------------------------------------------------------

collapse <- function (x) paste (x, collapse=" ", sep="")

#-----------------------------------------------------------------------------
#  to report github commit at startup, preprocess this source with:
#    sed s/XXXBUILDXXX/$commit/g matR/R/init.R > init.Rtemp
#    mv init.Rtemp matR/R/init.R
#-----------------------------------------------------------------------------

.onAttach <- function (libname, pkgname) { 
	if ("package:biom" %in% search())
		warning ("package \"biom\" is not compatible")
	ss <- " dbcb27"
	if (substr (ss, 2, 9) == "XXXBUILD") ss <- ""
	packageStartupMessage(pkgname, " (", packageVersion(pkgname), ss, ")")
	}
