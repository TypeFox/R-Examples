
#-----------------------------------------------------------------------------------------
#  methods extending "biom" class from package BIOM.utils.
#-----------------------------------------------------------------------------------------


"dimnames<-.biom" <- function (x, value) {
#-----------------------------------------------------------------------------------------
#  rownames, colnames are equivalent to the biom "id" field
#  
#  rownames(), colnames() are not generic.
#  but BIOM.utils effectively defines them for "biom" by virtue of dimnames(), which is.
#  here we do the same for rownames<-(), colnames()<- with dimnames<-()
#
#  maybe add safety check that none are duplicated
#-----------------------------------------------------------------------------------------
	within (x, {
		if (exists ("sparse", inherits=FALSE)) {
			sparse$dimnames <- value
			names (sparse$dimnames) <- c ("rows", "columns")
		} else
			dimnames(data) <- value } )
	x
	}


#-----------------------------------------------------------------------------------------
#  rows,columns --- based on "metadata" fields
#
#  in each case, the returned data.frame has
#     IDs in rownames
#     metadata in colnames
#
#  -->should return vector rather than data.frame of one column??			prob not.
#  -->these will need to change slightly when biom is reimplemented
#  -->could rewrite with list2df()
#
#  produce a list of character vectors
#  produce a corresponding list of logical index vectors
#  select the matching elements; this is essentially the desired data, but its shape may be ragged...
#  ...because every row may be structured differently
#  so gather all names of matching fields
#  and make the data rectangular, using indexing to add NA's where needed.
#  sapply() returned matrix, or vector in case of a single metadata field
#  in latter case, must recover its (lost) name,
#  before constructing the data.frame
#-----------------------------------------------------------------------------------------

rows <- function (x, pattern="*") {
	ll <- lapply (x$rows, unlist)
	ii <- lapply (ll, function (vv, p) grepl (p, names(vv)), pattern)
	yy <- mapply (`[`, ll, ii, SIMPLIFY=FALSE)
	rr <- sort (Reduce (union, lapply (yy, names)))
	ss <- sapply (yy, function (y) { y <- y [rr]; names(y) <- rr; y })
	if (is.vector (ss))
		ss <- matrix (ss, nrow=1, dimnames = list (rr, NULL))
	as.data.frame (t(ss), row.names = rownames (x))
	}


columns <- function (x, pattern="*") {
	ll <- lapply (x$columns, unlist)
	ii <- lapply (ll, function (vv, p) grepl (p, names(vv)), pattern)
	yy <- mapply (`[`, ll, ii, SIMPLIFY=FALSE)
	rr <- sort (Reduce (union, lapply (yy, names)))

	ss <- sapply (yy, function (y) { y <- y [rr]; names(y) <- rr; y })
	if (is.vector (ss))
		ss <- matrix (ss, nrow=1, dimnames = list (rr, NULL))

	as.data.frame (t(ss), row.names = dimnames (x) [[2]])
	}


#-----------------------------------------------------------------------------------------
#  metadata is "append-only" so the replacement functions are not typical:
#    rows (xx, "rating") <- rating.list
#-----------------------------------------------------------------------------------------

insertHelper <- function (xx, yy, name) {
	xx <- as.list (xx)
	xx [[name]] <- yy
	xx
	}

"rows<-" <- function (x, name, value) {
	x$rows <- mapply (insertHelper, x$rows, as.list (value), MoreArgs = list (name), SIMPLIFY=FALSE)
	x
	}

"columns<-" <- function (x, name, value) {
	x$columns <- mapply (insertHelper, x$columns, as.list (value), MoreArgs = list (name), SIMPLIFY=FALSE)
	x
	}


`[.biom` <- function (x, i, j, ...) {
#-----------------------------------------------------------------------------------------
#  subsetting allows indexing by logical, character (i.e., dimnames), numeric.
#  note that indexing can be used to reorder as usual (needs testing).
#-----------------------------------------------------------------------------------------
	m <- as.matrix (x, expand=TRUE) [i, j, drop=FALSE]
	x$rows <- x$rows [match (rownames(m), rownames(x))]
	x$columns <- x$columns [match (colnames(m), colnames(x))]
	x$date <- strftime (Sys.time())
	x$id <- paste0 ("derived with ", deparse (match.call(), width.cutoff=500))
	x$generated_by <- tagline()

	if (is.null (x$sparse)) {
		x$data <- m
	} else {
		x$sparse <- list(
			dim = dim(m),
			dimnames = list(
				rows = rownames(m),
				columns = colnames(m)))
		x$data <- dense2sparse(m)
		x$data [,1:2] <- x$data [,1:2] - 1
		}
	x
	}


merge.biom <- function (x, y, ...) {
#-----------------------------------------------------------------------------------------
#  requires here that columns are distinct; rows may not be
#
#  we require unique column names (biom ids)
#  we identify rows with matching names (biom ids)
#  merge the matrices on rownames
#  ensure column order (merge might actually guarantee this)
#  ...or sort columns lex?
#  sort rows lexicographically; seems like a good idea
#  ---> replace new NA's with ...anything?
#
#  construct new rows
#    list elements with ID matching IDs of rows of A <- copy from A
#    remainder: copy from B
#
#  whereas column metadata from the two objects can be simply combined (as below)
#  row order may be is not guaranteed
#-----------------------------------------------------------------------------------------
	new.column.names <- c (colnames(x), colnames(y))
	if (anyDuplicated (new.column.names))
		stop("merge prevented by duplicated columns")
	if (anyDuplicated (c(rownames(x), rownames(y))))
		warning ("merging with rows in common takes row metadata from \'x\'")
	if (x$type != y$type)
		warning ("merging different \'type\'s forces common \'type\'")

	nc <- length (new.column.names)
	new.columns <- append (x$columns, y$columns)

	new.row.names <- union (rownames (x), rownames (y))
	nr <- length (new.row.names)
	new.rows <- vector ("list", nr)
	names (new.rows) <- new.row.names					# names are convenient here,
	new.rows [rownames (y)] <- y$rows					# not required by biom format
	new.rows [rownames (x)] <- x$rows

####  unsparse merge
####    create a new matrix of the right size
####    and fill in data from x and y, both expanded

	if (is.null (x$sparse) || is.null (y$sparse)) {
		if (!is.null (x$sparse) || !is.null (y$sparse))
			warning ("\'sparse\' expanded to merge with \'dense\'")
		sparse <- NULL
		mm <- matrix(
			0, nrow=nr, ncol=nc,
			dimnames=list(
				new.row.names,
				new.column.names))
		mm [rownames(x), colnames(x)] <- as.matrix (x$data, TRUE)
		mm [rownames(y), colnames(y)] <- as.matrix (y$data, TRUE)
				
####  sparse merge:
####    update row indices (nonzero entries from both x and y)
####    update col indices (nonzero entries from y only)

		} else {
			sparse <- list (new.row.names, new.column.names)
			x$data [,1] <- match (rownames (x), new.row.names) [1 + x$data [,1]] - 1
			y$data [,1] <- match (rownames (y), new.row.names) [1 + y$data [,1]] - 1
			y$data [,2] <- ncol (x) + y$data [,2]
			mm <- rbind (x$data, y$data)
			}
	zz <- biom (mm, x$type, sparse)
	zz$rows <- new.rows
	zz$columns <- new.columns
	zz$id <- paste0 ("derived with ", deparse (match.call(), width.cutoff=500))
	zz$generated_by <- tagline()
	zz
	}
