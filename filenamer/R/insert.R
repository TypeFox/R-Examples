#' Insert tag or extension into a file name
#' 
#' This function inserts a tag or extension into a file name.
#' It can also replace an element of a file name.
#'
#' By default, tags are inserted at the ultimate position and extensions
#' at the penultimate position, if possible.
#' (That is, the final file extension will not change, unless the insertion
#' position is specified otherwise or the orginal file name had no extension.)
#' If \code{replace} is \code{TRUE}, the tag at the indicated position is 
#' replaced by the new tag instead.
#' 
#' @param x    file name (\code{character} or \code{filename})
#' @param ...  unused arguments
#' @return modified object of the original type
#' @export
#'
#' @examples
#' f <- as.filename("data_expr_2014-05-01.tsv")
#'
#' # new file name with inserted tags for saving normalized data
#' g <- insert(f, tag=c("mod", "norm"))
#' print(as.character(g))
#'
#' # new file name with inserted extension for saving sorted data
#' h <- insert(f, ext="sorted")
#' print(as.character(h))
#'
#' # new file name with different extension for saving in different format
#' i <- insert(f, ext="csv", replace=TRUE)
#' print(as.character(i))
#'
#' # insert another tag
#' j <- insert(g, tag="qc", tag.pos=2)
#' print(as.character(j))
#'
insert <- function(x, ...) UseMethod("insert");

#' @rdname insert
#' @param tag      one or more file name tags to insert
#' @param tag.pos  position at which to insert tag
#'                 (\code{NULL}: append at the end or replace tag)
#' @param ext      one or more file extension tags to insert
#' @param ext.pos  position at which to insert extension
#'                 (\code{NULL}: insert at penultimate position)
#' @param replace  if \code{TRUE}, tag or extension is replaced
#'                 (default: replace last tag)
#' @export
insert.filename <- function(
	x, tag=NULL, tag.pos=NULL, ext=NULL, ext.pos=NULL, replace=FALSE, ...
) {

	if (!is.null(tag)) {
		if (is.null(x$tag)) {
			# set tag
			x$tag <- tag;
		} else {
			n <- length(x$tag);
			if (replace) {
				# replace tag
				if (is.null(tag.pos)) {
					# replace last tag
					tag.pos <- n;
				} else if (tag.pos > n) {
					stop("Cannot replace nonexistent file name tag")
				}
				x$tag[tag.pos] <- tag;
			} else {
				if (is.null(tag.pos)) {
					# insert at end
					x$tag <- c(x$tag, tag);
				} else if (tag.pos == 1) {
					x$tag <- c(tag, x$tag);
				} else if (tag.pos > n) {
					stop("Cannot add file name tag at illegal position")
				} else {
					x$tag <- c(x$tag[1:(tag.pos-1)], tag, x$tag[tag.pos:length(x$tag)]);
				}
			}
		}
	}

	if (!is.null(ext)) {
		if (is.null(x$ext)) {
			# set extension
			x$ext <- ext;
		} else {
			n <- length(x$ext);
			if (replace) {
				if (is.null(ext.pos)) {
					# replace last extension
					ext.pos <- n;
				} else if (ext.pos > n) {
					stop("Cannot replace nonexistent extension tag")
				}
				x$ext[ext.pos] <- ext;
			} else {
				if (is.null(ext.pos)) {
					# insert at ultimate position
					ext.pos <- n;
				}
				if (ext.pos == 1) {
					x$ext <- c(ext, x$ext);
				} else if (ext.pos > n) {
					stop("Cannot add extension tag at illegal position")
				} else {
					x$ext <- c(x$ext[1:(ext.pos-1)],
							ext, x$ext[ext.pos:length(x$ext)]);
				}
			}
		}
	}

	x
}

#' @rdname insert
#' @export
insert.character <- function(x, ...) {
	as.character(insert.filename(as.filename(x), ...))
}
