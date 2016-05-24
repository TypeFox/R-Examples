.mime_part_finalizer <- function(x) {
  if (!is.null(x$file))
    file.remove(x$file)
}

.mime_part <- function(headers, file=NULL, text=NULL) {
  if (!is.null(file) && !is.null(text))
    stop("Can only provide file or text for mime part.")

  e <- environment()
  reg.finalizer(e, .mime_part_finalizer, onexit=TRUE)
  class(e) <- "mime_part"
  e
}

.write_mime_part <- function(mp, con=stdout()) {
  writeLines(paste(names(mp$headers), unlist(mp$headers), sep=": "),
             con, sep="\r\n")
  writeLines("", con, sep="\r\n")
  if (is.null(mp$file))
    writeLines(mp$text, con)
  else
    writeLines(readLines(mp$file), con, sep="\r\n")
}

#' @importFrom base64enc base64encode
.file_attachment <- function(fn, name,
                             type="application/octet-stream",
                             disposition="attachment") {
  if (missing(name))
    name <- basename(fn)

  text <- base64encode(fn, linewidth=72, newline="\n")
  headers <- list("Content-Type"=type,
                  "Content-Disposition"=sprintf("%s; filename=%s",
                    disposition, name),
                  "Content-Transfer-Encoding"="base64")
  
  .mime_part(headers=headers, text=text)
}

.plot_attachment <- function(plt, name=deparse(substitute(plt)), device, ...) {
  fn <- tempfile()
  device(file=fn, ...)
  print(plt)
  dev.off()
  ## FIXME: Guess content type from device!
  res <- .file_attachment(fn, name, type="application/pdf")
  file.remove(fn)
  res
}

##' Create a MIME part
##'
##' @param x Object to include
##' @param name Name of mime part. Usually the filename of the
##'   attachment as displayed by the e-mail client.
##' @param ... Possible further arguments for \code{mime_part}
##'   implementations.
##' @return An S3 \code{mime_part} object.
##' @export
mime_part <- function(x, name, ...)
  UseMethod("mime_part", x)

##' Default MIME part method
##'
##' Creates a string representation of the object \code{x} using
##' \code{dput}. This representation is then turned into a file
##' attachment.
##'
##' @param x R object
##' @param name Filename used for attachment (sans the .R extension)
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##'
##' @method mime_part default
##' @export
mime_part.default <- function(x, name, ...) {
  str <- dput(x)
  .mime_part(headers=list(
               "Content-Type"="text/plain",
               "Content-Disposition"=sprintf("attachment; file=%s.R", name)),
             text=str)
}

##' Creates a MIME part from a trellis plot object
##'
##' Writes a PDF file of the plot defined by \code{x} and turns this
##' PDF file into a file attachment.
##'
##' @param x A \code{trellis} (lattice) object
##' @param name Name of attachment (sans .pdf extension).
##' @param device Graphics device used to render the plot. Defaults to
##'   \code{pdf}.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##' 
##' @method mime_part trellis
##' @export
mime_part.trellis <- function(x, name=deparse(substitute(x)), device=pdf, ...)
  .plot_attachment(x, name=name, device=device, ...)

##' Creates a MIME part from a ggplot2 plot object
##'
##' Writes a PDF file of the plot defined by \code{x} and turns this
##' PDF file into a file attachment.
##'
##' @param x A \code{ggplot} object
##' @param name Name of attachment (sans .pdf extension).
##' @param device Graphics device used to render the plot. Defaults to
##'   \code{pdf}.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##' 
##' @method mime_part ggplot
##' @export
mime_part.ggplot <- function(x, name=deparse(substitute(x)), device=pdf, ...)
  .plot_attachment(x, name=name, device=device, ...)

##' Create a MIME part from a matrix.
##'
##' @param x Matrix
##' @param name Basename of file attachment that is generated.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object
##' 
##' @method mime_part matrix
##' @export
mime_part.matrix <- function(x, name=deparse(substitute(x)), ...) {
  f <- tempfile()
  on.exit(file.remove(f))
  write.table(x, file=f, ...)
  .file_attachment(f, name=sprintf("%s.txt", name), type="text/plain")
}

##' Create a MIME part from a \code{data.frame}.
##' 
##' @param x A \code{data.frame}.
##' @param name Basename of file attachment that is generated.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##'
##' @method mime_part data.frame
##' @export
mime_part.data.frame <- function(x, name=deparse(substitute(x)), ...) {
  f <- tempfile()
  on.exit(file.remove(f))
  write.table(x, file=f, ...)
  .file_attachment(f, name=sprintf("%s.txt", name), type="text/plain")
}

##' Create a MIME part from a character string. If the string matches
##' a filename, a MIME part containing that file is returned instead.
##' 
##' @param x Character string, possibly a filename.
##' @param name Name of attachment.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##' 
##' @method mime_part character
##' @export
mime_part.character <- function(x, name, ...) {
  if (length(x) == 1 && file.exists(x)) {
    .file_attachment(x, name, ...)
  } else {
     .mime_part(headers=list(
                 "Content-Type"="text/plain",
                 "Content-Disposition"="inline"),
               text=paste(x, collapse="\r\n"))
  }
}
