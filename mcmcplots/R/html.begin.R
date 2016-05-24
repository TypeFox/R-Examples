.html.begin <- function(outdir = tempdir(), filename = "index", extension = "html", title, cssfile) {
    doctype <-
        '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
  <html xmlns="http://www.w3.org/1999/xhtml">
'
    cssfile <- cssfile[cssfile != ""]
    csslink <- if (0 < length(cssfile)) {
        paste('<link rel="stylesheet" type="text/css" href="', cssfile, '" />',
              sep="")
    }

    file <- file.path(outdir, paste(filename, extension, sep="."))
    out <- doctype
    out <- c(out, '<head>')
    out <- c(out, paste('<title>', title, '</title>', sep=""))
    out <- c(out, csslink)
    out <- c(out, '</head>\n<body>')
    out <- paste(out, collapse="\n", sep="")
    cat(out, file=file, append=FALSE)
    invisible(file)
}
