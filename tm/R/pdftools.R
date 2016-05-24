pdf_info_via_xpdf <-
function(file, options = NULL)
{
    outfile <- tempfile("pdfinfo")
    on.exit(unlink(outfile))
    status <- system2("pdfinfo", c(options, shQuote(normalizePath(file))),
                      stdout = outfile)
    ## Could check the status ...
    ## This does not work ...
    ##   info <- as.list(read.dcf(outfile)[1L, ])

    tags <- c("Title", "Subject", "Keywords", "Author", "Creator",
              "Producer", "CreationDate", "ModDate", "Tagged", "Form",
              "Pages", "Encrypted", "Page size", "File size",
              "Optimized", "PDF version" )
    re <- sprintf("^(%s)",
                  paste(sprintf("%-16s", sprintf("%s:", tags)),
                        collapse = "|"))
    lines <- readLines(outfile, warn = FALSE)
    ind <- grepl(re, lines)
    tags <- sub(": *", "", substring(lines[ind], 1L, 16L))
    info <- split(sub(re, "", lines), cumsum(ind))
    names(info) <- tags
    fmt <- "%a %b %d %X %Y"
    if(!is.null(d <- info$CreationDate))
        info$CreationDate <- strptime(d, fmt)
    if(!is.null(d <- info$ModDate))
        info$ModDate <- strptime(d, fmt)
    if(!is.null(p <- info$Pages))
        info$Pages <- as.integer(p)
    info
}

pdf_info_via_gs <-
function(file)
{
    file <- normalizePath(file)

    gs_cmd <- tools::find_gs_cmd()

    out <- system2(gs_cmd,
                   c("-dNODISPLAY -q",
                     sprintf("-sFile=%s", shQuote(file)),
                     system.file("ghostscript", "pdf_info.ps",
                                 package = "tm")),
                   stdout = TRUE)
    out <- out[cumsum(out == "") == 2L][-1L]
    val <- sub("^[^:]+:[[:space:]]*", "", out)
    names(val) <- sub(":.*", "", out)
    val <- as.list(val)
    if(!is.null(d <- val$CreationDate))
        val$CreationDate <- PDF_Date_to_POSIXt(d)
    if(!is.null(d <- val$ModDate))
        val$ModDate <- PDF_Date_to_POSIXt(d)

    val
}

PDF_Date_to_POSIXt <-
function(s)
{
    ## Strip optional 'D:' prefix.
    s <- sub("^D:", "", s)
    ## Strip apostrophes in offset spec.
    s <- gsub("'", "", s)
    if(nchar(s) <= 14L) {
        s <- sprintf("%s%s", s,
                     substring("    0101000000", nchar(s) + 1L, 14L))
        strptime(s, "%Y%m%d%H%M%S")
    } else if(substring(s, 15L, 15L) == "Z") {
        strptime(substring(s, 1L, 14L), "%Y%m%d%H%M%S")
    } else {
        strptime(s, "%Y%m%d%H%M%S%z")
    }
}

pdf_text_via_gs <-
function(file)
{
    files <- normalizePath(file)

    gs_cmd <- tools::find_gs_cmd()

    tf <- tempfile("pdf")
    on.exit(unlink(tf))

    ## The current mechanism is first converting PDF to Postscript using
    ## the ps2write device, and then extract text using the ps2ascii.ps
    ## program.  This fails for some files (e.g.,
    ## /data/rsync/PKGS/AlleleRetain/inst/doc/AlleleRetain_User_Guide.pdf
    ## which Ghostscript also fails to render.  Note that rendering via
    ## gv works "fine": but this uses the pswrite device which produces
    ## bitmap (from which no text can be extracted, of course).
    ## Using the txtwrite device is simply too unstable: e.g.,
    ##   gs -dBATCH -dNOPAUSE -sDEVICE=txtwrite -dQUIET -sOutputFile=- \
    ##     /data/rsync/PKGS/AlleleRetain/inst/doc/AlleleRetain_User_Guide.pdf
    ## keeps segfaulting.
    ## An additional nuisance is that there seems no simple way to
    ## detect a ps2ascii.ps failure.

    ## Finally, note that we currently use -DSIMPLE: without this, more
    ## information would be made available, but require post-processing.

    ## Step 1.  Convert PDF to Postscript.
    res <- system2(gs_cmd,
                   c("-q -dNOPAUSE -dBATCH -P- -dSAFER -sDEVICE=ps2write",
                     sprintf("-sOutputFile=%s", tf),
                     "-c save pop -f",
                     shQuote(file)))
    ## Step 2.  Extract text.
    txt <- system2(gs_cmd,
                   c("-q -dNODISPLAY -P- -dSAFER -dDELAYBIND -dWRITESYSTEMDICT -dSIMPLE",
                     "-c save -f ps2ascii.ps",
                     tf,
                     "-c quit"),
                   stdout = TRUE)
    ## Argh.  How can we catch errors?
    ## The return values are always 0 ...
    if(any(grepl("Error handled by opdfread.ps", txt))) {
        stop(paste(c("Ghostscript failed, with output:", txt),
                   collapse = "\n"))
    }

    strsplit(paste(txt, collapse = "\n"), "\f")[[1L]]
}
