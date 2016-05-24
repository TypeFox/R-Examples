PDF_info <-
function(file)
{
    x <- PDF_doc_from_file(file)
    
    y <- .Call("PDF_info", x$doc)

    ## Convert dates here (could also do so in the C code).
    d <- c("CreationDate", "ModDate")
    y[d] <- lapply(y[d], .POSIXct, tz = "GMT")

    ## Simplify and expand info on page sizes.
    sizes <- matrix(y$Sizes, nrow = 2L)
    afmts <- sprintf("%g x %g pts", sizes[1L, ], sizes[2L, ])
    ufmts <- unique(afmts)
    pos <- match(afmts, ufmts)
    sizes <- sizes[, unique(pos), drop = FALSE]
    unms <- vapply(split(sizes, col(sizes)), PDF_paper_size_name, "")
    ind <- !is.na(unms)
    ufmts[ind] <- sprintf("%s [%s]", ufmts[ind], unms[ind])
    y$Sizes <- ufmts[pos]

    ## PDF version string should start with 'PDF-'.
    y$PDF_Version <- substring(y$PDF_Version, 5L)

    class(y) <- "PDF_info"

    y
}

format.PDF_info <-
function(x, ...)
{
    x <- x[!is.na(x)]
    nms <- names(x)
    labels <- as.list(nms)
    names(labels) <- nms
    pos <- which(nms == "Sizes")
    if(length(pos)) {
        if(length(unique(x[[pos]])) == 1L) {
            labels[[pos]] <- "Page size"
            x[[pos]] <- x[[pos]][1L]
        } else {
            labels[[pos]] <-
                sprintf("Page %s size",
                        format(seq_along(x[[pos]]), justify = "right"))
        }
    }
    labels <- sprintf("%s:", unlist(labels))
    sprintf("%s %s",
            format(labels, justify = "left"),
            unlist(lapply(x, format)))
}

print.PDF_info <-
function(x, ...)
{
    writeLines(format(x, ...))
    if(length(errors <- attr(x, "errors")))
        writeLines(c("\nPDF problems:", sprintf("  %s", errors)))
    invisible(x)
}

PDF_page_sizes <-
do.call(rbind,
        list("A0" =        c(2384L, 3371L),
             "A1" =        c(1685L, 2384L),
             "A2" =        c(1190L, 1684L),
             "A3" =        c( 842L, 1190L),
             "A4" =        c( 595L,  842L),
             "A5" =        c( 420L,  595L),
             "B4" =        c( 729L, 1032L),
             "B5" =        c( 516L,  729L),
             "letter" =    c( 612L,  792L),
             "tabloid" =   c( 792L, 1224L),
             "ledger" =    c(1224L,  792L),
             "legal" =     c( 612L, 1008L),
             "statement" = c( 396L,  612L),
             "executive" = c( 540L,  720L),
             "folio" =     c( 612L,  936L),
             "quarto" =    c( 610L,  780L),
             "10x14" =     c( 720L, 1008L)))
    
PDF_paper_size_name <-
function(size)
{
    w <- size[1L]
    h <- size[2L]
    n <- NA_character_
    pos <- which((abs(PDF_page_sizes[, 1L] - w) < 1) &
                 (abs(PDF_page_sizes[, 2L] - h) < 1))
    if(!length(pos)) {
        pos <- which((abs(PDF_page_sizes[, 2L] - w) < 1) &
                     (abs(PDF_page_sizes[, 1L] - h) < 1))
    }
    ## See e.g. http://en.wikipedia.org/wiki/Paper_size for the exact
    ## tolerances defined in the respective standards.
    if(length(pos))
        n <- rownames(PDF_page_sizes)[pos]
    n
}
