AlcesteSource <- function(x, encoding = "auto") {
    if(encoding == "auto")
        encoding <- stringi::stri_enc_detect(readBin(x, "raw", 1024))[[1]]$Encoding[1]

    if(is.null(encoding))
        encoding <- ""

    lines <- iconv(readLines(x, warn=FALSE),
                   from=encoding, to="UTF-8", sub="byte")

    newdocs <- grepl("^(\\*\\*\\*\\*|[[:digit:]]+ \\*)", lines)
    content <- split(lines, cumsum(newdocs))

    SimpleSource(encoding, length(content),
                 content=content, uri=x,
                 reader=readAlceste, class="AlcesteSource")
}

getElem.AlcesteSource <- function(x) list(content = x$content[[x$position]], uri = x$URI)
