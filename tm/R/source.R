## Author: Ingo Feinerer
## Sources

getSources <-
function()
   c("DataframeSource", "DirSource", "URISource", "VectorSource", "XMLSource",
     "ZipSource")

SimpleSource <-
function(encoding = "",
         length = 0,
         position = 0,
         reader = readPlain,
         ...,
         class)
{
    if (!is.character(encoding))
        stop("invalid encoding")
    if (!is.numeric(length) || (length < 0))
        stop("invalid length entry denoting the number of elements")
    if (!is.numeric(position))
        stop("invalid position")
    if (!is.function(reader))
        stop("invalid default reader")

    s <- list(encoding = encoding, length = length,
              position = position, reader = reader, ...)
    class(s) <- unique(c(class, "SimpleSource", "Source"))
    s
}

# A data frame where each row is interpreted as document
DataframeSource <-
function(x)
    SimpleSource(length = nrow(x), content = x, class = "DataframeSource")

# A directory with files interpreted as documents
DirSource <-
function(directory = ".", encoding = "", pattern = NULL,
         recursive = FALSE, ignore.case = FALSE, mode = "text")
{
    if (!identical(mode, "text") &&
        !identical(mode, "binary") &&
        !identical(mode, ""))
        stop(sprintf("invalid mode '%s'", mode))

    d <- dir(directory, full.names = TRUE, pattern = pattern,
             recursive = recursive, ignore.case = ignore.case)

    if (!length(d))
        stop("empty directory")

    isfile <- !file.info(d)[["isdir"]]
    if (any(is.na(isfile)))
        stop("non-existent or non-readable file(s): ",
             paste(d[is.na(isfile)], collapse = " "))

    SimpleSource(encoding = encoding, length = sum(isfile),
                 mode = mode, filelist = d[isfile], class = "DirSource")
}

# Documents identified by a Uniform Resource Identifier
URISource <-
function(x, encoding = "", mode = "text")
{
    if (!identical(mode, "text") &&
        !identical(mode, "binary") &&
        !identical(mode, ""))
        stop(sprintf("invalid mode '%s'", mode))

    SimpleSource(encoding = encoding, length = length(x), mode = mode, uri = x,
                 class = "URISource")
}

# A vector where each component is interpreted as document
VectorSource <-
function(x)
    SimpleSource(length = length(x), content = x, class = "VectorSource")

XMLSource <-
function(x, parser, reader)
{
    tree <- XML::xmlParse(x)
    content <- parser(tree)
    XML::free(tree)

    SimpleSource(length = length(content), reader = reader, content = content,
                 uri = x, class = "XMLSource")
}

# A ZIP file with its compressed files interpreted as documents
ZipSource <-
function(zipfile, pattern = NULL, recursive = FALSE, ignore.case = FALSE,
         mode = "text")
{
    if (!identical(mode, "text") &&
        !identical(mode, "binary") &&
        !identical(mode, ""))
        stop(sprintf("invalid mode '%s'", mode))

    SimpleSource(exdir = NULL,
                 files = NULL,
                 mode = mode,
                 pattern = pattern,
                 recursive = recursive,
                 ignore.case = ignore.case,
                 zipfile = zipfile,
                 class = "ZipSource")
}

# tau:::read_all_bytes
read_all_bytes <-
function(con, chunksize = 2 ^ 16)
{
    if(is.character(con)) {
        return(readBin(con, raw(), file.info(con)$size))
    }

    if(!isOpen(con)) {
        open(con, "rb")
        on.exit(close(con))
    }

    bytes <- list()
    repeat {
        chunk <- readBin(con, raw(), chunksize)
        bytes <- c(bytes, list(chunk))
        if(length(chunk) < chunksize) break
    }

    unlist(bytes)
}

readContent <-
function(x, encoding, mode)
{
    if (identical(mode, "text"))
        iconv(readLines(x, warn = FALSE), encoding, "UTF-8", "byte")
    else if (identical(mode, "binary"))
        read_all_bytes(x)
    else if (identical(mode, ""))
        NULL
    else
        stop("invalid mode")
}

open.SimpleSource <-
close.SimpleSource <-
function(con, ...)
    con
open.ZipSource <-
function(con, ...)
{
    x <- con
    exdir <- tempfile("ZipSource")
    dir.create(exdir, mode = "0700")

    destfile <- x$zipfile

    if (!file.exists(destfile)) {
        destfile <- tempfile()
        download.file(x$zipfile, destfile)
        on.exit(file.remove(destfile))
    }

    files <- unzip(destfile, list = TRUE)
    ## Directories have length 0
    files <- files[files$Length > 0, "Name"]
    ## Idea: Subdirectories contain file separators
    if (!x$recursive)
        files <- files[!grepl(.Platform$file.sep, files, fixed = TRUE)]
    ## Idea: pattern and ignore.case refer to the file name (like basename)
    ## Cf. also ?dir
    if (!is.null(x$pattern))
        files <- files[grepl(x$pattern, files, ignore.case = x$ignore.case)]

    unzip(destfile, files, exdir = exdir)

    x$exdir <- exdir
    x$files <- files
    x$length <- length(files)
    x
}

close.ZipSource <-
function(con, ...)
{
    x <- con
    if (!is.null(x$exdir)) {
        unlink(x$exdir, recursive = TRUE)
        x$exdir <- NULL
        x$files <- NULL
        x$length <- 0
    }
    x
}

eoi <-
function(x)
    UseMethod("eoi", x)
eoi.SimpleSource <-
function(x)
    x$length <= x$position

getElem <-
function(x)
    UseMethod("getElem", x)
getElem.DataframeSource <-
function(x)
    list(content = x$content[x$position, ],
         uri = NULL)
getElem.DirSource <-
function(x)
{
    filename <- x$filelist[x$position]
    list(content = readContent(filename, x$encoding, x$mode),
         uri = paste0("file://", filename))
}
getElem.URISource <-
function(x)
    list(content = readContent(x$uri[x$position], x$encoding, x$mode),
         uri = x$uri[x$position])
getElem.VectorSource <-
function(x)
    list(content = x$content[x$position],
         uri = NULL)
getElem.XMLSource <-
function(x)
    list(content = XML::saveXML(x$content[[x$position]]),
         uri = x$uri)
getElem.ZipSource <-
function(x)
{
       path <- file.path(x$exdir, x$files[x$position])
       list(content = readContent(path, x$encoding, x$mode),
            uri = paste0("file://", path))
}


length.SimpleSource <-
function(x)
    x$length

pGetElem <-
function(x)
    UseMethod("pGetElem", x)
pGetElem.DataframeSource <-
function(x)
    lapply(seq_len(x$length),
           function(y) list(content = x$content[y,],
                            uri = NULL))
pGetElem.DirSource <-
function(x)
    lapply(x$filelist,
           function(f) list(content = readContent(f, x$encoding, x$mode),
                            uri = paste0("file://", f)))
pGetElem.URISource <-
function(x)
    lapply(x$uri,
           function(uri) list(content = readContent(uri, x$encoding, x$mode),
                              uri = uri))
pGetElem.VectorSource <-
function(x)
    lapply(x$content,
           function(y) list(content = y,
                            uri = NULL))
pGetElem.ZipSource <-
function(x)
    lapply(file.path(x$exdir, x$files),
           function(f) list(content = readContent(f, x$encoding, x$mode),
                            uri = paste0("file://", f)))

reader <-
function(x)
    UseMethod("reader", x)
reader.SimpleSource <-
function(x)
    x$reader

stepNext <-
function(x)
    UseMethod("stepNext", x)
stepNext.SimpleSource <-
function(x)
{
    x$position <- x$position + 1
    x
}
