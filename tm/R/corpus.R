# Author: Ingo Feinerer

PCorpus <-
function(x,
         readerControl = list(reader = reader(x), language = "en"),
         dbControl = list(dbName = "", dbType = "DB1"))
{
    stopifnot(inherits(x, "Source"))

    readerControl <- prepareReader(readerControl, reader(x))

    if (!filehash::dbCreate(dbControl$dbName, dbControl$dbType))
        stop("error in creating database")
    db <- filehash::dbInit(dbControl$dbName, dbControl$dbType)

    x <- open(x)
    tdl <- vector("list", length(x))
    counter <- 1
    while (!eoi(x)) {
        x <- stepNext(x)
        elem <- getElem(x)
        doc <- readerControl$reader(elem,
                                    readerControl$language,
                                    as.character(counter))
        filehash::dbInsert(db, meta(doc, "id"), doc)
        tdl[[counter]] <- meta(doc, "id")
        counter <- counter + 1
    }
    x <- close(x)

    p <- list(content = tdl,
              meta = CorpusMeta(),
              dmeta = data.frame(row.names = seq_along(tdl)),
              dbcontrol = dbControl)
    class(p) <- c("PCorpus", "Corpus")
    p
}

Corpus <-
VCorpus <-
function(x, readerControl = list(reader = reader(x), language = "en"))
{
    stopifnot(inherits(x, "Source"))

    readerControl <- prepareReader(readerControl, reader(x))

    x <- open(x)
    tdl <- vector("list", length(x))
    # Check for parallel element access
    if (is.function(getS3method("pGetElem", class(x), TRUE)))
        tdl <- mapply(function(elem, id)
                          readerControl$reader(elem, readerControl$language, id),
                      pGetElem(x),
                      id = as.character(seq_along(x)),
                      SIMPLIFY = FALSE)
    else {
        counter <- 1
        while (!eoi(x)) {
            x <- stepNext(x)
            elem <- getElem(x)
            doc <- readerControl$reader(elem,
                                        readerControl$language,
                                        as.character(counter))
            tdl[[counter]] <- doc
            counter <- counter + 1
        }
    }
    x <- close(x)

    as.VCorpus(tdl)
}

`[.PCorpus` <-
function(x, i)
{
    if (!missing(i)) {
        x$content <- x$content[i]
        x$dmeta <- x$dmeta[i, , drop = FALSE]
    }
    x
}

`[.VCorpus` <-
function(x, i)
{
    if (!missing(i)) {
        x$content <- x$content[i]
        x$dmeta <- x$dmeta[i, , drop = FALSE]
        if (!is.null(x$lazy))
            x$lazy$index <- x$lazy$index[i]
    }
    x
}

.map_name_index <-
function(x, i)
{
    if (is.character(i))
        match(i, meta(x, "id", "local"))
    else
        i
}

`[[.PCorpus` <-
function(x, i)
{
    i <- .map_name_index(x, i)
    db <- filehash::dbInit(x$dbcontrol[["dbName"]], x$dbcontrol[["dbType"]])
    filehash::dbFetch(db, x$content[[i]])
}
`[[.VCorpus` <-
function(x, i)
{
    i <- .map_name_index(x, i)
    if (!is.null(x$lazy))
        .Call(copyCorpus, x, materialize(x, i))
    x$content[[i]]
}

`[[<-.PCorpus` <-
function(x, i, value)
{
    i <- .map_name_index(x, i)
    db <- filehash::dbInit(x$dbcontrol[["dbName"]], x$dbcontrol[["dbType"]])
    db[[x$content[[i]]]] <- value
    x
}
`[[<-.VCorpus` <-
function(x, i, value)
{
    i <- .map_name_index(x, i)
    # Mark new objects as inactive for lazy mapping
    if (!is.null(x$lazy))
        x$lazy$index[i] <- FALSE
    x$content[[i]] <- value
    x
}

as.list.PCorpus <- as.list.VCorpus <-
function(x, ...)
    setNames(content(x), as.character(lapply(content(x), meta, "id")))

as.VCorpus <-
function(x)
    UseMethod("as.VCorpus")
as.VCorpus.VCorpus <- identity
as.VCorpus.list <-
function(x)
{
    v <- list(content = x,
              meta = CorpusMeta(),
              dmeta = data.frame(row.names = seq_along(x)))
    class(v) <- c("VCorpus", "Corpus")
    v
}

outer_union <-
function(x, y, ...)
{
    if (nrow(x) > 0L)
        x[, setdiff(names(y), names(x))] <- NA
    if (nrow(y) > 0L)
        y[, setdiff(names(x), names(y))] <- NA
    res <- rbind(x, y)
    if (ncol(res) == 0L)
        res <- data.frame(row.names = seq_len(nrow(x) + nrow(y)))
    res
}

c.VCorpus <-
function(..., recursive = FALSE)
{
    args <- list(...)
    x <- args[[1L]]

    if (length(args) == 1L)
        return(x)

    if (!all(unlist(lapply(args, inherits, class(x)))))
        stop("not all arguments are of the same corpus type")

    v <- list(content = do.call("c", lapply(args, content)),
              meta = CorpusMeta(meta = do.call("c",
                lapply(args, function(a) meta(a, type = "corpus")))),
              dmeta = Reduce(outer_union, lapply(args, meta)))
    class(v) <- c("VCorpus", "Corpus")
    v
}

content.VCorpus <-
function(x)
{
    if (!is.null(x$lazy))
        .Call(copyCorpus, x, materialize(x))
    x$content
}

content.PCorpus <-
function(x)
{
    db <- filehash::dbInit(x$dbcontrol[["dbName"]], x$dbcontrol[["dbType"]])
    filehash::dbMultiFetch(db, unlist(x$content))
}

inspect <-
function(x)
    UseMethod("inspect", x)
inspect.PCorpus <- inspect.VCorpus <-
function(x)
{
    print(x)
    cat("\n")
    print(noquote(content(x)))
    invisible(x)
}

length.PCorpus <- length.VCorpus <-
function(x)
    length(x$content)

names.PCorpus <- names.VCorpus <-
function(x)
    as.character(meta(x, "id", "local"))

`names<-.PCorpus` <- `names<-.VCorpus` <-
function(x, value)
{
    meta(x, "id", "local") <- as.character(value)
    x
}

format.PCorpus <- format.VCorpus <-
function(x, ...)
{
    c(sprintf("<<%s>>", class(x)[1L]),
      sprintf("Metadata:  corpus specific: %d, document level (indexed): %d",
              length(meta(x, type = "corpus")),
              ncol(meta(x, type = "indexed"))),
      sprintf("Content:  documents: %d", length(x)))
}

writeCorpus <-
function(x, path = ".", filenames = NULL)
{
    filenames <- file.path(path,
      if (is.null(filenames))
          sprintf("%s.txt", as.character(meta(x, "id", "local")))
      else filenames)

    stopifnot(length(x) == length(filenames))

    mapply(function(doc, f) writeLines(as.character(doc), f), x, filenames)

    invisible(x)
}
