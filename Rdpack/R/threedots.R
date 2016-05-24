# > tools:::.strip_whitespace
.strip_whitespace <-
function (x)
{
    x <- sub("^[[:space:]]+", "", x)
    x <- sub("[[:space:]]+$", "", x)
    x
}


## tools:::.Rd_drop_nodes_with_tags,      but tools:::RdTags(e)  => Rdo_tags(e)
toolsdotdotdot.Rd_drop_nodes_with_tags <-
function (x, tags)
{
    recurse <- function(e) {
        if (is.list(e))
            structure(lapply(e[is.na(match(Rdo_tags(e), tags))],
                recurse), Rd_tag = attr(e, "Rd_tag"))
        else e
    }
    recurse(x)
}


## > tools:::.Rd_get_metadata,      but tools:::RdTags()  => Rdo_tags()
toolsdotdotdot.Rd_get_metadata <-
function (x, kind)   # e.g. kind = "keyword", see help page of Rd_db()
{
    x <- x[Rdo_tags(x) == sprintf("\\%s", kind)]
    if (!length(x))
        character()
    else unique(.strip_whitespace(sapply(x, as.character)))
}


## > utils:::.getHelpFile
utilsdotdotdot.getHelpFile <-
function (file)
{
    path <- dirname(file)
    dirpath <- dirname(path)
    if (!file.exists(dirpath))
        stop(gettextf("invalid %s argument", sQuote("file")),
            domain = NA)
    pkgname <- basename(dirpath)
    RdDB <- file.path(path, pkgname)
    if (!file.exists(paste(RdDB, "rdx", sep = ".")))
        stop(gettextf("package %s exists but was not installed under R >= 2.10.0 so help cannot be accessed",
            sQuote(pkgname)), domain = NA)
    toolsdotdotdotfetchRdDB(RdDB, basename(file))
}


## > tools:::fetchRdDB
toolsdotdotdotfetchRdDB <-
function (filebase, key = NULL)
{
    fun <- function(db) {
        vals <- db$vals
        vars <- db$vars
        datafile <- db$datafile
        compressed <- db$compressed
        envhook <- db$envhook
        fetch <- function(key) lazyLoadDBfetch(vals[key][[1L]],
            datafile, compressed, envhook)
        if (length(key)) {
            if (!key %in% vars)
                stop(gettextf("No help on %s found in RdDB %s",
                  sQuote(key), sQuote(filebase)), domain = NA)
            fetch(key)
        }
        else {
            res <- lapply(vars, fetch)
            names(res) <- vars
            res
        }
    }
    res <- lazyLoadDBexec(filebase, fun)
    if (length(key))
        res
    else invisible(res)
}
