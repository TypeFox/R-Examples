# check.index.R

# Check that an index vector specified by the user is ok to index an object.
# We want to preclude confusing R messages or behaviour later.
# An example is when max(index) > length(object) which quietly
# returns NA and can cause confusing downstream behaviour.
# This returns a vector suitable for indexing into object (will
# be identical to index unless index is a character vector).
#
# If index is a character vector, then regex matching is used against
# the names in the object, and an integer vector is returned.

check.index <- function(index, index.name, object,
    colnames        = NULL,
    is.col.index    = FALSE,
    allow.empty     = FALSE, # if index is char will warn if necessary regardless of allow.empty
    allow.zeros     = FALSE,
    allow.negatives = TRUE,
    allow.dups      = FALSE,
    treat.NA.as.one = FALSE,
    is.degree.spec  = FALSE) # special handling for degree1 and degree2 specs
{
    index.name <- quotify.short(index.name, "index", quote="'")

    # check that the given index and object can be evaluated
    try <- try(eval(index))
    if(is.try.err(try))
        stop0("illegal ", index.name)
    try <- try(eval(object))
    if(is.try.err(try))
        stop0("illegal ", quotify.short(object, quote="'"))

    is.col.index    <- check.boolean(is.col.index)
    allow.empty     <- check.boolean(allow.empty)
    allow.zeros     <- check.boolean(allow.zeros)
    allow.negatives <- check.boolean(allow.negatives)
    allow.dups      <- check.boolean(allow.dups)
    treat.NA.as.one <- check.boolean(treat.NA.as.one)

    if(is.null(index)) {
        if(!allow.empty)
           stop0(index.name, " is NULL and cannot be used as an index")
        return(NULL)
    }
    if(treat.NA.as.one && (length(index) == 1 && is.na(index)[1]))
        index <- 1
    if(anyNA(index))
        stop0("NA in ", index.name)
    if(NROW(index) != 1 && NCOL(index) != 1)
        stop0(index.name, " must be a vector not a matrix (",
              index.name, " has dimensions ",
              NROW(index), " x ", NCOL(index), ")")

    len <- get.len(object, is.col.index)

    if(is.character(index)) # currently only works for column names of object
        check.character.index(index, index.name, object, colnames, len,
                              is.col.index, allow.empty, is.degree.spec)
    else if(is.logical(index))
        check.logical.index(index, index.name, len, allow.empty)
    else if(is.numeric(index))
        check.numeric.index(index, index.name, len, allow.empty,
            allow.negatives, allow.dups, allow.zeros, treat.NA.as.one)
    else
        stop0(index.name, " must be an index vector (numeric, logical, or character)")
}
get.len <- function(object, is.col.index)
{
    if(is.col.index)
        len <- NCOL(object) # index is for columns of object
    else if(is.null(dim(object)))
        len <- length(object)
    else
        len <- NROW(object) # index is for rows of object
                            # NROW also works for lists

    stopifnot(length(len) == 1)
    stopifnot(len > 0)

    len
}
matchmult <- function(x, tab) # like match but return multiple matches if present
{
     matches <- integer(0)
     for(i in seq_along(x)) {
        xi <- x[i]
        for(itab in 1:length(tab))
            if(xi == tab[itab])
                matches <- c(matches, itab)
     }
     matches
}
# This does regex matching of index and returns an integer vector

check.character.index <- function(index, index.name, object, names, len,
                                  is.col.index, allow.empty, is.degree.spec)
{
    stopifnot(is.character(index))
    if(!is.character(index))
    # certain regular expressions match everything, even if names not avail
    if(length(index) == 1 && index %in% c("", ".", ".*"))
        return(1:len)
    if(is.col.index && is.null(names))
        names <- colnames(object)
    if(length(names) == 0 || !is.character(names))
        stop0(index.name,
              " specifies names but the names are unavailable")
    matches <- integer(0)
    warning.names <- integer(0) # these regexs don't match any column names
    for(name in index) {
        igrep <- grep(name, names)
        if(length(igrep))
            matches <- c(matches, igrep)
        else
            warning.names <- c(warning.names, name)
    }
    if(is.degree.spec) {
        if(is.null(dim(object)))          # vector, degree1
            matches <- matchmult(matches, object)
        else if(length(dim(object)) == 2) # 2D matrix, degree2
            matches <- c(matchmult(matches, object[,1]), matchmult(matches, object[,2]))
        else
            stop0("that kind of object is not yet supported for ", index.name)
    }
    new.index <- unique(matches[!is.na(matches)])
    for(name in warning.names)
        warning0("\"", name, "\" in ", unquote(index.name), " does not match any names\n",
                 "         Available names are ", paste.trunc(quotify(names)))
    new.index
}
check.logical.index <- function(index, index.name, len, allow.empty)
{
    stopifnot(is.logical(index))
    if(!allow.empty) {
        if(length(index) == 0)
            stop0("length(", unquote(index.name), ") == 0")
        if(length(index[index == TRUE]) == 0)
            stop0(index.name, " is all FALSE")
    }
    # note that a single FALSE or TRUE is ok regardless of length(object)
    if(length(index) > len && length(index) != 1) {
        stop0("logical index ", index.name, " is too long.\n",
              "       Its length is ", length(index),
              " and the max allowed length is ", len)
    }
    index
}
check.numeric.index <- function(index, index.name, len, allow.empty,
                                allow.negatives, allow.dups,
                                allow.zeros, treat.NA.as.one)
{
    stopifnot(is.numeric(index))
    if(!allow.empty) {
        if(length(index) == 0)
            stop0(index.name, " is empty, (its length is 0)")
        else if(all(index == 0))
            if(length(index) == 1)
                stop0(index.name, " is 0")
            else
                stop0(index.name, " is all zeros")
    }
    if(!is.integral(index))
        stop0(index.name, " is not an integer")
    if(any(index < 0) && any(index > 0))
        stop0("mixed negative and positive values in ", index.name)
    if(!allow.zeros && any(index == 0) && length(index) != 1)
        warning0("zero in ", index.name)
    if(!allow.negatives && any(index < 0))
        stop0("negative value in ", index.name)
    if(!allow.dups && any(duplicated(index)))
        warning0("duplicates in ", index.name)
    if(any(abs(index) > len)) {
        if(length(index) == 1)
            prefix <- paste0(unquote(index.name), "=", index, " but ")
        else
            prefix <- paste0(index.name, " is out of range, ")
        if(len != 1)
            stop0(prefix, "allowed values are 1 to ",  len)
        else if(treat.NA.as.one)
            stop0(prefix, "the only allowed value is 1 (or NA)")
        else
            stop0(prefix, "the only allowed value is 1")
    }
    index
}
