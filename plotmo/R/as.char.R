# as.char.R: brief description of an object as a string e.g. "c(1,2)"
#            also includes print_summary for summarizing matrices and data.frames

as.char <- function(object, maxlen=20)
{
    stopifnot(is.numeric(maxlen), length(maxlen) == 1, maxlen > 0)

    if(is.null(object))
        "NULL"

    else if(is.name(object))
        paste.trunc(object, maxlen=maxlen) # e.g. "..3" for unforced dot args

    else if(is.environment(object))
        environment.as.char(object)

    else if(is.call(object)) { # e.g. x is a call object in foo(x=1:3)
        s <- strip.space(paste(format(object), collapse=""))
        if(nchar(s) > maxlen)
            s <- paste0(substr(s, 1, maxlen), "...)")
        s
    }
    else if(NCOL(object) == 1 && is.character(object))
        paste.c(paste0("\"", object, "\""))

    else if(NCOL(object) == 1 && is.logical(object))
        paste.c(object)

    else if(NCOL(object) == 1 && is.numeric(object)) {
        # digits=4 is arb but seems about right, and zapsmall means more can
        # be displayed in limited space if just one val is say 3.553e-15
        paste.c(signif(zapsmall(object, digits=4), digits=4))
    }
    else if(length(dim(object)) == 2)
        sprintf("%s[%g,%g]", class(object)[1], NROW(object), NCOL(object))

    else if(class(object)[1] == "list") # not is.list() because e.g. lm objects are lists
        paste0("list(", paste.trunc(list.as.char(object), maxlen=maxlen+12), ")")

    else
        paste0(class(object)[1], ".object")
}
# compact description of a list
# maxlen is max length of each list element (not of the entire list)

list.as.char <- function(object, maxlen=20)
{
    stopifnot(is.list(object) || is.pairlist(object))
    s <- ""
    names <- names(object)
    for(i in seq_along(object)) {
        if(i != 1)
            s <- sprintf("%s, ", s)
        name.ok <- length(names) >= i && !is.na(names[i]) && nzchar(names[i])
        if(name.ok && names[i] == "...")
            s <- sprintf("%s...", s) # print dots as ... not as ...=pairlist.object
        else {
            if(name.ok)
                s <- sprintf("%s%s=", s, names(object)[i])
            s <- sprintf("%s%s", s, as.char(object[[i]], maxlen=maxlen))
        }
    }
    s    # one element character vector e.g "x=1, 2"
}
environment.as.char <- function(env, maxlen=60) # compact description
{
    if(is.null(env)) # illegal, but we still want to format it
        return("env(NULL)")

    stopifnot(is.environment(env))

    # format(env) returns "<environment: xxx>"
    stripped.env <- gsub("<environment: |>", "", format(env)[1])

    # if it's a standard environment return the environment's name
    if(grepl("^namespace:|^R_[[:alnum:]]+Env", stripped.env))
        stripped.env # something like "namespace:stats" or "R_GlobalEnv"

    else # return the names of the objects in the environment
        sprintf("env(%s)",
                paste.trunc(paste0(ls(env, all.names=TRUE), collapse=", "),
                            maxlen=maxlen))
}
# The main purpose of this routine is to summarize matrices and data.frames,
# but it will also (semi)gracefully handle any object that is passed to it.
#
# Note that this only does anything if trace >= 2.
#
# If x is a matrix or dataframe or similar, print first few rows and last row.
# If trace >= 4, then print all rows and cols, up to 1000 rows and 100 cols.
#
# the details argument:
#   0=don't print data, print the colnames truncated to one line of output
#   1=don't print data, print all colnames
#   -1=like print data but don't prefix the output with spaces
#   2=print the data

print_summary <- function(x, xname=trunc.deparse(substitute(x)),
                          trace=2, msg="", prefix="", details=2)
{
    stopifnot(is.numeric(trace))
    if(trace < 2)
        return()
    if(is.null(x)) {
        printf("%s: NULL\n", xname)
        return()
    }
    if(length(x) == 0) {
        printf("%s: length zero\n", xname)
        return()
    }
    # try(data.frame(), silent=TRUE) is not actually silent
    # for language objects, so handle them specially
    if(is.language(x)) {
        x$na.action <- NULL # don't want to print the na.action if there is one
        s <- try(format(x))
        max <- if(trace <= 2) 8 else 1000
        if(length(s) > max) {
            s <- s[1:max]
            s[max] <- paste(s[max], "\n...")
        }
        s <- gsub("[ \t\n]", "", s) # remove white space
        s <- gsub(",", ", ", s)     # replace comma with comma space
        s <- paste(s, collapse="\n    ", sep="")
        printf("%s%s%s:\n%s\n", prefix, xname, msg, s)
        return()
    }
    if(is.list(x) && !is.data.frame(x)) { # data.frames are lists, hence must check both
        if(details < 2 && trace < 4) {
            printf("%s: list with elements %s\n", xname, quotify.trunc(paste(names(x))))
            return()
        }
        printf("%s ", xname)
        str(x)
        return()
    }
    df <- try(my.data.frame(x, trace, stringsAsFactors=FALSE), silent=TRUE)
    if(is.try.err(df)) { # be robust for whatever gets passed to this function
        printf("print_summary: could not convert class \"%s\" to a data.frame\n",
               class(x)[1])
        printf("%s%s%s:\n", prefix, xname, msg)
        if(length(dim(x)) == 2) { # it's a matrix or other 2D object?
            if(trace >= 4) {
                try(print.with.strings.quoted(x))
                try(print(summary(x)))
            } else {
                try(print.with.strings.quoted(head(x)))
                printf("...\n")
            }
        } else
            try(print.with.strings.quoted(x))
        return()
    }
    if(details < 2 && trace < 4) {
        # don't print the data, just the dimensions and colnames
        if(details != -1)
            printf("        ")
        printf("%s%s[%d,%d]%s ",
                prefix, xname, nrow(df), ncol(df), msg)
        print.colnames(x, full=details == 2)
        return()
    }
    colnames <- safe.colnames(x)
    printf("%s%s[%d,%d]%s%s:\n",
           prefix, xname, nrow(df), ncol(df), msg,
           if(is.null(colnames)) " with no column names" else "")
    df.short <- df
    maxrows <- if(trace >= 4) 1000 else 5
    if(maxrows < nrow(df)) {
        df.short <- df[c(1:(maxrows-1), nrow(df)), , drop=FALSE]
        if(is.null(rownames(df.short)))
            rownames(df.short) <- c(1:(maxrows-1), nrow(df))
        rownames(df.short)[maxrows-2+1] <- "..."
    }
    maxcols <- if(trace >= 4) 100 else 10
    if(maxcols < ncol(df)) {
        df.short[,maxcols] <- "..."
        df.short <- df.short[, 1:maxcols, drop=FALSE]
        if(!is.null(colnames))
            colnames(df.short)[maxcols] <- "..."
    }
    try(print.with.strings.quoted(df.short))
    is.fac <- sapply(df, is.factor)
    if(is.null(colnames))
        colnames(df) <- sprintf("[,%d]", seq_len(NCOL(x)))
    if(any(is.fac)) {
        names <- paste0(colnames(df),
                        ifelse(sapply(df, is.ordered), "(ordered)", ""))
        if(sum(is.fac) == 1) # only one fac, so enough space to print levels too
            printf("  %s is a factor with levels: %s\n",
                   paste.trunc(names[is.fac]),
                   paste.trunc(levels(df[,is.fac])))
        else
            printf("  factors: %s\n", paste.trunc(names[is.fac]))
    }
    if(trace >= 4)
        try(print(summary(df)))
}
print.colnames <- function(x, full=FALSE)
{
    colnames <- safe.colnames(x)
    if(is.null(colnames))
        printf("with no column names\n")
    else {
        colnames[which(colnames == "")] <- "\"\""
        if(full) # full colnames (up to 1000 characters)
            printf("with colname%s %s\n",
                if(length(colnames(x)) > 1) "s" else "",
                paste.trunc(colnames, maxlen=getOption("width")-20))
        else # short version of colnames
            printf.wrap("with colname%s %s\n",
                if(length(colnames(x)) > 1) "s" else "",
                paste.trunc(colnames))
    }
}
# Like print but puts quotes around strings.
# Useful for disambiguating strings from factors.
#
# "..." is not quoted because it is used as a
# "something was deleted" indicator in print_summary

print.with.strings.quoted <- function(x)
{
    if(length(dim(x)) == 2)
        for(i in seq_len(NCOL(x)))
            if(is.character(x[,i]) && x[,i] != "...")
                x[,i] <- paste0("\"", x[,i], "\"")
    print(x)
}
