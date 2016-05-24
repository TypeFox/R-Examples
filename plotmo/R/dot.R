# dot.R: functions to access dot arguments
# Stephen Milborrow Mar 2015 Durban
#
# TODO when match.call is fixed (R 3.2.0), remove the dots arg in all these funcs
#      i.e. use the parent's dots

dot <- function(ARGNAME, ..., DEF=NA, EX=TRUE, NEW=NA)
{
    dots <- drop.unnamed.dots(match.call(expand.dots=FALSE)$...)
    argname <- process.argname(ARGNAME)
    exact <- process.exact(argname, EX)
    new <- process.new(NEW, argname, deparse(substitute(DEF)))
    for(i in seq_along(argname))
        if(!is.na(idot <- dotindex.aux(argname[i], dots, exact[i]))) {
            argval <- try(eval(dots[[idot]], parent.frame(1)))
            if(is.try.err(argval))
                stop0("cannot evaluate '", argname[i], "'")
            dotname <- names(dots)[idot]

            # TODO following commented out until we want to start
            #      issuing deprecated messages for earth and plotmo
            #
            # # the grepl prevents a warning if user uses say a dot arg of
            # # plain 'col' when ARGNAME="pt.col col.pt col"
            # if(is.specified(new) && argname[i] != new && grepl("\\.", dotname))
            #     warning0("'", dotname,
            #              "' is deprecated,  please use '", new, "' instead")

            return(argval)
        }
    DEF
}
# Like dot() but default is existing value of ARGNAME.
# For example, dotd("xlab", ...) is equivalent to dot("xlab", DEF=xlab, ...).
# TODO add to test suite
dotd <- function(ARGNAME, ..., EX=TRUE)
{
    if(is.dot("DEF", ...))
        stop0("'DEF' cannot be used with dotd")
    if(is.dot(ARGNAME, ..., EX=EX))
        dot(ARGNAME, ..., EX=EX)
    else # use the current value of ARGNAME as the default
        eval(as.name(ARGNAME), parent.frame(1))
}
is.dot <- function(ARGNAME, ..., EX=TRUE)
{
    dots <- drop.unnamed.dots(match.call(expand.dots=FALSE)$...)
    argname <- process.argname(ARGNAME)
    exact <- process.exact(argname, EX)
    for(i in seq_along(argname))
        if(!is.na(dotindex.aux(argname[i], dots, exact[i])))
            return(TRUE)
    FALSE
}
dotindex <- function(ARGNAME, ..., EX=TRUE)
{
    dots <- drop.unnamed.dots(match.call(expand.dots=FALSE)$...)
    argname <- process.argname(ARGNAME)
    exact <- process.exact(argname, EX)
    for(i in seq_along(argname))
        if(!is.na(idot <- dotindex.aux(argname[i], dots, exact[i])))
            return(idot)
    NA
}
drop.unnamed.dots <- function(dots)
{
    dots[which(names(dots) == "")] <- NULL
    dots
}
# allow comma or space separated argnames
# e.g. convert c("a", "b,c d") to c("a", "b", "c", "d")
process.argname <- function(argname)
{
    stopifnot(is.character(argname))
    argname <- gsub(" +|,+", ",", argname)  # convert space or multi commas to comma
    argname <- gsub("^,+|,+$", "", argname) # drop leading and trailing commas
    if(any(!nzchar(argname)))
        stop0("empty string in ARGNAME")
    unlist(strsplit(argname, split=","))    # convert to a vector
}
process.exact <- function(argname, exact)
{
    stopifnot(is.numeric(exact) || is.logical(exact),
              all((exact == 0) | (exact == 1)))
    if(length(exact) > length(argname))
        stop0("length(EX)=", length(exact),
              " is greater than length(ARGNAME)=", length(argname))
    rep(exact, length.out=length(argname)) # recycle
}
process.new <- function(new, argname, defname) # returns NA or a string
{
    if(anyNA(new))
        return(NA)
    if(is.numeric(new)) {
        if(length(new) != 1)
            stop0("length(NEW) != 1")
        if(new < 0 || floor(new) != new)
            stop0("NEW=", new, " is not allowed")
        if(new == 0) {
            if(!grepl("^[[:alnum:]._]+$", defname))
                stop0("NEW=0 cannot be used when DEF=",
                      defname, " (not an identifier)")
            # following helps prevent mistakes when e.g. defname=NA or NULL
            if(grepl("^[A-Z]+$", defname)) # all upper case
                stop0("NEW=0 cannot be used when DEF=", defname)
            return(defname)
        }
        if(new > length(argname))
            stop0("NEW=", new, " but length(ARGNAME) is only ", length(argname))
        return(argname[new])
    }
    # new is a string
    stopifnot.identifier(new, "NEW")
    new
}
dotindex.aux <- function(argname, dots, exact=FALSE) # workhorse
{
    stopifnot.identifier(argname, "ARGNAME")

    if(length(dots) == 0)
        return(NA)

    # first look for an exact match

    caller <- callers.name(n=2)
    index <- which(argname == names(dots))
    if (length(index) > 1)  # multiple exact matches?
        stop0("argument '", argname, "' for ", caller, "() is duplicated")
    if(length(index) == 0)  # no exact match
        index <- NA
    if(!is.na(index) || exact)
        return(index)

    # look for a partial match

    index <- which(!is.na(charmatch(names(dots), argname)))

    if(length(index) == 0)  # no match
        return(NA)

    if(length(index) == 1)  # single match
        return(index)

    # length(index) > 1 multiple matches

    stopifnot(all(index >= 0))
    name1 <- names(dots)[index[1]]
    name2 <- names(dots)[index[2]]

    if(name1 == name2)      # e.g. foo("abc", a=1, a=2)
        stop0("argument '", name1, "' for ", caller, "() is duplicated")

    # e.g. arguments 'a' and 'ab' both match 'abc' in foo()
    stop0("arguments '", name1, "' and '", name2,
         "' both match '", argname, "' in ", caller)
}
