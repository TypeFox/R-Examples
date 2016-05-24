# dotlib.R: miscellaneous functions for the dots routines

is.try.err <- function(object) class(object)[1] == "try-error"

# collapse, and truncate if strings in ... are too long
paste.trunc <- function(..., sep=" ", collapse=" ", maxlen=60)
{
    s <- paste(..., sep=sep, collapse=collapse)
    if(nchar(s) > maxlen) {
        stopifnot(maxlen > 3)
        s <- paste0(substr(s, 1, maxlen-3),
                    if(substr(s, maxlen-3, maxlen-3) == ".") ".." # avoid 4 dots
                    else                                     "...")
    }
    s
}
# Add the elements of the extra list to the original list.  Elements of the
# original list that have the same names as extra elements get overwritten.
#
# Like utils::modifyList(keep.null=TRUE) except:
# (i)  input args can be NULL (NULL is treated as an empty list)
# (ii) unnamed elements in extra are added to original (modifyList drops them)

merge.list <- function(original, extra)
{
    if(is.null(original))
        original <- list()
    if(is.null(extra))
        return(original)
    stopifnot(is.list(original))
    stopifnot(is.list(extra)) # pairlist would probably be ok too
    for(i in seq_along(extra)) {
        e <- extra[[i]]
        name <- names(extra)[i]
        if(is.null(name) || !nzchar(name)) # extra element is unnamed?
            original <- c(original, if(is.null(e)) list(NULL) else e)
        else if(is.null(e))
            original[name] <- list(NULL) # avoid "assign deletes elem if rhs is null"
        else
            original[[name]] <- e
    }
    original
}
# Evaluate each element of the list dots in the environment specified by n.
# (This function can actually be used any list, but the evaluating
# environment and enclosure are set up for dot arg lists.)
#
# TODO "scalar" is ugly, it is for par() alone and prevents
# e.g. errmsg graphical parameter "lty" has the wrong length

eval.dotlist <- function(dots, n=1, scalar=FALSE)
{
    stopifnot(is.list(dots) || is.pairlist(dots))
    env <- parent.frame(n)
    dotnames <- names(dots)
    for(i in seq_along(dots)) {
        e <- try(eval(dots[[i]], envir=env, enclos=env), silent=TRUE)
        if(!is.try.err(e)) {
            if(is.null(e))
                dots[i] <- list(NULL) # avoid "assign deletes elem if rhs is null"
            else if(!scalar || (dotnames[i] %in% PAR.VEC) || length(e) == 1)
                dots[[i]] <- e
            else
                dots[[i]] <- e[[1]] # select first element of e only
                                    # TODO it would be better to drop the element entirely
        }
    }
    dots
}
# Is the string s a valid R lexigraphic identifier?
# If allow.specials=TRUE we allow special chars used in DROP and KEEP strings.
# The name argument is used only in error messages.

stopifnot.identifier <- function(s, name=short.deparse(substitute(s)),
                                 allow.empty=FALSE, allow.specials=FALSE)
{
    if(!is.character(s))
        stop0(name, " is not a character variable (class(",
              name, ") is \"", class(s), "\")")
    if(length(s) != 1)
        stop0(name, " has more than one element\n       ",
              name, " = c(", paste.trunc("\"", s, "\"", sep=""), ")")
    if(!allow.empty && !nzchar(s))
        stop0(name, " is an empty string")
    # Note that we allow : for e.g. graphics::abline
    # TODO the following allows integers (no alphabetic characters), it shouldn't
    start <- if(allow.specials) # include , * $
                 regexpr("[^._:[:alnum:],*$]", s)
             else
                 regexpr("[^._:[:alnum:]]", s)
    if(start > 0)
        stop0("illegal character \"", substr(s, start, start),
              "\" in ", name, " = \"", s, "\"")
}
