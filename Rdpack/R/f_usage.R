parse_pairlist <- function(x){
    is.missing.arg <- function(arg) typeof(arg) == "symbol" && deparse(arg) == ""
                                  # x == NULL corresponds to functions with no arguments (also
                                  # length(NULL) is 0) also, NULL is a pairlist with length 0.
                                           # Is this function used with x other than pairlist?
    if(is.null(x) || length(x) == 0)       # If not, the test of length(x) is redundant.
        return(list(argnames = character(0), defaults = character(0)))

    nonmis <- x[ !sapply(x, is.missing.arg) ]
    wrk <- character(length(nonmis))
    names(wrk) <- names(nonmis)
    for(s in names(nonmis)){
        wrk[[s]] <- paste(deparse(nonmis[[s]], backtick = TRUE, width.cutoff = 500L)
                          , collapse = "\n")
    }
    list(argnames = names(x), defaults = wrk )
}
                                                                   # 2012-10-03 new arg. infix
pairlist2f_usage1 <- function(x, name, S3class = "", S4sig = "", infix = FALSE, fu = TRUE){
    structure(c(list(name=name, S3class=S3class, S4sig=S4sig, infix=infix, fu = fu),
                parse_pairlist(x)), class="f_usage")
}

format_funusage <- function(x, name = "", width = 72, realname){
    res <- paste(name,  "(", paste(x, collapse = ", "),  ")", sep="")

    if(is.numeric(width)  &&  nchar(res, type="width") > width){
        delim <- c("(", rep(", ", length(x) - 1), ")")
        wrk <- paste(c(name, x), delim, sep="")
        lens <- nchar(wrk, type="width")
        if(!missing(realname))
            lens[1] <- nchar(realname, type="width") + 1
        indent <- paste(rep(" ", lens[1]), collapse="")
        res <- character(0)
        curlen <- 0
        for(i in seq_along(wrk)){
            if(curlen + lens[i] > width){
                res <- paste(res, "\n", indent,  sep="")
                curlen <- lens[1]   #  = number of chars in `indent'
            }
            res <- paste(res, wrk[i], sep="")
            curlen <- curlen + lens[i]
        }
    }
    res
}

deparse_usage1 <- function(x, width = 72){
    if(!x$fu) # a variable, not function
        return( structure( x$name, names = x$name ) )
          # todo: maybe x$name tryabva da e character, as.character here should not be needed.
    if(as.character(x$name) %in% c("[","[[", "$", "@",  "[<-", "[[<-", "$<-",  "@<-", "!"))
        "dummy"
    else if(x$infix){  # infix operator
        if(grepl(".+<-$", x$name)){ # name end is "<-" but is not equal to it
            name2 <- sub("^(.+)<-$", "\\1", x$name)
            m <- length(x$argnames)
            res <- paste(name2, "(", paste(x$argnames[-m], collapse=", "), ")",
                         "<-", x$argnames[m])
        }else                               # todo: make sure  that the name is not in quotes!
            res <- paste(x$argnames, collapse = paste0(" ", x$name, " "))

        return(res)
    }

    res <- x$argnames
    names(res) <- x$argnames

    nams <- names(x$defaults)
    res[nams] <- paste(res[nams], "=", x$defaults)

    assop <- grepl(".+<-$", x$name) # name end is "<-" but is not equal to it
    name <- x$name
    if(assop){
        name <- sub("<-$", "", x$name)
        value <- res[length(res)]
        res <- res[-length(res)]
    }

    res <- if(!identical(x$S3class, ""))
               format_funusage(res, paste("\\method{", name, "}{", x$S3class, "}", sep=""),
                               realname = name )
           else if(!identical(x$S4sig, ""))
               format_funusage(res, paste("\\S4method{", name, "}{",
                                          paste0(x$S4sig, collapse = ", "),
                                          "}", sep=""), realname = name )
           else
               switch(name,
                      "$" =, "@" = paste0(res[1], name, res[2]),
                      "[" =, "[[" = paste0(res[1], name, paste0(res[-1], collapse = ", "),
                                                   .closeOp[name]),
                      "!" = paste0("!", res[1]),
                      ## default
                      format_funusage(res, name)
                      )

    if(assop)           # if assignment, add to the last line, usually the only one
        res[length(res)] <- paste0(res[length(res)], " <- ", value)
                   # "[<-"  = paste0(res[1], "[", paste0(res[c(-1,-length(res))],
                   #                        collapse = ", "), "]", " <- ", res[length(res)]),
                   # "[[<-" = paste0(res[1], "[[", paste0(res[c(-1,-length(res))],
                   #                        collapse = ", "), "]]", " <- ", res[length(res)]),
                   # "$<-"  = paste0(res[1], "$", res[2], " <- ", res[3]),
                   # "@<-"  = paste0(res[1], "@", res[2], " <- ", res[3]),

    res <- gsub("...", "\\dots", res, fixed=TRUE)
    structure( paste(res, collapse = ""), names=x$name )
}

as.character.f_usage <- function(x,...){
    deparse_usage1(x)
}

deparse_usage <- function (x){
    if(class(x) == "f_usage")
        return(deparse_usage1(x))

    nams <- names(x)
    if(!is.null(nams))            # remove names since sapply merges them with the names of
        names(x) <- NULL          # the list obtained by lapply()

    res <- sapply(x, deparse_usage1)
    if(is.null(names(res)))            # in most cases names(res) will be the same as nams
        names(res) <- nams             # but give preference to the ones returned by
                                       # deparse_usage1 which takes names from the objects.
                                       # This `if' will hardly ever kick in...
    res
}

.closeOp <- list("[" = "]", "[[" = "]]", "(" = ")", "{" = "}")
