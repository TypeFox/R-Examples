inspect_usage <- function(rdo){
    ut <- get_usage_text(rdo)
    if(grepl("=[ ]*[,)]", ut)){        # TODO: this is a patch; should be sorted out properly!
        cat("Encountered formal arguments in the form 'arg ='\n")
        cat("\twithout righ-hand side, dropping the equal sign to proceed.\n")
        ut <- gsub("=[ ]*([,)])", "\\1", ut)
    }

    urdo <- parse_usage_text(ut)
    ucur <- lapply(seq_along(urdo), function(x) get_usage(name = urdo[[x]]$name,
                                                          S3class = urdo[[x]]$S3class,
                                                          S4sig   = urdo[[x]]$S4sig,
                                                          infix   = urdo[[x]]$infix,
                                                          fu = urdo[[x]]$fu,
                                                          out.format="list"))

    mapply(compare_usage1, urdo, ucur, SIMPLIFY=FALSE)# ucur is of the same length as urdo but
}                                                     # it may contain NULLs in some positions

inspect_args <- function(rdo, i_usage){
    argnames_rdo <- Rdo_get_argument_names(rdo)

                                                                 # TODO: this is for testing!
    # stopifnot(all( tools:::.Rd_get_argument_names(rdo) ==
    #                Rdo_get_argument_names(rdo) ))

    if(missing(i_usage))
        i_usage <- inspect_usage(rdo)

    argnames_cur <- unique(unlist(lapply(i_usage,
                                         function(x) attr(x,"details")$cur_usage$argnames)))

    structure(identical(sort(argnames_rdo), sort(argnames_cur)),
              details = list( rdo_argnames = argnames_rdo
                          , cur_argnames = argnames_cur
                          , added_argnames = argnames_cur[!(argnames_cur %in% argnames_rdo)]
                          , removed_argnames = argnames_rdo[!(argnames_rdo %in% argnames_cur)]
              ))
}
                                                                              # new 2012-09-22
parse_text <- function(text, ...,  keep = TRUE){            # see comments in parse_usage_text
   ks <- options("keep.source")
   if(identical(keep, ks)){     # comment on 2014-03-28 : this condition will never be TRUE
                                #                         since ks is named!
       res <- parse(text=text, ...)
   }else{
       options(keep.source = keep)
       res <- parse(text=text, ...)
       options(ks) # restore previous value
   }
   res
}

parse_usage_text <- function(text){
    text <- gsub("\\\\dots", "...", text)                                    # deal with \dots

                        # deal with \method and \S3method macros for desciption of S3 methods.
    text <- gsub("\\\\(S3)?method\\{([^}]*)\\}\\{([^}]*)\\}", "S3_usage(\\2,\\3)", text)

                                            # todo: tova sa krapki za S3 metodi.  Need to pass
                                            #       syntactically correct text to parse below.
    text <- gsub("S3_usage\\((\\[+),", "S3_usage(`\\1`,", text) #  [, [[
    text <- gsub("S3_usage\\((\\$+),", "S3_usage(`\\1`,", text)   # $
    if(any(grepl("S3_usage\\([^,]*,function[ ]*\\)", text))){ # 2012-10-12 dobavyam "any"
        fun_flag <- TRUE
        text <- gsub("(S3_usage\\([^,]*,)function[ ]*\\)", "\\1fufufu_function)", text)
    }else
        fun_flag <- FALSE
                                                                   # deal with \S4method macro
    text <- gsub("\\\\S4method\\{([^}]*)\\}\\{([^}]*)\\}", "S4_usage(\\1,\\2)", text)

                               # parse_text ensures that the srcref attribute will be set.
                               # parse() normally does that but 'R CMD check' for example sets
    tmp <- parse_text(text)    # keep.source to FALSE, leading to errors in the examples.
                                                        # todo: check for errors after parse()
                       # a character vector with one element for each usage expression in text
    usages <- sapply(attr(tmp, "srcref"),
                     function(x) paste(as.character(x), sep="\n", collapse="\n"))

    if(fun_flag)
        usages <- gsub("fufufu_function", "function", usages)

    lapply(usages, parse_1usage_text)
}

                     # TODO: x %% y i podobni tyabva da se opravyat! vzh Arithmetic.Rd ot base
parse_1usage_text <- function(text){
                # modify usage texts to become function bodies embedded in calls to formals()
    f <- function(x) paste0(sub("^[^(]+\\(", "formals(function(", x), " NULL )")

    fu <- TRUE # in  most cases we are dealing with functions

    if(grepl("^[[:alpha:]._][[:alnum:]._]*$", text)){  # a variable, not function call
        fu <- FALSE
        name <- text
        S3class <- ""
        S4sig = ""
        infix <- FALSE
        res <- NULL

    }else if(grepl("^S3_usage",text)){
        S3class <- gsub("[^,]*,([^)]*)\\).*", "\\1", text)

                                        # patch, undo the quotes maybe put by parse_usage_text
        text <- gsub("S3_usage\\([`\"']([^,`]+)[`\"']", "S3_usage(\\1", text)
        name    <- gsub("[^(]*\\(([^,]*),.*", "\\1", text)

        if(grepl("S3_usage\\(.*\\)[ ]* <-[ ]*value[ ]*$", text)){   # S3 assignment method
            text <- gsub("(S3_usage\\()([^,]*)(.*)\\)[ ]* <-[ ]*value[ ]*$",
                         "\\1`\\2<-`\\3,value)", text)
            name <- paste0(name, "<-")
        }else if(!grepl("^[[:alpha:]._:][[:alnum:]._:]+$", name)){        # non-syntactic name
            text <- gsub("S3_usage\\(([^,]*,)", "S3_usage(`\\1`,", text)
        }

        text <- gsub("[^)]*\\)(.*)", paste(name, ".", S3class, "\\1", sep=""), text)
        S4sig = ""
        infix <- FALSE

        wrk <- f(text)
        res <- eval(parse(text = wrk))
    }else if(grepl("^S4_usage",text)){
        name  <- gsub("[^(]*\\(([^,]*),.*", "\\1", text)   # same as for S3 name above
        S3class <- ""
        S4sig <- gsub("[^,]*,([^)]*)\\).*", "\\1", text)   # same as for S3class above
                                             # but S4 signatures may have more than 1 element,
                                             # make a character vector from the single string
        s <- paste("names(formals(function(", S4sig, ") NULL))", sep="")

        S4sig <- eval(parse(text = s))
        text <- gsub("[^)]*\\)(.*)", paste0(name, "\\1"), text)
        infix <- FALSE

        wrk <- f(text)
        res <- eval(parse(text = wrk))
    }else{                                          # extract the function names from `usages'
        S3class <- ""
        S4sig = ""

        e <- parse(text = text)
        name <- as.character(e[[1]][[1]])
                             # !grepl( paste0("^",name)) - s regexp tryabva da se vnimava
                             #         ponezhe ako name e "+", tova e drugo.
                             # ne mozhe i fixed = TRUE, ponezhe iskam da match-na ot nachaloto
        infix <- !(substr(text,1,nchar(name)) == name)

        ec <- lapply(1:length(e[[1]]), function(x) as.character(e[[1]][[x]]))

        res <- NULL
        if(name == "<-"){                                        # todo: more care needed here
            e2 <- e[[1]][[2]]
            if(is.call(e2)){
                e2_named <- .make_named(e2)
                name <- paste0(e2[[1]], "<-")
                wrk <- paste0("formals(function(",
                              paste(e2_named[-1], collapse=","),
                              ",", e[[1]][[3]],   # rhs (i.e., value) should have 1 elem only
                              ") NULL )" )
            }else{ # simple assignment (todo: can it be anything else here?)
                wrk <- paste0("formals(function(", paste(ec[-1], collapse=","), ") NULL )" )
            }
        }else if(name == "!"){
            wrk <- sub("!", "", text)
            wrk <- paste0("formals(function(", wrk, ") NULL)")

        }else if(name == ":"){
            wrk <- "formals(function(from,to) NULL)"
        }else if(name %in% c("if", "for", "while", "repeat", "break", "next")){ # control
             res <- as.character(sapply(e[[1]], identity))
        }else if(infix){
            e <- parse(text = paste0("quote(",text,")"))
            ec <- .make_named(lapply(e[[1]][[2]], identity))
            wrk <- paste0("formals(function(", paste(ec[-1], collapse=","), ") NULL )" )
        }else{
            wrk <- f(text)
        }

        if(is.null(res))
            res <- eval(parse(text = wrk))
    }

    pairlist2f_usage1(res, name, S3class, S4sig, infix, fu) # convert pairlists obtained from
}                                                  # `formals()' into named "f_usage" objects.

.make_named <- function(v, sep = " = "){
    if(length(v) == 0)
        return(character(0))

    nams <- allNames(v)
    wrk <- sapply(seq_along(v),
                  function(x) if(nams[x]=="") paste(v[x]) else paste(nams[x], v[x], sep=sep))
    wrk
}
                                     # generate f_usage object for a function (needs clean up)
                                                              # todo: argument `...' not used?
get_usage <- function(object, name = NULL, force.function = FALSE, ...,
                      S3class = "", S4sig = "", infix = FALSE, fu = TRUE,
                      out.format = "text"){
    if (missing(name))                                  # based on a chunk from  utils::prompt
        if (is.character(object)){
            name <- object
            object <- NULL
        }else {
            name <- substitute(object)
            if (is.name(name))
                name <- as.character(name)
            else if (is.call(name) &&
                     (as.character(name[[1L]]) %in%  c("::", ":::", "getAnywhere"))) {
                name <- as.character(name)
                name <- name[length(name)]
            }
            else stop("cannot determine a usable name")
        }
                                                   # get(name, envir = asNamespace("mixAR"))
                                                   # do.call(getAnywhere,list(x))
    x <- if(!missing(object) && !is.null(object))
             object
         else if(!fu){
             x0 <- try(get(name), silent=TRUE)
             if(inherits(x0,"try-error")){
                 spec_values <- c("NULL", "TRUE", "FALSE", "NA",
                                  sapply(c(Inf, NaN, -Inf), as.character))
                 if(name %in% spec_values)
                     x0 <- spec_values[ name == spec_values ][1]
                 else
                     x0 <- NULL
             }
             "OK: variable"       # this value is not used further, it must not be NULL though
         }else if(!identical(S3class, "")){         # 2012-10-16 dobavyam getAnywhere, etc.
             x0 <- try(getS3method(name, S3class), silent=TRUE)
             if(inherits(x0,"try-error")){
                 x0 <- do.call(getAnywhere, list(paste0(name, ".", S3class)))
                 if(length(x0$objs) > 0)
                     x0$objs[[1]]
                 else
                     NULL
             }else
                 x0
         }else if(!identical(S4sig, "")) # todo: coordinate with the rest!
             getMethod(name, S4sig)     # transform S4sig here if it is not convenient to
                                        # keep it ready for use entry in f_usage objects.
         else{
             name0 <- if(grepl('^".*"$', name))  # non-syntactic name
                          sub('^"(.*)"$', "\\1", name)
                      else
                          name

             x0 <- try(get(name0, envir = parent.frame()), silent=TRUE)
             if(inherits(x0,"try-error")){
                 x0 <- do.call(getAnywhere, list(name0))
                 if(length(x0$objs) > 0) # todo: needs more work here. IN particular, there
                     x0$objs[[1]]        #       should be a package argument to avoid taking
                                         #       blindly whatever comes up.
                 else
                     NULL
             }else
                 x0
         }

    if(is.null(x))
        return(x)

    if (fu && !(is.function(x) || force.function)){
        warning("The object is not a function.")         # The return value may be appropriate
        return(name)                                     # for data objects.  # todo: Rethink!
    }

    argls <- if(fu){
                 if(!identical(S4sig, ""))
                     S4formals(x)
                 else if(!identical(S3class, ""))
                     formals(x)
                 else{
                     spec_args <- .special_args[name]
                     if(is.na(spec_args)){
                         wrk <- formals(x)
                         if(is.null(wrk))  # takes care for primitive functions, napr. seq.int
                             wrk <- formals(args(x))
                             # argls will still be NULL if x is a function with no arguments.
                             # Note that formals returns pairlist and the pairlist with
                             # zero elements is NULL (unlike list()). Taka che
                             # pairlist2f_usage() needs to know how to deal with this case.
                         wrk
                     }else
                         eval(parse(text = paste0("formals(function(", spec_args, ") NULL)")))
                 }
             }else
                 NULL
                                    # 2012-10-11 smenyam pairlist2f_usage na pairlist2f_usage1
    res <- pairlist2f_usage1(argls, name, S3class = S3class, S4sig = S4sig,
                            infix = infix, fu = fu)

    if(out.format != "list")
        res <- as.character(res)

    res
}

.special_args <- c("~"    = "y, model",
                   "@"    = "object, name",
                   "$"    = "x, name",
                   "||"   = "x, y",
                   "&&"   = "x, y",
                   "["    = "x, i, j, ..., drop = TRUE",
                   "[["   = "x, i, j, ..., exact = TRUE",
                   "@<-"  = "object, name, value",             # formals("@<-") actually works
                   "$<-"  = "x, name, value",
                   "[<-"  = "x, i, j, ..., value",
                   "[[<-" = "x, i, value",
                   ":"    = "from, to"
                   )

    # the comparison is symmetric but the interpretation assumes that ucur may be more recent.
compare_usage1 <- function(urdo, ucur){   # urdo - usage from Rdo file/object;
                                          # ucur -       generated from actual object
    obj_removed <- is.null(ucur) || is.na(ucur)
    obj_added   <- is.null(urdo) || is.na(urdo)

    if(!obj_added && !obj_removed  && urdo$S3class != ""){
        fn <- paste(urdo$name, ".", urdo$S3class, sep="")
        if(ucur$name == fn){
            ucur$name <- urdo$name
            ucur$S3class <- urdo$S3class
        }
    }

    status <- identical(urdo, ucur)

    alias <- if(obj_removed)             ""
             else if(ucur$S3class != "")
                 paste(ucur$name, ".", ucur$S3class, sep="")
             else if(!identical(ucur$S4sig, ""))
                 paste0(ucur$name, ",", paste(ucur$S4sig, collapse=","), "-method")
             else
                 ucur$name

    if(grepl('^".*"$', alias)){  # non-syntactic name, drop the quotes
        alias <- sub('^"(.*)"$', "\\1", alias)
    }


    if(obj_removed || obj_added)
        return( structure( status, details = list( obj_removed = obj_removed
                                                 , obj_added   = obj_added
                                                 , rdo_usage   = urdo
                                                 , cur_usage   = ucur
                                                 , alias       = alias
                                   )) )

    identical_names <- urdo$name == ucur$name

    identical_argnames <- identical(urdo$argnames, ucur$argnames)
    identical_defaults <- identical(urdo$defaults, ucur$defaults)
    identical_formals  <- identical_argnames & identical_defaults

    added_argnames <- ucur$argnames[ !(ucur$argnames %in% urdo$argnames) ]
    removed_argnames <- urdo$argnames[ !(urdo$argnames %in% ucur$argnames) ]

                                                   # note: !!! intersect() is binary operation
    s <- intersect( intersect(names(urdo$argnames), names(ucur$argnames)),
                    intersect(names(urdo$defaults), names(ucur$defaults)) )

    unchanged_defaults <- urdo$defaults[ ucur$defaults[s] == urdo$defaults[s] ]

    names_unchanged_defaults <- names(unchanged_defaults)[unchanged_defaults]

               # todo: more details for the case when !identical, e.g. equal up to reordering,
               #       added/removed defaults

    structure( status, details = list( identical_names          = identical_names
                                     , obj_removed              = obj_removed
                                     , obj_added                = obj_added
                                     , identical_argnames       = identical_argnames
                                     , identical_defaults       = identical_defaults
                                     , identical_formals        = identical_formals
                                     , added_argnames           = added_argnames
                                     , removed_argnames         = removed_argnames
                                     , names_unchanged_defaults = names_unchanged_defaults
                                     , rdo_usage                = urdo
                                     , cur_usage                = ucur
                                     , alias = alias
                       ))
}
