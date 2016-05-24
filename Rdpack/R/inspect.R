.parse_long_name <- function(fnam){           # 2012-11-04 new, code from .capture_promptAny()
    namreg <- "^(.+)-([^-]+)$"
    if(grepl(namreg, fnam)){                    # fnam is of the form xxxx-yyy (non-empty rhs)
        fname  <- gsub(namreg, "\\1", fnam)
        type   <- gsub(namreg, "\\2", fnam)   # without "-"
        ## suffix <- gsub("^([^-]+)(-.*)", "\\2", fnam)   # with "-"
    }else{
        fname <- fnam
        type = ""
    }

    c(name = fname, type = type)
}

parse_Rdname <- function(rdo){  # 2012-10-01 todo: Rdpack-internal ste dade type="internal"
                           # 2012-10-01 otkomentiram tova, to raboteshe vyarno po sluchaynost!
                                        # (ponezhe gsub() po-dolu go prevrasta v character()!)
                                    # nam <- rdo[[ which( tools:::RdTags(rdo) == "\\name" ) ]]
    ind <- Rdo_which_tag_eq(rdo, "\\name")
    nam <- c( rdo[[ c(ind, 1) ]])  # assumes length(ind) == 1
            # 2012-11-04 changed with the code after the commented out section
            #
            # if(grepl("^([^-]+)-.*", nam)){                    # nam is of the form xxxx-yyy
            #     fname  <- gsub("^([^-]+)-.*", "\\1", nam)
            #     type   <- gsub("^([^-]+)-(.*)", "\\2", nam)   # suffix without the '-'
            #     # suffix <- gsub("^([^-]+)(-.*)", "\\2", nam) #        with the "-"
            # }else{
            #     fname <- nam
            #     type = ""
            # }
    wrknam <- .parse_long_name(nam)
    fname <- wrknam["name"]
    type <- wrknam["type"]
                                                                           # malko kato krapka
    doctype <- Rdo_which_tag_eq(rdo, "\\docType")

    if(length(doctype) > 0  && type == "" ){  # todo: sravni s gornoto!
        type <- c( rdo[[ c(doctype[1],1) ]] ) # wrap n c() to remove attr.
    }

    list(fname = fname, type = type)
}


inspect_Rd <- function(rdo, package = NULL){                     # rdo: Rd object or file name
    if(is.character(rdo) && length(rdo)==1)
        rdo <- parse_Rd(rdo)

    if(!inherits(rdo, "Rd")  && is.null(attr(rdo,"Rd_Tag")))
        return(structure(FALSE,
                        details = "inspect_Rd: 1st arg. must be an Rd object or filename."))

    type <- toolsdotdotdot.Rd_get_metadata(rdo,"docType")
    if(length(type) == 0){
        type <- parse_Rdname(rdo)$type     # todo: clean up
        if(!(type %in% c("package", "")))   # for now; there is, e.g. ts-methods which
           type <- ""                       #          describes S3 methods.
    }

    switch(type,
           methods = inspect_Rdmethods(rdo, package = package),
           class   = inspect_Rdclass(rdo),
           package = {cat("Currently reprompt for packages is not implemented.",
                          "\tHowever, you could use the function promptPackageSexpr()",
                          "\tto create a self-updating shell. 'promptPackageSexpr'",
                          "\tis like 'promptPackage' but uses \\Sexpr commands rather",
                          "\tthan strings for automatically generated information.",
                          "Returning the Rd object unchanged.\n",
                          sep = "\n")
                      rdo
                     },
           data    = {cat("Currently reprompt for data is not available.",
                          "Returning the Rd object unchanged.\n",
                          sep = "\n")
                      rdo
                     },
           ## default
           inspect_Rdfun(rdo)
           )
}

inspect_Rdclass <- function(rdo){                                             # rdo: Rd object
    rdo <- inspect_slots(rdo)          # methods: items have the form \item{fname}{signature},
    rdo <- inspect_clmethods(rdo)      # compare with those from promptClass

    rdo                                                      # todo: inspect other things too?
}

inspect_Rdmethods <- function(rdo, package = NULL){                           # rdo: Rd object
    i_usage <- inspect_signatures(rdo, package = package)                # i_usage$added_sig
                                                                         # i_usage$removed_sig
    if(i_usage$changed)
        cat("\tMethods in the documentation are not the same as the current ones.\n")

             # 2014-08-23 new 'if' to tell the user to removed the doc. of non-existen methods
    if(length(i_usage$removed_sig) > 0){
        cat("\tMethods for the following signatures ",
            "where not found", if(is.null(package)) ""
                                 else paste0(" in package ", package)
            , ":\n", sep = "")
        for( sigus in i_usage$removed_sig){
            txt <- paste(sigus$argnames, sigus$defaults, collapse = ",\t", sep = " = ")
            cat("\t  ", txt, "\n")
        }
                                    # 2014-08-23 todo: maybe give a more cautious advice here?
        cat("\tPlease remove these from the documentation.\n\n")
    }

    if(length(i_usage$added_sig) > 0){
        cat("\tAppending new signatures to section \"Methods\"\n")

        ## 2014-08-23 TODO:
        ##       the code below now is executed only when  length(i_usage$added_sig) > 0;
        ##       before it was executed whenever i_usage$changed was TRUE.

        newsigs <- sapply(i_usage$added_sig, function(x) deparse_usage1(x))
        wrk <- lapply(c(newsigs), Rdo_sigitem)   # todo: slozhi blank line predi vseki item,
        names(wrk) <- NULL                       #       no da ne e chast ot itemite,
                                                 #       a mezhdu tyach.
        idescr <- .locate_top_tag(rdo, "Methods")
        if(length(idescr) == 0){  # no descibe enviroment there.
            wrk <- Rdo_macro(wrk, "\\describe")
            idescr <- .locate_sec_content(rdo, "Methods")
            rdo[[idescr]] <- c(rdo[[idescr]], list(wrk))
        }else
            rdo <- append_to_Rd_list(rdo, wrk, pos = idescr)
    }
                            # 2012-09-07 new; not under the if above,
                            # since aliases may be mismatched by manual editing (or otherwise)
    rdo <- update_aliases_tmp(rdo, package = package)

                                                            # if "\\usage" section exists.
    if(length(Rdo_which_tag_eq(rdo, "\\usage")) > 0)
        rdo <- inspect_Rdfun(rdo, alias_update = FALSE) # alias_update was added to Rd_fun
                                                        # to be used here
    rdo
}
                              # todo: more systematic messages about what has changed (or not)
inspect_Rdfun <- function(rdo, alias_update = TRUE){                          # rdo: Rd object
    i_usage <- inspect_usage(rdo)

    if(length(i_usage) == 0)   # todo: processing of arguments may still be needed here
        return(rdo)

    i_args <- inspect_args(rdo, i_usage)

                                   # 2011-12-08 new; corresponding additions to compare_usage1
                  # 2012-10-05 slagam as.character, ponezhe nyakoi mozhe da sa ot class `name`
                  #            todo: make sure that they are character in the first place?
    aliases <- sapply(i_usage, function(x) as.character(attr(x, "details")$alias))
    aliases <- aliases[aliases != ""]

    rdoaliases <- Rdo_collect_aliases(rdo)        # 2012-10-05 promenyam za da otchitam #ifdef

    new_aliases <- aliases[ !(aliases %in% rdoaliases) ]
    if(alias_update  && length(new_aliases) > 0)                          # update the aliases
        for(alias in unique(new_aliases))
            rdo <- Rdo_insert(rdo, char2Rdpiece(alias, "alias"))

    flag_removed <- sapply(i_usage, function(x) attr(x, "details")$obj_removed)
    if(any(flag_removed)){
        nams_removed <- sapply(i_usage, function(x) attr(x, "details")$rdo_usage$name)
        nams_removed <- nams_removed[flag_removed]
        obs_aliases <- nams_removed[ nams_removed %in% rdoaliases ]

        cat("The following objects are no longer described in this Rd object\n")
        cat("but it contains alias entries for them:\n")
        cat("\t", paste(obs_aliases), "\n")
        cat("You may wish to remove these \\alias{} entries from the Rd file.\n\n")
    }

    usage_changed <- !all( sapply(i_usage, identity))    # usage_changed <- TRUE # for testing

    if(usage_changed){                                               # update usage if changed
        f <- function(x){
            a <- as.character(attr(x,"details")$cur_usage)
            parse_Rdtext(a, section = "\\usage")
        }
        u <- lapply(i_usage, f)
        u <- .Rdinter(u, before_first = TRUE)
        wrk <- structure(do.call("c", u),  Rd_tag = "\\usage")

        rdo <- Rdo_replace_section(rdo, wrk)
    }

                   # todo: Rdo_append_argument probably should check length(newargs)>0 anyway?
    if(!i_args){                                                 # update arguments if changed
        newargs <- attr(i_args,"details")$added_argnames
        if(length(newargs)>0){
            cat("\nnewargs is:",  newargs, "\n")
            rdo <- Rdo_append_argument(rdo, newargs)
        }

        remargs <- attr(i_args,"details")$removed_argnames
        if(length(remargs)>0){
            cat("Argument(s): ",  remargs, "\n")
            cat("\tare no longer in any usage statements in this Rd object.\n")
            cat("Please remove the corresponding \\item's from section '\\arguments'.\n")
            cat("\n")
        }
    }

    rdo
}
