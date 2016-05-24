                           # todo: need to setup error processing. E.g. in 2.14-0 signatures
                           #       for 'initialize-methods' in package 'methods' have an entry
                           #       with a duplicated signature which causes error.
.get_signature <- function(rdo, pos){
    sigtxt <- rapply(rdo[[pos]],
                     function(x){
                         if(length(x)>1){  #when Rd_tag = "USERMACRO"
                             ""
                         }else if(grepl("^[ ]*signature\\(", x)){
                             wrk <- gsub("^[ ]*(signature\\([^)]*\\)).*", "\\1", x)

                             # in case manual formatting has introduced more spaces. However,
                             # this needs to be done more carefully. It would be best to
                             # parse the signature rather than to scrape the string. On the
                             # other hand, signatures are produced (usually) programmatically.
                             gsub("[[:space:]]+", " ", wrk)
                         }else
                             ""
                     }, classes = "character", deflt="")
    sigtxt[ sigtxt != ""]
}

.get_item_labsig <- function(rdo, pos, lab=TRUE){
    label <- if(lab) .get_item_label(rdo, pos)
             else    as.character(NA)

    sig <- .get_signature(rdo, pos)

    if(inherits(try(parse(text=sig),silent=TRUE), "try-error")){
        txt <- parse_Rdpiece(rdo[[pos]],result="text")

        sig <- gsub("^.*(signature\\([^)]*\\)).*", "\\1", paste(txt,collapse=" "))
        sig <- gsub("[[:space:]]+", " ", sig)
        sig <- gsub("'", "", sig)   # krapka
    }

    c(name = label, signature = sig)
}

.get_top_signature <- function(rdo, pos = NULL, sec = "Methods", lab = TRUE){
    if(is.null(pos))
        pos <- .locate_top_items(rdo, sec)
    lapply(pos, function(x) .get_item_labsig(rdo, x, lab = lab) )
}

                                         # this inspect signatures in documentation of a class
inspect_clmethods <- function(rdo, final = TRUE){
    fullname <- .get.name_content(rdo)$name     # name of the class, including suffix `-class'

    cur <- .capture_promptAny(fullname, final=final)    # gather info about the actual methods
    curnames <- .get_top_labels(cur, "Methods")
    curitems <- .locate_top_items(cur, "Methods")
    cursigs <- .get_top_signature(cur, pos = curitems, sec = "Methods")

    rdonames <- .get_top_labels(rdo, "Methods")                        # gather info about rdo
    items <- .locate_top_items(rdo, "Methods")
    rdosigs <- .get_top_signature(rdo, pos = items, sec = "Methods")

    cmp <- .asym_compare(rdosigs, cursigs)                    # compare methods in rdo and cur
    indxnew <- cmp$i_new
    indxdrp <- cmp$i_removed

    if(length(indxnew)>0){                                  # Now make inference based on cmp.
        cat("Undocumented methods found.\n")        # 2012-10-16 cat("Undocumented methods: ")
        cat("\tAdding items for them.\n")           #            print(cursigs[indxnew])
                                                    #            print(cursigs[indxnew])
        curnamesnew <- curnames[indxnew]
                                                      # todo: some clean up of the code below.
                                          # todo: sort by function name (first arg. na \\item)
        newitems <- .get_subset(cur, curitems[indxnew], rdtag = "\\describe")
        newitems <- .nl_and_indent(newitems)

                                                        # 2012-10-16 proverkata za empty items
        if(length(items) == 0){                 # no \describe environment
            idescr <- .locate_sec_content(rdo, "Methods")
            if(length(idescr)==1 && is.na(idescr)){        # NA returned if no section Methods
                                                           # todo: NA is inconvenient to check
                                 # not needed: methpos <- Rdo_get_insert_pos(rdo,"\\section")
                rdo <- Rdo_insert(rdo, char2Rdpiece("Functions with methods for this class:",
                                                    "Methods", force.sec = TRUE ),
                                  newline = FALSE)
                idescr <- .locate_sec_content(rdo, "Methods")
            }
            rdo[[idescr]] <- c(rdo[[idescr]], list(newitems))
        }else{
            dindx <- .locate_enclosing_tag(rdo, items[[1]], "\\describe")
            wrk <- c(rdo[[dindx]], newitems )

                               # ne mozhe tolkova prosto ponezhe v rdo[[indx]] mozhe da ima i
                               #    drugi elementi osven \item (napr \n
                               # rdo[[dindx]] <- structure(wrk[ order(c(curnames,rdonames)) ],
                               #                           Rd_tag = "\\describe")
            wrk2 <- wrk
            allnams <- c(rdonames,curnamesnew)
            rdo[[dindx]] <- if(length(wrk2) == length(allnams))    # only only \item's present
                                structure(wrk2[order(allnams)], Rd_tag = "\\describe")
                            else
                                structure(wrk2                , Rd_tag = "\\describe")
        }

    }

    if(length(indxdrp)>0){                   # todo: maybe put this note in a section in rdo?
        cat("The following methods are no longer present:\n")
        print(rdosigs[indxdrp])
        cat("\tPlease remove their descriptions manually.\n")
    }

    rdo
}

  ## This inspects signatures in documentation of methods. Signatures in documentation of
  ## classes are stored somewhat differently.  this was written before inspect_clmethods() and
  ## was geared towards using existing code for ordinary functions (mainly parse_usage_text()
inspect_signatures <- function(rdo, package = NULL, sec = "Methods"){
    rdosigs <- .get_top_signature(rdo, sec = "Methods", lab = FALSE)      # process signatures
    sigtxt <- sapply(rdosigs, function(x) x["signature"])                 # in `rdo'
    urdo <- parse_usage_text(sigtxt)

    curtxt <- get_sig_text(rdo, package=package)                   # process actual signatures
    ucur <- parse_usage_text(curtxt)

                     # 2012-09-07 - may tryabva da varne samite elementi, ne technite indeksi.
    icomp <- .asym_compare(urdo, ucur)

    removed_sig <- urdo[ icomp$i_removed ]
    added_sig   <- ucur[ icomp$i_new ]

    if(length(added_sig)>0){
        cat(length(added_sig), "new signatures found.\n")     # print(added_sig)
        cat("\tAdding items for them.\n")
    }                                                  # else  cat("\tNo new signatures.\n")

    if(length(removed_sig)>0){                                   # for testing only;
        cat(length(removed_sig), "removed signatures:\n")        # remove eventually
        print(removed_sig)
    }else
        cat("\tNo removed signatures.\n")

    list( changed =  length(added_sig) != 0  | length(removed_sig) != 0
         , added_sig = added_sig, removed_sig = removed_sig )
}

    # todo: It would be better to call promptMethods() to get the signatures but in version
    # R-2.13.x I had trouble with argument `where' (could not figure out how to use it to
    # restrict to funtions from a package; also, promptMethods() seems to call the deprecated
    # function getMethods()). Check how these things stand in R-2.14, there may be no problem
    # any more (checked, in 2.14-0 it is the same).
                          # curtxt0 <- showMethods(fname, printTo=FALSE)
                          # curtxt <- curtxt0[-1] # drop the description string
                          # curtxt <- gsub("^[ ]*\\(inherited from[^)]*\\)[ ]*", ""  , curtxt)
                          # curtxt <- curtxt[ curtxt != "" ]   # todo: improve!
                          # curtxt <- paste("signature(", curtxt, ")", sep="")

get_sig_text <- function(rdo, package = NULL){           # finds the actual signatures without
    fname <- .get.name_content(rdo)$short                # using promptXXX functions
    meths <- if(is.character(package))
                 findMethodSignatures(fname,where = asNamespace(package))
             else
                 findMethodSignatures(fname)

    nammeths <- colnames(meths)

                        # insert quotes to make this consistent with parse_usage_text();
                        # do it defensively (hardly necessary), in case they are already there
    sig <- apply(meths, c(1,2),
                 function(x) if(!all( grepl("^\"",x))) paste("\"",x,"\"",sep="") else x)

                                                           # get one column for each signature
    sig <- apply(sig, 1, function(x) paste(nammeths, x, sep=" = ") )
                                     # adjust `sig' so that `paste()' would do the right thing
    sig <- if(is.matrix(sig)) apply( t(sig), 1, function(x) paste(x, collapse=", "))
           else               sig                  # a vector if one

    paste("signature(", sig, ")", sep="")
}

S4formals <- function(fun, ...){                                    # could be made S4 generic
    if(!is(fun, "MethodDefinition"))
        fun <- getMethod(fun, ...)
                                                  # todo: check that this is ok
    formals(body(fun@.Data)[[c(2,3)]])            #       if not, fall back to formals(m1)
}

# showMethods("plot")
# showMethods("plot",includeDefs=TRUE)
# showMethods("plot",includeDefs=TRUE,print=FALSE)
# s <- showMethods("plot",includeDefs=TRUE, print=FALSE)
#
#
# m1 <- getMethod("plot",c(x="profile.mle", y="missing"))
#
# class(m1)
# showClass(class(m1))
#
# m1
# m1@defined
# m1@target
# m1@.Data
#
# str(m1)
#
# formals(m1)
#
# m1@.Data
# m1body <- body(m1@.Data)
# m1body
# m1body[[1]]         # `{`
# m1body[[2]]         # the expression: .local <- function(...) ...
# class(m1body[[2]])  # "<-"
# m1body[[c(2,1)]]    # `<-`
# m1body[[c(2,2)]]    # .local
# m1body[[c(2,3)]]    # the function, what I need
# formals(m1body[[c(2,3)]])  # the formals !
