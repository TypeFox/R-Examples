reprompt <- function(object, infile = NULL, Rdtext = NULL, final = TRUE,
                     type=NULL, package=NULL, methods = NULL, #  for the call to promptMethods
                        verbose = TRUE, filename = NULL, sec_copy = TRUE, ...){
    objmis <- missing(object)
    tidyflag <- from_infile <- FALSE
                                     # If 'object' is a string ending in ".Rd" and containing
                                     # at least one "/", it is taken to be "infile"; a
                                     # (somewhat dubious) convenience feature for the common
                                     # mistake of omitting the name of the "infile" argument.
    if(is.null(infile)  &&  length(object) == 1  &&  is.character(object)
                        && grepl("/.*[.]Rd$", object) )
        infile <- object

    if(!is.null(Rdtext)){                                         # process Rdtext, if present
        if(is.null(infile)){
            infile <- tempfile()
            cat(Rdtext, file = infile, sep = "\n")            # save parsed Rdtext to 'infile'
            on.exit(unlink(infile))
        }else
            cat("both 'infile' and 'Rdtext' are given, ignoring Rdtext\n")
    }

    if(!objmis && inherits(object,"Rd")){
        if(verbose) cat("Processing the Rd object...\n")
        if(!is.null(infile))
            cat("ignoring 'infile' and/or 'Rdtext' since 'object' is of class 'Rd'\n")
        rdo <- object
    }else if(!is.null(infile)){
        if(verbose) cat("\nParsing the Rd documentation in file", infile, "...\n")
        else cat("\n", basename(infile), ": ")
        rdo <- parse_Rd(infile)
        from_infile <- TRUE
    }else{
        if(verbose) cat("Rd source not supplied, looking for installed documentation.\n")

        fnam <- if(is.character(object)) object else deparse(substitute(object))

        rdo <- .capture_installed_help(fnam, type = type, package=package)
        if(inherits(rdo,"try-error"))
            cat("Rd source not supplied and installed documentation not found.\n")
        else{
            if(verbose) cat("Installed documentation found, processing it...\n")

            rdo <- .order_sections(rdo) # the sections may not be in canonical order in
            tidyflag <- TRUE            # instaled help
        }
    }

    if(inherits(rdo,"Rd")){                # do the main job: inspect the documentation object
        res <- inspect_Rd(rdo, package = package)
    }else{                                # documentation not found, try to generate fresh one
        if(verbose)
            cat("Trying a 'prompt' function to generate documentation for the object.\n")

                                                    # 2012-11-04 arg. type, package
        res <- .capture_promptAny(fnam, type = type, package = package,
                                  final=final, methods=methods)

        if(inherits(res,"try-error"))
            stop("unsuccessful attempt to create Rd doc. using a 'prompt' function.")
        else if (verbose)
            cat("\tsuccess: documentation generated using a 'prompt' function.\n")
    }

    if(tidyflag)
        res <- .Rd_tidy(res)   # tidy() could do more,

    if(is.null(filename)){   # generate appropriate file name; todo: may need some mangling?
        filename <- if(is.null(infile))
                         paste(res[[ Rdo_which_tag_eq(res, "\\name") ]],
                               ".Rd", sep="")
                    else basename(infile)         # do not overwrite unless in the current dir
    }
                                                                        # todo: error checking
    if(is.character(filename) || identical(filename,FALSE)){               # convert to Rdtext
        res <- Rdo2Rdf(res, ex_restore = TRUE,
                       file = if(is.character(filename)) filename else NULL,
                       srcfile = if(from_infile && sec_copy) infile else NULL )
        if(is.character(filename))
            res <- invisible(filename) # return only the file name in this case
    }
    res
}
                           # (promptMethods) todo: filename = FALSE is a useful
                           # alternative. In that case the text is returned in a named list
                           # containing one element for each Rd section (multiple occurences
                           # of sections like '\alias' are grouped together.
                                                    # 2012-11-04 new arg. type, package
.capture_promptAny <- function(fnam, type, package, final, ..., methods){
              # 2012-11-04 promenyam za da raboti is replacement methods, e.g. "[<--methods"
              #
              # if(grepl("^([^-]+)-.*", fnam)){                 # fnam is of the form xxxx-yyy
              #     fname  <- gsub("^([^-]+)-.*", "\\1", fnam)
              #     type   <- gsub("^([^-]+)-(.*)", "\\2", fnam)   # without "-"
              #     ## suffix <- gsub("^([^-]+)(-.*)", "\\2", fnam)   # with "-"
              # }else{
              #     fname <- fnam
              #     type = ""
              # }

           # 2012-11-04 replacing with the code after the comments
           #
           # namreg <- "^(.+)-([^-]+)$"
           # if(grepl(namreg, fnam)){           # fnam is of the form xxxx-yyy (non-empty rhs)
           #     fname  <- gsub(namreg, "\\1", fnam)
           #     namtype   <- gsub(namreg, "\\2", fnam)   # without "-"
           #     ## suffix <- gsub("^([^-]+)(-.*)", "\\2", fnam)   # with "-"
           # }else{
           #     fname <- fnam
           #     namtype = ""
           # }

    wrknam <- .parse_long_name(fnam)
    fname <- wrknam["name"]
    namtype <- wrknam["type"]

    if(missing(type) || is.null(type))
        type <- namtype
    else if(namtype != ""  &&  !identical(type, namtype)){
        cat("The name and type arguments give conflicting 'type' information.\n")
        cat("\tUsing argument 'type'.\n")
    }# else 'type'  has the value needed.

    wrk <- try(switch(type,
                      methods = {
                          if(is.null(methods) && !is.null(package))
                              methods <- findMethods(fname, where = asNamespace(package))

                          if(is.null(methods)) promptMethods(f=fname, filename = NA)
                          else          promptMethods(f=fname, filename = NA, methods=methods)
                      },
                      class   = promptClass(clName=fname, filename = NA),
                      package = promptPackageSexpr(fname, filename = NA),
                      ## default

                      # v tozi variant ima problemi za funktsii ot Rdpack, za koito parviyat
                      # "if" dava TRUE i sled tova stava greshka. Za funktsii ot drugi paketi
                      # tova ne e problem, ponezhe za tyach "if"-at dava FALSE, ako sa
                      # nevidimi.
                      #
                      # Tay kato tazi situatsiya mozhe da vaznikne po razlichni nachini,
                      # promenyam koda. Tryabva oste rabota za sluchaya kogato ima poveche ot
                      # edno ime...
                      ### if(exists(fname, envir = parent.frame())){
                      ###     prompt(object=fname, filename = NA, force.function=TRUE, ...)
                      ### }else{ # todo: needs more work here
                      ###     x0 <- do.call(getAnywhere,list(fname))
                      ###     browser()
                      ###     prompt(object=x0$objs[[1]], filename = NA, force.function=TRUE,
                      ###            name = fname, ...)
                      ### }
                      {
                          wrk0 <- try(prompt(name=fname, filename = NA, ...), silent=TRUE)
                          if(inherits(wrk0,"try-error")){
                              x0 <- do.call(getAnywhere,list(fname))
                              wrk0 <- prompt(object=x0$objs[[1]], filename=NA, name = fname,
                                             # force.function=TRUE,
                                             ...)
                          }
                          wrk0 # todo: needs more work here. IN particular, there should be a
                           #       package argument to avoid taking blindly whatever comes up.
                      }   )
               , silent = TRUE)

    if(inherits(wrk,"try-error"))
        res <- wrk
    else{
        res <- .parse_Rdlines(wrk)
                                          # if successful, res is not inspected here
                                          # since it is generated from the actual definitions.
        if(final && type != "package"){ # put dummy title and description
                                        # to avoid errors when installing a package
            wrk <- char2Rdpiece("~~ Dummy title ~~", "title")
            res <- Rdo_replace_section(res, wrk)

            wrk <- char2Rdpiece("~~ Dummy description ~~", "description")
            res <- Rdo_replace_section(res, wrk)

                                       # tidy a bit, e.g. to start each section on new line,
                                       # which may not be the case for installed documentation
            res <- .Rd_tidy(res)   # tidy() could do more; e.g. reorder sections
        }
    }

    res
}

.capture_installed_help <- function(fnam, type = NULL, package = NULL, suffix = NULL){
                         # this does not work, package seems not evaluated or deparsed
                         #      hlp <- help(paste(fnam, "-methods", sep=""), package=package)
                         # TODO: the last example in "help()" amy be helpful here.
                         #
    namreg <- "^(.+)-([^-]+)$"                             # 2012-11-04 new namreg and related
    fullname <- if(grepl(namreg, fnam))   # fnam is of the form xxxx-yyy
                    fnam
                else if(!is.null(type) && is.character(type) && type!="")
                    paste(fnam, "-", type, sep="")
                else if(!is.null(suffix))
                    paste(fnam, suffix, sep="")
                else
                    fnam

    hlp <- help(fullname)                                                   # todo: more care!
    hlpfile <- as.vector(hlp) # removes attributes
                                    # todo: but may be of length > 1,  e.g. for "initialize"
                                    #   cat("hlpfile has ", length(hlpfile), " element(s):\n")
                                    #   print(hlpfile)
    if(!is.null(package)){                  # try first "/package/" to avoid surprise matches;
                             # see what happens with package = "methods", without the slashes;
                             # also, 'package' may be part of the name of another package.
        indx <- which( grepl(paste("/", package, "/", sep=""), hlpfile))
        if(length(indx)==0)             # ... but if nothing matched, try without the slashes.
            indx <- which(grepl(package, hlpfile))
        hlpfile <- hlpfile[ indx ]
    }

    if(length(hlpfile) > 1){    # todo: mozhe da se napravi v loop to collect a bunch of sig's
        hlpfile <- hlpfile[1]
        cat("length(hlpfile)>1, taking the first element.\n")
    }

    try(utilsdotdotdot.getHelpFile(hlpfile), silent=TRUE)
}

                                 # 'usage' may be an "f_usage" object obtained e.g. by a
                                 # previous call to get_usage() or generated programmatically.
promptUsage <- function (..., usage){                          # todo: add formatting options?
    if(missing(usage)) get_usage(..., out.format="text")
    else               as.character(usage)
}
