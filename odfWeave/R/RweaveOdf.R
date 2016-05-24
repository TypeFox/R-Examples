RweaveOdf <- function()
{
    list(setup = RweaveOdfSetup,
         runcode = RweaveOdfRuncode,
         writedoc = RweaveOdfWritedoc,
         finish = RweaveOdfFinish,
         checkopts = RweaveOdfOptions)
}

RweaveOdfSetup <-
    function(file, syntax,
             output=NULL, quiet=FALSE, debug=FALSE, echo=TRUE,
             eval=TRUE, ...)
{
    theDots <- list(...)
    # pass the control object here
    
    if(is.null(output)){
        prefix.string <- basename(sub(syntax$extension, "", file))
        output <- paste(prefix.string, "xml", sep=".")
    }
    else{
        prefix.string <- basename(sub("\\.xml$", "", output))
    }
    if(!quiet) cat("  Writing to file ", output, "\n",
                   "  Processing code chunks ...\n", sep="")

    # Since we always pass encoding="UTF-8" to the Sweave function,
    # the encoding attribute should be set to "UTF-8" if any non-ASCII
    # characters were found in the input XML file, or to "ASCII" if not.
    # I'll issue a warning message if it is something else, and
    # use "UTF-8" encoding for the output file in any case.
    encoding <- attr(file, "encoding")
    if (! encoding %in% c("UTF-8", "ASCII"))
        warning("unexpected encoding attribute found on file: ", encoding)
    output <- file(output, open="w+", encoding="UTF-8")

    options <- list(prefix=TRUE, prefix.string=prefix.string,
                    engine="R", print=FALSE, eval=eval,
                    fig=FALSE, term=TRUE,
                    echo=echo, results="verbatim", 
                    strip.white="true")

    ## to be on the safe side: see if defaults pass the check
    options <- RweaveOdfOptions(options)

    # now attach the control arg passed by ... (if any) to options
    if("control" %in% names(theDots)) options$control <- theDots$control
    
    list(output=output, 
         debug=debug, quiet=quiet, syntax = syntax,
         options=options, chunkout=list())
}

RweaveOdfRuncode <- function(object, chunk, options, control)
{

    if(!(options$engine %in% c("R", "S"))){
        return(object)
    }

    if(!object$quiet){
        cat("  ", formatC(options$chunknr, width=2), ":")
        if(options$echo) cat(" echo")
        if(options$eval){
            if(options$print) cat(" print")
            if(options$term) cat(" term")
            cat(" ", options$results, sep="")
        }
        if(!is.null(options$label))
            cat("(label=", options$label, ")", sep="")
        cat("\n")
    }

    chunkprefix <- RweaveChunkPrefix(options)
    chunkout <- object$output
    SweaveHooks(options, run=TRUE)
    chunkexps <- try(parse(text=chunk), silent=TRUE)
    RweaveTryStop(chunkexps, options)

    if(length(chunkexps)==0)
        return(object)

    rCont <- odfTranslate(getOption("continue"), toR = FALSE)
    rPrompt <- odfTranslate(getOption("prompt"), toR = FALSE)
   
    # put this in a function
    codeMarkup <- RCodeTags()
    endTag <- "</text:p>" 

    for(nce in 1:length(chunkexps))
    {
        ce <- chunkexps[[nce]]
        dce <- deparse(ce, width.cutoff=0.75*getOption("width"))
        
        #convert some character, such as < or &   
        dceForXml <- odfTranslate(dce, toR = FALSE)

            
        if(object$debug)
            cat("\nRnw> ", paste(dce, collapse="\n+  "),"\n")
            
        # this block will print the R code (if echo = TRUE)            
        if(options$echo)
        {
            dceForXml2 <- paste(ifelse(seq(along = dceForXml) == 1, rPrompt, rCont), dceForXml)
            # now wrap this result in xml tags before printing
            # using style  names
            taggedDce <- paste(codeMarkup$input, dceForXml2, endTag, "\n", sep = "")               
            
            cat("\n", 
                paste(taggedDce, collapse=""),
                file=chunkout, append=TRUE, sep="")
        }

        # evaluate the R code and write to a file
        tmpcon <- file()
        sink(file=tmpcon)
        err <- NULL
        if(options$eval) err <- RweaveEvalWithOpt(ce, options)
#     cat("\n") # make sure final line is complete
        sink()
        
        # read them back in
        output <- readLines(tmpcon)
        close(tmpcon)
        ## delete empty output
        if(length(output)==1 & output[1]=="") output <- NULL

        RweaveTryStop(err, options)

        if(object$debug)
            cat(paste(output, collapse="\n"))
        # write the output to the file
        if(length(output)>0 & (options$results != "hide"))
        {
            if(options$results == "verbatim")
            {
               taggedOutput <- paste(codeMarkup$output, odfTranslate(output, toR = FALSE), endTag, "\n", sep = "")             
               output <- paste(taggedOutput,collapse="\n")
            }
# I'll have to find an example of when this matters            
            
#            if(options$strip.white %in% c("all", "true")){
#                output <- sub("^[[:space:]]*\n", "", output)
#                output <- sub("\n[[:space:]]*$", "", output)
#                if(options$strip.white=="all")
#                    output <- sub("\n[[:space:]]*\n", "\n", output)
#            }
            cat(output, file=chunkout, append=TRUE)
            remove(output)
        }
    }

    if(options$fig && options$eval && ! isTRUE(.odfEnv$fig.cancel))
    {
         deviceInfo <- getImageDefs()
         
         imageName <- paste(
            get("picPath", envir = .odfEnv),
            "/",
            chunkprefix,
            ".",
            deviceInfo$type,
            sep = "")

         figGen(plotName = imageName)    
            
         err <- try({SweaveHooks(options, run=TRUE);
                     eval(chunkexps, envir=.GlobalEnv)})
                     
         dev.off()
         
         if(inherits(err, "try-error")) stop(err) 
         
         plotMarkup <- odfInsertPlot(
            imageName, 
            name = gsub("-", "", chunkprefix),            
            height = deviceInfo$dispHeight, 
            width = deviceInfo$dispWidth,
            caption = .odfEnv$fig.caption)
         cat(plotMarkup, file=chunkout, append=TRUE)
    }

    # set fig.caption and fig.cancel to NULL, even if "fig" was false
    .odfEnv$fig.caption <- NULL            
    .odfEnv$fig.cancel <- NULL            

    return(object)
}

RweaveOdfWritedoc <- function(object, chunk)
{

    while(any(pos <- grep(object$syntax$docexpr, chunk)))
    {
        cmdloc <- regexpr(object$syntax$docexpr, chunk[pos[1]])
        cmd <- substr(chunk[pos[1]], cmdloc,
                      cmdloc+attr(cmdloc, "match.length")-1)
        cmd <- sub(object$syntax$docexpr, "\\1", cmd)
        if(object$options$eval){
            val <- as.character(eval(parse(text=cmd), envir=.GlobalEnv))
            ## protect against character(0), because sub() will fail
            if(length(val)==0) val <- ""
        }
        else {
            val <- cmd
         }
        chunk[pos[1]] <- sub(object$syntax$docexpr, val, chunk[pos[1]])
    }
    while(any(pos <- grep(object$syntax$docopt, chunk)))
    {
        opts <- sub(paste(".*", object$syntax$docopt, ".*", sep=""),
                    "\\1", chunk[pos[1]])
        object$options <- SweaveParseOptions(opts, object$options, RweaveOdfOptions)
        chunk[pos[1]] <- sub(object$syntax$docopt, "", chunk[pos[1]])
    }

    cat(chunk, sep="\n", file=object$output, append=TRUE)
    return(object)
}

RweaveOdfFinish <- function(object, error=FALSE)
{
    if(!object$quiet && !error)
        cat("\n",
            gettextf("  '%s' has been Sweaved",
                     summary(object$output)$description),
            "\n", sep = "")
    close(object$output)
    if(length(object$chunkout) > 0)
        for(con in object$chunkout) close(con)
}

RweaveOdfOptions <- function(options)
{

    ## ATTENTION: Changes in this function have to be reflected in the
    ## defaults in the init function!

    ## convert a character string to logical
    
    # Let's not check the control file
    if("control" %in% names(options))
    {
      controlData <- options$control
      options <- options[names(options) != "control"]
    } else {
       controlData <- NULL
    }

    c2l <- function(x){
        if(is.null(x)) return(FALSE)
        else return(as.logical(toupper(as.character(x))))
    }

    NUMOPTS <- c("width", "height")
    NOLOGOPTS <- c(NUMOPTS, "results", "prefix.string",
                   "engine", "label", "strip.white")

    for(opt in names(options)){
        if(! (opt %in% NOLOGOPTS)){
            oldval <- options[[opt]]
            if(!is.logical(options[[opt]])){
                options[[opt]] <- c2l(options[[opt]])
            }
            if(is.na(options[[opt]]))
                stop(gettextf("invalid value for '%s' : %s", opt, oldval),
                     domain = NA)
        }
        else if(opt %in% NUMOPTS){
            options[[opt]] <- as.numeric(options[[opt]])
        }
    }

    options$results <- tolower(as.character(options$results))
    options$results <- match.arg(options$results,
                                 c("verbatim", "xml", "hide"))

    options$strip.white <- tolower(as.character(options$strip.white))
    options$strip.white <- match.arg(options$strip.white,
                                     c("true", "false", "all"))

    # add the control info back in (if any)
    if(!is.null(controlData)) options$control <- controlData
    
    options
}

RCodeTags <- function()
{
   styleNames <-  getStyles()[c("input", "output")]
   styleTag <- ifelse(
      unlist(styleNames) != "",
      paste("text:style-name=\"", styleNames, "\"", sep = ""),
      "")
   startTags <- as.list(paste("<text:p ", styleTag, ">", sep = ""))
   names(startTags) <- names(styleNames)
   startTags    
}

