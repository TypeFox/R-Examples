"RweaveHTML" <- function()
{
    list(setup = RweaveHTMLSetup,
         runcode = RweaveHTMLRuncode,
         writedoc = RweaveHTMLWritedoc,
         finish = RweaveHTMLFinish,
         checkopts = RweaveHTMLOptions)
}

"RweaveHTMLSetup" <-
    function(file, syntax,
             output=NULL, quiet=FALSE, debug=FALSE, echo=TRUE,
             eval=TRUE, split=FALSE, cssfile="R2HTML.css",havecss=FALSE,width=500,height=500,border=1,png=TRUE)
{
    # This driver requires R2HTML package to work...
    #if(!require(R2HTML)) stop("R2HTML package is required.")
    if(is.null(output)){
        prefix.string <- basename(sub(syntax$extension, "", file))
        output <- paste(prefix.string, "html", sep=".")
    }
    else{
        prefix.string <- basename(sub("\\.html$", "", output))
    }
    if(!quiet) cat("Writing to file ", output, "\n",
                   "Processing code chunks ...\n", sep="")
    output <- file(output, open="w+")
    options <- list(prefix=TRUE, prefix.string=prefix.string,
                    engine="R", print=FALSE, eval=eval,
                    fig=FALSE, png=png,width=width, height=height, term=TRUE,
                    echo=echo, results="Robj", split=split,
                    strip.white=TRUE, include=TRUE,align="center",caption=NULL,bg="white",pointsize=12)
    list(output=output, debug=debug, quiet=quiet, syntax = syntax,
         options=options, chunkout=list(),cssfile=cssfile,havecss=havecss)
}

"RweaveHTMLRuncode" <- function(object, chunk, options)
{
    if(!(options$engine %in% c("R", "S"))) return(object)
    if(!object$quiet){
        cat(formatC(options$chunknr, width=2), ":")
        if(options$echo) cat(" echo")
        if(options$eval){
            if(options$print) cat(" print")
            if(options$term) cat(" term")
            cat("", options$results)
            if(options$fig){
                if(options$png) cat(" png")
            }
        }
        if(!is.null(options$label))
            cat(" (label=", options$label, ")", sep="")
        cat("\n")
    }


    #chunkprefix <- utils:::RweaveChunkPrefix(options)
    chunkprefix <- RweaveChunkPrefix(options)

    if(options$split){
        chunkout <- object$chunkout[[chunkprefix]]
        if(is.null(chunkout)){
            chunkout <- file(paste(chunkprefix, "html", sep="."), "w")
            if(!is.null(options$label))
                object$chunkout[[chunkprefix]] <- chunkout
        }
    }
    else
        chunkout <- object$output

    assign(".HTML.file",chunkout,pos=.HTMLEnv, immediate=TRUE)
    #utils:::SweaveHooks(options, run=TRUE)
    SweaveHooks(options, run=TRUE)

    chunkexps <- try(parse(text=chunk), silent=TRUE)
    #utils:::RweaveTryStop(chunkexps, options)
    RweaveTryStop(chunkexps, options)
    openSinput <- FALSE
    openSchunk <- FALSE

    if(length(chunkexps)==0)
        return(object)

    for(nce in 1:length(chunkexps))
    {
        ce <- chunkexps[[nce]]
        #dce <- deparse(ce, width.cutoff=0.75*getOption("width"))
        if(object$debug)
            cat("\nRnw> ", paste(ce, collapse="\n+  "),"\n")
        if(options$echo){
            if(!openSinput){
                if(!openSchunk){
                    cat("<!-- begin{Schunk} !-->\n",
                        file=chunkout, append=TRUE)
                    openSchunk <- TRUE
                }
                cat("<!-- begin{Sinput} !-->",
                    file=chunkout, append=TRUE)
                openSinput <- TRUE
            }
            cat("\n", paste(HTMLCommand(deparse(ce)),
                      collapse=paste("\n", getOption("continue"), sep="")),
                file=chunkout, append=TRUE, sep="")
        }

        # tmpcon <- textConnection("output", "w")
        # avoid the limitations (and overhead) of output text connections
         tmpcon <- file()
         sink(file=tmpcon)
        err <- NULL
        #if(options$eval) err <- utils:::RweaveEvalWithOpt(ce, options)
        if(options$eval) err <- RweaveEvalWithOpt(ce, options)
         cat("\n") # make sure final line is complete
         sink()
         output <- readLines(tmpcon)
         close(tmpcon)
        # delete empty output
        if(length(output)==1 & output[1]=="") output <- NULL

        #utils:::RweaveTryStop(err, options) #### !!!  err$value peut etre exporte via HTML(err.value)
        RweaveTryStop(err, options) #### !!!  err$value peut etre exporte via HTML(err.value)

        if(object$debug)
            cat(paste(output, collapse="\n"))

        if(length(output)>0 & (options$results!="hide")){
            if(!openSchunk){
                cat("<!-- begin{Schunk} !--> \n",
                    file=chunkout, append=TRUE)
                openSchunk <- TRUE
            }
            if(openSinput){
                cat("\n<!-- end{Sinput} !-->\n", file=chunkout, append=TRUE)
                openSinput <- FALSE
            }
            if (options$results=="Robj") HTML(err$value, file=chunkout, append=TRUE)
            if (options$results=="html") cat(err$value, file=chunkout, append=TRUE)
            remove(output)

        }
    }
    if(openSinput){
        cat("\n<!--\\end{Sinput}!-->\n", file=chunkout, append=TRUE)
    }
    if(openSchunk){
        cat("\n<!--\\end{Schunk}!-->\n", file=chunkout, append=TRUE)
    }

    if(is.null(options$label) & options$split)
        close(chunkout)

    if(options$fig && options$eval){
        if(options$png){
            png(filename=paste(chunkprefix, "png", sep="."),width=options$width,height=options$height,bg=options$bg,pointsize=options$pointsize)

            #err <- try({utils:::SweaveHooks(options, run=TRUE);
            err <- try({SweaveHooks(options, run=TRUE);
                        eval(chunkexps, envir=.GlobalEnv)})
            dev.off()
            if(inherits(err, "try-error")) stop(err)
        }
        if(options$include)
            cat("<p align='",options$align,"'><img height=",options$HTMLheight, " width=",options$HTMLwidth," src='", chunkprefix, ".png'",if (!is.null(options$border)) paste("border=",options$border,sep=""),">",if(!is.null(options$caption)) paste("<br><font class='caption='>",options$caption,"</font>",sep=""),"</p>", sep="",
                file=object$output, append=TRUE)
    }
    return(object)
}

"RweaveHTMLWritedoc" <- function(object, chunk)
{
    # Very temporary and ugly fix: importing function definition from
    # latest R source code (r45768)
    InternalSweaveParseOptions <-  function(text, defaults=list(), check=NULL)
    {
    x <- sub("^[[:space:]]*(.*)", "\\1", text)
    x <- sub("(.*[^[:space:]])[[:space:]]*$", "\\1", x)
    x <- unlist(strsplit(x, "[[:space:]]*,[[:space:]]*"))
    x <- strsplit(x, "[[:space:]]*=[[:space:]]*")

    ## only the first option may have no name: the chunk label
    if(length(x)>0){
        if(length(x[[1]])==1){
            x[[1]] <- c("label", x[[1]])
        }
    }
    else
        return(defaults)

    if(any(sapply(x, length)!=2))
        stop(gettextf("parse error or empty option in\n%s", text), domain = NA)

    options <- defaults

    for(k in 1:length(x))
        options[[ x[[k]][1] ]] <- x[[k]][2]

    if(!is.null(options[["label"]]) && !is.null(options[["engine"]]))
        options[["label"]] <- sub(paste("\\.", options[["engine"]], "$",
                                        sep=""),
                                  "", options[["label"]])

    if(!is.null(check))
        options <- check(options)

    options
    }



   if(any(grep("text/css", chunk)))
        object$havecss <- TRUE

    if(!object$havecss){
        if(any(grep("<body>", chunk, ignore.case = TRUE))) chunk <- gsub("<body>",paste("\n<link rel=stylesheet type=text/css href=",object$cssfile,"><body>",sep="") ,chunk,ignore.case=TRUE)
        else {
        	if(any(grep("</head>", chunk, ignore.case = TRUE))) chunk <- gsub("</head>",paste("\n<link rel=stylesheet type=text/css href=",object$cssfile,"></head>",sep="") ,chunk,ignore.case=TRUE)
        	else chunk <- gsub("<html>",paste("<html>","\n<link rel=stylesheet type=text/css href=",object$cssfile,">",sep="") ,chunk,ignore.case=TRUE)
        }
        object$havecss <- TRUE
    }
    while(any(pos <- grep(object$syntax$docexpr, chunk)))
    {
        cmdloc <- regexpr(object$syntax$docexpr, chunk[pos[1]])
        cmd <- substr(chunk[pos[1]], cmdloc,
                      cmdloc+attr(cmdloc, "match.length")-1)
        cmd <- sub(object$syntax$docexpr, "\\1", cmd)
        if(object$options$eval)
            val <- as.character(eval(parse(text=cmd), envir=.GlobalEnv))
        else
            val <- paste("<font class='Rcmd'>", cmd, "</font>", sep="")

        chunk[pos[1]] <- sub(object$syntax$docexpr, val, chunk[pos[1]])
    }
    while(any(pos <- grep(object$syntax$docopt, chunk)))
    {
        opts <- sub(paste(".*", object$syntax$docopt, ".*", sep=""),
                    "\\1", chunk[pos[1]])
        object$options <- InternalSweaveParseOptions(opts, object$options, RweaveHTMLOptions)
        chunk[pos[1]] <- sub(object$syntax$docopt, "", chunk[pos[1]])
    }
    cat(chunk, sep="\n", file=object$output, append=TRUE)
    return(object)
}

"RweaveHTMLFinish" <- function(object, error=FALSE)
{
    if(!object$quiet && !error)
        cat(paste("file ",summary(object$output)$description),"is completed", "\n")
    close(object$output)
    if(length(object$chunkout)>0){
        for(con in object$chunkout) close(con)
    }
}

"RweaveHTMLOptions" <- function(options)
{
    ## convert a character string to logical
    c2l <- function(x){
        if(is.null(x)) return(FALSE)
        else return(as.logical(toupper(as.character(x))))
    }
    NUMOPTS <- c("width", "height")
    NOLOGOPTS <- c(NUMOPTS, "results", "prefix.string",
                   "engine", "label","align","caption","border","height","width","HTMLheight","HTMLwidth","bg","pointsize")
    for(opt in names(options)){
        if(! (opt %in% NOLOGOPTS)){
            oldval <- options[[opt]]
            if(!is.logical(options[[opt]])){
                options[[opt]] <- c2l(options[[opt]])
            }
            if(is.na(options[[opt]]))
                stop(paste("invalid value for", opt, ":", oldval))
        }
        else if(opt %in% NUMOPTS){
            options[[opt]] <- as.numeric(options[[opt]])
        }
    }
    options$results <- match.arg(options$results,c("Robj","html", "hide"))
    options
}

#----------------------------------------------------------------------------------------------------#

SweaveSyntaxHTML <- SweaveSyntaxNoweb
SweaveSyntaxHTML$docexpr <- "<[/]?Sexpr([^>]*)>"
SweaveSyntaxHTML$syntaxname <- "<[/]?SweaveSyntax([^>]*)>"
SweaveSyntaxHTML$trans$docexpr <- "<[/]?Sexpr\\1>"
SweaveSyntaxHTML$trans$syntaxname <- "<!--SweaveSyntax{SweaveSyntaxHTML}!-->"
