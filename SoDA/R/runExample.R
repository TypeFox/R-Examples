runExample <- function(what,  where = "SoDA", run = TRUE, ..., echo = TRUE, prompt.echo, wd ) {
    file = exampleFiles(what, where, TRUE, TRUE) # one file, with path
    fileB = sub("^.*/","",file)
    if(length(grep("[.][RS]$", file)) >0 ) {
        if(run) {
            if(missing(wd))
              wd <- system.file(package=where)
            if(!nzchar(wd)) {
                oldwd <- getwd()
                setwd(wd)
                on.exit(setwd(oldwd))
            }
            message(gettextf("Running example file \"%s\"", file),
                    domain = NA)
            if(missing(prompt.echo))
              prompt.echo <- paste(">",abbreviate(substr(fileB, 1, nchar(file)-2),4)," ", sep="")
            source(file, ..., echo = echo, prompt.echo = prompt.echo)
        }
        else {
            value <- parse( file)
            if(length(value) == 1) #e.g., and often, a function definition
              value <- value[[1]]
            value
        }
    }
    else {  # should offer the user some choices?
       file
    }
}

.guessFile <- function(what, files) {
    pattern <- paste(what, ".R", sep="")
    which <- match(pattern, files)
    if(!is.na(which))
      return(files[which])
    pattern <- paste("^",what, sep="")
    which <- grep(pattern, files)
    if(length(which) > 0)
      return(files[which])
    ## should try for other languages, in case regexp failed?
    character()
}
        
        
