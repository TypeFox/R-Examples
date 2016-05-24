rsm.design <- function(design, fo){
    ## this is not a method for function design, but a mere function
    ## design contains the coding information
    
    ## the idea was to translate a user-specified formula into the format 
    ## needed by rsm, but this is harder than expected
    
    ## this does not work yet, unless fo is missing
    ## therefore not exported or documented
    if (!"design" %in% class(design))
        stop("function rsm.design works for class design only.")
    di <- design.info(design)
    if (is.null(di$response.names)) stop("design does not have a response variable")
    descoded <- design
    if (!is.null(di$coding)){
         coded <- sapply(di$coding, function(obj) strsplit(as.character(obj)," ")[[2]])
         uncoded <- sapply(di$coding, function(obj) strsplit(as.character(obj)," ")[[3]][1])
         descoded <- coded.data(design, formulas = di$coding)
         if (!missing(fo))
            fo <- as.formula(gsub(uncoded, coded, as.character(fo), fixed=TRUE))
         }
    if (missing(fo)) {
      if (is.null(di$coding)) fo <- paste(di$response.names[1],"~ SO(",
                  paste(names(di$factor.names),collapse=","), ")")
        else fo <- paste(di$response.names[1],"~ SO(",
                  paste(coded, collapse=","), ")")
          }
    ## trying to generate a valid call to rsm
    ## issue is that data cannot be understood within rsm,
    ##    because the function does not work with the "on the fly" data it gets
    ##    from here
    ## an analogous call to lm would work
    aus <- eval(parse(text=paste("rsm(",fo, ", data = descoded)")))
    aus
}
