###############################################################################
## Determine estimates by minimizing a given criterion
###############################################################################
MCEstimator <- function(x, ParamFamily, criterion, crit.name, 
                        startPar = NULL, 
                        Infos, trafo = NULL, penalty = 1e20, validity.check = TRUE,
                        asvar.fct, na.rm = TRUE, ..., .withEvalAsVar = TRUE){

    ## preparation: getting the matched call
    es.call <- match.call()
    dots <- match.call(expand.dots = FALSE)$"..."

    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

    ## some checking
    if(!is.numeric(x))
      stop(gettext("'x' has to be a numeric vector"))   
    if(!is(ParamFamily, "ParamFamily"))
      stop(gettext("'ParamFamily' has to be of class 'ParamFamily'"))
    if(!is.function(criterion))
      stop(gettext("'criterion' has to be a function"))

    ## manipulation of the arg list to method mceCalc
    argList <- c(list(x = x, PFam = ParamFamily, criterion = criterion, 
                   startPar = startPar, penalty = penalty))

    if(missing(validity.check)) validity.check <- TRUE
       argList$validity.check <- validity.check
    if(missing(Infos))      Infos <- NULL
    argList <- c(argList, Infos = Infos)
    if(missing(crit.name)) crit.name <- ""               
    argList <- c(argList, crit.name = crit.name)               
    if(!is.null(dots))      argList <- c(argList, dots)

    ## call to mceCalc
    res0 <- do.call(mceCalc, argList)
    
    asv <- if("FisherInfo" %in% slotNames(ParamFamily)){
              function(ParamFamily, param)
                                  solve(FisherInfo(ParamFamily, param = param))
           }else NULL
    
    argList <- c(list(res0, PFam = ParamFamily, 
                              trafo = trafo, 
                              res.name = paste("Minimum", crit.name, 
                                               "estimate", sep=" ", collapse=""), 
                              call = quote(es.call),
                              .withEvalAsVar=.withEvalAsVar))

    if(!is.null(asv))   argList <- c(argList, asvar.fct = asv)
    if(!is.null(dots))  argList <- c(argList, dots)

    ## digesting the results of mceCalc
    res <- do.call(.process.meCalcRes, argList)
    res@completecases <- completecases
    
    return(res)
}
