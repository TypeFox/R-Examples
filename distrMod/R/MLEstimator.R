###############################################################################
## Implementation of Maximum Likelihood estimation in i.i.d. setup
###############################################################################


## Maximum-Likelihood estimator
MLEstimator <- function(x, ParamFamily, startPar = NULL, 
                        Infos, trafo = NULL, penalty = 1e20,
                        validity.check = TRUE, na.rm = TRUE,
                        ..., .withEvalAsVar = TRUE){

    ## preparation: getting the matched call
    es.call <- match.call()
    dots <- match.call(expand.dots = FALSE)$"..."

    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

    ## some checking
    if(!is.numeric(x))
      stop(gettext("'x' has to be a numeric vector"))   
    if(is.null(startPar)) startPar <- startPar(ParamFamily)(x,...)


    ## manipulation of the arg list to method mceCalc
    argList <- c(list(x = x, PFam = ParamFamily, startPar = startPar, 
                      penalty = penalty))

    if(missing(validity.check)) validity.check <- TRUE
       argList$validity.check <- validity.check
    if(missing(Infos))      Infos <- NULL
        argList <- c(argList, Infos = Infos)
    if(!is.null(dots))      argList <- c(argList, dots)

    ## call to mleCalc
    res0 <- do.call(mleCalc, argList)

    asv <- if("FisherInfo" %in% slotNames(ParamFamily)){
              function(PFam = ParamFamily, param, ...)
                                  solve(FisherInfo(PFam, param = param))
           }else NULL
    
    argList <- list(res0, PFam = ParamFamily, trafo = trafo,
                      res.name = "Maximum likelihood estimate",
                      call = quote(es.call), .withEvalAsVar=.withEvalAsVar)

    if(!is.null(asv))   argList <- c(argList, asvar.fct = asv)
    if(!is.null(dots))  argList <- c(argList, dots)

    ## digesting the results of mceCalc
    res <- do.call(what = ".process.meCalcRes", args = argList)
    
    names(res@criterion) <- "negative log-likelihood"
    res@estimate.call <- es.call
    res@name <- "Maximum likelihood estimate"
    res@completecases <- completecases
    
    return(res)
}
