BsProb.design <- function(design, mFac = NULL, response=NULL, select=NULL, mInt = 2, p = 0.25, g = 2,
    ng = 1, nMod = 10){
 ## accessor function for BsProb from package BsMD
    if (!"design" %in% class(design)) stop("This function is for class design objects only.")
    di <- design.info(design)
    if (!(length(grep("FrF2",di$type))>0 |
           length(grep("pb",di$type))>0)) {
           if (!(di$type=="full factorial" & all(di$nlevels==2)))
        stop("The design obj must be of a type containing FrF2 or pb.")
        }
   if (is.null(di$response.names)) stop("The design must have at least one response.")
   if (is.null(mFac)) mFac <- di$nfactors   ## modified later, if select is used
    
   if (!is.null(select)){
       ## take care of selecting only part of the factors
       if (!is.numeric(select)) stop("select must be numeric")
       if (!all(floor(select)==select)) stop("select must contain integer numbers")
       if (any(select<1 | select>di$nfactors)) stop("select must contain numbers betweeen 1 and ", di$nfactors, " only")
       select <- unique(select)
       if (length(select)<2) stop("at least 2 factors must be selected")
       faclab <- names(di$factor.names)[select]
       if (mFac>length(select)) mFac <- length(select)
       }
    else faclab <- names(di$factor.names)
    
    X <- desnum(design)[,faclab]
    if (!is.null(response)){
      if (!response %in% di$response.names)
        stop("Requested response is not a response variable in design.")
        }
    else response <- di$response.names[1] ## use first response variable per default
    y <- unlist(design[,response])
    
    BsMD::BsProb(X, y, mFac=mFac, mInt=mInt, p=p, g=g, ng=ng, nMod=nMod)
}