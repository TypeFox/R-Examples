
### Parse and evaluate a formula, return the data as object of class 
### `ModelEnv'

ModelEnvFormula <- function(formula, data = list(), subset = NULL, 
                            na.action = NULL, frame = NULL, 
                            enclos = sys.frame(sys.nframe()),
                            other = list(), designMatrix = TRUE,
                            responseMatrix = TRUE,
                            setHook = NULL, ...)
{
  
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")

    ### NA-handling will for the ModelFrame objects later on...
    mf$na.action <- stats::na.pass

    MEF <- new("ModelEnvFormula")
    MEF@formula <- c(ParseFormula(formula, data=data)@formula, other)
    MEF@hooks$set <- setHook
    
    if (is.null(frame)) frame <- parent.frame()

    mf$subset <- try(subset)
    if (inherits(mf$subset, "try-error")) mf$subset <- NULL

    MEF@get <- function(which, data=NULL, frame=parent.frame(), envir = MEF@env)
    {
        if(is.null(data))
            RET <- get(which, envir = envir, inherits=FALSE)
        else{
            oldData <- get(which, envir = envir, inherits=FALSE)
            if (!use.subset) mf$subset <- NULL
            mf$data <- data
            mf$formula <- MEF@formula[[which]]
            RET <- eval(mf, frame, enclos = enclos)
            checkData(oldData, RET)
        }
        return(RET)
    }
    
    MEF@set <- function(which = NULL, data = NULL, frame = parent.frame(),
                        envir = MEF@env)
    {
        if (is.null(which)) which <- names(MEF@formula)
        if (any(duplicated(which)))
            stop("Some model terms used more than once")
        
        for (name in which){
            
            if (length(MEF@formula[[name]]) != 2)
                stop("Invalid formula for ", sQuote(name))
            
            mf$data <- data
            mf$formula <- MEF@formula[[name]]

            if (!use.subset) mf$subset <- NULL
            MF <- eval(mf, frame, enclos = enclos)
            if (exists(name, envir = envir, inherits = FALSE))
                checkData(get(name, envir = envir, inherits = FALSE), MF)
            assign(name, MF, envir = envir)
            mt <- attr(MF, "terms")
            
            ## <FIXME>
            ## maybe we don't want to save input and response
            ## in the cases below?
            ## </FIXME>
            if (name == "input" && designMatrix) {
                assign("designMatrix",
                       model.matrix(mt, data = MF, ...),
                       envir = envir)
            }

            if (name == "response" && responseMatrix) {
                attr(mt, "intercept") <- 0
                assign("responseMatrix",
                       model.matrix(mt, data=MF, ...),
                       envir = envir)
            }
        }
        MEapply(MEF, MEF@hooks$set, clone=FALSE)
    }

    use.subset <- TRUE
    MEF@set(which = NULL, data = data, frame = frame)
    use.subset <- FALSE
    
    ### handle NA's
    if (!is.null(na.action))
        MEF <- na.action(MEF)
    MEF
}

### compare basic properties of two data.frames

checkData <- function(old, new) {

    if (!is.null(old)){
        
        if(!all(names(old) %in% names(new)))
            stop("New data must contain the same columns as the original data")
        
        if (!identical(lapply(old, class), lapply(new[names(old)], class)))
            stop("Classes of new data do not match original data")

        if (!identical(lapply(old, levels), lapply(new[names(old)], levels)))
            stop("Levels in factors of new data do not match original data")
    }
}
  
### parse a formula and return the different pieces as `FormulaParts'
### object

ParseFormula <- function(formula, data = list()) {

    formula <- terms(formula, data = data)
    attributes(formula) <- NULL

    if (length(formula) == 3) {
        fresponse <- formula[c(1,2)]
        frhs <- formula[c(1,3)]
        ### if (frhs[[2]] == "1") frhs <- NULL
    }
  
    if (length(formula) == 2) {
        fresponse <- NULL   
        frhs <- formula
    }
  
    finput <- frhs
    fblocks <- frhs

    ### <FIXME>
    ### will fail for `y ~ . | blocks' constructs
    ### </FIXME>

    if (!is.null(frhs) && length(frhs[[2]]) > 1) {
        if (deparse(frhs[[2]][[1]]) == "|") {
            finput[[2]] <- frhs[[2]][[2]]
            fblocks[[2]] <- frhs[[2]][[3]]
        } else {
            fblocks <- NULL
        }
    } else {
        fblocks <- NULL
    }
  
    RET = new("FormulaParts")
  
    RET@formula$response <- fresponse
    RET@formula$input <- finput
    RET@formula$blocks <- fblocks

    return(RET)
}


###**********************************************************

## A simple model environment where designMatrix and responseMatrix
## are directly specified. Usefull for models without a formula
## interface. This is much more limited than ModelEnvFormula, but can
## be faster because no formula parsing is necessary. The subset
## argument needs to be a indexing vector into the design and response
## matrix, respectively. Funny things may happen if the matrices have
## no column names and the @[gs]et slots are used in combination with
## new data <FIXME>is proper handling of that case possible?</FIXME>

ModelEnvMatrix <- function(designMatrix=NULL, responseMatrix=NULL,
                           subset = NULL, na.action = NULL,
                           other=list(), ...)
{    
    MEM <- new("ModelEnv")

    N <- max(nrow(designMatrix), nrow(responseMatrix))
    
    if(is.null(subset) && N>0) subset <- 1:N
    
    if(!is.null(designMatrix))
        assign("designMatrix",
               as.matrix(designMatrix)[subset,,drop=FALSE],
               envir = MEM@env)
    
    if(!is.null(responseMatrix))
        assign("responseMatrix",
               as.matrix(responseMatrix)[subset,,drop=FALSE],
               envir = MEM@env)

    for(n in names(other)){
        if(is.matrix(other[[n]]))
            assign(n,
                   other[[n]][subset,,drop=FALSE],
                   envir = MEM@env)
        else
            assign(n,
                   other[[n]][subset],
                   envir = MEM@env)
    }
    
    MEM@get <- function(which, data=NULL, frame=NULL, envir = MEM@env)
    {
        if(is.null(data))
            RET <- get(which, envir = envir, inherits=FALSE)
        else
        {    
            if(is.null(colnames(data)))
                colnames(data) <- createColnames(data)
        
            oldNames <- colnames(get(which, envir = envir,
                                     inherits=FALSE))
            RET <- data[,oldNames,drop=FALSE]
        }
        return(RET)
    }
    
    MEM@set <- function(which = NULL, data = NULL, frame=NULL,
                        envir = MEM@env)
    {
        if(is.null(which))
            which <- c("designMatrix", "responseMatrix")
        
        if(is.null(data))
            stop("No data specified")
        
        if (any(duplicated(which)))
            stop("Some model terms used more than once")
        
        if(is.null(colnames(data)))
            colnames(data) <- createColnames(data)
        
        for (name in which){
            
            oldNames <- colnames(get(name, envir = envir,
                                     inherits=FALSE))
            
            assign(name, as.matrix(data[,oldNames,drop=FALSE]),
                   envir = envir)            
        }
    }

    ## handle NA's
    if (!is.null(na.action))
        MEM <- na.action(MEM)
    MEM
}

## Make sure that every matrix has column names
createColnames <- function(data)
{
    paste("V",1:ncol(data),sep=".")
}
