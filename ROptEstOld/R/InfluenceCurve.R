## generating function
InfluenceCurve <- function(name, Curve = EuclRandVarList(EuclRandVariable(Domain = Reals())), 
                           Risks, Infos){
    if(missing(name))
        name <- "influence curve"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))

    if(!is(Domain(Curve[[1]]), "EuclideanSpace"))
        stop("The domain of 'Curve' has to be a Euclidean space")
    if(!is.character(Infos))
        stop("'Infos' contains no matrix of characters")
    for(char in names(Risks))
        if(!extends(char, "RiskType"))
            stop(paste(char, "is no valid 'RiskType'"))
    if(ncol(Infos)!=2)
        stop("'Infos' must have two columns")
    
    IC1 <- new("InfluenceCurve")
    IC1@name <- name
    IC1@Curve <- Curve
    IC1@Risks <- Risks
    IC1@Infos <- Infos
    
    return(IC1)
}

## access methods
setMethod("name", "InfluenceCurve", function(object) object@name)
setMethod("Curve", "InfluenceCurve", function(object) object@Curve)
setMethod("Risks", "InfluenceCurve", function(object) object@Risks)
setMethod("Infos", "InfluenceCurve", function(object) object@Infos)

## add risk or information
setMethod("addRisk<-", "InfluenceCurve", 
    function(object, value){ 
        if(!is.list(value))
            stop("'value' has to be a list")
        for(char in names(value))
            if(!extends(char, "RiskType"))
                stop(paste(char, "is no valid 'RiskType'"))

        if(length(value) == 1){
            if(names(value) %in% names(object@Risks)){
                ind <- match(names(value), names(object@Risks))
                object@Risks[[ind]] <- c(object@Risks[[ind]], value)
            }else
                object@Risks <- c(object@Risks, value)
        }else{
            for(i in 1:length(value))
                if(names(value[[i]]) %in% names(object@Risks)){
                    ind <- match(names(value[[i]]), names(object@Risks))
                    object@Risks[[ind]] <- c(object@Risks[[ind]], value[[i]])
                }else
                    object@Risks <- c(object@Risks, value[[i]])
        }

        object 
    })
setMethod("addInfo<-", "InfluenceCurve", 
    function(object, value){ 
        object@Infos <- rbind(object@Infos, " " = value) 
        if(length(value)!=2)
            stop("length of 'value' is != 2")
        if(!is.character(value))
            stop("'value' is no vector of characters")
        object 
    })

## replace methods
setReplaceMethod("name", "InfluenceCurve", 
    function(object, value){ 
        object@name <- value 
        addInfo(object) <- c("name<-", "The slot 'name' has been changed")                
        object
    })
setReplaceMethod("Risks", "InfluenceCurve", 
    function(object, value){ 
        object@Risks <- value 
        for(char in names(value))
            if(!extends(char, "RiskType"))
                stop(paste(char, "is no valid 'RiskType'"))
        addInfo(object) <- c("Risks<-", "The slot 'Risks' has been changed")                
        object 
    })
setReplaceMethod("Infos", "InfluenceCurve", 
    function(object, value){ 
        object@Infos <- value 
        if(!is.character(value))
            stop("'value' is no matrix of characters")
        if(ncol(value)!=2)
            stop("'value' has to be a matrix with two columns")
        object
    })

## wrapped access methods
setMethod("Map", "InfluenceCurve", function(f, ...) Map(f@Curve))
setMethod("Domain", "InfluenceCurve", function(object) Domain(object@Curve))
setMethod("Range", "InfluenceCurve", function(object) Range(object@Curve))
