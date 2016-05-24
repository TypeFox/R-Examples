### defining the ics S4-class

### class definition
###
###
setClass("ics",representation(gKurt="numeric", UnMix="matrix", S1="matrix", S2="matrix", S1name="character", S2name="character",
         Scores="data.frame", DataNames="character", StandardizeB="character", StandardizegKurt="logical"))

### checking the validity of an ics object
###
###
setValidity("ics",function(object){
    if(!is(object@gKurt, "numeric")) return("Generalized kurtosis values of ics objects must be numeric")
    if(!(is(object@UnMix, "matrix") )) return("slot 'UnMix' of a ics object must be a numeric matrix")
    if(!(is.numeric(object@UnMix) )) return("slot 'UnMix' of a ics object must be a numeric matrix")
    if(!is(object@S1name, "character") || length(object@S1name)!=1) return("slot 'S1' of a ics object must be the name of a scatter function")

    if(!is(object@S2name, "character") || length(object@S2name)!=1) return("slot 'S2' of a ics object must be the name of a scatter function")
    if(!(is(object@Scores, "data.frame"))) return("slot 'Scores' of a ics object must be a numeric data frame")
    if(!all(sapply(object@Scores, is.numeric))) return("slot 'Scores' of a ics object must be a numeric data frame")
    if(!is(object@DataNames, "character")) return("slot 'DataNames' of a ics object must give the column names of the data matrix")
    if(!is(object@StandardizeB, "character") || length(object@StandardizeB)!=1) return("slot 'StandardizeB' of a ics object must be the name of a standardization method of 'UnMix'")

    if(length(object@gKurt)!=dim(object@UnMix)[2]) return("length of 'gKurt' must correspond to the number of columns of 'UnMix'")
    if(length(object@gKurt)!=dim(object@Scores)[2]) return("length of 'gKurt' must correspond to the number of columns of 'Scores'")
    if(length(object@gKurt)!=length(object@DataNames)) return("length of 'gKurt' must be the same as length of 'DataNames'")
    if(length(object@gKurt)<2) return("at least bivariate data is needed")
    
    if(!is(object@StandardizegKurt, "logical") || length(object@StandardizegKurt)!=1) return("slot 'StandardizegKurt' of a ics object must be 'TRUE' or 'FALSE'")
    
    if(!(is(object@S1, "matrix") )) return("slot 'S1' of a ics object must be a numeric matrix")
    if(!(is.numeric(object@S1) )) return("slot 'S1' of a ics object must be a numeric matrix")
    
    if(!(is(object@S2, "matrix") )) return("slot 'S2' of a ics object must be a numeric matrix")
    if(!(is.numeric(object@S2) )) return("slot 'S2' of a ics object must be a numeric matrix")
    return(TRUE)
})

### show / print method for an ics object
### only general kurtosis and unmixing matrix as output
###


setMethod("show",signature(object="ics"),
function(object)
    {
    tmp <- list(gKurt=object@gKurt,
                UnMix=object@UnMix)
    print(tmp,quote=F)
    invisible(tmp)  
    }
)


### plot method for an ics object
### -> scatterplot matrix, for larger matrices only first and last 3 coordinates are plotted
###

setMethod("plot",signature(x="ics",y="missing"),
function(x,index=NULL,...)
    {
    p<-ncol(x@UnMix)
    if (is.null(index) & p<=6) pairs(x@Scores,...)
    if (is.null(index) & p>6) pairs(x@Scores[,c(1:3,p-2:0)],...)
    if (length(index)==1) stop("index must be NULL or at least a vector of length 2")
    if (length(index)>1) pairs(x@Scores[,index],...)
    }
)

### summary method for an ics object
### -> more detailed, nicer output than print
###

setMethod("summary",signature(object="ics"),
function(object,digits=4)
    {   
    cat("\nICS based on two scatter matrices \n")
    cat("S1: ", object@S1name) 
    cat("\nS2: ",object@S2name)
    cat("\n")
    cat("\nThe generalized kurtosis measures of the components are:\n")
    print(format(round(object@gKurt,digits)),quote=F)
    cat("\n")
    cat("\nThe Unmixing matrix is:\n")
    print(round(object@UnMix,digits))
    invisible(object)
    }
)

### fitted method for an ics object
### only of interest if index is used.
### Otherwise returns just the original data

setMethod("fitted",signature(object="ics"),
function(object, index=NULL)
    {
    p<-ncol(object@Scores)
    if (is.null(index)==FALSE && max(index)>p) stop("undefined columns selected")
    Mix<-solve(object@UnMix)
    if (is.null(index)) index=1:p
    fits<-tcrossprod(as.matrix(object@Scores[,index]), Mix[,index])
    return(as.data.frame(fits))
    }
)

### coef method for an ics object
### extracts the unmixing matrix
### 

setMethod("coef",signature(object="ics"),
function(object)
    {
    return(object@UnMix)
    }
)
