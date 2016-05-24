# make data plotable

#' @title buildModelCI
#' 
#' @description Construct Confidence Interval Values
#' 
#' @details Takes a model and builds a data.frame holding the coefficient value and the confidence interval values.
#'
#' @author Jared P. Lander
#' @aliases buildModelCI
#' @export buildModelCI
#' @import plyr
#' @param model A Fitted model such as from lm, glm
#' @param \dots Arguments passed on onto other methods
#' @return A \code{\link{data.frame}} listing coeffcients and confidence bands.
#' @seealso \code{\link{coefplot}} \code{\link{multiplot}}
#' @examples
#'
#' data(diamonds)
#' model1 <- lm(price ~ carat + cut, data=diamonds)
#' coefplot:::buildModelCI(model1)
#' coefplot(model1)
#'
buildModelCI <- function(model, ...)
{
    UseMethod(generic="buildModelCI")
}

#' @title buildModelCI.default
#' 
#' @description Construct Confidence Interval Values
#' 
#' @details Takes a model and builds a data.frame holding the coefficient value and the confidence interval values.
#'
#' @author Jared P. Lander
#' @aliases buildModelCI.default
#' @export buildModelCI.default
#' @export
#' @method buildModelCI default
#' @import plyr
#' @param model A Fitted model such as from lm, glm
#' @param innerCI How wide the inner confidence interval should be, normally 1 standard deviation.  If 0, then there will be no inner confidence interval.
#' @param outerCI How wide the outer confidence interval should be, normally 2 standard deviations.  If 0, then there will be no outer confidence interval.
#' @param sort Determines the sort order of the coefficients.  Possible values are c("natural", "magnitude", "alphabetical")
#' @param decreasing logical; Whether the coefficients should be ascending or descending
#' @param predictors A character vector specifying which variables to keep.  Each individual variable has to be specfied, so individual levels of factors must be specified.  We are working on making this easier to implement, but this is the only option for now.
#' @param coefficients A character vector specifying which factor variables to keep.  It will keep all levels and any interactions, even if those are not listed.
#' @param strict If TRUE then predictors will only be matched to its own coefficients, not its interactions
#' @param newNames Named character vector of new names for coefficients
#' @param numeric logical; If true and factors has exactly one value, then it is displayed in a horizontal graph with constinuous confidence bounds.; not used for now.
#' @param intercept logical; Whether the Intercept coefficient should be plotted
#' @param interceptName Specifies name of intercept it case it is not the default of "(Intercept").
#' @param \dots See Details for information on \code{factors}, \code{only} and \code{shorten}
#' @param name A name for the model, if NULL the call will be used
#' @return A \code{\link{data.frame}} listing coeffcients and confidence bands.
#' @seealso \code{\link{coefplot}} \code{\link{multiplot}}
#' @examples
#'
#' data(diamonds, package='ggplot2')
#' model1 <- lm(price ~ carat + cut, data=diamonds)
#' coefplot:::buildModelCI(model1)
#' coefplot(model1)
#'
buildModelCI.default <- function(model, outerCI=2, innerCI=1, intercept=TRUE, numeric=FALSE, 
                         sort=c("natural", "magnitude", "alphabetical"), predictors=NULL, strict=FALSE, coefficients=NULL, 
                         newNames=NULL,
                         decreasing=TRUE, name=NULL, interceptName="(Intercept)", ...)
{
    sort <- match.arg(sort)

    #print(structure(as.list(match.call()[-1]), class = "uneval")$model)
    # get model information
    modelCI <- extract.coef(model, ...)
    
    # if the user has specified predictors calculate which coefficient they go with
    keptCoefsFromPredictors <- getCoefsFromPredictors(model=model, predictors=predictors, strict=strict)
    
    # if individual coefficients were specified use them, if not it will be null
    keptCoefsFromCoefficients <- coefficients
    
    if(!is.null(predictors) || !is.null(coefficients))
    {
        modelCI <- modelCI[modelCI$Coefficient %in% unique(c(keptCoefsFromPredictors, keptCoefsFromCoefficients)), ]
    }
    
    if(NROW(modelCI) == 0)
    {
        return(NULL)
    }
    
    if(!is.null(newNames))
    {
        modelCI$Coefficient <- revalue(x=modelCI$Coefficient, replace=newNames, warn_missing=FALSE)
    }
    
    # build confidence bounds columns
#     modelCI <- within(modelCI, {LowOuter <- Value - outerCI*SE;
#                   HighOuter <- Value + outerCI*SE;
#                   LowInner <- Value - innerCI*SE;
#                   HighInner <- Value + innerCI*SE})
    
    modelCI[, c("HighInner", "LowInner", "HighOuter", "LowOuter")] <- modelCI$Value + modelCI$SE %*% matrix(c(innerCI, -innerCI, outerCI, -outerCI), nrow=1)
    
    # get rid of SE column
    modelCI$SE <- NULL
    
    # if no intercept is desired, remove it
    if(!intercept)
    {
        modelCI <- modelCI[rownames(modelCI) != interceptName, ]
    }

    # make column for model name
    # if a name for the model is provided, use it, otherwise use the call
    if(is.null(name))
    {
        #modelCI$Model <- as.character(paste(model$call, collapse="_"))
        modelCI$Model <- paste(as.character(structure(as.list(match.call()[-1]), class = "uneval")$model), collapse="")
    }else
    {
        modelCI$Model <- name
    }

    ## possible orderings of the coefficients
    ordering <- switch(sort,
                       natural=order(1:nrow(modelCI), decreasing=decreasing), 	# the way the data came in
                       magnitude=order(modelCI$Value, decreasing=decreasing), 		#  size order
                       alphabetical=order(modelCI$Coefficient, decreasing=decreasing), 	# alphabetical order
                       order(1:nrow(modelCI))		# default, the way it came in
    )
    
    # implement the ordering
    modelCI <- modelCI[ordering, ]
    modelCI$Coefficient <- factor(modelCI$Coefficient, levels=modelCI$Coefficient)
    
    return(modelCI)
}


# #' Melt the modelCI
# #'
# #' Melt a modelCI into a form suitable for plotting
# #'
# #' \code{\link{buildModelCI}} builds a data.frame for plotting.  This function melts it into plottable form and seperates the coefficient data from the SE data into seprate data.frames
# #'
# #' @author Jared P. Lander www.jaredlander.com
# #' @aliases meltModelCI
# #' @seealso \code{\link{coefplot}} \code{\link{buildModelCI}}
# #' @param modelCI A \code{\link{data.frame}} as built by \code{\link{buildModelCI}}
# #' @param keepCols The columns in modelCI that should be kept as there can be extras
# #' @param id.vars The columns to use as ID variables in \code{\link{melt}}
# #' @param variable.name Used in \code{\link{melt}} for naming the column that stores the melted variables
# #' @param value.name Used in \code{\link{melt}} for naming the column that stores the melted values
# #' @param innerCols The columns to be included in the \code{\link{data.frame}} of inner standard errors
# #' @param outerCols The columns to be included in the \code{\link{data.frame}} of outer standard errors
# #' @return A list consisting of
# #' \item{modelMelt}{Melted modelCI with all values}
# #' \item{modelMeltOuter}{modelMelt with only values associated with the outer standard errors}
# #' \item{modelMeltInner}{modelMelt with only values associated with the inner standard errors}
# #' @examples
# #'
# #' data(diamonds)
# #' model1 <- lm(price ~ carat + cut, data=diamonds)
# #' \dontrun{modeled <- coefplot:::buildModelCI(model1)
# #' coefplot:::meltModelCI(modeled)}
# #' coefplot(model1)
# #'
# meltModelCI <- function(modelCI, 
#                         keepCols=c("LowOuter", "HighOuter", "LowInner", "HighInner", "Coefficient", "Value", "Model"), 
#                         id.vars=c("Coefficient", "Model"), variable.name="Type", 
#                         value.name="Value", outerCols=c("LowOuter", "HighOuter"), 
#                         innerCols=c("LowInner", "HighInner"))
# {
#     # melt the data frame so it is suitable for ggplot
#     modelMelt <- melt(data=modelCI[, keepCols], id.vars=id.vars, variable.name=variable.name, value.name=value.name)
#     
#     # just the outerCI info
#     modelMeltOuter <- modelMelt[modelMelt$Type %in% outerCols, ]	# pull out the outer (95% default) CI
#     
#     # just the innerCI info
#     modelMeltInner <- modelMelt[modelMelt$Type %in% innerCols, ]	# pull out the inner (68% default) CI
#     
#     # return the data.frames
#     return(list(modelMelt=modelMelt, modelMeltOuter=modelMeltOuter, modelMeltInner=modelMeltInner))
# }
