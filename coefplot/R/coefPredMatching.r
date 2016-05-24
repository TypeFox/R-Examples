# 
# baseball <- data.frame(Bat=sample(1:100, 20, replace=T), Batter=sample(c("David", "Batley", "Bob", "Ace"), 20, replace=T), Hits=sample(1:20, 20, replace=T))
# baseball <- data.frame(Batter=sample(1:100, 10000, replace=T), Bat=sample(c("David", "Batley", "Bob", "Ace"), 10000, replace=T), Hits=sample(1:100, 10000, replace=T), Team=sample(c("Braves", "Red Sox", "Yankees", "Phillies"), size=10000, replace=T))
# bMod <- lm(Hits ~ Bat*Batter, baseball)
# bLinMod <- rxLinMod(Hits ~ Bat*Batter + Bat*Team, baseball)
# matchCoefs(bMod)

#' @title get.assign
#' @description The assignment vector for a model
#' @details Gets relative positions of predictors
#' @aliases get.assign
#' @author Jared P. Lander
#' @param model Fitted model
#' @param \dots Further arguments
#' @return The assignment vector
get.assign <- function(model, ...)
{
    UseMethod("get.assign")
}

#' @title get.assign.lm
#' @description The assignment vector for an lm model
#' @details Gets relative positions of predictors
#' @aliases get.assign.lm
#' @author Jared P. Lander
#' @param model Fitted model
#' @param \dots Further arguments
#' @return The assignment vector
get.assign.lm <- function(model, ...)
{
    model$assign
}

#' @title get.assign.glm
#' @description The assignment vector for a glm model
#' @details Gets relative positions of predictors
#' @aliases get.assign.glm
#' @author Jared P. Lander
#' @param model Fitted model
#' @param \dots Further arguments
#' @return The assignment vector
get.assign.glm <- function(model, ...)
{
    # build model.matrix
    theMat <- stats::model.matrix(object=model$formula, data=model$data)
    # get assignment
    attr(theMat, "assign")
}

#' @title matchCoefs
#' @description Match coefficients to predictors
#' @details Matches coefficients to predictors using information from model matrices
#' @author Jared P. Lander
#' @aliases matchCoefs
#' @param model Fitted model
#' @param \dots Further arguments
#' @return a data.frame matching predictors to coefficients
#' @examples
#' \dontrun{
#' require(reshape2)
#' require(plyr)
#' data("tips", package="reshape2")
#' mod1 <- lm(tip ~ total_bill * sex + day, tips)
#' mod2 <- lm(tip ~ total_bill * sex + day - 1, tips)
#' mod3 <- glm(tip ~ total_bill * sex + day, tips, family=gaussian(link="identity"))
#' mod4 <- lm(tip ~ (total_bill + sex + day)^3, tips)
#' mod5 <- lm(tip ~ total_bill * sex + day + I(total_bill^2), tips)
#' coefplot:::matchCoefs(mod1)
#' coefplot:::matchCoefs(mod2)
#' coefplot:::matchCoefs(mod3)
#' coefplot:::matchCoefs(mod4)
#' coefplot:::matchCoefs(mod5)
#' }
# 
matchCoefs <- function(model, ...)
{
    UseMethod(generic="matchCoefs")
}

#' @title matchCoefs.default
#' @description Match coefficients to predictors
#' @details Matches coefficients to predictors using information from model matrices
#' @author Jared P. Lander
#' @aliases matchCoefs.default
#' @import reshape2
#' @param model Fitted model
#' @param \dots Further arguments
#' @return a data.frame matching predictors to coefficients
matchCoefs.default <- function(model, ...)
{
    # get the terms
    theTerms <- model$terms
    # get the assignment position
    #thePos <- model$assign
    thePos <- get.assign(model)
    # get intercept indicator
    inter <- attr(theTerms, "intercept")
    # get coef names
    coefNames <- names(stats::coef(model))
    # get pred names
    predNames <- attr(theTerms, "term.labels")
    # expand out pred names to match coefficient names
    predNames <- predNames[thePos]
    # if there's an intercept term add it to the pred names
    if(inter == 1)
    {
        predNames <- c("(Intercept)", predNames)
    }
    
    # build data.frame linking term to coefficient name
    matching <- data.frame(Term=predNames, Coefficient=coefNames, stringsAsFactors=FALSE)
    
    ## now match individual predictor to term
    # get matrix as data.frame
    factorMat <- as.data.frame(attr(theTerms, "factor"))
    # add column from rownames as identifier
    factorMat$.Pred <- rownames(factorMat)
    factorMat$.Type <- attr(theTerms, "dataClasses")
    
    # melt it down for comparison
    factorMelt <- melt(factorMat, id.vars=c(".Pred", ".Type"), variable.name="Term")
    factorMelt$Term <- as.character(factorMelt$Term)
    
    # only keep rows where there's a match
    factorMelt <- factorMelt[factorMelt$value == 1, ]
    
    # again, bring in coefficient if needed
    if(inter == 1)
    {
        factorMelt <- rbind(data.frame(.Pred="(Intercept)", .Type="(Intercept)", Term="(Intercept)", value=1, stringsAsFactors=FALSE), factorMelt)
    }
    
    # join into the matching data.frame
    matching <- join(matching, factorMelt, by="Term")
    
    return(matching)
}

# theTerms <- mod1$terms
# thePos <- mod1$assign
# inter <- attr(theTerms, "intercept")
# coefNames <- names(coef(mod1))
# coefNames
# predNames <- attr(theTerms, "term.labels")
# predNames <- predNames[thePos]
# factorMat <- as.data.frame(attr(theTerms, "factor"))
# factorMat$.Pred <- rownames(factorMat)
# factorMelt <- melt(factorMat, id.vars=".Pred")
# factorMelt <- factorMelt[factorMelt$value == 1, ]
# 
# 

#' @title getCoefsFromPredictors
#' 
#' @description Generic function for finding which coefficients go with which predictors
#' 
#' @details The user specifies predictors whose coefficients should be included in the coefplot.
#' @aliases getCoefsFromPredictors
#' @author Jared P. Lander
#' @return A character vector of coefficients listing the coefficients that match the predictor
#' @param model A fitted model
#' @param predictors A character vector of predictors to match against
#' @param \dots further arguments
getCoefsFromPredictors <- function(model, predictors, ...)
{
    UseMethod(generic="getCoefsFromPredictors", object=model)
}


#' @title getCoefsFromPredictors.default
#' 
#' @description Default function (lm, glm) for matching coefficients with predictors
#' 
#' @details The user specifies predictors whose coefficients should be included in the coefplot.
#' @aliases getCoefsFromPredictors.default
#' @author Jared P. Lander
#' @return A character vector of coefficients listing the coefficients that match the predictor
#' @param model A fitted model
#' @param predictors A character vector of predictors to match against.  Interactions can be explicitely specified by VariableA:VariableB.
#' @param strict Logical specifying if interactions terms should be included (\code{FALSE}) or just the main terms (\code{TRUE}).
#' @param \dots further arguments
getCoefsFromPredictors.default <- function(model, predictors=NULL, strict=FALSE, ...)
{
    # if no predictors indicated just return NULL
    if(is.null(predictors))
    {
        return(NULL)
    }

    # build data.frame matching predictors with coefficients
    matchPredsCoefs <- matchCoefs(model)
    
    # find out which coefficients we'll be keeping
    # if strict, it will only match the singleton term
    # if not strict it will match interactions
    if(!strict)
    {
        toKeepNotStrict <- which(matchPredsCoefs$.Pred %in% predictors)
    }else
    {
        toKeepNotStrict <- NULL
    }
    toKeepStrict <- which(matchPredsCoefs$Term %in% predictors)
    
    keptCoefsFromPredictors <- unique(matchPredsCoefs$Coefficient[unique(c(toKeepNotStrict, toKeepStrict))])

    return(keptCoefsFromPredictors)
}


#' @title getCoefsFromPredictorsRevo
#' 
#' @description Function thhat does the work for Revo models for matching coefficients with predictors
#' 
#' @details The user specifies predictors whose coefficients should be included in the coefplot.
#' @aliases getCoefsFromPredictorsRevo
#' @author Jared P. Lander
#' @return A character vector of coefficients listing the coefficients that match the predictor.  As of now interactions cannot be explicitely specified.
#' @param model A fitted model
#' @param predictors A character vector of predictors to match against
#' @param strict Logical specifying if interactions terms should be included (\code{FALSE}) or just the main terms (\code{TRUE}).
#' @param \dots further arguments
getCoefsFromPredictorsRevo <- function(model, predictors=NULL, strict=FALSE, ...)
{
    # if no predictors indicated just return NULL
    if(is.null(predictors))
    {
        return(NULL)
    }
    
    predMatcher <- subSpecials(predictors)[[1]]
    names(predMatcher) <- predictors
    
    # get coefficient names
    coefNames <- names(stats::coef(model))
    
    # if strict do one search
    if(strict)
    {
        # strict
        toKeep <- lapply(predMatcher, FUN=doRegex, matchAgainst=coefNames, pattern="^%s(=+[^,]*)*$")
    }else
    {
        # not strict
        toKeep <- lapply(predMatcher, FUN=doRegex, matchAgainst=coefNames, pattern="(^| )%s($|,|=| for)")
    }
    
    keepFrame <- ldply(toKeep, function(x) data.frame(Target=x))
    
    keptCoefsFromPredictors <- coefNames[keepFrame$Target]
    
    return(keptCoefsFromPredictors)
}


#' @title doRegex
#' 
#' @description Helper function for matching coefficients
#' 
#' @details Only used by \code{\link{getCoefsFromPredictorsRevo}} for finding matches between predictors and coefficients
#' 
#' @aliases doRegex
#' @author Jared P. Lander
#' @param x Root pattern to search for
#' @param matchAgainst Text to search through
#' @param pattern Regex pattern to build x into
#' @return A list of indices of matchAgainst that is matched
doRegex <- function(x, matchAgainst, pattern="(^| )%s($|,|=)")
{
    grep(pattern=sprintf(pattern, x), x=matchAgainst)
}


#' @title getCoefsFromPredictors.rxLinMod
#' 
#' @description Function for matching coefficients with predictors for rxLinMod
#' 
#' @details The user specifies predictors whose coefficients should be included in the coefplot.
#' @aliases getCoefsFromPredictors.rxLinMod
#' @author Jared P. Lander
#' @return A character vector of coefficients listing the coefficients that match the predictor
#' @param model A fitted model
#' @param predictors A character vector of predictors to match against
#' @param strict Logical specifying if interactions terms should be included (\code{FALSE}) or just the main terms (\code{TRUE}).
#' @param \dots further arguments
getCoefsFromPredictors.rxLinMod <- function(model, predictors=NULL, strict=FALSE, ...)
{
    getCoefsFromPredictorsRevo(model, predictors, strict)
}


#' @title getCoefsFromPredictors.rxLogit
#' 
#' @description Function for matching coefficients with predictors for rxLogit
#' 
#' @details The user specifies predictors whose coefficients should be included in the coefplot.
#' @aliases getCoefsFromPredictors.rxLogit
#' @author Jared P. Lander
#' @return A character vector of coefficients listing the coefficients that match the predictor
#' @param model A fitted model
#' @param predictors A character vector of predictors to match against
#' @param strict Logical specifying if interactions terms should be included (\code{FALSE}) or just the main terms (\code{TRUE}).
#' @param \dots further arguments
getCoefsFromPredictors.rxLogit <- function(model, predictors=NULL, strict=FALSE, ...)
{
    getCoefsFromPredictorsRevo(model, predictors, strict)
}


#' @title getCoefsFromPredictors.rxGlm
#' 
#' @description Function for matching coefficients with predictors for rxGlm
#' 
#' @details The user specifies predictors whose coefficients should be included in the coefplot.
#' @aliases getCoefsFromPredictors.rxGlm
#' @author Jared P. Lander
#' @return A character vector of coefficients listing the coefficients that match the predictor
#' @param model A fitted model
#' @param predictors A character vector of predictors to match against
#' @param strict Logical specifying if interactions terms should be included (\code{FALSE}) or just the main terms (\code{TRUE}).
#' @param \dots further arguments
getCoefsFromPredictors.rxGlm <- function(model, predictors=NULL, strict=FALSE, ...)
{
    getCoefsFromPredictorsRevo(model, predictors, strict)
}
