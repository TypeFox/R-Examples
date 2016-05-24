# this file extracts coefficients from varying models
#' extract.coef.default
#' 
#' Extract Coefficient Information from Models
#' 
#' Gets the coefficient values and standard errors, and variable names from a model.
#' 
#' @author Jared P. Lander
#' @method extract.coef default
#' @aliases extract.coef.default
#' @param model Model object to extract information from.
#' @param \dots Further arguments
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' require(ggplot2)
#' library(coefplot)
#' data(diamonds)
#' mod1 <- lm(price ~ carat + cut + x, data=diamonds)
#' extract.coef(mod1)
#' }
#' 
extract.coef.default <- function(model, ...)
{
    # get summary of model
    theSumm <- summary(model)
    # get coef and standard error
#     print(theSumm)
#     print(head(model))
    info <- as.data.frame(theSumm$coefficients[, 1:2])
    names(info) <- c("Value", "SE")
    # make a variable tracking the name
    info$Coefficient <- rownames(info)
    
    return(info)
}

#' extract.coef.lm
#' 
#' Extract Coefficient Information from lm Models
#' 
#' Gets the coefficient values and standard errors, and variable names from an lm model.
#' 
#' @author Jared P. Lander
#' @aliases extract.coef.lm
#' @method extract.coef lm
#' @inheritParams extract.coef.default
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data(diamonds)
#' library(coefplot)
#' mod1 <- lm(price ~ carat + cut + x, data=diamonds)
#' extract.coef(mod1)
#' }
#' 
extract.coef.lm <- function(model, ...)
{
    extract.coef.default(model=model, ...)
}


#' extract.coef.glm
#' 
#' Extract Coefficient Information from glm Models
#' 
#' Gets the coefficient values and standard errors, and variable names from a glm model.
#' 
#' @author Jared P. Lander
#' @aliases extract.coef.glm
#' @inheritParams extract.coef.default
#' @method extract.coef glm
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data(diamonds)
#' library(coefplot)
#' mod2 <- glm(price > 10000 ~ carat + cut + x, data=diamonds, family=binomial(link="logit"))
#' extract.coef(mod2)
#' }
#' 
extract.coef.glm <- function(model, ...)
{
    extract.coef.default(model=model, ...)
}


#' extract.coef
#' 
#' Extract Coefficient Information from glm Models
#' 
#' Gets the coefficient values and standard errors, and variable names from a glm model.
#' 
#' @author Jared P. Lander
#' @aliases extract.coef
#' @export extract.coef
#' @param model Model object to extract information from.
#' @param \dots Further arguments
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data(diamonds)
#' library(coefplot)
#' mod1 <- lm(price ~ carat + cut + x, data=diamonds)
#' mod2 <- glm(price > 10000 ~ carat + cut + x, data=diamonds, family=binomial(link="logit"))
#' mod3 <- lm(price ~ carat*cut + x, data=diamonds)
#' extract.coef(mod1)
#' extract.coef(mod2)
#' extract.coef(mod3)
#' 
#' mod4 <- rxLinMod(price ~ carat*cut + x, diamonds)
#' }
#' 
extract.coef <- function(model, ...)
{
    UseMethod(generic="extract.coef", object=model)
}


#' extract.coef.rxLinMod
#' 
#' Extract Coefficient Information from rxLinMod Models
#' 
#' Gets the coefficient values and standard errors, and variable names from an rxLinMod model.
#' 
#' @author Jared P. Lander
#' @aliases extract.coef.rxLinMod
#' @inheritParams extract.coef.default
#' @method extract.coef rxLinMod
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data(diamonds)
#' mod3 <- rxLinMod(price ~ carat + cut + x, data=diamonds)
#' extract.coef(mod3)
#' }
#' 
extract.coef.rxLinMod <- function(model, ...)
{
    # get summary
    theSumm <- summary(model)[[1]]
    # get coefs
    info <- as.data.frame(theSumm$coefficients[, 1:2])
    # give good names
    names(info) <- c("Value", "SE")
    # get variable names
    info$Coefficient <- rownames(info)
    
    return(info)
}


#' extract.coef.rxGlm
#' 
#' Extract Coefficient Information from rxGlm Models
#' 
#' Gets the coefficient values and standard errors, and variable names from an rxGlm model.
#' 
#' @author Jared P. Lander
#' @aliases extract.coef.rxGlm
#' @inheritParams extract.coef.default
#' @method extract.coef rxGlm
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data(diamonds)
#' mod4 <- rxGlm(price ~ carat + cut + x, data=diamonds)
#' mod5 <- rxGlm(price > 10000 ~ carat + cut + x, data=diamonds, fmaily="binomial")
#' extract.coef(mod4)
#' extract.coef(mod5)
#' }
#' 
extract.coef.rxGlm <- function(model, ...)
{
    extract.coef.default(model=model, ...)
}


#' extract.coef.rxLogit
#' 
#' Extract Coefficient Information from rxLogit Models
#' 
#' Gets the coefficient values and standard errors, and variable names from an rxLogit model.
#' 
#' @author Jared P. Lander
#' @aliases extract.coef.rxLogit
#' @inheritParams extract.coef.default
#' @method extract.coef rxLogit
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data(diamonds)
#' mod6 <- rxLogit(price > 10000 ~ carat + cut + x, data=diamonds)
#' extract.coef(mod6)
#' }
#' 
extract.coef.rxLogit <- function(model, ...)
{
    extract.coef.default(model=model, ...)
}

#' @title extract.coef.glmnet
#' @description Extract Coefficient Information from Models
#' @details Gets the coefficient values and variable names from a model.  Since glmnet does not have standard errors, those will just be NA.
#' @author Jared P. Lander
#' @method extract.coef glmnet
#' @aliases extract.coef.glmnet
#' @param model Model object from which to extract information.
#' @param lambda Value of penalty parameter
#' @param \dots Further arguments
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' library(glmnet)
#' library(ggplot2)
#' library(useful)
#' data(diamonds)
#' diaX <- build.x(price ~ carat + cut + x - 1, data=diamonds, contrasts = TRUE)
#' diaY <- build.y(price ~ carat + cut + x - 1, data=diamonds)
#' modG1 <- glmnet(x=diaX, y=diaY)
#' extract.coef(modG1)
#' }
#' 
extract.coef.glmnet <- function(model, lambda=stats::median(model$lambda), ...)
{
    # get coefs at given s
    theCoef <- as.matrix(stats::coef(model, s=lambda))
    coefDF <- data.frame(Value=theCoef, SE=NA_real_, Coefficient=rownames(theCoef))
    coefDF <- coefDF[theCoef != 0, ]
    names(coefDF)[1] <- "Value"

    return(coefDF)
}

#' @title extract.coef.cv.glmnet
#' @description Extract Coefficient Information from Models
#' @details Gets the coefficient values and variable names from a model.  Since glmnet does not have standard errors, those will just be NA.
#' @author Jared P. Lander
#' @method extract.coef cv.glmnet
#' @aliases extract.coef.cv.glmnet
#' @param model Model object from which to extract information.
#' @param lambda Value of penalty parameter.  Can be either a numeric value or one of "lambda.min" or "lambda.1se"
#' @param \dots Further arguments
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' library(glmnet)
#' library(ggplot2)
#' library(useful)
#' data(diamonds)
#' diaX <- build.x(price ~ carat + cut + x - 1, data=diamonds, contrasts = TRUE)
#' diaY <- build.y(price ~ carat + cut + x - 1, data=diamonds)
#' modG1 <- cv.glmnet(x=diaX, y=diaY, k=5)
#' extract.coef(modG1)
#' }
#' 
extract.coef.cv.glmnet <- function(model, lambda="lambda.min", ...)
{
    extract.coef.glmnet(model, lambda=lambda, ...)
}


#' @title extract.coef.maxLik
#' @description Extract Coefficient Information from Models
#' @details Gets the coefficient values and variable names from a model.
#' @author Jared P. Lander
#' @method extract.coef maxLik
#' @export extract.coef.maxLik
#' @export
#' @aliases extract.coef.maxLik
#' @param model Model object from which to extract information.
#' @param \dots Further arguments
#' @return A \code{\link{data.frame}} containing the coefficient, the standard error and the variable name.
#' @examples
#' \dontrun{
#' library(maxLik)
#' loglik <- function(param) {
#'  mu <- param[1]
#'  sigma <- param[2]
#'  ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
#' ll
#' }
#' x <- rnorm(1000, 1, 2) # use mean=1, stdd=2
#' N <- length(x)
#' res <- maxLik(loglik, start=c(0,1)) # use 'wrong' start values
#' extract.coef(res)
#' }
#' 
extract.coef.maxLik <- function(model, ...)
{
    # get coefficients
    theCoef <- stats::coef(model)
    # get coef names
    coefNames <- names(theCoef)
    
    ## if names do not exist make the names the number of the coef
    if(is.null(coefNames))
    {
        coefNames <- seq(along=theCoef)
    }
    
    data.frame(Value=theCoef, SE=summary(model)$estimate[, 'Std. error'], Coefficient=coefNames)
}