
library("multcomp")

###<FIXME> compare results of mmod and glht.mlf </FIXME>

### code by Christian Ritz
"mmod" <- function(modelList, varName, seType = "san")
{
    require(multcomp, quietly = TRUE)
    require(sandwich, quietly = TRUE)

    if (length(seType) == 1) {seType <- rep(seType, length(modelList))}
    if (length(varName) == 1) {varName <- rep(varName, length(modelList))}    

    ## Extracting score contributions from the individual model fits    
    makeIIDdecomp <- function(modelObject, varName) 
    {
        numObsUsed <- ifelse(inherits(modelObject, "coxph"), modelObject$n, nrow(modelObject$model))
        iidVec0 <- bread(modelObject)[varName, , drop = FALSE] %*% t(estfun(modelObject))
        moNAac <- modelObject$na.action
        numObs <- numObsUsed + length(moNAac)
        iidVec <- rep(0, numObs)
        if (!is.null(moNAac))
        {
            iidVec[-moNAac] <- sqrt(numObs/numObsUsed) * iidVec0
        }
        else {
            iidVec <- iidVec0
        }
        list(iidVec = iidVec, numObsUsed = numObsUsed, numObs = numObs)
    }
    numModels <- length(modelList)
    if (identical(length(varName), 1))
    {
        varName <- rep(varName, numModels)
    }
    iidList <- mapply(makeIIDdecomp, modelList, varName, SIMPLIFY = FALSE)
    iidresp <- matrix(as.vector(unlist(lapply(iidList, function(listElt) {listElt[[1]]}))), nrow = numModels, byrow = TRUE)
    pickFct <- function(modelObject, varName, matchStrings) 
    {
        as.vector(na.omit((coef(summary(modelObject))[varName, ])[matchStrings]))
    }
          
    ## Retrieving parameter estimates from the individual fits
    estVec <- as.vector(unlist(mapply(pickFct, modelList, varName, MoreArgs = list(matchStrings = c("Estimate", "coef")))))
    # "Estimate" or "coef" used in glm(), lm() and coxph() summary output, respectively
    
    ## Calculating the estimated variance-covariance matrix of the parameter estimates
    numObs <- iidList[[1]]$numObs
    covar <- (iidresp %*% t(iidresp)) / numObs
    vcMat <- covar / numObs  # Defining the finite-sample variance-covariance matrix

    ## Replacing sandwich estimates by model-based standard errors
    modbas <- seType == "mod"
    if (any(modbas))
    {  
        corMat <- cov2cor(vcMat)
        ## Retrieving standard errors for the specified estimate from the individual fits
        modSE <- as.vector(unlist(mapply(pickFct, modelList, varName, MoreArgs = list(matchStrings = c("Std. Error", "se(coef)")))))
        
        sanSE <- sqrt(diag(vcMat))
        sanSE[modbas] <- modSE[modbas]
        vcMat <- diag(sanSE) %*% corMat %*% diag(sanSE)
    } 

    ## Naming the parameter vector (easier way to extract the names of the model fits provided as a list in the first argument?)
    names1 <- sub("list", "", deparse(substitute(modelList)), fixed = TRUE)
    names2 <- sub("(", "", names1, fixed = TRUE)
    names3 <- sub(")", "", names2, fixed = TRUE)
    names4 <- sub(" ", "", names3, fixed = TRUE)
    names(estVec) <- unlist(strsplit(names4, ","))
             
    return(parm(coef = estVec, vcov = vcMat, df = 0))
}
    




set.seed(29)
## Combining linear regression and logistic regression
y1 <- rnorm(100)
y2 <- factor(y1 + rnorm(100, sd = .1) > 0)
x1 <- gl(4, 25) 
x2 <- runif(100, 0, 10)

m1 <- lm(y1 ~ x1 + x2)
m2 <- glm(y2 ~ x1 + x2, family = binomial())
## Note that the same explanatory variables are considered in both models
##  but the resulting parameter estimates are on 2 different scales (original and log-odds scales)

## Simultaneous inference for the same parameter in the 2 model fits
simult.x12 <- mmod(list(m1, m2), c("x12", "x12"))
summary(glht(simult.x12))

## Simultaneous inference for different parameters in the 2 model fits
simult.x12.x13 <- mmod(list(m1, m2), c("x12", "x13"))
summary(glht(simult.x12.x13))

## Simultaneous inference for different and identical parameters in the 2 model fits
simult.x12x2.x13 <- mmod(list(m1, m1, m2), c("x12", "x13", "x13"))
summary(glht(simult.x12x2.x13))
confint(glht(simult.x12x2.x13))


## Examples for binomial data

## Two independent outcomes
y1.1 <- rbinom(100, 1, 0.5)
y1.2 <- rbinom(100, 1, 0.5)
group <- factor(rep(c("A", "B"), 50))

modely1.1 <- glm(y1.1 ~ group, family = binomial)
modely1.2 <- glm(y1.2 ~ group, family = binomial)

mmObj.y1 <- mmod(list(modely1.1, modely1.2), "groupB")
simult.y1 <- glht(mmObj.y1)
summary(simult.y1)

## Two perfectly correlated outcomes
y2.1 <- rbinom(100, 1, 0.5)
y2.2 <- y2.1
group <- factor(rep(c("A", "B"), 50))

modely2.1 <- glm(y2.1 ~ group, family = binomial)
modely2.2 <- glm(y2.2 ~ group, family = binomial)

mmObj.y2 <- mmod(list(modely2.1, modely2.2), "groupB")
simult.y2 <- glht(mmObj.y2)
summary(simult.y2)

