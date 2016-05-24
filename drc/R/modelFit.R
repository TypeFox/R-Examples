"modelFit" <- function(object, test = NULL, method = c("gof", "cum"))
{
    method <- match.arg(method)

    ## Fitting ANOVA model
    testList <- switch(object$"type", 
    "continuous" = drmLOFls(), 
    "binomial" = drmLOFbinomial(), 
    "Poisson" = drmLOFPoisson())
  
    switch(object$"type", 
    binomial = gofTest(object, testList$"gofTest"),  
    continuous = lofTest(object, testList$"anovaTest"))
}


"lofTest" <- function(object, anovaTest)
{
    if (!is.null(anovaTest))
    {
        dose <- object$"dataList"$"dose"
        resp <- object$"dataList"$"resp"
        curveid <- object$"dataList"$"curveid"
        if (is.null(object$"boxcox"))
        {
            bcAdd <- 0
        } else {
            bcAdd <- object$"boxcox"$"bcAdd"
        }
        afList <- anovaFormula(dose, resp, curveid, bcAdd)
        anovaForm <- afList$"anovaFormula"
        anovaData <- afList$"anovaData"        
        anovaModel <- anovaTest(anovaForm, anovaData)
        if (is.null(anovaModel))
        {
            return(returnFct())
        }

        anovaDF <- df.residual(anovaModel$"anovaFit")
        nlsDF <- df.residual(object)
        dfModel <- c(anovaDF, nlsDF)
        dfDiff <- c(NA, (nlsDF - anovaDF))
                                
        if (identical(anovaModel$"test", "F"))
        {
            anovaSS <- deviance(anovaModel$"anovaFit")
            anovaDF <- df.residual(anovaModel$"anovaFit")
            dfModel <- c(anovaDF, nlsDF)
            nlsSS <- object$"fit"$"value"
            loglik <- c(anovaSS, nlsSS)

            testStat <- (nlsSS - anovaSS)/dfDiff[2]/(anovaSS/anovaDF)
            pVal <- c(NA, pf(testStat, dfDiff[2], anovaDF, lower.tail = FALSE))
            testStat <- c(NA, testStat)

            headName<-"Lack-of-fit test\n"
            rowNames<-c("ANOVA", "DRC model")
            colNames<-c("ModelDf", "RSS", "Df", "F value", "p value")    
        }
            
        if (identical(anovaModel$"test", "lr"))
        {
            anovaDF <- anovaDF + (object$"sumList"$"lenData" - dim(anovaModel$"anovaFit"$"data")[1])
            dfModel <- c(anovaDF, nlsDF)
            dfDiff <- c(NA, (nlsDF - anovaDF))    
            loglik <- c(logLik(anovaModel$"anovaFit"), logLik(object))

            testStat <- 2*(loglik[1] - loglik[2])
            pVal <- c(NA, 1 - pchisq(testStat, dfDiff[2]))
            testStat <- c(NA, testStat)

            headName <- "Goodness-of-fit test\n"
            rowNames <- c("ANOVA", "DRC model")
            colNames <- c("ModelDf", "Log lik", "Df", "Chisq value", "p value")                
        }
        return(returnFct(dfModel, loglik, dfDiff, testStat, pVal, headName, colNames, rowNames))
         
    } else {
        return(returnFct())
    }
}

"gofTest" <- function(object, gofTest)
{       
    gofTest <- gofTest(object$"dataList"$"resp", weights(object), fitted(object), df.residual(object))
    
    if (!is.null(gofTest))
    {
        returnFct(c(NA, NA), c(NA, NA), c(NA, gofTest[2]), c(NA, gofTest[1]), 
        c(NA, 1 - pchisq(gofTest[1], gofTest[2])), "Goodness-of-fit test\n",
        c("", "", "Df", "Chisq value", "p value"), c("", "DRC model"))
    } else {  
        returnFct()
    }
}


"returnFct" <- function(dfModel = c(NA, NA), loglik = c(NA, NA), dfDiff = c(NA, NA), testStat = c(NA, NA), 
pVal = c(NA, NA), headName = "No test available\n", colNames = c("ModelDf", "Log lik", "Df", "Chisq value", "p value"), 
rowNames = c("", "DRC model"))
{
    dataFra <- data.frame(dfModel, loglik, dfDiff, testStat, pVal)
    dimnames(dataFra) <- list(rowNames, colNames)
    structure(dataFra, heading = headName, class = c("anova", "data.frame"))  
}