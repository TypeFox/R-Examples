"anova.drc" <-
function(object, ..., details = TRUE, test = NULL)
{
    if (length(list(object, ...)) > 1) 
    {
        return(anova.drclist(object, ..., details = details, test = test))
    } else {
#        obj1 <- object
        stop("Use the function modelFit()")
#        modelFit(object, test = test)
    }

#    if (inherits(object, "bindrc"))
#    {
#        loglik <- object$loglik[c(1, 3)]
#        dfModel <- object$loglik[c(2, 4)]
#        dfDiff <- c(NA, (dfModel[2] - dfModel[1]))
#
#        testStat <- 2*(loglik[1] - loglik[2])
#        pVal <- c(NA, 1 - pchisq(testStat, dfDiff[2]))
#        testStat <- c(NA, testStat)
#
#        headName<-"Analysis of deviance table\n"
#        rowNames<-c("Saturated model", "DRC model")
#        colNames<-c("ModelDf", "Loglik", "Df", "LR value", "p value")
#
#        dataFra <- data.frame(dfModel, loglik, dfDiff, testStat, pVal)
#
#        dimnames(dataFra) <- list(rowNames, colNames)
#        structure(dataFra, heading=headName, class=c("anova", "data.frame"))
#        
#    } else {

        ## Testing against ANOVA
#        listElt <- obj1[[4]]
#
#        wayStr <- obj1$"sumList"$"alternative"  
#        if (wayStr == 1) {wayStr <- "One-way"}
#        if (wayStr == 2) {wayStr <- "Two-way"} 


#        anovaSS <- listElt[3]
#        anovaDF <- listElt[4]
#        nlsSS <- listElt[5]
#        nlsDF <- listElt[6]
        
        
#        anovaModel <- obj1$"anova"
#        if ( (!is.null(anovaModel)) && (is.null(test)) ) 
#        {
#            anovaDF <- df.residual(anovaModel$"anovaFit")
#            nlsDF <- df.residual(obj1)
#            dfModel <- c(anovaDF, nlsDF)
#            dfDiff <- c(NA, (nlsDF - anovaDF))
#                                
#            if (anovaModel$"test"=="F")
#            {
#                anovaSS <- deviance(anovaModel$"anovaFit")
#                anovaDF <- df.residual(anovaModel$"anovaFit")
#                dfModel <- c(anovaDF, nlsDF)
#                nlsSS <- obj1$"fit"$"value"
##                nlsDF <- df.residual(obj1)
#
#                loglik <- c(anovaSS, nlsSS)
##                dfModel <- c(anovaDF, nlsDF)
##                dfDiff <- c(NA, (nlsDF-anovaDF))
#
#                testStat <- (nlsSS-anovaSS)/dfDiff[2]/(anovaSS/anovaDF)
#                pVal <- c(NA, pf(testStat,dfDiff[2], anovaDF, lower.tail=FALSE))
#                testStat <- c(NA, testStat)
#
#                headName<-"Lack-of-fit test\n"
#                rowNames<-c(paste(wayStr, "ANOVA"), "DRC model")
#                colNames<-c("ModelDf", "RSS", "Df", "F value", "p value")    
#            }
#            
#            if (anovaModel$"test"=="lr")
#            {
#                anovaDF <- anovaDF + (obj1$"sumList"$"lenData" - dim(anovaModel$"anovaFit"$"data")[1])
#                dfModel <- c(anovaDF, nlsDF)
#                dfDiff <- c(NA, (nlsDF - anovaDF))
#                
#                loglik <- c(logLik(anovaModel$"anovaFit"), logLik(obj1))
#
#                testStat <- 2*(loglik[1]-loglik[2])
#                pVal <- c(NA, 1 - pchisq(testStat, dfDiff[2]))
#                testStat <- c(NA, testStat)
#
#                headName <- "Goodness-of-fit test\n"
#                rowNames <- c(paste(wayStr, "ANOVA"), "DRC model")
#                colNames <- c("ModelDf", "Log lik", "Df", "Chisq value", "p value")                
#            }           
#        } 
#        if ( (is.null(anovaModel)) || ( (!is.null(test)) && (test == "od") ) )
#        {        
##        else {
##            gof <- (obj1$"fct"$"gofTest")(obj1)
#            gofTest <- obj1$"gofTest"
#            if (!is.null(gofTest)) 
#            {
##                gof <- gofTest(obj1)
#            
#                lenData <- obj1$"sumList"$"lenData"
#                dfModel <- c(NA, NA)
#                loglik <- c(NA, NA)
#                dfDiff <- c(NA, gofTest[2])
#                testStat <- c(NA, gofTest[1])
#                pVal <- c(NA, 1 - pchisq(testStat[2], dfDiff[2]))
#
#                headName <- "Goodness-of-fit test\n"
#                rowNames <- c("", "DRC model")
#                colNames <- c("", "", "Df", "Chisq value", "p value")                            
#            } else {  # in case no test is available
#                dfModel <- c(NA, NA)
#                loglik <- c(NA, NA)
#                dfDiff <- c(NA, NA)
#                testStat <- c(NA, NA)
#                pVal <- c(NA, NA)
#                headName <- "No test available\n"
#                rowNames <- c("", "DRC model")
#                colNames <- c("ModelDf", "Log lik", "Df", "Chisq value", "p value")   
#            }
#            
##            headName<-"Goodness-of-fit test\n"
##            rowNames<-c("", "DRC model")
##            colNames<-c("ModelDf", "", "Df", "Chisq value", "p value")                            
#        }
#        if ( (!is.null(anovaModel)) && (!is.null(test)) && (test == "Chisq") )  # overruling any previous calculations!
#        {
#            lv1 <- logLik(anovaModel$"anovaFit")
#            lv2 <- logLik(obj1)            
#        
#            dfModel <- c(attr(lv1, "df"), attr(lv2, "df"))
#            loglik <- c(lv1, lv2)
#            dfDiff <- c(NA, diff(dfModel))  # dfModel[1] - dfModel[2])
#        
#            testStat <- -2*(lv2 - lv1)
#            pVal <- c(NA, pchisq(testStat, dfDiff[2], lower.tail = FALSE))
#            testStat <- c(NA, testStat)
#
#            headName <- "Lack-of-fit test\n"
#            rowNames <- c(paste(wayStr, "ANOVA"), "DRC model")
#            colNames <- c("ModelDf", "Log lik", "Df", "Chisq value", "p value")   
#        }


#        loglik <- c(anovaSS, nlsSS)
#        dfModel <- c(anovaDF, nlsDF)
#        dfDiff <- c(NA, (nlsDF-anovaDF))

#        testStat <- (nlsSS-anovaSS)/dfDiff[2]/(anovaSS/anovaDF)
#        pVal <- c(NA, pf(testStat,dfDiff[2], anovaDF, lower.tail=FALSE))
#        testStat <- c(NA, testStat)

#        headName<-"ANOVA table\n"
#        rowNames<-c(paste(wayStr, "ANOVA"), "DRC model")
#        colNames<-c("ModelDf", "RSS", "Df", "F value", "p value")

#        dataFra <- data.frame(dfModel, loglik, dfDiff, testStat, pVal)
#
#        dimnames(dataFra) <- list(rowNames, colNames)
#        structure(dataFra, heading = headName, class = c("anova", "data.frame"))
#    }
 
#    ## Testing two models against each other
#    if (!anovaTest)
#    {
#        rowNames <- c(deparse(substitute(obj1)), deparse(substitute(obj2)))
#
#        sumObj1 <- summary(obj1)
#        sumObj2 <- summary(obj2)
#
#        if (!Ftest)
#        {
#            loglik <- c(sumObj1[[4]][1],sumObj2[[4]][1])
#            dfModel <- c(sumObj1[[4]][2],sumObj2[[4]][2])
#            testStat <- (2*abs(loglik[1]-loglik[2]))
#            dfDiff <- c(NA,abs(dfModel[1]-dfModel[2]))
#
#            pVal <- c(NA,1-pchisq(testStat,dfDiff[2]))
#            testStat <- c(NA,testStat)
#
#            headName <- "ANOVA-like table\n"
#            colNames <- c("ModelDf", "Loglik", "Df", "LR value", "p value")
#        } else {
#            
#            sumVec1 <- obj1[[4]]
#            sumVec2 <- obj2[[4]]
#
#            if (sumVec2[6]>sumVec1[6]) {sumTemp <- sumVec1; sumVec1 <- sumVec2; sumVec2 <- sumTemp; rowNames <- rowNames[c(2,1)]}
#
#            loglik <- c(sumVec1[5],sumVec2[5])
#            dfModel <- c(sumVec1[6],sumVec2[6])
#            dfDiff <- c((loglik[1]-loglik[2])/(dfModel[1]-dfModel[2]), loglik[2]/dfModel[2])
#            testStat <- dfDiff[1]/dfDiff[2]
#
#            pVal <- c(NA,1-pf(testStat, dfModel[1]-dfModel[2], dfModel[2]))
#            testStat <- c(NA,testStat)
#
#            headName <- "ANOVA table\n"
#            colNames <- c("Df", "Sum Sq", "Mean Sq", "F value", "p value")
#        }
#    }

}
