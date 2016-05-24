###############################################################################
###############################################################################
## Class "SurveyData"
##
## Contains the specifications and the data for a survey
## to substantiate freedom from disease 
## using "individual sampling" or "limited sampling".
##
## Package: FFD
##
## Ian Kopacka
## 2011-06-30
###############################################################################
###############################################################################

## CLASS DEFINITION:
###############################################################################
###############################################################################
setClass(
        Class = "SurveyData", 
        representation = representation(
                nAnimalVec = "numeric",
                riskGroupVec = "character",
                riskValueData = "data.frame",
                populationData = "data.frame",
                designPrevalence = "numeric",
                alpha = "numeric",
                intraHerdPrevalence = "numeric",
                diagSensitivity = "numeric",
                costHerd = "numeric", 
                costAnimal = "numeric"        
        ),
        validity = function(object)
        { 
            if (any(object@nAnimalVec < 0))
            {
                stop ("[SurveyData validation]: Slot 'nAnimalVec' must not contain negative values.")
            }
            if (any((object@nAnimalVec  - as.integer(object@nAnimalVec)) != 0))
            {
                stop ("[SurveyData validation]: Slot 'nAnimalVec' must contain integer values.")
            }
            
            if ((length(object@nAnimalVec)*dim(object@populationData)[1] > 0)&(length(object@nAnimalVec) != dim(object@populationData)[1])){
                stop ("[SurveyData validation]: Columns of the slot 'populationData' must have the same length as slot 'nAnimalVec'.")
            }    
            if ((length(object@riskGroupVec) > 0)&(length(object@nAnimalVec) > 0)&
                    (length(object@nAnimalVec) != length(object@riskGroupVec))){
                stop ("[SurveyData validation]: Slot 'riskGroupVec' must have the same length as slot 'nAnimalVec'.")
            }  
            if (dim(object@riskValueData)[1] > 0){
                if (dim(object@riskValueData)[2] < 2){
                    stop ("[SurveyData validation]: Data frame in slot 'riskValueData' must have 2 columns: 'riskGroup', 'riskValue'.")             
                }
                if (any(!(object@riskGroupVec %in% object@riskValueData[,1]))){
                    stop ("[SurveyData validation]: All values from slot 'riskGroupVec' must be contained in the first column of slot 'riskValueData'.")            
                }
                if (!(class(object@riskValueData[,2]) %in% c("numeric", "integer"))){
                    stop ("[SurveyData validation]: The second column of slot 'riskValueData' must contain numeric values.")
                }
                if (min(object@riskValueData[,2]) <= 0){
                    stop ("[SurveyData validation]: The second column of slot 'riskValueData' must contain positive values.")
                }
            }       
            if ((length(object@designPrevalence) > 0)){ 
                if((object@designPrevalence <= 0)|(object@designPrevalence >= 1))
                {
                    stop ("[SurveyData validation]: Slot 'designPrevalence' must contain values in (0,1).")
                }
            }
            if ((length(object@alpha) > 0)){ 
                if((object@alpha <= 0)|(object@alpha >= 1))
                {
                    stop ("[SurveyData validation]: Slot 'alpha' must contain values in (0,1).")
                }
            }
            if ((length(object@intraHerdPrevalence) > 0)){ 
                if((object@intraHerdPrevalence <= 0)|(object@intraHerdPrevalence >= 1))
                {
                    stop ("[SurveyData validation]: Slot 'intraHerdPrevalence' must contain values in (0,1).")
                }
            }
            if ((length(object@diagSensitivity) > 0)){ 
                if((object@diagSensitivity <= 0)|(object@diagSensitivity >= 1))
                {
                    stop ("[SurveyData validation]: Slot 'diagSensitivity' must contain values in (0,1).")
                }   
            }     
            return(TRUE)
        }    
)

## METHODS:
###############################################################################
###############################################################################


#################################################### SHOW:
##########################################################

setMethod("show", signature(object="SurveyData"),
        function(object){
            displayLimit <- 20
            cat("Object of class 'SurveyData':\n")
            cat("Slots:\n")        
            if (length(object@designPrevalence) > 0){
                cat("@designPrevalence:      ", sprintf("%.3f", object@designPrevalence), "\n")
            } else {
                cat("@designPrevalence:       NO DATA\n")
            }
            if (length(object@alpha) > 0){
                cat("@alpha:                 ", sprintf("%.2f", object@alpha), "\n")
            } else {
                cat("@alpha:                  NO DATA\n")
            }
            if (length(object@intraHerdPrevalence) > 0){
                cat("@intraHerdPrevalence:   ", sprintf("%.2f", object@intraHerdPrevalence), "\n")
            } else {
                cat("@intraHerdPrevalence:    NO DATA\n")
            }
            if (length(object@diagSensitivity) > 0){
                cat("@diagSensitivity:       ", sprintf("%.2f", object@diagSensitivity), "\n")
            } else {
                cat("@diagSensitivity:        NO DATA\n")
            }
            if (length(object@costHerd) > 0){
                cat("@costHerd:              ", sprintf("%.2f", object@costHerd), "\n")
            } else {
                cat("@costHerd:               NO DATA\n")
            }
            if (length(object@costAnimal) > 0){
                cat("@costAnimal:            ", sprintf("%.2f", object@costAnimal), "\n")
            } else {
                cat("@costAnimal:             NO DATA\n")
            }       
            cat("@nAnimalVec:")
            if (length(object@nAnimalVec) > 0){
                cat("\n")
                if (length(object@nAnimalVec) > displayLimit){
                    show(object@nAnimalVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@nAnimalVec), "displayed)\n")
                } else {
                    show(object@nAnimalVec)
                }
            } else {        
                cat("             NO DATA\n")
            } 
            cat("@riskGroupVec:")
            if (length(object@riskGroupVec) > 0){
                cat("\n")
                if (length(object@riskGroupVec) > displayLimit){
                    show(object@riskGroupVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@riskGroupVec), "displayed)\n")
                } else {
                    show(object@riskGroupVec)
                }
            } else {        
                cat("           NO DATA\n")
            } 
            cat("@riskValueData:")
            if (dim(object@riskValueData)[1] > 0){
                cat("\n")
                if (dim(object@riskValueData)[1] > displayLimit){
                    show(head(object@riskValueData, n = displayLimit))
                    cat("   (only first",displayLimit,"rows of",dim(object@riskValueData)[1], "displayed)\n")
                } else {
                    show(object@riskValueData)
                }
            } else {        
                cat("          NO DATA\n")
            }  
            cat("@populationData:")
            if (dim(object@populationData)[1] > 0){
                cat("\n")
                if (dim(object@populationData)[1] > displayLimit){
                    show(head(object@populationData, n = displayLimit))
                    cat("   (only first",displayLimit,"rows of",dim(object@populationData)[1], "displayed)\n")
                } else {
                    show(object@populationData)
                }
            } else {        
                cat("        NO DATA\n")
            }         
        }
)

################################################# SUMMARY:
##########################################################

setMethod("summary", signature(object = "SurveyData"),
        function(object, output = c("proportion", "percent")){
            output <- match.arg(output)
            outVec <- c("Survey Parameters:", "------------------")
            if (length(object@designPrevalence) > 0){
                if (output == "proportion"){            
                    outChar <- sprintf("%.3f", object@designPrevalence)
                } else {
                    outChar <- paste(sprintf("%.2f", object@designPrevalence*100), "%")
                } 
                outVec <- c(outVec, paste("Design Prevalence:              ", 
                                outChar))
            } else {
                outVec <- c(outVec, "Design Prevalence:               NO DATA")         
            }
            if (length(object@alpha) > 0){
                outVec <- c(outVec, paste("Significance level:             ", 
                                sprintf("%.2f", object@alpha)))         
            } else {
                outVec <- c(outVec, "Significance level:              NO DATA") 
            }
            if (length(object@intraHerdPrevalence) > 0){
                if (output == "proportion"){            
                    outChar <- sprintf("%.3f", object@intraHerdPrevalence)
                } else {
                    outChar <- paste(sprintf("%.2f", object@intraHerdPrevalence*100), "%")
                }   
                outVec <- c(outVec, paste("Intra herd prevalence:          ", 
                                outChar))
            } else {
                outVec <- c(outVec, "Intra herd prevalence:           NO DATA")            
            }
            if (length(object@diagSensitivity) > 0){
                if (output == "proportion"){            
                    outChar <- sprintf("%.3f", object@diagSensitivity)
                } else {
                    outChar <- paste(sprintf("%.2f", object@diagSensitivity*100), "%")
                }  
                outVec <- c(outVec, paste("Sensitivity of diagnostic test: ", 
                                outChar))
            } else {
                outVec <- c(outVec, "Sensitivity of diagnostic test:  NO DATA")
            }
            if (length(object@costHerd) > 0){
                outVec <- c(outVec, paste("Cost per tested herd:           ", 
                                sprintf("%.2f", object@costHerd)))
            }
            if (length(object@costAnimal) > 0){
                outVec <- c(outVec, paste("Cost per tested animal:         ", 
                                sprintf("%.2f", object@costAnimal)))
            }
            outVec <- c(outVec,"","Survey Data:", "------------")
            if (length(object@nAnimalVec) > 0){
                outVec <- c(outVec, paste("Number of herds:            ", length(object@nAnimalVec)))
                outVec <- c(outVec, paste("Total number of animals:    ", sum(object@nAnimalVec)))
                outVec <- c(outVec, "Number of animals per herd:")
                tempSummary <- rbind(names(summary(object@nAnimalVec)),
                        summary(object@nAnimalVec))
                tempSummary <- format(tempSummary, justify = "centre")
                tempSummaryMat <- as.vector(apply(X = tempSummary, MARGIN = 1, FUN = function(X){
                                    Reduce(function(x,y) paste(x,y), X)
                                }))
                outVec <- c(outVec, tempSummaryMat)        
            } else { 
                outVec <- c(outVec, "No animal data.") 
            }       
            ## Risk group vector:
            if (length(object@riskGroupVec) > 0){
                outVec <- c(outVec, paste("Number of Risk groups:      ", length(unique(object@riskGroupVec))))
                outVec <- c(outVec, "Distribution of Risk groups:")
                tempSummary <- rbind(names(table(object@riskGroupVec)),
                        table(object@riskGroupVec))
                tempSummary <- format(tempSummary, justify = "centre")
                tempSummaryMat <- as.vector(apply(X = tempSummary, MARGIN = 1, FUN = function(X){
                                    Reduce(function(x,y) paste(x,y), X)
                                }))
                outVec <- c(outVec, tempSummaryMat)        
            } else { 
                outVec <- c(outVec, "No risk group data.") 
            }
            ## Risk value data:
            if (dim(object@riskValueData)[1] > 0){
                outVec <- c(outVec, "Risk values (relative):")
                tempVar <- capture.output(print(object@riskValueData, row.names = FALSE))
                tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                outVec <- c(outVec, tempVar)        
            } else { 
                outVec <- c(outVec, "No risk value data.") 
            }
            if (dim(object@populationData)[1] > 0){
                outVec <- c(outVec, "Additional population data: ")            
                tempStr <- capture.output(str(object@populationData))   
                tempStr <- gsub(pattern = "\t", replacement = "  ", x = tempStr)
                outVec <- c(outVec, tempStr)   
            } else {   
                outVec <- c(outVec, " no additional population data.")
            } 
            outTab <- as.table(matrix(outVec,ncol = 1))
            dimnames(outTab) <- list(rep("", length(outVec)), "")
            outTab
        }
)


###############################################################################
###############################################################################
## Class "LtdSampling"
##
## Contains the parameters and the data for a survey
## to substantiate freedom from disease 
## using "limited sampling". Additionally to the survey
## parameters (design prevalence, overall significance,
## intra-herd prevalence, sensitivity of the diagnostic test, 
## cost per tested animal and cost per tested herd) the object 
## contains the mean herd sensitivity, the number of herds to 
## be tested, the mean overall number of animals to be tested
## and the expected costs.
##
## Package: FFD
##
## Ian Kopacka
## 2010-07-13
###############################################################################
###############################################################################

## CLASS DEFINITION:
###############################################################################
###############################################################################
setClass(
        Class = "LtdSampling", 
        representation = representation(
                surveyData = "SurveyData",
                sampleSizeLtd = "numeric",
                meanHerdSensitivity = "numeric",
                meanHerdSensPerRG = "numeric",
                nHerds = "numeric",
                nHerdsPerRiskGroup = "numeric",
                nSampleFixVec = "numeric",
                probVec = "numeric",
                nAnimalsMean = "numeric",
                expectedCost = "numeric"     
        ),
        validity = function(object)
        {      
            if (any(object@sampleSizeLtd <= 0)){
                stop ("[LtdSampling validation]: Slot 'sampleSizeLtd' must be positive.")
            }
            if (any((object@sampleSizeLtd - as.integer(object@sampleSizeLtd)) != 0)){
                stop ("[LtdSampling validation]: Slot 'sampleSizeLtd' must be an integer.")
            }        
            if ((length(object@meanHerdSensitivity) > 0)){ 
                if((object@meanHerdSensitivity <= 0)|(object@meanHerdSensitivity >= 1))
                {
                    stop ("[LtdSampling validation]: Slot 'meanHerdSensitivity' must contain values in (0,1).")
                }
            }
            if (length(object@meanHerdSensPerRG) > 0){
                if (length(object@meanHerdSensPerRG) != dim(object@surveyData@riskValueData)[1]){
                    stop ("[LtdSampling validation]: Slot 'meanHerdSensPerRG' must have the same length as the number of rows in 'surveyData@riskValueData'.")              
                }
                if (!all(names(object@meanHerdSensPerRG) %in% object@surveyData@riskValueData[,1])){
                    stop ("[LtdSampling validation]: The names of slot 'meanHerdSensPerRG' must be contained in the first column of 'surveyData@riskValueData'.")               
                }   
                if((min(object@meanHerdSensPerRG) <= 0)|(max(object@meanHerdSensPerRG) >= 1))
                {
                    stop ("[LtdSampling validation]: Slot 'meanHerdSensPerRG' must contain values in (0,1).")
                }
            }       
            if (any(object@nHerds <= 0)){
                stop ("[LtdSampling validation]: Slot 'nHerds' must be positive.")
            }
            if (any((object@nHerds - as.integer(object@nHerds)) != 0)){
                stop ("[LtdSampling validation]: Slot 'nHerds' must be an integer.")
            } 
            if (length(object@nHerdsPerRiskGroup) > 0){
                if (length(object@nHerdsPerRiskGroup) != dim(object@surveyData@riskValueData)[1]){
                    stop ("[LtdSampling validation]: Slot 'nHerdsPerRiskGroup' must have the same length as the number of rows in 'surveyData@riskValueData'.")             
                }
                if (!all(names(object@nHerdsPerRiskGroup) %in% object@surveyData@riskValueData[,1])){
                    stop ("[LtdSampling validation]: The names of slot 'nHerdsPerRiskGroup' must be contained in the first column of 'surveyData@riskValueData'.")              
                }           
            }       
            if (length(object@nAnimalsMean) > 0){
                if (object@nAnimalsMean <= 0){
                    stop ("[LtdSampling validation]: Slot 'nAnimalsMean' must be a positive numeric")
                }
            }
            return(TRUE)
        }    
)

## METHODS:
###############################################################################
###############################################################################


#################################################### SHOW:
##########################################################

setMethod("show", signature(object="LtdSampling"),
        function(object){
            displayLimit <- 20
            cat("Object of class 'LtdSampling':\n")
            cat("Slots:\n")  
            cat("@surveyData:\n")
            cat("----------------------------\n")
            show(object@surveyData)
            cat("----------------------------\n")
            
            if (length(object@sampleSizeLtd) > 0){
                cat("@sampleSizeLtd:         ", sprintf("%d", object@sampleSizeLtd), "\n")
            } else {
                cat("@sampleSizeLtd:          NO DATA\n")
            }
            if (length(object@meanHerdSensitivity) > 0){
                cat("@meanHerdSensitivity:   ", sprintf("%.2f", object@meanHerdSensitivity), "\n")
            } else {
                cat("@meanHerdSensitivity:    NO DATA\n")
            }
            if (length(object@meanHerdSensPerRG) > 0){
                tempVar <- capture.output(print(object@meanHerdSensPerRG))
                tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                cat("@meanHerdSensPerRG:     ", tempVar[1], "\n")
                tempVar <- tempVar[-1]          
                for (ii in seq(along = tempVar)) cat("                        ", 
                            tempVar[ii], "\n")                      
            } else {
                cat("@meanHerdSensPerRG:      No risk groups used\n")
            }
            if (length(object@nHerds) > 0){
                cat("@nHerds:                ", sprintf("%d", object@nHerds), "\n")
            } else {
                cat("@nHerds:                 NO DATA\n")
            }
            if (length(object@nHerdsPerRiskGroup) > 0){
                tempVar <- capture.output(print(object@nHerdsPerRiskGroup))
                tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                cat("@nHerdsPerRiskGroup:    ", tempVar[1], "\n")
                tempVar <- tempVar[-1]          
                for (ii in seq(along = tempVar)) cat("                        ", 
                            tempVar[ii], "\n")                      
            } else {
                cat("@nHerdsPerRiskGroup:     No risk groups used\n")
            }   
            if (length(object@nSampleFixVec) > 0){
                cat(paste("@nSampleFixVec:          ", 
                                Reduce(function(x,y) paste(x,y, sep = " "), 
                                        as.character(object@nSampleFixVec)), "\n", sep = ""))                   
            } else {
                cat("@nSampleFixVec:          No risk groups used\n")
            }
            if (length(object@probVec) > 0){
                cat(paste("@probVec:                ", 
                                Reduce(function(x,y) paste(x,y, sep = " "), 
                                        as.character(object@probVec)), "\n", sep = ""))                 
            } else {
                cat("@probVec:                No risk groups used\n")
            }
            if (length(object@nAnimalsMean) > 0){
                cat("@nAnimalsMean:          ", sprintf("%.2f", object@nAnimalsMean), "\n")
            } else {
                cat("@nAnimalsMean:           NO DATA\n")
            }
            if (length(object@expectedCost) > 0){
                cat("@expectedCost:          ", sprintf("%.2f", object@expectedCost), "\n")
            } else {
                cat("@expectedCost:           NO DATA\n")
            }        
        }
)

################################################# SUMMARY:
##########################################################

setMethod("summary", signature(object="LtdSampling"),
        function(object, output = c("proportion", "percent")){
            output <- match.arg(output)
            displayLimit <- 20
            outVec <- c("LIMITED SAMPLING:","", 
                    as.vector(summary(object@surveyData, output)),"",
                    "Sampling strategy:",
                    "------------------")
            if (length(object@sampleSizeLtd) > 0){
                outVec <- c(outVec, paste("Fixed number of animals to test per herd: ", 
                                sprintf("%d",object@sampleSizeLtd)))
            } else {
                outVec <- c(outVec, "No data on number of animals to test per herd.")
            }
            if (length(object@meanHerdSensitivity) > 0){
                if (output == "proportion"){            
                    outChar <- sprintf("%.3f",object@meanHerdSensitivity)
                } else {
                    outChar <- paste(sprintf("%.2f",object@meanHerdSensitivity*100), "%")
                }   
                outVec <- c(outVec, paste("Mean herd sensitivity:                    ", 
                                outChar))
            } else {
                outVec <- c(outVec, "No data on mean herd sensitivity.")            
            }
            if (length(object@nHerds) > 0){
                outVec <- c(outVec, paste("Number of herds to test:                  ", 
                                sprintf("%d",object@nHerds)))
            } else {
                outVec <- c(outVec,"No data on number of herds to test.")
            }
            if (length(object@nHerdsPerRiskGroup) > 0){
                tempVar <- capture.output(print(object@nHerdsPerRiskGroup))
                tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                outVec <- c(outVec, paste("Number of herds to test per risk group:    ",
                                tempVar[1], sep = "")) 
                tempVar <- tempVar[-1]
                tempVar <- paste("                                          ", tempVar)
                outVec <- c(outVec, tempVar)  
            }   
            if (length(object@nAnimalsMean) > 0){
                outVec <- c(outVec, paste("Expected total number of animals to test: ", 
                                sprintf("%.2f",object@nAnimalsMean)))
            } else {
                outVec <- c(outVec, "No data on expected total number of animals to test.")
            }
            if (length(object@expectedCost) > 0){
                outVec <- c(outVec, paste("Expected total costs of the survey:       ", 
                                sprintf("%.2f",object@expectedCost)))
            } else {
                outVec <- c(outVec,"No data on expected total costs of the survey.")
            }  
            outTab <- as.table(matrix(outVec,ncol = 1))
            dimnames(outTab) <- list(rep("", length(outVec)), "")
            outTab
        }
)

#################################################### HTML:
##########################################################

setMethod("HTML", signature(x = "LtdSampling"),
        function(x, filename = "LtdSampling", outdir = getwd(), CSSFile="ffd.css", 
                Title = "Limited Sampling", append = TRUE,...){ 
            
            ## Create html-file:
            target <- HTMLInitFile(outdir = outdir , filename = filename, 
                    CSSFile = CSSFile, Title = Title,...)
            
            ## Write content:
            ###########################################################################
            ###########################################################################
            cat("\n\n<h1>Survey to substantiate Freedom From Disease</h1>\n", file = target,
                    append = TRUE)
            cat("<h2>Limited Sampling</h2>\n<br>\n", file = target, append = TRUE)
            
            ## Survey Parameters:
            cat("<h3>Survey Parameters:</h3>\n", file = target, append = TRUE)  
            surveyDat <- x@surveyData      
            theTable <- matrix(c(as.character(surveyDat@designPrevalence), 
                            as.character(surveyDat@alpha), 
                            as.character(surveyDat@intraHerdPrevalence), 
                            as.character(surveyDat@diagSensitivity), 
                            as.character(surveyDat@costHerd), 
                            as.character(surveyDat@costAnimal)), ncol = 1)
            indVec <- c(length(surveyDat@designPrevalence)>0, 
                    length(surveyDat@alpha)>0, 
                    length(surveyDat@intraHerdPrevalence)>0,
                    length(surveyDat@diagSensitivity)>0,
                    length(surveyDat@costHerd)>0,
                    length(surveyDat@costAnimal)>0)
            rownames(theTable) <- c("Design Prevalence ", "Overall Significance (alpha) ",
                    "Intra herd prevalence ", "Sensitivity of diagnostic test ",
                    "Cost per herd ", "Cost per animal ")[indVec]
            HTML(theTable, align = "left")
            cat("<br>\n\n", file = target, append = TRUE)
            
            if (length(surveyDat@nAnimalVec) > 0){
                cat("<h3>Data Description:</h3>\n", file = target, append = TRUE)  
                theTable <- matrix(c(length(surveyDat@nAnimalVec), 
                                sum(surveyDat@nAnimalVec), 
                                min(surveyDat@nAnimalVec),
                                median(surveyDat@nAnimalVec),
                                max(surveyDat@nAnimalVec)), ncol = 1)
                rownames(theTable) <- c("Number of herds ", 
                        "Total number of animals ",
                        "Minimal herd size ",
                        "Median herd size ",
                        "Maximal herd size ")
                HTML(theTable, align = "left")
                cat("<br>\n\n", file = target, append = TRUE)
            }
            
            ## Sample Parameters:
            cat("<h3>Sampling strategy:</h3>\n", file = target, append = TRUE)  
            theTable <- matrix(c(as.character(x@sampleSizeLtd),
                            as.character(round(x@meanHerdSensitivity,2)),
                            as.character(x@nHerds),
                            as.character(round(x@nAnimalsMean,2)), 
                            as.character(round(x@expectedCost,2))), ncol = 1)
            indVec <- c(length(x@sampleSizeLtd)>0, 
                    length(x@meanHerdSensitivity)>0, 
                    length(x@nHerds)>0,
                    length(x@nAnimalsMean)>0,
                    length(x@expectedCost)>0)
            rownames(theTable) <- c("Fixed number of animals to test per herd ", 
                    "Mean herd sensitivity ",
                    "Number of herds to test ", 
                    "Expected total number of animals to test ",
                    "Expected total costs of the survey ")[indVec]
            HTML(theTable, align = "left")
            cat("<br>\n\n", file = target, append = TRUE)   
            
            ## Close body and html tag
            cat("\n\n<br>\n<hr size=1>\n<font size=-1>\nGenerated on: <i>", format(Sys.time(), "%b %d %X %Y"), 
                    "</i> - <b>FFD</b>\n<hr size=1>\n\n", file = target, append = TRUE)
            
            cat("\n</body>\n</html>\n", file = target, append = TRUE)        
            
            ## Check if css.file exists:
            css.file <- file.path(outdir,CSSFile, fsep = .Platform$file.sep)
            ## If file does not exist create one:
            createStyleFile(css.file)                    
            
            ## Return value:
            invisible(x) 
        })

################################################## SAMPLE:
########################################################## 

setMethod("sample", signature(x = "LtdSampling"),
        function(x, size = c("fixed", "dynamic"), replace, prob){
            ## Set default value for flag 'size':
            if (missing(size)) size <- "fixed"
            size <- match.arg(size)         
			
            ## Sample with fixed sample size:
            if (size == "fixed"){ 
                ## No risk groups:
                ##################
                if (length(x@nHerdsPerRiskGroup) == 0){         
                    ## Sample x@nHerds herds using simple random sampling:
                    indexSample <- sort(sample(x = seq(along = x@surveyData@nAnimalVec),
                                    size = x@nHerds, replace = FALSE))
                    ## Compute the a-posteriori alpha error:
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec[indexSample], 
                            method = "limited", sampleSizeLtd = x@sampleSizeLtd, 
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, diagSpecificity = 1)            
                    aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaErrorVector, 
                            nPopulation = length(x@surveyData@nAnimalVec), 
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")  
                    nSampleArgVec <- numeric()       
                } else {
                    ## Sampling with risk groups:
                    #############################
                    indexVec <- seq(along = x@surveyData@nAnimalVec)
                    indexList <- split(x = indexVec, f = x@surveyData@riskGroupVec)
                    namesRiskGroupVec <- names(x@nHerdsPerRiskGroup)
                    nSampleRiskGroupVec <- as.vector(x@nHerdsPerRiskGroup)
                    indexSampleList <- lapply(seq(along = indexList), function(ii){
                                riskGroup <- names(indexList)[[ii]]     
                                nSample <- nSampleRiskGroupVec[namesRiskGroupVec == riskGroup]
                                if (length(indexList[[ii]]) > 1){
                                    indexOutVec <- sample(x = indexList[[ii]],
                                            size = nSample, replace = FALSE)    
                                } else {
                                    if (nSample == 1){
                                        indexOutVec <- indexList[[ii]]
                                    } else {
                                        stop ("Sample size larger than Population size.")   
                                    }                       
                                }               
                            })
                    ## Concatenate the sample indices for the risk groups:
                    indexSample <- sort(Reduce(function(x,y) c(x,y), indexSampleList),
                            decreasing = FALSE)
                    ## Compute the a-posteriori alpha error:
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec[indexSample], 
                            method = "limited", sampleSizeLtd = x@sampleSizeLtd, 
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, diagSpecificity = 1)            
                    ## Determine the number of herds in each risk group:
                    riskValueDf <- x@surveyData@riskValueData[,1:2]
                    names(riskValueDf) <- c("riskGroup", "riskValues")
                    riskValueDf$riskGroup <- as.character(riskValueDf$riskGroup)
                    riskValueDf$id <- seq(along = riskValueDf[,1])
                    riskGroupTab <- table(x@surveyData@riskGroupVec)
                    riskGroupDf <- data.frame(riskGroup = as.character(names(riskGroupTab)), 
                            nPopulation = as.vector(riskGroupTab))
                    riskValueDf <- merge(x = riskValueDf, y = riskGroupDf, by = "riskGroup",
                            sort = FALSE)   
                    riskValueDf <- riskValueDf[order(riskValueDf$id, decreasing = FALSE),]
                    aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                            groupVec = x@surveyData@riskGroupVec[indexSample], 
                            groupLevels = riskValueDf$riskGroup,
                            nPopulationVec = riskValueDf$nPopulation,
                            nRelRiskVec = riskValueDf$riskValues,
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")  
					nSampleArgVec <- nSampleRiskGroupVec
                }           
            } else {
                ## No risk groups:
                ##################
                if (length(x@nHerdsPerRiskGroup) == 0){
                    ## Data frame with herd-based alpha-errors for each herd
                    ## (=1-herd sensitivity):
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec, 
                            method = "limited", sampleSizeLtd = x@sampleSizeLtd, 
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, diagSpecificity = 1)
                    alphaDataFrame <- data.frame(size = x@surveyData@nAnimalVec, 
                            alpha = alphaErrorVector, id = seq(along = x@surveyData@nAnimalVec))
                    ## Permutate data frame:
                    alphaDataFrame <- alphaDataFrame[sample(x = alphaDataFrame$id, 
                                    size = length(alphaDataFrame$id), replace = FALSE),]
                    ## Dynamic sampling: find the smallest integer k such that the 
                    ## sample containing the first k rows of the data frame have
                    ## an a-posteriori error <=  x@surveyData@alpha    
                    k <- x@nHerds  
                    aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                            nPopulation = length(x@surveyData@nAnimalVec), 
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")
                    if (aPostAlpha > x@surveyData@alpha){
                        while (aPostAlpha > x@surveyData@alpha){
                            k <- k + 1
                            aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                                    nPopulation = length(x@surveyData@nAnimalVec), 
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")               
                        }            
                    } else {
                        while (aPostAlpha <= x@surveyData@alpha){
                            k <- k - 1
                            aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                                    nPopulation = length(x@surveyData@nAnimalVec), 
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")               
                        }         
                        k <- k + 1
                        aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                                nPopulation = length(x@surveyData@nAnimalVec), 
                                nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                method = "approx")                               
                    }
                    indexSample <- sort(alphaDataFrame$id[1:k])  
                    nSampleArgVec <- numeric()
                } else {
                    ## Sampling with risk groups:
                    #############################
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec, 
                            method = "limited", sampleSizeLtd = x@sampleSizeLtd, 
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, diagSpecificity = 1)
                    alphaDataFrame <- data.frame(size = x@surveyData@nAnimalVec, 
                            alpha = alphaErrorVector, id = seq(along = x@surveyData@nAnimalVec),
                            riskGroupVec = as.character(x@surveyData@riskGroupVec))
                    ## Permutate data frame:
                    alphaDataFrame <- alphaDataFrame[sample(x = alphaDataFrame$id, 
                                    size = length(alphaDataFrame$id), replace = FALSE),]    
                    alphaList <- split(x = alphaDataFrame, f = alphaDataFrame$riskGroupVec)
                    ## Dynamic sampling: find the smallest integer k such that the 
                    ## sample containing the first k rows of the data frame have
                    ## an a-posteriori error <=  x@surveyData@alpha    
                    k <- x@nHerds   
                    nSampleFixVec <- x@nSampleFixVec
                    names(nSampleFixVec) <- x@surveyData@riskValueData$riskGroup                
                    nSamplePropVec <- as.vector(x@probVec*
                                    table(x@surveyData@riskGroupVec)[x@surveyData@riskValueData$riskGroup]
                                            [is.na(nSampleFixVec)])
                    ## Proportion for sample sizes:
                    if (sum(is.na(nSampleFixVec)) > 1){
                        nSamplePropVec <- nSamplePropVec/sum(nSamplePropVec)
                    } else {
                        nSamplePropVec <- 1
                    }
                    ## Population size of the risk groups:
                    nPopulationVec <- table(x@surveyData@riskGroupVec)[x@surveyData@riskValueData$riskGroup]
                    ## Names of the risk groups:
                    groupLevels <- as.character(x@surveyData@riskValueData$riskGroup)
                    ## Sample size for the risk groups:
                    nSampleArgVec <- nSampleFixVec
                    nSampleArgVec[is.na(nSampleFixVec)] <- 
                            roundConstantSum(k*nSamplePropVec, output = 0)
                    alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                                riskGroup <- names(alphaList)[[ii]]     
                                nSample <- nSampleArgVec[riskGroup] 
                                if (nSample > 0){
                                    alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                                } else {
                                    alphaOut <- numeric()
                                }                           
                            }) 
                    alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                    groupVec <- names(alphaList)    
                    groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                    aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                            groupVec = groupVec, 
                            groupLevels = groupLevels,
                            nPopulationVec = nPopulationVec[groupLevels],
                            nRelRiskVec = x@surveyData@riskValueData$riskValues,
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")
                    if (aPostAlpha > x@surveyData@alpha){
                        while (aPostAlpha > x@surveyData@alpha){
                            k <- k + 1
                            nSampleArgVec <- nSampleFixVec
                            nSampleArgVec[is.na(nSampleFixVec)] <- 
                                    roundConstantSum(k*nSamplePropVec, output = 0)
                            alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                                        riskGroup <- names(alphaList)[[ii]]     
                                        nSample <- nSampleArgVec[riskGroup] 
                                        if (nSample > 0){
                                            alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                                        } else {
                                            alphaOut <- numeric()
                                        }                           
                                    }) 
                            alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                            groupVec <- names(alphaList)    
                            groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                            aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                                    groupVec = groupVec, 
                                    groupLevels = groupLevels,
                                    nPopulationVec = nPopulationVec[groupLevels],
                                    nRelRiskVec = x@surveyData@riskValueData$riskValues,
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")               
                        }            
                    } else {
                        while (aPostAlpha <= x@surveyData@alpha){
                            k <- k - 1
                            nSampleArgVec <- nSampleFixVec
                            nSampleArgVec[is.na(nSampleFixVec)] <- 
                                    roundConstantSum(k*nSamplePropVec, output = 0)
                            alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                                        riskGroup <- names(alphaList)[[ii]]     
                                        nSample <- nSampleArgVec[riskGroup] 
                                        if (nSample > 0){
                                            alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                                        } else {
                                            alphaOut <- numeric()
                                        }                           
                                    }) 
                            alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                            groupVec <- names(alphaList)    
                            groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                            aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                                    groupVec = groupVec, 
                                    groupLevels = groupLevels,
                                    nPopulationVec = nPopulationVec[groupLevels],
                                    nRelRiskVec = x@surveyData@riskValueData$riskValues,
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")              
                        }         
                        k <- k + 1
                        nSampleArgVec <- nSampleFixVec
                        nSampleArgVec[is.na(nSampleFixVec)] <- 
                                roundConstantSum(k*nSamplePropVec, output = 0)
                        alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                                    riskGroup <- names(alphaList)[[ii]]     
                                    nSample <- nSampleArgVec[riskGroup] 
                                    if (nSample > 0){
                                        alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                                    } else {
                                        alphaOut <- numeric()
                                    }                           
                                }) 
                        alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                        groupVec <- names(alphaList)    
                        groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                        aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                                groupVec = groupVec, 
                                groupLevels = groupLevels,
                                nPopulationVec = nPopulationVec[groupLevels],
                                nRelRiskVec = x@surveyData@riskValueData$riskValues,
                                nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                method = "approx")                              
                    }
                    indexSampleList <- lapply(seq(along = alphaList), function(ii){
                        riskGroup <- names(alphaList)[[ii]]     
                        nSample <- nSampleArgVec[riskGroup] 
                        if (nSample > 0){
                            indexOut <- alphaList[[ii]]$id[1:nSample]
                        } else{
                            indexOut <- numeric()
                        }
                        return(indexOut)
                    })
                    indexSample <- sort(Reduce(function(x,y) c(x,y), 
                        indexSampleList), decreasing = FALSE)                   
                }
            }
            
            #######################################################################
            ## Return value:
            out <- list(indexSample = indexSample, aPostAlpha = aPostAlpha)         
            ## Add population data of sample:    
            if (dim(x@surveyData@populationData)[1] > 0){        
                out$sample <- x@surveyData@populationData[indexSample,]
            }
            ## Add sample size of the risk groups:
            if (length(nSampleArgVec) > 0){
                out$nSamplePerRiskGroup <- nSampleArgVec                
            }           
            return(out)              
        }
)

###############################################################################
###############################################################################
## Class "LtdSamplingSummary"
##
## Contains the parameters and the data for a survey
## to substantiate freedom from disease 
## using "limited sampling". Additionally to the survey
## parameters (design prevalence, overall significance,
## intra-herd prevalence, sensitivity of the diagnostic test, 
## cost per tested animal and cost per tested herd) the object 
## contains the mean herd sensitivity, the number of herds to 
## be tested, the mean overall number of animals to be tested
## and the expected costs for a range of possible sample limits
## (= fixed number of animals to test per herd).
##
## Package: FFD
##
## Ian Kopacka
## 2010-07-13
###############################################################################
###############################################################################
setClass(
        Class = "LtdSamplingSummary", 
        representation = representation(
                surveyData = "SurveyData",
                sampleSizeLtdVec = "numeric",
                meanHerdSensVec = "numeric",
                meanHerdSensPerRGMx = "matrix",
                nHerdsVec = "numeric",
                nHerdsPerRiskGroupMx = "matrix",
                nSampleFixVec = "numeric",
                probVec = "numeric",
                nAnimalsMeanVec = "numeric",
                expectedCostVec = "numeric"     
        ),
        validity = function(object)
        {      
            if (any(object@sampleSizeLtdVec <= 0)){
                stop ("[LtdSamplingSummary validation]: Slot 'sampleSizeLtdVec' must contain positive values.")
            }
            if (any((object@sampleSizeLtdVec - as.integer(object@sampleSizeLtdVec)) != 0)){
                stop ("[LtdSamplingSummary validation]: Slot 'sampleSizeLtdVec' must contain an integer vector.")
            }      
            if (length(object@meanHerdSensVec) != length(object@sampleSizeLtdVec)){
                stop ("[LtdSamplingSummary validation]: Slots 'meanHerdSensVec' and 'sampleSizeLtdVec' must have the same length.")         
            }       
            if (any((object@meanHerdSensVec <= 0)|(object@meanHerdSensVec >= 1))){
                stop ("[LtdSamplingSummary validation]: Slot 'meanHerdSensVec' must contain values in (0,1).")
            }  
            if (dim(object@meanHerdSensPerRGMx)[1] > 0){
                if (dim(object@meanHerdSensPerRGMx)[2] != dim(object@surveyData@riskValueData)[1]){
                    stop ("[LtdSamplingSummary validation]: Slot 'meanHerdSensPerRGMx' have the same number of columns as the number of rows in 'surveyData@riskValueData'.")              
                }
                if (!all(colnames(object@meanHerdSensPerRGMx) %in% object@surveyData@riskValueData[,1])){
                    stop ("[LtdSamplingSummary validation]: The names of slot 'meanHerdSensPerRGMx' must be contained in the first column of 'surveyData@riskValueData'.")              
                } 
                if (dim(object@meanHerdSensPerRGMx)[1] != length(object@sampleSizeLtdVec)){
                    stop ("[LtdSamplingSummary validation]: The number of rows of slot 'meanHerdSensPerRGMx' must equal the length of 'sampleSizeLtdVec'.")
                }
                if (any((object@meanHerdSensPerRGMx <= 0)|(object@meanHerdSensPerRGMx >= 1))){
                    stop ("[LtdSamplingSummary validation]: Slot 'meanHerdSensPerRGMx' must contain values in (0,1).")
                }  
            }
            if (length(object@nHerdsVec) != length(object@sampleSizeLtdVec)){
                stop ("[LtdSamplingSummary validation]: Slots 'nHerdsVec' and 'sampleSizeLtdVec' must have the same length.")           
            }
            if (any(object@nHerdsVec <= 0)){
                stop ("[LtdSamplingSummary validation]: Slot 'nHerdsVec' must contain positive values.")
            }
            if (any((object@nHerdsVec - as.integer(object@nHerdsVec)) != 0)){
                stop ("[LtdSamplingSummary validation]: Slot 'nHerdsVec' must contain an integer vector.")
            } 
            if (dim(object@nHerdsPerRiskGroupMx)[1] > 0){
                if (dim(object@nHerdsPerRiskGroupMx)[2] != dim(object@surveyData@riskValueData)[1]){
                    stop ("[LtdSamplingSummary validation]: Slot 'nHerdsPerRiskGroupMx' have the same number of columns as the number of rows in 'surveyData@riskValueData'.")              
                }
                if (!all(colnames(object@nHerdsPerRiskGroupMx) %in% object@surveyData@riskValueData[,1])){
                    stop ("[LtdSamplingSummary validation]: The names of slot 'nHerdsPerRiskGroupMx' must be contained in the first column of 'surveyData@riskValueData'.")              
                } 
                if (dim(object@nHerdsPerRiskGroupMx)[1] != length(object@sampleSizeLtdVec)){
                    stop ("[LtdSamplingSummary validation]: The number of rows of slot 'nHerdsPerRiskGroupMx' must equal the length of 'sampleSizeLtdVec'.")
                }
            }
            if (length(object@nAnimalsMeanVec) != length(object@sampleSizeLtdVec)){
                stop ("[LtdSamplingSummary validation]: Slots 'nAnimalsMeanVec' and 'sampleSizeLtdVec' must have the same length.")         
            }
            if (any(object@nAnimalsMeanVec <= 0)){
                stop ("[LtdSamplingSummary validation]: Slot 'nAnimalsMeanVec' must contain positive values")
            }        
            return(TRUE)
        }    
)

## METHODS:
###############################################################################
###############################################################################


#################################################### SHOW:
##########################################################

setMethod("show", signature(object="LtdSamplingSummary"),
        function(object){
            displayLimit <- 20
            cat("Object of class 'LtdSamplingSummary':\n")
            cat("Slots:\n")  
            cat("@surveyData:\n")
            cat("----------------------------\n")
            show(object@surveyData)
            cat("----------------------------\n")
            
            cat("@sampleSizeLtdVec:")
            if (length(object@sampleSizeLtdVec) > 0){
                cat("\n")
                if (length(object@sampleSizeLtdVec) > displayLimit){
                    show(object@sampleSizeLtdVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@sampleSizeLtdVec), "displayed)\n")
                } else {
                    show(object@sampleSizeLtdVec)
                }
            } else {        
                cat("       NO DATA\n")
            }
            cat("@meanHerdSensVec:")
            if (length(object@meanHerdSensVec) > 0){
                cat("\n")
                if (length(object@meanHerdSensVec) > displayLimit){
                    show(object@meanHerdSensVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@meanHerdSensVec), "displayed)\n")
                } else {
                    show(object@meanHerdSensVec)
                }
            } else {        
                cat("        NO DATA\n")
            }   
            cat("@meanHerdSensPerRGMx:")
            if (dim(object@meanHerdSensPerRGMx)[1] > 0){
                cat("\n")
                if (dim(object@meanHerdSensPerRGMx)[1] > displayLimit){
                    print(as.data.frame(object@meanHerdSensPerRGMx[1:displayLimit,]), row.names = FALSE)
                    cat("   (only first",displayLimit,"rows of",dim(object@meanHerdSensPerRGMx)[1], "displayed)\n")
                } else {
                    print(as.data.frame(object@meanHerdSensPerRGMx), row.names = FALSE)
                }
            } else {        
                cat("   NO DATA\n")
            }
            cat("@nHerdsVec:")
            if (length(object@nHerdsVec) > 0){
                cat("\n")
                if (length(object@nHerdsVec) > displayLimit){
                    show(object@nHerdsVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@nHerdsVec), "displayed)\n")
                } else {
                    show(object@nHerdsVec)
                }
            } else {        
                cat("              NO DATA\n")
            }
            cat("@nHerdsPerRiskGroupMx:")
            if (dim(object@nHerdsPerRiskGroupMx)[1] > 0){
                cat("\n")
                if (dim(object@nHerdsPerRiskGroupMx)[1] > displayLimit){
                    print(as.data.frame(object@nHerdsPerRiskGroupMx[1:displayLimit,]), row.names = FALSE)
                    cat("   (only first",displayLimit,"rows of",dim(object@nHerdsPerRiskGroupMx)[1], "displayed)\n")
                } else {
                    print(as.data.frame(object@nHerdsPerRiskGroupMx), row.names = FALSE)
                }
            } else {        
                cat("   NO DATA\n")
            }
            cat("@nSampleFixVec:")
            if (length(object@nSampleFixVec) > 0){
                cat("\n")
                if (length(object@nSampleFixVec) > displayLimit){
                    show(object@nSampleFixVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@nSampleFixVec), "displayed)\n")
                } else {
                    show(object@nSampleFixVec)
                }
            } else {        
                cat("          NO DATA\n")
            }
            cat("@probVec:")
            if (length(object@probVec) > 0){
                cat("\n")
                if (length(object@probVec) > displayLimit){
                    show(object@probVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@probVec), "displayed)\n")
                } else {
                    show(object@probVec)
                }
            } else {        
                cat("                NO DATA\n")
            }       
            cat("@nAnimalsMeanVec:")
            if (length(object@nAnimalsMeanVec) > 0){
                cat("\n")
                if (length(object@nAnimalsMeanVec) > displayLimit){
                    show(object@nAnimalsMeanVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@nAnimalsMeanVec), "displayed)\n")
                } else {
                    show(object@nAnimalsMeanVec)
                }
            } else {        
                cat("        NO DATA\n")
            }
            cat("@expectedCostVec:")
            if (length(object@expectedCostVec) > 0){
                cat("\n")
                if (length(object@expectedCostVec) > displayLimit){
                    show(object@expectedCostVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@expectedCostVec), "displayed)\n")
                } else {
                    show(object@expectedCostVec)
                }
            } else {        
                cat("        NO DATA\n")
            }               
        }
)

################################################# SUMMARY:
##########################################################

#setMethod("summary", signature(object="LtdSamplingSummary"),
#    function(object){
#        displayLimit <- 20
#        cat("LIMITED SAMPLING:\n\n")
#        summary(object@surveyData)
#        cat("\n")
#        cat("Cost optimal sampling strategy:\n")
#        cat("-------------------------------\n")
#        if(length(object@expectedCostVec) > 0){
#            indCostOpt <- which.min(object@expectedCostVec)
#            cat("Fixed number of animals to test per herd: ", sprintf("%d",object@sampleSizeLtdVec[indCostOpt]), "\n")
#            cat("Mean herd sensitivity:                    ", sprintf("%.2f",object@meanHerdSensVec[indCostOpt]), "\n")
#            cat("Number of herds to test:                  ", sprintf("%d",object@nHerdsVec[indCostOpt]), "\n")
#            cat("Expected total number of animals to test: ", sprintf("%.2f",object@nAnimalsMeanVec[indCostOpt]), "\n")
#            cat("Expected total costs of the survey:       ", sprintf("%.2f",object@expectedCostVec[indCostOpt]), "\n")
#        } else {
#            cat("No data regarding costs.\n")
#        }        
#    }
#)
setMethod("summary", signature(object = "LtdSamplingSummary"),
        function(object, output = c("proportion", "percent")){
            output <- match.arg(output)
            displayLimit <- 20
            outVec <- c("LIMITED SAMPLING DIAGNOSTICS:", "", 
                    as.vector(summary(object@surveyData, output)), "",
                    "Cost optimal sampling strategy:",
                    "-------------------------------")        
            if(length(object@expectedCostVec) > 0){
                indCostOpt <- which.min(object@expectedCostVec)
                if (output == "proportion"){
                    meanHerdSensChar <- sprintf("%.3f",object@meanHerdSensVec[indCostOpt])              
                } else {
                    meanHerdSensChar <- paste(sprintf("%.2f",object@meanHerdSensVec[indCostOpt]*100),
                            "%")
                }           
                outVec <- c(outVec,
                        paste("Fixed number of animals to test per herd: ", 
                                sprintf("%d",object@sampleSizeLtdVec[indCostOpt])),
                        paste("Mean herd sensitivity:                    ", 
                                meanHerdSensChar))
                if (dim(object@meanHerdSensPerRGMx)[1] > 0){
                    if (output == "proportion"){
                        meanHerdSensRG <- round(object@meanHerdSensPerRGMx[indCostOpt,],3)              
                    } else {
                        meanHerdSensRG <- round(object@meanHerdSensPerRGMx[indCostOpt,]*100,2)
                    }               
                    tempVar <- capture.output(print(meanHerdSensRG))
                    tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                    outVec <- c(outVec, paste("Mean herd sensitivity per risk group:      ",
                                    tempVar[1], sep = "")) 
                    tempVar <- tempVar[-1]
                    tempVar <- paste("                                          ", tempVar)
                    outVec <- c(outVec, tempVar)  
                }
                outVec <- c(outVec,
                        paste("Number of herds to test:                  ", 
                                sprintf("%d",object@nHerdsVec[indCostOpt])))    
                if (dim(object@nHerdsPerRiskGroupMx)[1] > 0){
                    tempVar <- capture.output(print(object@nHerdsPerRiskGroupMx[indCostOpt,]))
                    tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                    outVec <- c(outVec, paste("Number of herds to test per risk group:    ",
                                    tempVar[1], sep = "")) 
                    tempVar <- tempVar[-1]
                    tempVar <- paste("                                          ", tempVar)
                    outVec <- c(outVec, tempVar)  
                }
                
                outVec <- c(outVec,
                        paste("Expected total number of animals to test: ", 
                                sprintf("%.2f",object@nAnimalsMeanVec[indCostOpt])),
                        paste("Expected total costs of the survey:       ", 
                                sprintf("%.2f",object@expectedCostVec[indCostOpt])))
            } else {
                outVec <- c(outVec,"No data regarding costs.")
            } 
            outTab <- as.table(matrix(outVec,ncol = 1))
            dimnames(outTab) <- list(rep("", length(outVec)), "")
            outTab
        }
)

#################################################### HTML:
##########################################################

setMethod("HTML", signature(x = "LtdSamplingSummary"),
        function(x, filename = "LtdSamplingSummary", outdir = getwd(), CSSFile = "ffd.css", 
                Title = "Limited Sampling Diagnostics", append = TRUE,...){ 
            
            ## Create html-file:
            target <- HTMLInitFile(outdir = outdir , filename = filename, 
                    CSSFile = CSSFile, Title = Title,...)
            
            ## Write content:
            ###########################################################################
            ###########################################################################
            cat("\n\n<h1>Survey to substantiate Freedom From Disease</h1>\n", file = target,
                    append = TRUE)
            cat("<h2>Limited Sampling</h2>\n<br>\n", file = target, append = TRUE)
            
            ## Survey Parameters:
            cat("<h3>Survey Parameters:</h3>\n", file = target, append = TRUE)  
            surveyDat <- x@surveyData      
            surveyDat <- x@surveyData      
            theTable <- matrix(c(as.character(surveyDat@designPrevalence), 
                            as.character(surveyDat@alpha), 
                            as.character(surveyDat@intraHerdPrevalence), 
                            as.character(surveyDat@diagSensitivity), 
                            as.character(surveyDat@costHerd), 
                            as.character(surveyDat@costAnimal)), ncol = 1)
            indVec <- c(length(surveyDat@designPrevalence)>0, 
                    length(surveyDat@alpha)>0, 
                    length(surveyDat@intraHerdPrevalence)>0,
                    length(surveyDat@diagSensitivity)>0,
                    length(surveyDat@costHerd)>0,
                    length(surveyDat@costAnimal)>0)
            rownames(theTable) <- c("Design Prevalence ", "Overall Significance (alpha) ",
                    "Intra herd prevalence ", "Sensitivity of diagnostic test ",
                    "Cost per herd ", "Cost per animal ")[indVec]
            HTML(theTable, align = "left")
            cat("<br>\n\n", file = target, append = TRUE)
            
            ## Data description:
            if (length(surveyDat@nAnimalVec) > 0){
                cat("<h3>Data Description:</h3>\n", file = target, append = TRUE)  
                theTable <- matrix(c(length(surveyDat@nAnimalVec), 
                                sum(surveyDat@nAnimalVec), 
                                min(surveyDat@nAnimalVec),
                                median(surveyDat@nAnimalVec),
                                max(surveyDat@nAnimalVec)), ncol = 1)
                rownames(theTable) <- c("Number of herds ", 
                        "Total number of animals ",
                        "Minimal herd size ",
                        "Median herd size ",
                        "Maximal herd size ")
                HTML(theTable, align = "left")
                cat("<br>\n\n", file = target, append = TRUE)
            }
            
            ## Cost optimal sample Parameters:
            cat("<h3>Cost optimal sampling strategy:</h3>\n", file = target, append = TRUE)
            if(length(x@expectedCostVec) > 0){
                indCostOpt <- which.min(x@expectedCostVec)
                theTable <- matrix(c(as.character(x@sampleSizeLtdVec[indCostOpt]),
                                as.character(round(x@meanHerdSensVec[indCostOpt],2)),
                                as.character(x@nHerdsVec[indCostOpt]),
                                as.character(round(x@nAnimalsMeanVec[indCostOpt],2)), 
                                as.character(round(x@expectedCostVec[indCostOpt],2))), ncol = 1)
                indVec <- c(length(x@sampleSizeLtdVec)>0, 
                        length(x@meanHerdSensVec)>0, 
                        length(x@nHerdsVec)>0,
                        length(x@nAnimalsMeanVec)>0,
                        length(x@expectedCostVec)>0)
                rownames(theTable) <- c("Fixed number of animals to test per herd ", 
                        "Mean herd sensitivity ",
                        "Number of herds to test ", 
                        "Expected total number of animals to test ",
                        "Expected total costs of the survey ")[indVec]
                HTML(theTable, align = "left")                   
            } else {
                HTML("No data regarding costs.")
            }
            cat("<br>\n\n", file = target, append = TRUE)
            
            ## Diagnostic plots:
            cat("<h3>Diagnostic plots:</h3>\n", file = target, append = TRUE)
            survey.Data <- x@surveyData
            
            if(length(survey.Data@nAnimalVec) > 0){   
                par(mfrow = c(1, 1))     
                ## Mean herd sensitivity:
                plot(x@sampleSizeLtdVec, x@meanHerdSensVec, type = "l",
                        xlab = "Sample limit", ylab = "Mean herd sensitivity")
                HTMLplot(Width = 500, Height = 500, Caption = "Herd sensitivity",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"1", sep = "_"))
                ## Number of herds to be tested:
                plot(x@sampleSizeLtdVec, x@nHerdsVec, type = "l",
                        xlab = "Sample limit", ylab = "No. of herds to be tested")
                HTMLplot(Width = 500, Height = 500, Caption = "Sample size",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"2", sep = "_"))
                ## Total number of animals to be tested:
                plot(x@sampleSizeLtdVec, x@nAnimalsMeanVec, type = "l",
                        xlab = "Sample limit", ylab = "Expected total no. of animals to be tested")
                HTMLplot(Width = 500, Height = 500, Caption = "Total number of animals to test",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"3", sep = "_"))
                ## Expected cost:
                plot(x@sampleSizeLtdVec, x@expectedCostVec, type = "l",
                        xlab = "Sample limit", ylab = "Expected cost") 
                HTMLplot(Width = 500, Height = 500, Caption = "Expected costs of the survey",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"4", sep = "_"))                         
            } else {
                HTML("No data to produce plots.")
            }        
            
            ## Close body and html tag
            cat("\n\n<br>\n<hr size=1>\n<font size=-1>\nGenerated on: <i>", format(Sys.time(), "%b %d %X %Y"), 
                    "</i> - <b>FFD</b>\n<hr size=1>\n\n", file = target, append = TRUE)
            
            cat("\n</body>\n</html>\n", file = target, append = TRUE)        
            
            ## Check if css.file exists:
            css.file <- file.path(outdir,CSSFile, fsep = .Platform$file.sep)
            ## If file does not exist create one:
            createStyleFile(css.file)                    
            
            ## Return value:
            invisible(x) 
        })


#################################################### PLOT:
##########################################################

setMethod("plot", signature(x = "LtdSamplingSummary"),
        function(x,y,...){
            survey.Data <- x@surveyData
            if(length(survey.Data@nAnimalVec > 0)){        
                par(mfrow = c(2, 2))      
                # ## Distribution of the herd sizes:
#            hist(survey.Data@nAnimalVec, xlab = "Herd size", ylab = "Frequency",
#                main = "", col = "#DBA400")               
                ## Mean herd sensitivity:
                plot(x@sampleSizeLtdVec, x@meanHerdSensVec, type = "l",
                        xlab = "Sample limit", ylab = "Mean herd sensitivity")
                ## Number of herds to be tested:
                plot(x@sampleSizeLtdVec, x@nHerdsVec, type = "l",
                        xlab = "Sample limit", ylab = "No. of herds to be tested")
                ## Total number of animals to be tested:
                plot(x@sampleSizeLtdVec, x@nAnimalsMeanVec, type = "l",
                        xlab = "Sample limit", ylab = "Expected total no. of animals to be tested")
                ## Expected cost:
                plot(x@sampleSizeLtdVec, x@expectedCostVec, type = "l",
                        xlab = "Sample limit", ylab = "Expected cost")  
                ## Titel:
                par(oma = c(2,1,3,1)) 
                title("Analysis limited sampling", outer = TRUE)                                  
            } else {
                cat("Object of class 'LtdSamplingSummary' contains no data.\n")
            }
        }
)


###############################################################################
###############################################################################
## Class "IndSampling"
##
## Contains the parameters and the data for a survey
## to substantiate freedom from disease 
## using "individual sampling". Additionally to the survey
## parameters (design prevalence, overall significance,
## intra-herd prevalence, sensitivity of the diagnostic test, 
## cost per tested animal and cost per tested herd) the object 
## contains the number of herds to be tested, the mean overall number of 
## animals to be tested, the expected costs and a lookup table containing the
## number of animals to test depending on the herd size.
##
## Package: FFD
##
## Ian Kopacka
## 2010-07-16
###############################################################################
###############################################################################

## CLASS DEFINITION:
###############################################################################
###############################################################################
setClass(
        Class = "IndSampling", 
        representation = representation(
                surveyData = "SurveyData",
                herdSensitivity = "numeric",
                nHerds = "numeric",
                nHerdsPerRiskGroup = "numeric",
                nSampleFixVec = "numeric",
                probVec = "numeric",
                nAnimalsMean = "numeric",
                expectedCost = "numeric",
                lookupTable = "matrix"     
        ),
        validity = function(object)
        {      
            if (any((object@herdSensitivity <= 0)|(object@herdSensitivity >= 1))){
                stop ("[IndSampling validation]: Slot 'herdSensitivity' must contain values in (0,1).")
            }
            if (any(object@nHerds <= 0)){
                stop ("[IndSampling validation]: Slot 'nHerds' must contain positive values.")
            }
            if (any((object@nHerds - as.integer(object@nHerds)) != 0)){
                stop ("[IndSampling validation]: Slot 'nHerds' must contain an integer vector.")
            }
            if (length(object@nHerdsPerRiskGroup) > 0){
                if (length(object@nHerdsPerRiskGroup) != dim(object@surveyData@riskValueData)[1]){
                    stop ("[IndSampling validation]: Slot 'nHerdsPerRiskGroup' have the same length as the number of rows in 'surveyData@riskValueData'.")              
                }
                if (!all(names(object@nHerdsPerRiskGroup) %in% object@surveyData@riskValueData[,1])){
                    stop ("[IndSampling validation]: The names of slot 'nHerdsPerRiskGroup' must be contained in the first column of 'surveyData@riskValueData'.")              
                }           
            } 
            if (any(object@nAnimalsMean <= 0)){
                stop ("[IndSampling validation]: Slot 'nAnimalsMean' must contain positive values")
            }
            return(TRUE)
        }    
)

## METHODS:
###############################################################################
###############################################################################


#################################################### SHOW:
##########################################################

setMethod("show", signature(object="IndSampling"),
        function(object){
            displayLimit <- 20
            cat("Object of class 'IndSampling':\n")
            cat("Slots:\n")  
            cat("@surveyData:\n")
            cat("----------------------------\n")
            show(object@surveyData)
            cat("----------------------------\n")
            
            if (length(object@herdSensitivity) > 0){
                cat("@herdSensitivity:       ", sprintf("%.2f", object@herdSensitivity), "\n")
            } else {
                cat("@herdSensitivity:        NO DATA\n")
            }
            if (length(object@nHerds) > 0){
                cat("@nHerds:                ", sprintf("%d", object@nHerds), "\n")
            } else {
                cat("@nHerds:                 NO DATA\n")
            }
            if (length(object@nHerdsPerRiskGroup) > 0){
                tempVar <- capture.output(print(object@nHerdsPerRiskGroup))
                tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                cat("@nHerdsPerRiskGroup:    ", tempVar[1], "\n")
                tempVar <- tempVar[-1]          
                for (ii in seq(along = tempVar)) cat("                        ", 
                            tempVar[ii], "\n")                      
            } else {
                cat("@nHerdsPerRiskGroup:     No risk groups used\n")
            }   
            if (length(object@nSampleFixVec) > 0){
                cat(paste("@nSampleFixVec:          ", 
                                Reduce(function(x,y) paste(x,y, sep = " "), 
                                        as.character(object@nSampleFixVec)), "\n", sep = ""))                   
            } else {
                cat("@nSampleFixVec:          No risk groups used\n")
            }
            if (length(object@probVec) > 0){
                cat(paste("@probVec:                ", 
                                Reduce(function(x,y) paste(x,y, sep = " "), 
                                        as.character(object@probVec)), "\n", sep = ""))                 
            } else {
                cat("@probVec:                No risk groups used\n")
            }
            if (length(object@nAnimalsMean) > 0){
                cat("@nAnimalsMean:          ", sprintf("%.2f", object@nAnimalsMean), "\n")
            } else {
                cat("@nAnimalsMean:           NO DATA\n")
            }
            if (length(object@expectedCost) > 0){
                cat("@expectedCost:          ", sprintf("%.2f", object@expectedCost), "\n")
            } else {
                cat("@expectedCost:           NO DATA\n")
            }  
            cat("@lookupTable:")  
            if (dim(object@lookupTable)[1] > 0){
                cat("\n")
                show(object@lookupTable)
            } else {        
                cat("            NO DATA\n")
            }    
        }
)

################################################# SUMMARY:
##########################################################

setMethod("summary", signature(object="IndSampling"),
        function(object, output = c("proportion", "percent")){
            output <- match.arg(output)
            displayLimit <- 20
            outVec <- c("INDIVIDUAL SAMPLING:", "",
                    as.vector(summary(object@surveyData, output)), "",
                    "Sampling strategy:",
                    "------------------")
            if (length(object@herdSensitivity) > 0){
                if (output == "proportion"){            
                    outChar <- sprintf("%.3f",object@herdSensitivity)
                } else {
                    outChar <- paste(sprintf("%.2f",object@herdSensitivity*100), "%")
                }           
                outVec <- c(outVec, paste("Herd sensitivity:                         ", 
                                outChar))
            } else {
                outVec <- c(outVec,"No herd sensitivity specified.")
            }
            if (length(object@nHerds) > 0){
                outVec <- c(outVec, paste("Number of herds to test:                  ", 
                                sprintf("%d",object@nHerds)))
            } else {
                outVec <- c(outVec, "No data on number of herds to test.")
            }
            if (length(object@nHerdsPerRiskGroup) > 0){
                tempVar <- capture.output(print(object@nHerdsPerRiskGroup))
                tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                outVec <- c(outVec, paste("Number of herds to test per risk group:    ",
                                tempVar[1], sep = "")) 
                tempVar <- tempVar[-1]
                tempVar <- paste("                                          ", tempVar)
                outVec <- c(outVec, tempVar)  
            }
            if (length(object@nAnimalsMean) > 0){
                outVec <- c(outVec, paste("Expected total number of animals to test: ", 
                                sprintf("%.2f",object@nAnimalsMean)))
            } else {
                outVec <- c(outVec, 
                        "No data on expected total number of animals to test.")
            }
            if (length(object@expectedCost) > 0){
                outVec <- c(outVec, paste("Expected total costs of the survey:       ", 
                                sprintf("%.2f",object@expectedCost)))
            } else {
                outVec <- c(outVec,"No data on expected total costs of the survey.")
            } 
            if (dim(object@lookupTable)[1] > 0){            
                ## Format lookup Table for pretty summary:
                ## Index for herd sizes where entire herd is tested:
                if (any(object@lookupTable[,"N_upper"] == object@lookupTable[,"sampleSize"])){
				    indAll <- range(which(object@lookupTable[,"N_upper"] == object@lookupTable[,"sampleSize"]))
                    nRows <- dim(object@lookupTable)[1]
                    ## Herd sizes where entire herd is tested:
                    tempDf <- data.frame(herdSize = paste(object@lookupTable[indAll[1],"N_lower"], "-", 
                        object@lookupTable[indAll[2],"N_upper"]), sampleSize = "entire herd",
                        stringsAsFactors = FALSE)
			    } else {
					tempDf <- data.frame(herdSize = character(), sampleSize = character())
					nRows <- dim(object@lookupTable)[1]
					indAll <- c(0,0)
				}
                ## Rest of the herd sizes:
                if (nRows > indAll[2]){
                    tempDf2 <- data.frame(herdSize = sapply((indAll[2]+1):nRows, function (ii){
                                        N_lower <- object@lookupTable[ii,"N_lower"]
                                        N_upper <- object@lookupTable[ii,"N_upper"]
                                        if (N_lower == N_upper){
                                            out <- as.character(N_upper)
                                        } else {
                                            out <- paste(N_lower, "-", N_upper)
                                        }
                                        return(out)
                                    }), sampleSize = object@lookupTable[(indAll[2]+1):nRows,"sampleSize"],
                            stringsAsFactors = FALSE)
                    ## Paste the two data frames together:
                    tempDf <- rbind(tempDf,tempDf2)
                }
                ## Format strings for displaying:
                herdSizeString <- format(c("Herd size", tempDf$herdSize), justify = "centre")
                sampleSizeString <- format(c("No. of animals to test", tempDf$sampleSize), justify = "centre")
                outString <- paste(herdSizeString, sampleSizeString, sep = "  |  ")
                outString <- c(paste(" ",outString[1], sep = ""), 
                        Reduce(function(x,y) paste(x,y,sep = ""),rep("-",nchar(outString[1]))),
                        outString[-1])
                outString <- paste(rep(" ", length(outString)), outString)
                outVec <- c(outVec, "Lookup table for the number of animals to test per herd:",
                        "", outString)
            } else {
                outVec <- c(outVec,"No lookup table specified.")
            } 
            outTab <- as.table(matrix(outVec,ncol = 1))
            dimnames(outTab) <- list(rep("", length(outVec)), "")
            outTab
        }
)

#################################################### HTML:
##########################################################

setMethod("HTML", signature(x = "IndSampling"),
        function(x, filename = "IndSampling", outdir = getwd(), CSSFile = "ffd.css", 
                Title = "Individual Sampling", append = TRUE,...){ 
            
            ## Create html-file:
            target <- HTMLInitFile(outdir = outdir , filename = filename, 
                    CSSFile = CSSFile, Title = Title,...)
            
            ## Write content:
            ###########################################################################
            ###########################################################################
            cat("\n\n<h1>Survey to substantiate Freedom From Disease</h1>\n", file = target,
                    append = TRUE)
            cat("<h2>Individual Sampling</h2>\n<br>\n", file = target, append = TRUE)
            
            ## Survey Parameters:
            cat("<h3>Survey Parameters:</h3>\n", file = target, append = TRUE)  
            surveyDat <- x@surveyData      
            surveyDat <- x@surveyData      
            theTable <- matrix(c(as.character(surveyDat@designPrevalence), 
                            as.character(surveyDat@alpha), 
                            as.character(surveyDat@intraHerdPrevalence), 
                            as.character(surveyDat@diagSensitivity), 
                            as.character(surveyDat@costHerd), 
                            as.character(surveyDat@costAnimal)), ncol = 1)
            indVec <- c(length(surveyDat@designPrevalence)>0, 
                    length(surveyDat@alpha)>0, 
                    length(surveyDat@intraHerdPrevalence)>0,
                    length(surveyDat@diagSensitivity)>0,
                    length(surveyDat@costHerd)>0,
                    length(surveyDat@costAnimal)>0)
            rownames(theTable) <- c("Design Prevalence ", "Overall Significance (alpha) ",
                    "Intra herd prevalence ", "Sensitivity of diagnostic test ",
                    "Cost per herd ", "Cost per animal ")[indVec]
            HTML(theTable, align = "left")
            cat("<br>\n\n", file = target, append = TRUE)
            
            ## Data description:
            if (length(surveyDat@nAnimalVec) > 0){
                cat("<h3>Data Description:</h3>\n", file = target, append = TRUE)  
                theTable <- matrix(c(length(surveyDat@nAnimalVec), 
                                sum(surveyDat@nAnimalVec), 
                                min(surveyDat@nAnimalVec),
                                median(surveyDat@nAnimalVec),
                                max(surveyDat@nAnimalVec)), ncol = 1)
                rownames(theTable) <- c("Number of herds ", 
                        "Total number of animals ",
                        "Minimal herd size ",
                        "Median herd size ",
                        "Maximal herd size ")
                HTML(theTable, align = "left")
                cat("<br>\n\n", file = target, append = TRUE)
            }
            
            ## Sample Parameters:
            cat("<h3>Sampling strategy:</h3>\n", file = target, append = TRUE)  
            theTable <- matrix(c(as.character(round(x@herdSensitivity,2)),
                            as.character(x@nHerds),
                            as.character(round(x@nAnimalsMean,2)), 
                            as.character(round(x@expectedCost,2))), ncol = 1)
            indVec <- c(length(x@herdSensitivity)>0, 
                    length(x@nHerds)>0,
                    length(x@nAnimalsMean)>0,
                    length(x@expectedCost)>0)
            rownames(theTable) <- c("Herd sensitivity ",
                    "Number of herds to test ", 
                    "Expected total number of animals to test ",
                    "Expected total costs of the survey ")[indVec]
            HTML(theTable, align = "left")
            if (dim(x@lookupTable)[1] > 0){
                HTML("<br>")       
                cat("<br>\n\n<h4>Number of animals to test per herd:</h4>\n", file = target, 
                        append = TRUE)        
                ## Format lookup Table for pretty summary:
                ## Index for herd sizes where entire herd is tested:
                indAll <- range(which(x@lookupTable[,"N_upper"] == x@lookupTable[,"sampleSize"]))
                nRows <- dim(x@lookupTable)[1]
                ## Herd sizes where entire herd is tested:
                tempDf <- data.frame(herdSize = paste(x@lookupTable[indAll[1],"N_lower"], "-", 
                                x@lookupTable[indAll[2],"N_upper"]), sampleSize = "entire herd",
                        stringsAsFactors = FALSE)
                ## Rest of the herd sizes:
                if (nRows > indAll[2]){
                    tempDf2 <- data.frame(herdSize = sapply((indAll[2]+1):nRows, function (ii){
                                        N_lower <- x@lookupTable[ii,"N_lower"]
                                        N_upper <- x@lookupTable[ii,"N_upper"]
                                        if (N_lower == N_upper){
                                            out <- as.character(N_upper)
                                        } else {
                                            out <- paste(N_lower, "-", N_upper)
                                        }
                                        return(out)
                                    }), sampleSize = x@lookupTable[(indAll[2]+1):nRows,"sampleSize"],
                            stringsAsFactors = FALSE)
                    ## Paste the two data frames together:
                    tempDf <- rbind(tempDf,tempDf2)            
                }
                names(tempDf) <- c("Herd size"," No. of animals to test")
                HTML(tempDf, align = "left", row.names = FALSE)
            }
            
            ## Close body and html tag
            cat("\n\n<br>\n<hr size=1>\n<font size=-1>\nGenerated on: <i>", format(Sys.time(), "%b %d %X %Y"), 
                    "</i> - <b>FFD</b>\n<hr size=1>\n\n", file = target, append = TRUE)
            
            cat("\n</body>\n</html>\n", file = target, append = TRUE)        
            
            ## Check if css.file exists:
            css.file <- file.path(outdir,CSSFile, fsep = .Platform$file.sep)
            ## If file does not exist create one:
            createStyleFile(css.file)                    
            
            ## Return value:
            invisible(x)         
        })

################################################## SAMPLE:
##########################################################

setMethod("sample", signature(x = "IndSampling"),
        function(x, size = c("fixed", "dynamic"), replace, prob){
#    function(x, size, replace, prob){
#        ## Set default value for flag 'size':
            if (missing(size)) size <- "fixed"
            size <- match.arg(size)            
			
            ## Sample with fixed sample size:
            if (size == "fixed"){ 
                ## No risk groups:
                ##################
                if (length(x@nHerdsPerRiskGroup) == 0){ 
                    ## Sample x@nHerds herds using simple random sampling:
                    indexSample <- sort(sample(x = seq(along = x@surveyData@nAnimalVec),
                                    size = x@nHerds, replace = FALSE))
                    ## Compute the a-posteriori alpha error:
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec[indexSample], 
                            method = "individual", 
                            herdSensitivity = x@herdSensitivity, 
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, 
                            diagSpecificity = 1)            
                    aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaErrorVector, 
                            nPopulation = length(x@surveyData@nAnimalVec), 
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")  
                    nSampleArgVec <- numeric()
                } else {
                    ## Sampling with risk groups:
                    #############################
                    indexVec <- seq(along = x@surveyData@nAnimalVec)
                    indexList <- split(x = indexVec, f = x@surveyData@riskGroupVec)
                    namesRiskGroupVec <- names(x@nHerdsPerRiskGroup)
                    nSampleRiskGroupVec <- as.vector(x@nHerdsPerRiskGroup)
                    indexSampleList <- lapply(seq(along = indexList), function(ii){
                        riskGroup <- names(indexList)[[ii]]     
                        nSample <- nSampleRiskGroupVec[namesRiskGroupVec == riskGroup]
                        if (length(indexList[[ii]]) > 1){
                            indexOutVec <- sample(x = indexList[[ii]],
                                size = nSample, replace = FALSE)    
                        } else {
                            if (nSample == 1){
                                indexOutVec <- indexList[[ii]]
                            } else {
                                stop ("Sample size larger than Population size.")   
                            }                       
                        }               
                    })
                    ## Concatenate the sample indices for the risk groups:
                    indexSample <- sort(Reduce(function(x,y) c(x,y), indexSampleList),
                            decreasing = FALSE)
                    ## Compute the a-posteriori alpha error:
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec[indexSample], 
                            method = "individual", 
                            herdSensitivity = x@herdSensitivity,  
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, 
                            diagSpecificity = 1)            
                    ## Determine the number of herds in each risk group:
                    riskValueDf <- x@surveyData@riskValueData[,1:2]
                    names(riskValueDf) <- c("riskGroup", "riskValues")
                    riskValueDf$riskGroup <- as.character(riskValueDf$riskGroup)
                    riskValueDf$id <- seq(along = riskValueDf[,1])
                    riskGroupTab <- table(x@surveyData@riskGroupVec)
                    riskGroupDf <- data.frame(riskGroup = as.character(names(riskGroupTab)), 
                            nPopulation = as.vector(riskGroupTab))
                    riskValueDf <- merge(x = riskValueDf, y = riskGroupDf, by = "riskGroup",
                            sort = FALSE)   
                    riskValueDf <- riskValueDf[order(riskValueDf$id, decreasing = FALSE),]
                    aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                            groupVec = x@surveyData@riskGroupVec[indexSample], 
                            groupLevels = riskValueDf$riskGroup,
                            nPopulationVec = riskValueDf$nPopulation,
                            nRelRiskVec = riskValueDf$riskValues,
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")  
					nSampleArgVec <- nSampleRiskGroupVec
                }
            } else {
            ## Dynamic sampling:
            ####################
                ## No risk groups:
                ##################
                if (length(x@nHerdsPerRiskGroup) == 0){
                    ## Data frame with herd-based alpha-errors for each herd
                    ## (=1-herd sensitivity):
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec, 
                            method = "individual", herdSensitivity = x@herdSensitivity, 
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, diagSpecificity = 1)
                    alphaDataFrame <- data.frame(size = x@surveyData@nAnimalVec, 
                            alpha = alphaErrorVector, id = seq(along = x@surveyData@nAnimalVec))
                    ## Permutate data frame:
                    alphaDataFrame <- alphaDataFrame[sample(x = alphaDataFrame$id, 
                                    size = length(alphaDataFrame$id), replace = FALSE),]
                    ## Dynamic sampling: find the smallest integer k such that the 
                    ## sample containing the first k rows of the data frame have
                    ## an a-posteriori error <=  x@surveyData@alpha    
                    k <- x@nHerds  
                    aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                            nPopulation = length(x@surveyData@nAnimalVec), 
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")
                    if (aPostAlpha > x@surveyData@alpha){
                        while (aPostAlpha > x@surveyData@alpha){
                            k <- k + 1
                            aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                                    nPopulation = length(x@surveyData@nAnimalVec), 
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")               
                        }            
                    } else {
                        while (aPostAlpha <= x@surveyData@alpha){
                            k <- k - 1
                            aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                                    nPopulation = length(x@surveyData@nAnimalVec), 
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")               
                        }         
                        k <- k + 1
                        aPostAlpha <- computeAposterioriError(alphaErrorVector = alphaDataFrame$alpha[1:k], 
                                nPopulation = length(x@surveyData@nAnimalVec), 
                                nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                method = "approx")                               
                    }
                    indexSample <- sort(alphaDataFrame$id[1:k]) 
                    nSampleArgVec <- numeric()
                } else {
                    ## Sampling with risk groups:
                    #############################
                    ## Data frame with herd-based alpha-errors for each herd
                    ## (=1-herd sensitivity):
                    alphaErrorVector <- computeAlpha(nAnimalVec = x@surveyData@nAnimalVec, 
                            method = "individual", herdSensitivity = x@herdSensitivity, 
                            intraHerdPrevalence = x@surveyData@intraHerdPrevalence, 
                            diagSensitivity = x@surveyData@diagSensitivity, diagSpecificity = 1)
                    alphaDataFrame <- data.frame(size = x@surveyData@nAnimalVec, 
                            alpha = alphaErrorVector, id = seq(along = x@surveyData@nAnimalVec),
                            riskGroupVec = as.character(x@surveyData@riskGroupVec))
                    ## Permutate data frame:
                    alphaDataFrame <- alphaDataFrame[sample(x = alphaDataFrame$id, 
                            size = length(alphaDataFrame$id), replace = FALSE),]    
                    alphaList <- split(x = alphaDataFrame, f = alphaDataFrame$riskGroupVec)
                    ## Dynamic sampling: find the smallest integer k such that the 
                    ## sample containing the first k rows of the data frame have
                    ## an a-posteriori error <=  x@surveyData@alpha    
                    k <- x@nHerds   
                    nSampleFixVec <- x@nSampleFixVec
                    names(nSampleFixVec) <- x@surveyData@riskValueData$riskGroup                
                    nSamplePropVec <- as.vector(x@probVec*
                        table(x@surveyData@riskGroupVec)[x@surveyData@riskValueData$riskGroup]
                        [is.na(nSampleFixVec)])
                    ## Proportion for sample sizes:
                    if (sum(is.na(nSampleFixVec)) > 1){
                        nSamplePropVec <- nSamplePropVec/sum(nSamplePropVec)
                    } else {
                        nSamplePropVec <- 1
                    }
                    ## Population size of the risk groups:
                    nPopulationVec <- table(x@surveyData@riskGroupVec)[x@surveyData@riskValueData$riskGroup]
                    ## Names of the risk groups:
                    groupLevels <- as.character(x@surveyData@riskValueData$riskGroup)
                    
                    ## Compute baseline alpha error:
                    ################################
                    ## Sample size for the risk groups:
                    nSampleArgVec <- nSampleFixVec
                    nSampleArgVec[is.na(nSampleFixVec)] <- 
                            roundConstantSum(k*nSamplePropVec, output = 0)
                    alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                        riskGroup <- names(alphaList)[[ii]]     
                        nSample <- nSampleArgVec[riskGroup] 
                        if (nSample > 0){
                            alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                        } else {
                            alphaOut <- numeric()
                        }                           
                    }) 
                    alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                    groupVec <- names(alphaList)    
                    groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                    aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                            groupVec = groupVec, 
                            groupLevels = groupLevels,
                            nPopulationVec = nPopulationVec[groupLevels],
                            nRelRiskVec = x@surveyData@riskValueData$riskValues,
                            nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                            method = "approx")
                    ## Modify sample (add/remove farms):
                    ####################################
                    if (aPostAlpha > x@surveyData@alpha){
                        while (aPostAlpha > x@surveyData@alpha){
                            k <- k + 1
                            nSampleArgVec <- nSampleFixVec
                            nSampleArgVec[is.na(nSampleFixVec)] <- 
                                    roundConstantSum(k*nSamplePropVec, output = 0)
                            alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                                        riskGroup <- names(alphaList)[[ii]]     
                                        nSample <- nSampleArgVec[riskGroup] 
                                        if (nSample > 0){
                                            alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                                        } else {
                                            alphaOut <- numeric()
                                        }                           
                                    }) 
                            alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                            groupVec <- names(alphaList)    
                            groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                            aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                                    groupVec = groupVec, 
                                    groupLevels = groupLevels,
                                    nPopulationVec = nPopulationVec[groupLevels],
                                    nRelRiskVec = x@surveyData@riskValueData$riskValues,
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")               
                        }            
                    } else {
                        while (aPostAlpha <= x@surveyData@alpha){
                            k <- k - 1
                            nSampleArgVec <- nSampleFixVec
                            nSampleArgVec[is.na(nSampleFixVec)] <- 
                                    roundConstantSum(k*nSamplePropVec, output = 0)
                            alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                                        riskGroup <- names(alphaList)[[ii]]     
                                        nSample <- nSampleArgVec[riskGroup] 
                                        if (nSample > 0){
                                            alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                                        } else {
                                            alphaOut <- numeric()
                                        }                           
                                    }) 
                            alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                            groupVec <- names(alphaList)    
                            groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                            aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                                    groupVec = groupVec, 
                                    groupLevels = groupLevels,
                                    nPopulationVec = nPopulationVec[groupLevels],
                                    nRelRiskVec = x@surveyData@riskValueData$riskValues,
                                    nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                    method = "approx")              
                        }         
                        k <- k + 1
                        nSampleArgVec <- nSampleFixVec
                        nSampleArgVec[is.na(nSampleFixVec)] <- 
                                roundConstantSum(k*nSamplePropVec, output = 0)
                        alphaErrorList <- lapply(seq(along = alphaList), function(ii){
                                    riskGroup <- names(alphaList)[[ii]]     
                                    nSample <- nSampleArgVec[riskGroup] 
                                    if (nSample > 0){
                                        alphaOut <- alphaList[[ii]]$alpha[1:nSample]
                                    } else {
                                        alphaOut <- numeric()
                                    }                           
                                }) 
                        alphaErrorVector <- Reduce(function(x,y) c(x,y), alphaErrorList)
                        groupVec <- names(alphaList)    
                        groupVec <- rep(groupVec, nSampleArgVec[groupVec])              
                        aPostAlpha <- computeAposterioriErrorRiskGroups(alphaErrorVector = alphaErrorVector, 
                                groupVec = groupVec, 
                                groupLevels = groupLevels,
                                nPopulationVec = nPopulationVec[groupLevels],
                                nRelRiskVec = x@surveyData@riskValueData$riskValues,
                                nDiseased = max(round(length(x@surveyData@nAnimalVec)*x@surveyData@designPrevalence),1), 
                                method = "approx")                              
                    }
                    indexSampleList <- lapply(seq(along = alphaList), function(ii){
                        riskGroup <- names(alphaList)[[ii]]     
                        nSample <- nSampleArgVec[riskGroup] 
                        if (nSample > 0){
                            indexOut <- alphaList[[ii]]$id[1:nSample]
                        } else{
                            indexOut <- numeric()
                        }
                        return(indexOut)
                    })
                    indexSample <- sort(Reduce(function(x,y) c(x,y), 
                        indexSampleList), decreasing = FALSE)   
                }
            }
            
            #######################################################################
            ## Return value:
            out <- list(indexSample = indexSample, aPostAlpha = aPostAlpha)         
            ## Add population data of sample:    
            if (dim(x@surveyData@populationData)[1] > 0){        
                out$sample <- x@surveyData@populationData[indexSample,]
            }
            ## Add sample size of the risk groups:
            if (length(nSampleArgVec) > 0){
                out$nSamplePerRiskGroup <- nSampleArgVec                
            }           
            return(out)        
        }
)


###############################################################################
###############################################################################
## Class "IndSamplingSummary"
##
## Contains the parameters and the data for a survey
## to substantiate freedom from disease 
## using "individual sampling". Additionally to the survey
## parameters (design prevalence, overall significance,
## intra-herd prevalence, sensitivity of the diagnostic test, 
## cost per tested animal and cost per tested herd) the object 
## contains the number of herds to 
## be tested, the mean overall number of animals to be tested
## and the expected costs for a range of possible herd
## sensitivities.
##
## Package: FFD
##
## Ian Kopacka
## 2010-07-15
###############################################################################
###############################################################################
setClass(
        Class = "IndSamplingSummary", 
        representation = representation(
                surveyData = "SurveyData",
                herdSensVec = "numeric",
                nHerdsVec = "numeric",
                nHerdsPerRiskGroupMx = "matrix",
                nSampleFixVec = "numeric",
                probVec = "numeric",
                nAnimalsMeanVec = "numeric",
                expectedCostVec = "numeric"     
        ),
        validity = function(object)
        {      
            if (any((object@herdSensVec <= 0)|(object@herdSensVec >= 1))){
                stop ("[IndSamplingSummary validation]: Slot 'herdSensVec' must contain values in (0,1).")
            }
            if (any(object@nHerdsVec <= 0)){
                stop ("[IndSamplingSummary validation]: Slot 'nHerdsVec' must contain positive values.")
            }
            if (any((object@nHerdsVec - as.integer(object@nHerdsVec)) != 0)){
                stop ("[IndSamplingSummary validation]: Slot 'nHerdsVec' must contain an integer vector.")
            }
            if (any(object@nAnimalsMeanVec <= 0)){
                stop ("[IndSamplingSummary validation]: Slot 'nAnimalsMeanVec' must contain positive values")
            }     
            if (length(object@herdSensVec) != length(object@nHerdsVec)){
                stop ("[IndSamplingSummary validation]: Slots 'herdSensVec' and 'nHerdsVec' must have the same length.")         
            } 
            if (length(object@herdSensVec) != length(object@nAnimalsMeanVec)){
                stop ("[IndSamplingSummary validation]: Slots 'herdSensVec' and 'nAnimalsMeanVec' must have the same length.")         
            } 
            if (dim(object@nHerdsPerRiskGroupMx)[1] > 0){
                if (dim(object@nHerdsPerRiskGroupMx)[2] != dim(object@surveyData@riskValueData)[1]){
                    stop ("[IndSamplingSummary validation]: Slot 'nHerdsPerRiskGroupMx' have the same number of columns as the number of rows in 'surveyData@riskValueData'.")              
                }
                if (!all(colnames(object@nHerdsPerRiskGroupMx) %in% object@surveyData@riskValueData[,1])){
                    stop ("[IndSamplingSummary validation]: The names of slot 'nHerdsPerRiskGroupMx' must be contained in the first column of 'surveyData@riskValueData'.")              
                } 
                if (dim(object@nHerdsPerRiskGroupMx)[1] != length(object@herdSensVec)){
                    stop ("[IndSamplingSummary validation]: The number of rows of slot 'nHerdsPerRiskGroupMx' must equal the length of 'herdSensVec'.")
                }
            }
            return(TRUE)
        }    
)

## METHODS:
###############################################################################
###############################################################################


#################################################### SHOW:
##########################################################
setMethod("show", signature(object="IndSamplingSummary"),
        function(object){
            displayLimit <- 20
            cat("Object of class 'IndSamplingSummary':\n")
            cat("Slots:\n")  
            cat("@surveyData:\n")
            cat("----------------------------\n")
            show(object@surveyData)
            cat("----------------------------\n")
            cat("@herdSensVec:")
            if (length(object@herdSensVec) > 0){
                cat("\n")
                if (length(object@herdSensVec) > displayLimit){
                    show(object@herdSensVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@herdSensVec), "displayed)\n")
                } else {
                    show(object@herdSensVec)
                }
            } else {        
                cat("        NO DATA\n")
            }
            cat("@nHerdsVec:")
            if (length(object@nHerdsVec) > 0){
                cat("\n")
                if (length(object@nHerdsVec) > displayLimit){
                    show(object@nHerdsVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@nHerdsVec), "displayed)\n")
                } else {
                    show(object@nHerdsVec)
                }
            } else {        
                cat("              NO DATA\n")
            }
            cat("@nHerdsPerRiskGroupMx:")
            if (dim(object@nHerdsPerRiskGroupMx)[1] > 0){
                cat("\n")
                if (dim(object@nHerdsPerRiskGroupMx)[1] > displayLimit){
                    print(as.data.frame(object@nHerdsPerRiskGroupMx[1:displayLimit,]), row.names = FALSE)
                    cat("   (only first",displayLimit,"rows of",dim(object@nHerdsPerRiskGroupMx)[1], "displayed)\n")
                } else {
                    print(as.data.frame(object@nHerdsPerRiskGroupMx), row.names = FALSE)
                }
            } else {        
                cat("   NO DATA\n")
            }
            cat("@nSampleFixVec:")
            if (length(object@nSampleFixVec) > 0){
                cat("\n")
                if (length(object@nSampleFixVec) > displayLimit){
                    show(object@nSampleFixVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@nSampleFixVec), "displayed)\n")
                } else {
                    show(object@nSampleFixVec)
                }
            } else {        
                cat("          NO DATA\n")
            }
            cat("@probVec:")
            if (length(object@probVec) > 0){
                cat("\n")
                if (length(object@probVec) > displayLimit){
                    show(object@probVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@probVec), "displayed)\n")
                } else {
                    show(object@probVec)
                }
            } else {        
                cat("                NO DATA\n")
            }       
            cat("@nAnimalsMeanVec:")
            if (length(object@nAnimalsMeanVec) > 0){
                cat("\n")
                if (length(object@nAnimalsMeanVec) > displayLimit){
                    show(object@nAnimalsMeanVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@nAnimalsMeanVec), "displayed)\n")
                } else {
                    show(object@nAnimalsMeanVec)
                }
            } else {        
                cat("        NO DATA\n")
            }
            cat("@expectedCostVec:")
            if (length(object@expectedCostVec) > 0){
                cat("\n")
                if (length(object@expectedCostVec) > displayLimit){
                    show(object@expectedCostVec[1:displayLimit])
                    cat("   (only first",displayLimit,"elements of",length(object@expectedCostVec), "displayed)\n")
                } else {
                    show(object@expectedCostVec)
                }
            } else {        
                cat("        NO DATA\n")
            }                   
        }
)

################################################# SUMMARY:
##########################################################

#setMethod("summary", signature(object="IndSamplingSummary"),
#    function(object){
#        displayLimit <- 20
#        cat("INDIVIDUAL SAMPLING:\n\n")
#        summary(object@surveyData)
#        cat("\n")
#        cat("Cost optimal sampling strategy:\n")
#        cat("-------------------------------\n")
#        if(length(object@expectedCostVec) > 0){
#            indCostOpt <- which.min(object@expectedCostVec)
#            cat("Herd sensitivity:                         ", sprintf("%.2f",object@herdSensVec[indCostOpt]), "\n")
#            cat("Number of herds to test:                  ", sprintf("%d",object@nHerdsVec[indCostOpt]), "\n")
#            cat("Expected total number of animals to test: ", sprintf("%.2f",object@nAnimalsMeanVec[indCostOpt]), "\n")
#            cat("Expected total costs of the survey:       ", sprintf("%.2f",object@expectedCostVec[indCostOpt]), "\n")            
#        } else {
#            cat("No data regarding costs.\n")
#        }        
#    }
#)
setMethod("summary", signature(object="IndSamplingSummary"),
        function(object, output = c("proportion", "percent")){
            output <- match.arg(output)
            displayLimit <- 20
            outVec <- c("INDIVIDUAL SAMPLING DIAGNOSTICS:", "", 
                    as.vector(summary(object@surveyData, output)), "",
                    "Cost optimal sampling strategy:",
                    "-------------------------------")
            if(length(object@expectedCostVec) > 0){
                indCostOpt <- which.min(object@expectedCostVec)
                if (output == "proportion"){
                    herdSensChar <- sprintf("%.3f",object@herdSensVec[indCostOpt])              
                } else {
                    herdSensChar <- paste(sprintf("%.2f",object@herdSensVec[indCostOpt]*100),
                            "%")
                }
                outVec <- c(outVec,
                        paste("Herd sensitivity:                         ", 
                                herdSensChar),
                        paste("Number of herds to test:                  ", 
                                sprintf("%d",object@nHerdsVec[indCostOpt])))
                if (dim(object@nHerdsPerRiskGroupMx)[1] > 0){
                    tempVar <- capture.output(print(object@nHerdsPerRiskGroupMx[indCostOpt,]))
                    tempVar <- gsub(pattern = "\t", replacement = "  ", x = tempVar)
                    outVec <- c(outVec, paste("Number of herds to test per risk group:    ",
                                    tempVar[1], sep = "")) 
                    tempVar <- tempVar[-1]
                    tempVar <- paste("                                          ", tempVar)
                    outVec <- c(outVec, tempVar)  
                }
                
                outVec <- c(outVec,
                        paste("Expected total number of animals to test: ", 
                                sprintf("%.2f",object@nAnimalsMeanVec[indCostOpt])),
                        paste("Expected total costs of the survey:       ", 
                                sprintf("%.2f",object@expectedCostVec[indCostOpt])))            
            } else {
                outVec <- c(outVec,"No data regarding costs.")
            } 
            outTab <- as.table(matrix(outVec,ncol = 1))
            dimnames(outTab) <- list(rep("", length(outVec)), "")
            outTab
        }
)

#################################################### HTML:
##########################################################

setMethod("HTML", signature(x = "IndSamplingSummary"),
        function(x, filename = "IndSamplingSummary", outdir = getwd(), CSSFile = "ffd.css", 
                Title = "Individual Sampling Diagnostics", append = TRUE,...){ 
            
            ## Create html-file:
            target <- HTMLInitFile(outdir = outdir , filename = filename, 
                    CSSFile = CSSFile, Title = Title,...)
            
            ## Write content:
            ###########################################################################
            ###########################################################################
            cat("\n\n<h1>Survey to substantiate Freedom From Disease</h1>\n", file = target,
                    append = TRUE)
            cat("<h2>Individual Sampling</h2>\n<br>\n", file = target, append = TRUE)
            
            ## Survey Parameters:
            cat("<h3>Survey Parameters:</h3>\n", file = target, append = TRUE)  
            surveyDat <- x@surveyData      
            theTable <- matrix(c(as.character(surveyDat@designPrevalence), 
                            as.character(surveyDat@alpha), 
                            as.character(surveyDat@intraHerdPrevalence), 
                            as.character(surveyDat@diagSensitivity), 
                            as.character(surveyDat@costHerd), 
                            as.character(surveyDat@costAnimal)), ncol = 1)
            indVec <- c(length(surveyDat@designPrevalence)>0, 
                    length(surveyDat@alpha)>0, 
                    length(surveyDat@intraHerdPrevalence)>0,
                    length(surveyDat@diagSensitivity)>0,
                    length(surveyDat@costHerd)>0,
                    length(surveyDat@costAnimal)>0)
            rownames(theTable) <- c("Design Prevalence ", "Overall Significance (alpha) ",
                    "Intra herd prevalence ", "Sensitivity of diagnostic test ",
                    "Cost per herd ", "Cost per animal ")[indVec]
            HTML(theTable, align = "left")
            cat("<br>\n\n", file = target, append = TRUE)
            
            ## Data description:
            if (length(surveyDat@nAnimalVec) > 0){
                cat("<h3>Data Description:</h3>\n", file = target, append = TRUE)  
                theTable <- matrix(c(length(surveyDat@nAnimalVec), 
                                sum(surveyDat@nAnimalVec), 
                                min(surveyDat@nAnimalVec),
                                median(surveyDat@nAnimalVec),
                                max(surveyDat@nAnimalVec)), ncol = 1)
                rownames(theTable) <- c("Number of herds ", 
                        "Total number of animals ",
                        "Minimal herd size ",
                        "Median herd size ",
                        "Maximal herd size ")
                HTML(theTable, align = "left")
                cat("<br>\n\n", file = target, append = TRUE)
            }
            
            ## Cost optimal sample Parameters:
            cat("<h3>Cost optimal sampling strategy:</h3>\n", file = target, append = TRUE)
            if(length(x@expectedCostVec) > 0){
                indCostOpt <- which.min(x@expectedCostVec)
                theTable <- matrix(c(as.character(x@herdSensVec[indCostOpt]),
                                as.character(x@nHerdsVec[indCostOpt]),
                                as.character(round(x@nAnimalsMeanVec[indCostOpt],2)), 
                                as.character(round(x@expectedCostVec[indCostOpt],2))), ncol = 1)
                indVec <- c(length(x@herdSensVec)>0, 
                        length(x@nHerdsVec)>0,
                        length(x@nAnimalsMeanVec)>0,
                        length(x@expectedCostVec)>0)
                rownames(theTable) <- c("Herd sensitivity ", 
                        "Number of herds to test ", 
                        "Expected total number of animals to test ",
                        "Expected total costs of the survey ")[indVec]            
                HTML(theTable, align = "left")                    
            } else {
                HTML("No data regarding costs.")
            }
            cat("<br>\n\n", file = target, append = TRUE)
            
            ## Diagnostic plots:
            cat("<h3>Diagnostic plots:</h3>\n", file = target, append = TRUE)
            survey.Data <- x@surveyData
            
            if(length(survey.Data@nAnimalVec) > 0){        
                ## Mean number of animals per herd:
                #cat("<h4>Mean number of animals to test per herd</h4>\n", file = target, append = TRUE)
                par(mfrow = c(1, 1))
                plot(x@herdSensVec, x@nAnimalsMeanVec/x@nHerdsVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "Mean no. of animals per herd to be tested")
                HTMLplot(Width = 500, Height = 500, Caption = "Number of animals to test per herd",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"1", sep = "_"))
                ## Number of herds to be tested:
                #cat("<h4>Sample size</h4>\n", file = target, append = TRUE)
                plot(x@herdSensVec, x@nHerdsVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "No. of herds to be tested")
                HTMLplot(Width = 500, Height = 500, Caption = "Sample size",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"2", sep = "_"))
                ## Total number of animals to be tested:
                #cat("<h4>Number of animals to test</h4>\n", file = target, append = TRUE)
                plot(x@herdSensVec, x@nAnimalsMeanVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "Expected total no. of animals to be tested")
                HTMLplot(Width = 500, Height = 500, Caption = "Total number of animals to test",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"3", sep = "_"))
                ## Expected cost:
                #cat("<h4>Expected costs</h4>\n", file = target, append = TRUE)
                plot(x@herdSensVec, x@expectedCostVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "Expected cost")  
                HTMLplot(Width = 500, Height = 500, Caption = "Expected costs of the survey",
                        GraphFileName = paste("GRAPH",format(Sys.time(), "%b%d_%H%M%S"),"4", sep = "_"))                         
            } else {
                HTML("No data to produce plots.")
            }        
            
            ## Close body and html tag
            cat("\n\n<br>\n<hr size=1>\n<font size=-1>\nGenerated on: <i>", format(Sys.time(), "%b %d %X %Y"), 
                    "</i> - <b>FFD</b>\n<hr size=1>\n\n", file = target, append = TRUE)
            
            cat("\n</body>\n</html>\n", file = target, append = TRUE)        
            
            ## Check if css.file exists:
            css.file <- file.path(outdir,CSSFile, fsep = .Platform$file.sep)
            ## If file does not exist create one:
            createStyleFile(css.file)                    
            
            ## Return value:
            invisible(x) 
        })

#################################################### PLOT:
##########################################################

setMethod("plot", signature(x = "IndSamplingSummary"),
        function(x,y,...){
            survey.Data <- x@surveyData
            if(length(survey.Data@nAnimalVec > 0)){        
                par(mfrow = c(2, 2))      
                # ## Distribution of the herd sizes:
#            hist(survey.Data@nAnimalVec, xlab = "Herd size", ylab = "Frequency",
#                main = "", col = "#DBA400")               
                ## Mean number of animals per herd:
                plot(x@herdSensVec, x@nAnimalsMeanVec/x@nHerdsVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "Mean no. of animals per herd to be tested")
                ## Number of herds to be tested:
                plot(x@herdSensVec, x@nHerdsVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "No. of herds to be tested")
                ## Total number of animals to be tested:
                plot(x@herdSensVec, x@nAnimalsMeanVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "Expected total no. of animals to be tested")
                ## Expected cost:
                plot(x@herdSensVec, x@expectedCostVec, type = "l",
                        xlab = "Herd sensitivity", ylab = "Expected cost")           
                ## Titel:
                par(oma = c(2,1,3,1)) 
                title("Analysis individual sampling", outer = TRUE)                  
            } else {
                cat("Object of class 'IndSamplingSummary' contains no data.\n")
            }
        }
)
