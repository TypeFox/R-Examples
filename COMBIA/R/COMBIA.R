#' @import hash
#' @import oro.nifti
#' @import lattice
#' @import gdata
#' @import latticeExtra
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom stats optimize quantile sd update
#' @importFrom utils combn read.table write.csv
### {}


# Package gdata  is used to read and write excel file format
# Package hash   is used for implementation of hash data structure
# Package lattice and latticeExtra are used for  lattice plots
# Package oro.nifti is used for  colors in the plots 

# Function to calculate predicted survival values
calculateSiPredicted <- function( conc, a, h){
    predictedSIs<- (1 / (1 + ((conc/a)^h)))
    return(predictedSIs )
  }

# Function to calculate sum of sequare error
J <- function(sis, conc, a, h){
    siPredicted <- apply( matrix(conc, nrow=1 , ncol=length(conc)), 2, calculateSiPredicted, a, h )
    j <- (1/(2*length(conc))) * sum((sis-(siPredicted))^2) 
    return(j)
  }

# Function to calculate gradiant of change of IC50 values
gradJ_a <- function(sis, conc, a, h){
    siPredicted <- apply( matrix(conc, nrow=1, ncol=length(conc)), 2, calculateSiPredicted, a, h )
    errors <- sis - siPredicted
    gradj_a <- -1 * (1/length(conc)) *   sum((errors * ((siPredicted)^2) ) * ((conc/a)^h) * (h/a) ) 
    return(gradj_a)
  }

# Function to calculate gradiant of change of hill coefficient(h) values
gradJ_h <- function(sis, conc, a, h){
    siPredicted <- apply( matrix(conc, nrow=1, ncol=length(conc)), 2, calculateSiPredicted, a, h )
    errors <- sis - siPredicted
    gradj_h <- (1/length(conc)) *  sum( (errors* ((siPredicted)^2) ) * ((conc/a)^h) * log(conc/a) )  
    return(gradj_h)
  }

# Function to extract maximum error
Jmax <- function(sis, conc, a, h){
    siPredicted <- apply( matrix(conc, nrow=1, ncol=length(conc)), 2, calculateSiPredicted, a, h )
    return( max(abs(sis-siPredicted)) )
  }


# Implementation of nonlinear least sequare (nls) regression
nlsComIntAct <- function(sis, conc, a, h) {
  
  if (a<=0 ||h<=0 ){
    stop("Reasonable prediction of IC50 or Hill Coefficient can not be made on this data ")
  }
  
  
  hstartX <- 0    # Variable to store h values
  c50startX <-  0 # Variable to store IC50 values
  EX <- 0         # Variable to store mean of error sequare
  EXmax <- 0      # Variable to store maximum of mean of error sequare
  gXa <- 0        # Variable to store gradiant  of IC50
  gXh <- 0        # Variable to store gradiant of h
  ErrorValues <- 0 # For diagnostic purposes
  ErrorValuesMax <- 0 # For diagnostic purposes
  alpha_aXStore <- 0 # For diagnostic purposes
  alpha_hStore <- 0 # For diagnostic purposes
  
  # Starting guess values of h and IC50
  i <- 1
  hstartX[i] <- h
  c50startX[i] <- a
    
  gXa[i] <- gradJ_a(sis, conc, c50startX[i], hstartX[i])
  gXh[i] <- gradJ_h(sis, conc, c50startX[i], hstartX[i])

  EX[i] <- J(sis, conc, c50startX[i], hstartX[i])
  EXmax[i] <- Jmax(sis, conc, c50startX[i], hstartX[i])
  ErrorValues[i] <- EX[i] # For diagnostics
  ErrorValuesMax[i] <- EXmax[i] # For diagnostics
  
  # Step size
  alpha_aX <- c50startX[i]/ (1000 * abs(gXa[i]) )
  alpha_hX <-   hstartX[i]/(1000 * abs(gXh[i]))
  alpha_aXStore[i] <- alpha_aX # For diagnostics
  alpha_hStore[i] <- alpha_hX # For diagnostics
  
  flagContinue <- TRUE
  cStore<- 0
  hStore <- 0
  cStore[i] <- c50startX[i]
  hStore[i] <- hstartX[i]
  
  maxItr <- 50
  
  while(flagContinue==TRUE){
    i <- i+1
    c50startX[i] <- c50startX[i-1] - alpha_aX * gXa[i-1]
    hstartX[i] <- hstartX[i-1] - alpha_hX * gXh[i-1]
    
    EX[i] <- J(sis, conc, c50startX[i], hstartX[i])
    EXmax[i] <- Jmax(sis, conc, c50startX[i], hstartX[i])
    ErrorValues[i] <- EX[i] # For diagnostics
    ErrorValuesMax[i] <- EXmax[i] #For diagnostics
    
    # Stop Criteria
    if (EXmax[i] < 0.05){
      flagContinue <- FALSE
    }
  
    if ( (i-1) >= maxItr ){
      flagContinue <- FALSE
    }
  
  # Adaptive selection of step size
  cStore[i] <- c50startX[i]
  hStore[i] <- hstartX[i]
  
    if (EX[i] < EX[i-1]){
      
      alpha_aX <- 1.05 * alpha_aX
      alpha_hX <- 1.05 * alpha_hX 
      alpha_aXStore[i] <- alpha_aX # For diagnostics
      alpha_hStore[i] <- alpha_hX  # For diagnostics
      
      gXa[i] <- gradJ_a(sis, conc, c50startX[i], hstartX[i])
      gXh[i] <- gradJ_h(sis, conc, c50startX[i], hstartX[i])
      
    } else {
    
      c50startX[i] <- c50startX[i-1] 
      hstartX[i] <- hstartX[i-1] 
      EX[i] <- EX[i-1]
      alpha_aX <- 0.95 * alpha_aX
      alpha_hX <- 0.95 * alpha_hX 
      alpha_aXStore[i] <- alpha_aX  # For diagnostics
      alpha_hStore[i] <- alpha_hX  # For diagnostics
      gXa[i] <- gXa[i-1]
      gXh[i] <- gXh[i-1]
    }
  
  } # End of while
return(list(c(c50startX[i], hstartX[i])) )
}

# Function to select satrting values  for  parameters IC50 and h from the data 
startingGuessIC50nH <- function(drugObs_Mean,ConcentrationCleaned )
  {
  
    if (any(drugObs_Mean < 0.5) )
    {
      drugObs_Mean2Index <-  which(max( drugObs_Mean[which(drugObs_Mean < 0.50)]) == drugObs_Mean)      
      drugObs_Mean2 <- drugObs_Mean[drugObs_Mean2Index ] 
      drugObs_MeanIC2 <- ConcentrationCleaned[drugObs_Mean2Index ]
      if (drugObs_Mean2 <= 0 ) {drugObs_Mean2 <- 0.01}
      } else{
        drugObs_Mean2 <- min(drugObs_Mean) 
        if (drugObs_Mean2 <= 0 ){drugObs_Mean2 <- 0.01}
        drugObs_MeanIC2 <- ConcentrationCleaned[ which(min(drugObs_Mean) == drugObs_Mean) ]
    }
  
  
  if ((drugObs_Mean2 < 0.5)  )
  {
    searchAbleValues <- drugObs_Mean[1: (drugObs_Mean2Index-1) ]
    ConcentrationCleanedSearchable <- ConcentrationCleaned[1: (drugObs_Mean2Index-1)]
    minInd80  <- order(abs(searchAbleValues- 0.8))[1:2]
    drugObs_Mean1 <-   mean(searchAbleValues[ minInd80[!is.na(minInd80)==TRUE ] ])
    if (drugObs_Mean1 >= 1 ){drugObs_Mean1 <- 0.99}
    drugObs_MeanIC1 <- mean(ConcentrationCleanedSearchable[  minInd80[!is.na(minInd80)==TRUE ]  ] )
    } else {
      drugObs_Mean1 <- mean(drugObs_Mean[3:4] )
      if (drugObs_Mean1 >= 1 ){drugObs_Mean1 <- 0.99}
      drugObs_MeanIC1 <- mean(ConcentrationCleaned[3:4])
    }
  
    return(c(drugObs_Mean2,drugObs_MeanIC2, drugObs_Mean1,drugObs_MeanIC1  ))
  }


##' This function applies Loewe model.
##' @param xConcentration X drug concentrations tested in an experiment.
##' @param yConcentration Y drug concentrations tested in an experiment.
##' @param drugYObs_Mean Concentration wise mean survival of y drug treatments.  
##' @param drugXObs_Mean Concentration wise mean survival of x drug treatments. 
##' @return Loewe model values. 
##' @examples 
##' xConcentration <- c(0.00,0.20, 0.39,  0.78,  1.56, 3.12, 6.25, 12.50, 25.00, 50.0) 
##' yConcentration <- c(128,  64,  32,  16,   8,   4,   2,   0)
##' drugXObs_Mean <- c(0.9747255, 0.9197924, 0.9520692, 0.9517162, 0.9032701, 0.7892114,
##'                      0.6768190, 0.6524227, 0.4561164)
##' drugYObs_Mean <- rev( c( 0.93, 0.89, 0.73, 0.42, 0.24, 0.21, 0.11) )
##' rslt <- loeweModel( xConcentration, yConcentration, drugYObs_Mean, drugXObs_Mean)
##' @author Muhammad kashif 
##' @export
loeweModel<- function ( xConcentration, yConcentration,
                        drugYObs_Mean, drugXObs_Mean)
  {
    noOfRows <- length(yConcentration)
    noOfCols <- length(xConcentration)
  
    # Variables to store calculated parameters
    cx50 <- 0 
    hx   <- 0
    cy50 <- 0
    hy   <- 0
  
    # Indices of concentrations to be search for nls implementation
    xConcentrationCleaned <- xConcentration[2: length(xConcentration)] 
    yConcentrationRevCleaned <- rev(yConcentration[1: (length(yConcentration) - 1) ] )
    drugYObs_Mean_Rev <- rev(drugYObs_Mean)
    
   drugXObs_StartingGuess <- startingGuessIC50nH(drugXObs_Mean,xConcentrationCleaned )

   drugXObs_Mean2 <- drugXObs_StartingGuess[1] # Survival value that is just below IC50 or minimum
   drugXObs_MeanIC2 <- drugXObs_StartingGuess[2] # Concentration that's treatment resulted in  above mentioned survival value
   drugXObs_Mean1 <- drugXObs_StartingGuess[3]  # Survival values that is close to 80
   drugXObs_MeanIC1 <- drugXObs_StartingGuess[4] # Concentration that's treatment resulted in  above mentioned survival value close to 80
   
   drugYObs_StartingGuess <- startingGuessIC50nH(drugYObs_Mean_Rev, yConcentrationRevCleaned )

   drugYObs_Mean_Rev2 <- drugYObs_StartingGuess[1]
   drugYObs_Mean_RevIC2 <-drugYObs_StartingGuess[2]
   drugYObs_Mean_Rev1 <-drugYObs_StartingGuess[3]
   drugYObs_Mean_RevIC1 <-drugYObs_StartingGuess[4]
  
  
  # Starting guess values of h and IC50
  hstX <- log10( (drugXObs_Mean1/drugXObs_Mean2) * (( 1- drugXObs_Mean2)/(1-drugXObs_Mean1 ) )  ) / 
          log10(drugXObs_MeanIC2/drugXObs_MeanIC1)

  hstY <- log10((drugYObs_Mean_Rev1/drugYObs_Mean_Rev2) * (( 1- drugYObs_Mean_Rev2)/(1-drugYObs_Mean_Rev1 ) )  ) / 
          log10(drugYObs_Mean_RevIC2/drugYObs_Mean_RevIC1)
  
  c50stX <- drugXObs_MeanIC1/ (( (1-drugXObs_Mean1)/drugXObs_Mean1 ) ^ (1/hstX) )
    
  c50stY <- drugYObs_Mean_RevIC1/ (( (1-drugYObs_Mean_Rev1)/drugYObs_Mean_Rev1 ) ^ (1/hstY) )

  paramValuesX <- 0
  paramValuesY <- 0
  
  paramValuesX <- nlsComIntAct(drugXObs_Mean, xConcentrationCleaned, c50stX, hstX)
  paramValuesY <- nlsComIntAct(drugYObs_Mean_Rev, yConcentrationRevCleaned, c50stY, hstY) 
  
  cx50 <- paramValuesX[[1]][1]
  hx   <- paramValuesX[[1]][2]
  cy50 <- paramValuesY[[1]][1]
  hy   <- paramValuesY[[1]][2]
  
  
  # length(yConcentration):1 will make y drug moving from higher to lower cocnetration as in combination case
  # so variables concentration.product and  combXYObs_Mean have same dimensions
  concentration.product <- expand.grid(xConcentration[2:length(xConcentration)], yConcentration[1: (length(yConcentration)-1)] )  # X and Y cons with untreated well removed 
  
  # Inverse survival function
  Inv_Survival <- function(s, h, c50,alpha) {
      ( ( (((1- ( s ^  (1/ alpha))    )/  (s ^ (1/ alpha))     )) ^ (1/h  ) )* c50 ) 
  }
  
  # Function SurvivalAB implements the equation 6 in the manuscript corresponding to loewe response surface  
  SurvivalAB <- function(s, da, db, cah, ca50, cbh, cb50,alpha){ 
    abs(((da/ Inv_Survival(s, cah, ca50,alpha)) + 
           (db/ Inv_Survival(s, cbh, cb50, alpha)))- 1)
  }
  
  survivalindex.normal <- 0
  for (i in 1: nrow(concentration.product)){  
    survivalindex.normal.intermediate <- optimize(SurvivalAB, c(0.0001,1), concentration.product[i,1], concentration.product[i,2],
                                                  hx, cx50, hy, cy50, alpha=1, tol= 1e-16)
    
    survivalindex.normal[i] <- as.numeric(survivalindex.normal.intermediate[1])
  }
  
  # Variables concentration.product and  combXYObs_Mean have same dimensions
  # and by default matrix will create matrix bycol and therefore it will be 
  # in wrong order  and following code will correct the dimensions.
  loeweSynObs_Model <- matrix(survivalindex.normal, ncol=noOfCols-1, nrow=noOfRows-1, byrow=TRUE)
  
  return(loeweSynObs_Model)
  
} # End of loeweModel function



##' This function calculates Loewe synergy/antagonism and associated BIs.
##' @param rawDataPreProcessed Raw preprocessed experimental data. 
##' @param xConcentration X drug concentrations tested in an experiment.
##' @param yConcentration Y drug concentrations tested in an experiment.
##' @param nBoot Number of times to bootstrap in order  to calculate BIs.
##' @return Three lists, first containing Loewe Synergy/Antagonism, lower and upper bound of corresponding BI.  
##' Second list consists of global BI for maximum synergy observed in the experiment and third 
##' contains global BI of maximum antagonism.
##' @examples
##' library(gdata)
##' dataFile <- system.file("/raw/rawDataPreProcessed.csv", package="COMBIA")
##' dataSample <- read.csv(dataFile, header=FALSE )
##' xConc <- c(0.00,  0.20,  0.39,  0.78,  1.56,3.12,  6.25, 12.50, 25.00, 50) 
##' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
##' noOFBoot <- 50
##' rslt <- applyLoewe(as.matrix(dataSample), xConc, yConc, noOFBoot)
##' @author Muhammad kashif
##' @export
applyLoewe <- function( rawDataPreProcessed, xConcentration, yConcentration, nBoot)
  {
    # calculate total number of rows and colmuns
    noOfRows <- length(yConcentration); 
    noOfCols <- length(xConcentration); 
 
    # calculate number of replicates 
    rawDataPreProcessed_NA <- rawDataPreProcessed;
    rawDataPreProcessed_NA[which(rawDataPreProcessed == 0, arr.ind=TRUE)] <- NA
    replicateCount_individual <-   apply(rawDataPreProcessed_NA, 2, function(x) length (which(!is.na(x)) ))
    
    # Prepare data for Loewe analysis
    
    # Mean of data
    totalNumberofReplicates <- nrow(rawDataPreProcessed)
    rawDataPreProcessedMean <- apply(rawDataPreProcessed_NA,2,mean, na.rm=TRUE) 
    
    # Exract data of single constituent drugs and combinations
    rawDataPreProcessed_mat_temp <-matrix(rawDataPreProcessedMean, noOfRows, noOfCols)
    drugYObsMean <- as.vector(rawDataPreProcessed_mat_temp[1:noOfRows-1 ,1])  # all rows of first column, descending order
    drugXObsMean <- as.vector(rawDataPreProcessed_mat_temp[noOfRows,2:noOfCols]) # all columns of last row ascending oredr
    combXYObsMean <- rawDataPreProcessed_mat_temp[1:noOfRows-1, 2:noOfCols ]# x ascending and y descending
  
    # Call function to applying Loewe model  
    loeweSynObsModel <- loeweModel( xConcentration, yConcentration,
                                    drugYObsMean,drugXObsMean )
     # Note: Data in variable loeweSynObs_Model is formated simialr to combXYObs_Mean 
    
    # Calculate Loewe synergy
    # Data is in variable loeweSynergy is orderd as x drug descending and y drug ascending in each row 
    loeweSynergy <-  loeweSynObsModel - combXYObsMean # positive is synergy
    loeweSynergy_vector <- as.vector(loeweSynergy)
    # End of loewe Application

    # Calculate BIs 

    # Extract all residues
    myResidualMat <- matrix(rep(NA, totalNumberofReplicates*noOfRows* noOfCols),
                      nrow=totalNumberofReplicates, ncol=noOfRows* noOfCols)
    for(k in (1:(noOfRows * noOfCols) ))
      { 
        myResidualMat[ ,k] <- rawDataPreProcessed_NA[,k]- 
        rawDataPreProcessedMean[k]
      }

    # Clean for NAs
    myResidualMat_cleaned <- myResidualMat[!is.na(myResidualMat)]
    # Flip signs
    myResidualMat_cleanedFlip <- -1* myResidualMat_cleaned
    # All residues
    residuesAll <- c(as.vector(myResidualMat_cleaned), as.vector(myResidualMat_cleanedFlip) );

    # Bootstrap
    synergy_Loewe_pred_Boot <- matrix(rep(NA, ((noOfCols-1) * (noOfRows-1))  * nBoot ),
                           nrow=(noOfCols-1) * (noOfRows-1) , ncol=nBoot)
    bootCntr <- 1;
    while (bootCntr <= nBoot) 
      { 
        loeweSynObs_Model_Boot<-0  #Initialization 
        meanBars<-0
          for( i in (1:(noOfCols*noOfRows) ))
            {
              meanBars[i]<- mean(rawDataPreProcessedMean[i]+ sample(residuesAll , replicateCount_individual[i] , replace=TRUE))
            }
          meanBars_temp <- matrix(meanBars, noOfRows, noOfCols)
          drugYObsMeanBar <- as.vector(meanBars_temp[1:noOfRows-1 ,1])  # all rows of first column, descending order
          drugXObsMeanBar <- as.vector(meanBars_temp[noOfRows,2:noOfCols]) # all columns of last row ascending oredr
          combXYObsMeanBar <- meanBars_temp[1:noOfRows-1, 2:noOfCols ]# x ascending and y descending
        

        loeweSynObs_Model_Boot <- loeweModel(xConcentration, yConcentration,
                                                      drugYObsMeanBar, drugXObsMeanBar)
                  loeweSynergy_Boot <-  loeweSynObs_Model_Boot - combXYObsMeanBar# positive is synergy
                  # of values whose order is same of xy-combo-mean synergy_Loewe_pred_Boot[,1]
                  synergy_Loewe_pred_Boot[,bootCntr] <- as.vector(loeweSynergy_Boot )  ;
                  bootCntr<-bootCntr+1
      } # End of bootstrap 

    quantilesCI_ConcComb<- matrix(rep(0, 3*((noOfCols-1) * (noOfRows-1)) ), 3, ((noOfCols-1) * (noOfRows-1)) )
      for (pv in (1:((noOfCols-1) * (noOfRows-1))) )
        {
          # 0.025= Lower Bound and 0.975 is Upper bound
          quantilesCI_ConcComb[2:3,pv]<- apply(matrix(synergy_Loewe_pred_Boot[pv,], 1, nBoot) , 1, quantile, c(.025, 0.975  ))
        }  

    quantilesCI_ConcComb[1,]<-loeweSynergy_vector
    quantilesCI_Max<-    apply(matrix(apply(synergy_Loewe_pred_Boot,2, max), 1, nBoot) , 1, quantile, c(.025, 0.975  ))
    quantilesCI_Min<-    apply(matrix(apply(synergy_Loewe_pred_Boot,2, min), 1, nBoot) , 1, quantile, c(.025, 0.975  ))

  return(list(quantilesCI_ConcComb, quantilesCI_Max, quantilesCI_Min))

}# end of applyLoewe function


##' Given a well/wells alphanumeric name e.g. "B2", this function translates it to 
##' numeric values.
##' @param range Well name e-g B2.
##' @param excelFormate TRUE if well name is in spreadsheet format, otherwise cell culture
##' plate format.
##' @return Index of starting row, ending row, starting column and ending column.
##' @examples
##' rng <-c("B2")
##' exclF <- TRUE
##' rslt <-  extractValuesFromRange(rng, exclF)
##' @author Muhammad kashif
##' @export
extractValuesFromRange <- function(range, excelFormate)
{
  if (excelFormate==FALSE){
      columnstartindex <- as.integer(substr(range[1], 2, 4)) # Extract strating column
      rowstartindex <- which( letters[1 : 26] == tolower(substr(range[1], 0, 1))) # Extract strating row
      columnendindex <- as.integer(substr(range[2], 2, 4))  # Extract ending column
      rowendindex<- which( letters[1 : 26] ==tolower( substr(range[2], 0, 1))) # Extract ending row
  
      # Update variables for reformating data
      noOfRows<- (rowendindex- rowstartindex)+1 # Calculate no of rows in a single data
      noOfCols<- (columnendindex- columnstartindex)+1 # Calculate no of cols in a single data
      return(c(rowstartindex,rowendindex,columnstartindex,columnendindex  ))  
  
  } else { #  if input of range is similar to excel format
      rowstartindex <- as.integer(substr(range[1], 2, 4)) # Extract strating column
      columnstartindex <- which( letters[1 : 26] ==tolower( substr(range[1], 0, 1)) ) # Extract strating row
      rowendindex <- as.integer(substr(range[2], 2, 4))  # Extract ending column
      columnendindex <- which( letters[1 : 26] ==tolower( substr(range[2], 0, 1)) ) #Extract ending row
      
      # Update variables for reformating data
      noOfRows <- (rowendindex- rowstartindex)+1 # Calculate no of rows in a single data
      noOfCols <- (columnendindex- columnstartindex)+1 # Calculate no of cols in a single data
      return(c(rowstartindex,rowendindex,columnstartindex,columnendindex  ))    
  }
  
}



########################################
##' This function prints and saves 2D and 3D plots of synergy/antagonism.
##' @param processedData A matrix to plot.
##' @param xConcentration X drug concentrations.
##' @param yConcentration Y drug concentrations.
##' @param xDrug X drug name.
##' @param yDrug Y drug name.
##' @param cellLine Cell line or experiment name.
##' @return Plot the values. 
##' @examples
##' dataFile <- system.file("/raw/processedData.csv", package="COMBIA")
##' procData <- read.csv( dataFile, header=FALSE)
##' xConc <- c(0.00,  0.20, 0.39, 0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50) 
##' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
##' xD <- "X_Drug"
##' yD <- "Y_Drug"
##' clN <- "myCell"
##' rslt <- synAntPlot(as.matrix(procData),xConc,yConc, xD, yD, clN)  
##' @author Muhammad kashif
##' @export
synAntPlot<- function(processedData, xConcentration, yConcentration, xDrug, yDrug
                      , cellLine) 
  {
    cellLine <- cellLine;
    noOfRows <- length(yConcentration); 
    noOfCols <- length(xConcentration); 
  
    # Color Key
    ColorFun <- colorRampPalette(tim.colors(255))
    #Levelplot plots the columns and first column at botto so change arrangement of data as:
    #Plot only those are significantat = do.breaks( c(-60 , 60), 255)
    xConcentration_label <- xConcentration
    yConcentration_label <- yConcentration
  
    objectSynAntPlot <- levelplot( 100* (t(processedData[nrow(processedData):1, ])) , col.regions= ColorFun(255), 
                                    at = do.breaks( c(-50 , 50), 255),
                                    scales=list( x= list( at=c(1: (noOfCols-1) ), labels=xConcentration_label[2:length(xConcentration_label)] ,
                                                cex=1.3 ),  y= list( at=c(1:(noOfRows-1) ), 
                                                labels=rev(yConcentration_label[1:(length(yConcentration_label)-1)]), cex=1.3)
                                                ),
                                    xlab= as.vector(xDrug) , ylab= as.vector(yDrug), main=""
                                    ,colorkey=list(labels= list(cex=1.3,at=c(-50,-25,0,25,50), 
                                                          labels= as.character(c(-50,-25,0,25,50)))
                                             ) )
   
   
    # Update the font size                           
    objectSynAntPlotUpdated <- update(objectSynAntPlot,par.settings = list(fontsize = list(text = 14, points = 14),
                                                                            par.ylab.text = list(cex = 1.5) ,
                                                                            par.xlab.text = list(cex = 1.5) 
                                                                           ))
    # Print Formated plot
    print(objectSynAntPlotUpdated)
  
      
    myPathRoot<- system.file( package="COMBIA")
    #Add prefiX Combined_AnalyzedData
    myPath <- paste(myPathRoot,   "AnalyzedData_", sep="")
    myPath <- paste(myPath,   xDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   yDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   cellLine, sep="")
    myPath2 <- paste(myPath,   "_3D.png", sep="")
    myPath <- paste(myPath,   ".png", sep="")
    png(filename=myPath, width=800, height=600)
    print(objectSynAntPlotUpdated)
    dev.off()

    #For 3D graphs
      processedDataMat<-as.matrix(processedData)
      colnames(processedDataMat) <- c(1:ncol(processedDataMat))
      rownames(processedDataMat) <- c(1:nrow(processedDataMat))
      h <- cloud( processedDataMat , 
              col.facet= level.colors(  as.matrix(processedDataMat) , at = do.breaks( c(-0.5 , 0.5), 255), colors=TRUE,
                        col.regions=ColorFun(255)), 
              zlim=c(-0.5, 0.5),
              colorkey=list(col=ColorFun(255), at = do.breaks( c(-0.5 , 0.5), 255) ,
                         labels= list(cex=1.3,at=c(-0.5,-0.25,0,0.25,0.5), 
                                      labels= as.character(c(-50,-25,0,25,50)))
              ),
              xlab= paste(as.vector(yDrug),"          ") , ylab= paste("                      ", as.vector(xDrug)), main="", zlab="",
              scales=list( y= list(arrows=FALSE, at=c(1: (noOfCols-1) ), lab= xConcentration_label[2:length(xConcentration_label)], cex=1.1),  
                           x= list(arrows=FALSE, at=c(1:(noOfRows-1) ), lab=  
                                     c(yConcentration_label[1:(length(yConcentration_label)-2)], paste(yConcentration_label[(length(yConcentration_label)-1)], "   ", sep="") ), 
                                   cex=1.1),
                           z= list(arrows=FALSE, at=c(-0.5,-0.25,0,0.25,0.5),lab=c(-50,-25,0,25,50), cex=1.2)
              ),
              panel.3d.cloud=panel.3dbars, type="h",reference=TRUE,
              screen = list(z = -40, x = -25)
            )
  

  hUpdated <- update(h, par.settings = list(fontsize = list(text = 14, points = 14),
                                         par.ylab.text = list(cex = 1.5) ,
                                         par.xlab.text = list(cex = 1.5) 
                    ))
    
  print(hUpdated)
  
  png(myPath2, width=850, height=600)
  print(hUpdated)
  dev.off()
  
  
  # End of data saving
  }

##' Function "synergySignificant" calculates significant synergy/antagonism concentration combinations and saves synergy/antagonism values.
##' @param synergyCalculationLists List of synergy/antagonism calculations.
##' @param noOfRows Number of rows in the experiment.
##' @param noOfCols Number of columns in the experiment.
##' @param xDrug Name of drug at x-axis.
##' @param yDrug Name of drug at y-axis.
##' @param cellLine Cell line or experiment name.
##' @return Processed data. 
##' @examples
##' dataFile <- system.file("/raw/rawDataPreProcessed.csv", package="COMBIA")
##' dataSample <- read.csv(dataFile, header=FALSE)
##' nR <- 8
##' nC <- 10
##' rslt <- applyBliss(nR, nC,  as.matrix(dataSample ), 100) 
##' synergySignificant(rslt, nR, nC,"A", "B", "Cell")
##' @author Muhammad kashif
##' @export
synergySignificant <- function(synergyCalculationLists, noOfRows, noOfCols, xDrug, yDrug, cellLine )
  {
    #Separate the lists synergyCalculationLists
    synergy_Calculation_Matrix <- matrix(unlist(synergyCalculationLists[1]), nrow=3, ncol=(noOfRows-1)* (noOfCols-1))
    CIForBestSynergy <- unlist(synergyCalculationLists[2])
    CIForBestSynergy[3] <- max(synergy_Calculation_Matrix[1,])
    CIForBestAntagonism <- unlist(synergyCalculationLists[3]) 
    CIForBestAntagonism[3] <- min(synergy_Calculation_Matrix[1,])
  
    # Formating analyzed data same way as of the input data
    # Note: While calclating Synergy  first row counterpart in synergy_Calculation_Matrix
    # should be postive and should be nagative for antagonism.
  
    # Combinations that are significantlly (95% BI) synergistic or antagonistic
    levelOut <-  rep(0, ( (noOfRows-1)* (noOfCols-1) ) )   # Variable to hold the data to print out
    indecesOfSynergy  <- which (synergy_Calculation_Matrix[1,]>0)  # Indices of synergistic combinations
    indecesOfSynergy_Significant95  <- which( synergy_Calculation_Matrix[2,indecesOfSynergy] > 0) # Only those synergistic indices that are significant 95% 
    indecesOfAntagonism <- which (synergy_Calculation_Matrix[1,]<0) # Indices of antagonsitic combinations
    indecesOfAntagonism_Significant95  <- which(synergy_Calculation_Matrix[3,indecesOfAntagonism] < 0) # Only those antagonistic indices that are significant 95% 
    
    # Update output variable
    levelOut[indecesOfSynergy[indecesOfSynergy_Significant95] ] <- 
            synergy_Calculation_Matrix[1,indecesOfSynergy[indecesOfSynergy_Significant95] ]
  
    levelOut[indecesOfAntagonism[indecesOfAntagonism_Significant95] ] <- 
            synergy_Calculation_Matrix[1,indecesOfAntagonism[indecesOfAntagonism_Significant95] ]
  
    processedData <- matrix(levelOut, nrow=noOfRows-1, ncol=noOfCols-1 )
    
    myPathRoot=system.file("/raw", "", package="COMBIA")
    
    myPath <- paste(myPathRoot,   "AnalyzedData_", sep="")
    myPath <- paste(myPath,   xDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   yDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   cellLine, sep="")
    myPath <- paste(myPath,   ".csv", sep="")
    write.csv(100*processedData, file=myPath )
  
  
  #Save CIS for Synergy and Antagonism
    myPath <- 0
    myPath <- paste(myPathRoot,   "AnalyzedData_", sep="")
    myPath <- paste(myPath,   xDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   yDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   cellLine, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,"CIs", sep="")
    myPath <- paste(myPath,   ".csv", sep="")
    write.csv(100*c(CIForBestSynergy, CIForBestAntagonism), file=myPath )
  
  return(processedData)
  }# End of synergySignificant


##' On basis of input well ranges, the function extracts values of these ranges from data.
##' @param rawDataUnProcessed A data matrix.
##' @param wellRanges Ranges of wells.
##' @param wellplace Place of treated (case) well range.
##' @param simple TRUE if survival values are already calculated otherwise it is FALSE.
##' @param excelFormate True if well ranges are in excel format.
##' @return Replicate values.
##' @examples 
##' dataFile <- system.file("/raw/testData.xls", package="COMBIA")
##' rData <- read.xls( dataFile, sheet=1, skip=0, sep=",", nrows=41, 
##'                     fill=TRUE, header=FALSE,
##'                     blank.lines.skip = FALSE)[,1:13]
##' wellR= c( "l3:l10","m3:m10","b3:k10",  "l13:l20","m13:m20","b13:k20", 
##'             "l23:l30","m23:m30","b23:k30",  "l33:l40","m33:m40","b33:k40")
##' rslt <-  extractReplicateValues(rData, wellR, excelFormate=TRUE )
##' @author Muhammad kashif
##' @import gdata
##' @export
extractReplicateValues <- function(rawDataUnProcessed, wellRanges, wellplace=3, simple=FALSE, excelFormate=FALSE)
  { 
    if (simple==FALSE){
        cnt <- 0
        replicateCounter <- 0
        # List to hold the replicated data
        replicatedData <- vector("list", length(wellRanges)/3)
          while (cnt < (length(wellRanges)) ){
            range <- unlist(strsplit( wellRanges[ cnt + wellplace ] ,":")) # Split the range i-e from A2:D6 in to A2 and D6
            numericValues <- extractValuesFromRange(range, excelFormate)
            replicateCounter <- replicateCounter +1;
            replicatedData[[replicateCounter]] <- rawDataUnProcessed[c(numericValues[1]:numericValues[2] ),
                                                      c(numericValues[3]:numericValues[4] )]
          cnt <- cnt+3;
          }
  return(replicatedData)
  }
  else
  {
    cnt <- 0
    replicateCounter <- 0
    # List to hold the replicated data
    replicatedData<- vector("list", length(wellRanges)/3)
    while (cnt < (length(wellRanges)) ){
      range <- unlist(strsplit( wellRanges[ cnt + wellplace ] ,":")) # Split the range i-e from A2:D6 in to A2 and D6
      numericValues <- extractValuesFromRange(range, excelFormate)
      replicateCounter <- replicateCounter + 1;
      replicatedData[[replicateCounter]] <- rawDataUnProcessed[c(numericValues[1]:numericValues[2] ),
                                                        c(numericValues[3]:numericValues[4] )]
      cnt=cnt+3;
    }
    return(replicatedData)
  }
  
  }



##' The function calculates coefficient of variance (CV).
##' @param vals Data.
##' @return CV of input values.
##' @examples 
##' mData<- matrix(1:10, 2,5)
##' rslt<-  cVCal(mData)
##' @author Muhammad kashif
##' @export
cVCal<- function(vals){
        # Function to calculate CV[sd/mean ] of a vector
         cv <- sd(vals,na.rm=TRUE)/mean(vals,na.rm=TRUE);
         return(cv)
        }



##' Function to calculate possible unique  perturbations of replicates.
##' @param totalNumberofReplicates Total replicate number.
##' @return unique possible perturbations.
##' @examples
##' rslt<- createUniquePertbs(5) 
##' @author Muhammad kashif
##' @export
createUniquePertbs <- function(totalNumberofReplicates)
  {
    # Code to make unique perturbations of replicates that will be used incase if CV 
    # will be greater than 30%. 
    n <- totalNumberofReplicates; # Number of replicates will be decreasing
    perturblist <- vector("list",n) 
    for(i in (1:n) ){ 
        perturblist[[i]] <- combn(n,i,unique )
    }
  return(perturblist)
}  



##' This function removes outliers from data. Outliers are removed if data have more than two replicates.
##' Outlier is defined as an observation that has maximum contribution to data CV when CV is higher than a 
##' user defined threshold.
##' @param arrangeReplicates A data matrix.
##' @param minThersholdForCVCal CV threshold for outlier definition.
##' @param minThersholdForCV  Threshold of survival values to be excluded from CV calculation.
##' @return Data with out outliers.
##' @examples
##' dataFile <- system.file("/raw/rawDataPreProcessed.csv", package="COMBIA")
##' dataSample <- read.csv(dataFile, header=FALSE )
##' minThersholdForCV <- 0.3
##' minThersholdForCVCal <- 0.1
##' removeOutliers( as.matrix(dataSample ) ,minThersholdForCV,
##'        minThersholdForCVCal) 
##' @author Muhammad kashif
##' @export
removeOutliers <- function(arrangeReplicates, minThersholdForCVCal, minThersholdForCV)
  {
    # Add NA to data that was empty 
    arrangeReplicates_NA <- arrangeReplicates;
    arrangeReplicates_NA[which(arrangeReplicates == 0)] <- NA
    # Count replicates for each data value
    replicateCount_individual <-   apply(arrangeReplicates_NA, 2, function(x) length (which(!is.na(x)) ))
    # Count which replicates should be used, each list represents each row
    if (class(apply(arrangeReplicates_NA, 2, function(x) which(!is.na(x)) ))=="matrix")
        {
          replicateCount_indices <-   as.list(data.frame(apply(arrangeReplicates_NA, 2, function(x) which(!is.na(x)) )))
        } else {
          replicateCount_indices <-   apply(arrangeReplicates_NA, 2, function(x) which(!is.na(x)) )  
        }
    totalCombinationsnSingle <- (ncol(arrangeReplicates_NA))  
    for(i in  1:totalCombinationsnSingle  ){ 
        # Calculate CV
        # If survial values are less than a thresholds CV is not calculated
        if (length(which( arrangeReplicates_NA[,i] > minThersholdForCVCal) ) > 0)
            {
              currCV <- cVCal(c(arrangeReplicates_NA[,i]) )
              # if CV is greater than threshold then apply function of removal of maximum CV contributing values 
              if ( (currCV > minThersholdForCV)==TRUE ){
                  cvfound <- 0
                  currentNoOfPerts <- 0
                  currentPertList <- 0
        
                  if (replicateCount_individual[i] > 2){
                      perturblist <- createUniquePertbs(replicateCount_individual[i])
                      rangeOfPerturbation <- 2:(length(perturblist)-1)
                      cntPertb <- length(rangeOfPerturbation)
                      indecesOfReplicatesIn_arrangeReplicates <- replicateCount_indices[[i]]
                      while(cntPertb>0){
                            currentPertList <- perturblist[[rangeOfPerturbation[cntPertb]]]
                            pertSIval <- matrix(arrangeReplicates_NA[  indecesOfReplicatesIn_arrangeReplicates[currentPertList] ,i ], 
                                                      nrow=nrow(currentPertList), ncol=ncol(currentPertList) )
                            if(cntPertb > 1)
                              {
                                if(min(apply(pertSIval, 2, cVCal)) < minThersholdForCV)
                                  {
                                    indPert <- which(min(apply(pertSIval, 2, cVCal))==apply(pertSIval, 2, cVCal))
                                    indicesToBeZerod <- setdiff(c(replicateCount_indices[[i]]),
                                                            c(indecesOfReplicatesIn_arrangeReplicates[currentPertList[,indPert]]))
                                    arrangeReplicates[indicesToBeZerod,i] <- 0
                                    cvfound <- 1
                                break 
                              }
            
                              } else { 
                                    indPert <- which(min(apply(pertSIval, 2, cVCal))==apply(pertSIval, 2, cVCal))
                                    indicesToBeZerod <- setdiff(c(replicateCount_indices[[i]]),
                                                                c(indecesOfReplicatesIn_arrangeReplicates[currentPertList[,indPert]]))
                                    arrangeReplicates[indicesToBeZerod,i] <- 0
                              }
                      cntPertb <- cntPertb-1 
                      } # End of while
        
                    }# End of if replciates >2
                  } # End of if (length(which(currSIvals > 0.1) ) > 0)
                } # End of perturb loop
            } # End of for(i in 80)
   
    return(arrangeReplicates)
  }# End OF VARIABLITY CHECKING FUNCTION




##' This function performs synergy analysis based on Bliss model and calculates associated BIs and global BI.
##' @param noOfRows Number of data rows in an experiment.
##' @param noOfCols Number of data columns in an experiment.
##' @param rawDataPreProcessed Data matrix to be analyzed.
##' @param nBoot Number of bootstraps to calculates BIs.
##' @return Three lists are returned, first list contains Bliss synergy/antagonism values, lower and upper bounds of 
##' corresponding BI. Second list consists of global BI of maximum synergistic combination and third list 
##' consists of global BI of maximum antagonistic combination.
##' @examples 
##' dataFile <- system.file("/raw/rawDataPreProcessed.csv", package="COMBIA")
##' dataSample <- read.csv(dataFile, header=FALSE )
##' nR <- 8
##' nC <- 10
##' rslt <- applyBliss(nR, nC,  as.matrix(dataSample ), 500) 
##' @author Muhammad kashif
##' @export
applyBliss <- function(noOfRows, noOfCols, rawDataPreProcessed, nBoot)
    {
      
      # Number of replicates 
      rawDataPreProcessed_NA <- rawDataPreProcessed;
      rawDataPreProcessed_NA[which(rawDataPreProcessed == 0)] <- NA
      replicateCount_individual <-   apply(rawDataPreProcessed_NA, 2, function(x) length (which(!is.na(x)) ))
  
      totalNumberofReplicates <- nrow(rawDataPreProcessed)
      rawDataPreProcessedMean <- apply(rawDataPreProcessed_NA, 2, mean, na.rm=TRUE) 
  
      rawDataPreProcessed_mat_temp <- matrix(rawDataPreProcessedMean, noOfRows, noOfCols)
      drugYObsMean <- as.vector(rawDataPreProcessed_mat_temp[1:noOfRows-1 ,1])  # All rows of first column, descending order
      drugXObsMean <- as.vector(rawDataPreProcessed_mat_temp[noOfRows,2:noOfCols]) # All columns of last row ascending oredr
      combXYObsMean <- rawDataPreProcessed_mat_temp[1:noOfRows-1, 2:noOfCols ]# x ascending and y descending
      combXYObsMeanVector <- as.vector(combXYObsMean)
  
      # Apply Bliss on experimental data
      # Application of Bliss model only for mean values
      comb_xy_model_Bliss_inter <- expand.grid( drugXObsMean, drugYObsMean )  # Use all values , same formate as of  input data
      comb_xy_model_Bliss <-  (comb_xy_model_Bliss_inter[,1] * comb_xy_model_Bliss_inter[,2]) # Generate model data for normal cells
      # Variable comb_xy_model_Bliss_Mat has same formate as in the raw data fetch accordoring to wellRange parameter
      comb_xy_model_Bliss_Mat <- matrix(comb_xy_model_Bliss,noOfRows-1,noOfCols-1, byrow=TRUE) # Format data similar to input data format
  
      # Calculation of Bliss synergy values
      synergy_Bliss_obs <-   comb_xy_model_Bliss_Mat - combXYObsMean  # If positive then synergy
      synergy_Bliss_obsVector <- as.vector(synergy_Bliss_obs)
  
      # Extract all residues
        myResidualMat <-  matrix(rep(NA, totalNumberofReplicates*noOfRows* noOfCols),
                              nrow=totalNumberofReplicates, ncol=noOfRows* noOfCols)
        for(k in (1:(noOfRows * noOfCols) )){ #k=1
            myResidualMat[ ,k] <- rawDataPreProcessed_NA[,k]- 
            rawDataPreProcessedMean[k]
        }
      # Clean for NAs
      myResidualMat_cleaned <- myResidualMat[!is.na(myResidualMat)]
      # Flip signs
      myResidualMat_cleanedFlip <- -1* myResidualMat_cleaned
      # All residuels
      residuesAll <- c(as.vector(myResidualMat_cleaned), as.vector(myResidualMat_cleanedFlip) );
  
      
      replicateCount_individual_temp <- matrix(replicateCount_individual, noOfRows, noOfCols)
      drugYObsReplicates <- replicateCount_individual_temp[1:noOfRows-1 ,1]  # All rows of first column, descending order
      drugXObsReplicates <- replicateCount_individual_temp[noOfRows,2:noOfCols] # All columns of last row ascending oredr
      combXYObsReplicates <- replicateCount_individual_temp[1:noOfRows-1, 2:noOfCols ] # x ascending and y descending
      combXYObsReplicatesVector <- as.vector(combXYObsReplicates)
  
      # Bootstrap
        synergy_Best_Bar <- matrix(rep(NA, ((noOfCols-1) * (noOfRows-1))  * nBoot ),
                           nrow=(noOfCols-1) * (noOfRows-1) , ncol=nBoot)
        for (bootCntr in (1:nBoot) )
            {
              counter=1;
              for(i in  (1: length(drugXObsMean)) ) #xDrug mean
                  {
                  for (j in (1:length(drugYObsMean))) #yDrug mean
                    {
                    xBar <- 0
                    yBar <- 0
                    xyBar <- 0
        
                    xBar <- drugXObsMean[i]+ sample(residuesAll , drugXObsReplicates[i] , replace=TRUE)
                    yBar <- drugYObsMean[j]+ sample(residuesAll , drugYObsReplicates[j] , replace=TRUE)
                    xyBar <- combXYObsMeanVector[counter] + sample(residuesAll , combXYObsReplicatesVector[counter] , replace=TRUE)
        
                    xBarMean <- mean(xBar)
                    yBarMean <- mean(yBar)
                    xyBarMean <- mean(xyBar)
                    synergy_Best_Bar[counter, bootCntr] <-   (xBarMean*yBarMean)-xyBarMean
                    counter <- counter+1;
                    }
                  }
              }
  
        quantilesCI_ConcComb <- matrix(rep(0, 3*length(combXYObsMeanVector)), 3, length(combXYObsMeanVector) )
        for (pv in (1:length(combXYObsMeanVector)) )
            {
              #0.025= Lower Bound and 0.975 is Upper bound
              quantilesCI_ConcComb[2:3,pv] <- apply(matrix(synergy_Best_Bar[pv,], 1, nBoot) , 1, quantile, c(.025, 0.975  ))
            }  
  
        quantilesCI_ConcComb[1,] <- synergy_Bliss_obsVector
        # Maximum max(synergy_Bliss_obsVector)
        quantilesCI_Max <-    apply(matrix(apply(synergy_Best_Bar,2, max), 1, nBoot) , 1, quantile, c(.025, 0.975  ))
        # Minimum min(synergy_Bliss_obsVector)
        quantilesCI_Min <-    apply(matrix(apply(synergy_Best_Bar,2, min), 1, nBoot) , 1, quantile, c(.025, 0.975  ))
        return(list(quantilesCI_ConcComb, quantilesCI_Max, quantilesCI_Min))

    }  



##' It combines data from multiple files/experiments into a matrix.
##' Multiple experiments can be performed at the same or 
##' different concentration ranges. The function can also be used if 
##' there is only one file/experiments.
##' @param yConcentration Y drug concentrations. 
##' @param xConcentration X drug concentrations.
##' @param replNo Number of replicates in all files.
##' @param file  File name/names. 
##' @param totalNumberofReplicates Total number of replicates per files.
##' @param siReplicates Data.
##' @return Combined data of replicate survival indices from multiple experiments.
##' @examples 
##' xConc <- c(0.00,  0.20,  0.39,  0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50) 
##' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
##' rN <- 4
##' fN <- 1
##' trN <- 4
##' dataFile <- system.file("/raw/rawDataPreProcessed.csv", package="COMBIA")
##' dataSample <- read.csv(dataFile, header=FALSE )
##' replList <- list(vector, 4)
##' for( i in 1:4)
##' { replList[[i]] <- dataSample[i,] }
##' rslt <- combineDataFromMultipleFiles(list(yConc), 
##' list(xConc), rN,fN,trN, replList )
##' @author Muhammad kashif
##' @export
combineDataFromMultipleFiles <- function(yConcentration, xConcentration, replNo, 
                                        file,totalNumberofReplicates, siReplicates )
   {

  
    # Make a big matrix that contain all data
    allY  <-   unique(as.vector(sapply(yConcentration, function(x){as.numeric(x)})))
    allY  <-  allY[order(allY, decreasing=TRUE)] # Keep y concentration in proper order, large to small
    
    allX  <-  unique(as.vector(sapply(xConcentration, function(x){as.numeric(x)})))
    allX  <-  allX[order(allX, decreasing=FALSE)] # Keep x concentratuion in proper order, small to large
    bigDataList <- list("vector", replNo-1)
    bigMatrix <- matrix(rep(0,(length(allX)*length(allY))), nrow= length(allY), ncol=length(allX) )

    # Add data of all replicates from all files at proper location in bigMatrix
    replNoNew <- 0
    for (noOfFiles in (1:length(file) )){
      
        curr_xConcentration <- as.vector(xConcentration[[noOfFiles]])
        curr_yConcentration <- as.vector(yConcentration[[noOfFiles]])
     
         for (repl in (1:totalNumberofReplicates[noOfFiles] )){ 
          replNoNew<-replNoNew+1;
          # Convert data into matrix column wise
          currentSIData <- matrix(siReplicates[[replNoNew]], nrow=length( unlist(yConcentration[[noOfFiles]] ) ) , ncol=length(xConcentration[[noOfFiles]]) )
          for (i in (1:length(allX)) ){
            for (j in (1:length(allY)) ){
              
                xIndex <- which(curr_xConcentration==allX[i])
                yIndex <- which(curr_yConcentration==allY[j]) 
                if( (length(xIndex)==0) ||(length(yIndex)==0)   ) 
                    {
                      bigMatrix[j,i]<-0
                    }else{
                      bigMatrix[j,i]<-as.numeric(currentSIData[yIndex,xIndex])
                    }# end of descision structure
            }# end of akkY
          }#end of allX
        bigDataList[[replNoNew]] <- bigMatrix
        bigMatrix <- matrix(rep(0,(length(allX)*length(allY))), nrow= length(allY), ncol=length(allX) )
      }#end of replciates in a file
    }#end of files
  
    arrangeReplicates <- matrix(rep(0,sum(totalNumberofReplicates) * length(allX) *length(allY) ),
                            nrow=sum(totalNumberofReplicates), ncol=length(allY)*length(allX))
    for ( i in 1:sum(totalNumberofReplicates)){#i=1
        arrangeReplicates[i,] <- as.numeric(bigDataList[[i]])
        }
    
  
  return(arrangeReplicates)
  }# End of combineDataFromMultipleFiles



##' Read data from macSynergyII format and clean for outliers.
##' @param file Name of file containing data.
##' @param sheet Sheet number in a spread sheet workbook.
##' @param nrow Number of rows in that sheet to be read.
##' @param wellRangesExcel TRUE if input wells ranges are in excel format.
##' @param minThersholdForCVCal CV threshold for data outliers.
##' @param minThersholdForCV Threshold of survival values not used for CV calculations.
##' @return Matrix of outlier removed data replicates.
##' @examples
##' fl <- system.file("/raw/testData.xls", package="COMBIA")
##' sh <- 1
##' wellR <- list(c( "l3:l10","m3:m10","b3:k10",  "l13:l20","m13:m20","b13:k20", 
##'            "l23:l30","m23:m30","b23:k30",  "l33:l40","m33:m40","b33:k40"))
##' minThersholdForCV <- 0.3
##' minThersholdForCVCal <- 0.1
##' rslt <- readMacSynergyValues(fl, sh, nrow=41, wellR,  
##' minThersholdForCVCal, minThersholdForCV)
##' @author Muhammad kashif
##' @export
##' @import gdata
readMacSynergyValues<- function(file, sheet, nrow=41, wellRangesExcel,
                                minThersholdForCVCal, minThersholdForCV){
  # It will be vector that conatianing total number of replciates per file.
    totalNumberofReplicates <- as.numeric(lapply(wellRangesExcel, length))/3 
    siReplicates <- vector("list",sum(totalNumberofReplicates))
    replNo <- 1
    yConcentration <- vector("list", length(totalNumberofReplicates))
    xConcentration <- vector("list", length(totalNumberofReplicates))
    for(cn in (1:length(file)) )
      {
        plateData <- read.xls( file[cn], sheet=sheet, skip=0, sep=",", nrows=41, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[,1:13]
        # Extracct replicate data for every control, case and empty well
        # Create variables to hold data
        controlValues <- vector("list", totalNumberofReplicates[[cn]])
        emptyValues <-  vector("list", totalNumberofReplicates[[cn]])
        caseValues <-    vector("list", totalNumberofReplicates[[cn]]) 
        yConcentration[[cn]] <- round(as.numeric(as.matrix(plateData[3:10, 1])),3)
        xConcentration[[cn]] <- round(as.numeric(as.matrix(plateData[11, 2:11])),3)  
        controlValues <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=1, excelFormate=TRUE) # extract control wells
        emptyValues   <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=2, excelFormate=TRUE) # extract empty wells
        caseValues    <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=3, excelFormate=TRUE) # extract case wells
        sur <- function (x,y,z) {(x-z)/(y-z) }
        for (macI in (1:totalNumberofReplicates[cn]) )
            { 
              print(paste("Ratio between empty and control",  ( mean( as.integer(as.vector(controlValues[[macI]])), na.rm=TRUE  ) /mean( as.integer(as.vector(emptyValues[[macI]])), na.rm=TRUE)  ) ))
              print( paste(   paste(paste("CV for control:", macI  ), ":")  ,  100 * (sd( as.integer(as.vector(  controlValues[[macI]]) ), na.rm=TRUE )/mean( as.integer(as.vector(controlValues[[macI]] )), na.rm=TRUE) ) ) )
              # Calculate survival values
              # A matrix is decomposed by column wise
              siReplicates[[replNo]] <- lapply(as.numeric(as.matrix(caseValues[[macI]])),
                                      sur, y=mean( as.integer(as.vector(controlValues[[macI]]) )),
                                      z=mean( as.integer(as.vector(emptyValues[[macI]])))  )
              replNo <- replNo+1;
            }  
      } # loop no of files
  
    arrangeReplicates <- combineDataFromMultipleFiles(yConcentration, xConcentration,replNo, 
                                        file,totalNumberofReplicates, siReplicates )

    # Function to remove outliers
    rawDataPreProcessed <- removeOutliers(arrangeReplicates, minThersholdForCVCal, minThersholdForCV )
  return(rawDataPreProcessed);
  } # End of MacSynergy function




##' Read data from Fluorometric microculture cytotoxic assay (FMCA) data format and clean for outliers.
##' @param file Name of file to be read.
##' @param platetype 384 etc. 
##' @param keyposition Bar code position in the header row of each assayed plate.
##' @param selectionkey 65000.
##' @param platekey Barcode.
##' @param wells Well ranges.
##' @param minThersholdForCVCal CV threshold for data outliers.
##' @param minThersholdForCV Threshold of survival values not used for CV calculations.
##' @param yConcentration Concentrations of y drug.
##' @param xConcentration Concentrations of x drug.
##' @return Matrix of outlier removed replicated values. 
##' @examples 
##' fl <- system.file("/raw/FluoOptima_384_2014-03-28test.txt", package="COMBIA")
##' wls <- list(c("A11:H11", "A12:H12","A1:H10",   "I11:P11", "I12:P12","I1:P10", 
##'         "A23:H23", "A24:H24","A13:H22",   "I23:P23", "I24:P24","I13:P22")
##'                         )
##' pltype <- "384"
##' keypos <- 2     
##' seleckey <- "65000"
##' barCode <- 7049
##' minThersholdForCVCal <- 0.1 
##' minThersholdForCV <- 0.3
##' xConc <- c(0.00,  0.20,  0.39,0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50.00) 
##' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
##' readFMCAValues(fl, pltype, keypos, seleckey, barCode,
##'               wls, minThersholdForCVCal, minThersholdForCV, xConc, yConc   )
##' @author Muhammad kashif
##' @export
readFMCAValues <- function(file, platetype, keyposition,      
                          selectionkey, platekey, wells,
                          minThersholdForCVCal,
                          minThersholdForCV, 
                          yConcentration,
                          xConcentration
                          ){
  
  # Extracct no of replicates 
  # It will be vector that conatianing total number of replciates per file.
  totalNumberofReplicates <- as.numeric(lapply(wells, length))/3 
  siReplicates <- vector("list",sum(totalNumberofReplicates))
  replNo <- 1
  for(cn in (1:length(file)) )
    {
    # call to function that read raw FMCA data and convert it into corresponding survival values
     rawDataUnProcessed <- readFluostarPlates(filename= unlist(file[cn]), platetype=platetype[cn], keyposition =keyposition,      
                                                  selectionkey=selectionkey, platekey=platekey[cn], wells=wells[[cn]]                  
                                            )
     replicateValues <- extractReplicateValues(rawDataUnProcessed, wells[[cn]], wellplace=3, excelFormate=FALSE) # extract control wells
     for(cm in 1:length(replicateValues) )
        {
          siReplicates[[(replNo) ]] <- replicateValues[[cm]]
          replNo <- replNo+ 1;
        }
    } #Loop no of files

  arrangeReplicates <- combineDataFromMultipleFiles(list(yConcentration), list(xConcentration),replNo, 
                                                    file,totalNumberofReplicates, siReplicates
                                                    )
  # Function to remove outliers
  rawDataPreProcessed <- removeOutliers(arrangeReplicates,  minThersholdForCVCal, minThersholdForCV
                                        )
  
  return(rawDataPreProcessed);
  } # End of read FMCA function



##' Read data from other than MacSynergy and FMCA format and clean for outliers.
##' @param file Name of file to be read.
##' @param sheet Sheet number of a workbook. 
##' @param wellRangesExcel well ranges in excel format.
##' @param platetype 384.
##' @param minThersholdForCVCal CV threshold for data outliers.
##' @param minThersholdForCV Threshold of survival values not used for CV calculations.
##' @param survivalFunc A function to calculate survival values, 
##' by default survival is calculated as, Survival= treated- background/untreated control-background.
##' @param xConcentration Concentrations of drug at x-axis of data.
##' @param yConcentration Concentrations of drug at y-axis of data.
##' @return Matrix of outlier removed replicates.
##' @examples 
##' fl <- system.file("/raw/FluoOptima_384_2014-03-28test.xls", package="COMBIA")
##' wls <- list(  c(  "K1:K8", "L1:L8","A1:J8",     "K9:K16", "L9:L16","A9:J16", 
##'                   "W1:W8", "X1:X8","M1:V8",     "W9:W16", "X9:X16","M9:V16")
##'                   )
##' sh <- 1
##' pltype <- "384"
##' minThersholdForCVCal <- 0.1
##' minThersholdForCV<- 0.3
##' survivalFunc <- function (x,y,z) {(x-z)/(y-z)}
##' xConc <- c(0.00,  0.20,  0.39,  0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50.00) 
##' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
##' rslt <- readOtherValues(fl, sh, wls, pltype, minThersholdForCVCal, 
##'                 minThersholdForCV, survivalFunc, xConc, yConc )
##' @author Muhammad kashif
##' @export
readOtherValues<- function(file, sheet,  wellRangesExcel, platetype,
                            minThersholdForCVCal, minThersholdForCV,  survivalFunc,
                           xConcentration, yConcentration)
                          {
      if (platetype=="384"){
          rowsPerPlate <- 16; colsPerPlate<-24  
          } else{
          rowsPerPlate <- 8; colsPerPlate<-12 
        }
  
  
    # Extracct no of replicates
    # It will be vector that conatianing total number of replciates per file.
    totalNumberofReplicates <- as.numeric(lapply(wellRangesExcel, length))/3 
    siReplicates <- vector("list",sum(totalNumberofReplicates))
    replNo <- 1
    for(cn in (1:length(file)) )
        {
          # In raw formate
          plateData <- read.xls( unlist(file[cn]), sheet=sheet, skip=0, sep=",", nrows=rowsPerPlate, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[,1:colsPerPlate]
          # Extracct replicate data for every control, case and empty well
          # Create variables 
          controlValues <- vector("list", totalNumberofReplicates[[cn]])
          emptyValues <-   vector("list", totalNumberofReplicates[[cn]])
          caseValues <-    vector("list", totalNumberofReplicates[[cn]]) 
    
          controlValues <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=1, excelFormate=TRUE) # Extract control wells
          emptyValues   <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=2, excelFormate=TRUE) # Extract empty wells
          caseValues    <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=3, excelFormate=TRUE) # Extract case wells
    
          for (macI in (1:totalNumberofReplicates[cn]) )
              { #macI<-1
                
                print(paste("Ratio between empty and control",( mean( as.integer(as.vector(controlValues[[macI]])), na.rm=TRUE  ) / mean( as.integer(as.vector(emptyValues[[macI]])), na.rm=TRUE)   ) ))
                print( paste(   paste(paste("CV for control:", macI  ), ":")  ,  100 * (sd( as.integer(as.vector(  controlValues[[macI]]) ), na.rm=TRUE )/mean( as.integer(as.vector(controlValues[[macI]] )), na.rm=TRUE) ) ) )
      
                # Calculate Survival
                # A matrix is decomposed by column wise
                siReplicates[[replNo]] <- lapply(as.numeric(as.matrix(caseValues[[macI]])),
                                      survivalFunc, y=mean( as.integer(as.vector(controlValues[[macI]]) )),z=mean( as.integer(as.vector(emptyValues[[macI]])),na.rm=TRUE)  )
                
                replNo <- replNo+1;
              }  
    
        } # loop no of files
  
      arrangeReplicates <- combineDataFromMultipleFiles(list(yConcentration), list(xConcentration),replNo, 
                                                   file,totalNumberofReplicates, siReplicates )
    # Call function to remove outliers
      rawDataPreProcessed <- removeOutliers(arrangeReplicates,
                                             minThersholdForCVCal, minThersholdForCV
                                          )
      print(rawDataPreProcessed)
      return(rawDataPreProcessed);
    } # End of other format function






##' Function calculates significant synergy/antagonism in experimental data as per Bliss or Loewe model.
##' @param filename Name of data file.
##' @param sheet Optional, sheet number if spreadsheet file is used for input.
##' @param model  bliss or loewe.
##' @param inputFormates Any of the three c("fmca","macsynergy","others").
##' @param platetype Optional default is 384.
##' @param keyposition Optional default is 2.
##' @param selectionkey Optional default is 65000.
##' @param platekey Optional barcode.
##' @param minThersholdForCVCal Optional default is 0.15. 
##' @param minThersholdForCV Optional default is 0.3.
##' @param wells  Defines thes experiment layout in the well ranges, these ranges should  be in triplet form that is 
##' 1-control wells range, 2-empty wells range and 3-case wells range. An experiment
##' can have multiple replicates thus having multiple triplicates of well ranges.
##' @param yConcentration Y drug concentrations.
##' @param xConcentration X drug concentrations.
##' @param xDrug X drug name.
##' @param yDrug Y drug name.
##' @param cellLine Cell line/experiment name.
##' @param survivalFunc Optional default is function (x,y,z) {(x-z)/(y-z)} 
##' where x is combination/drug treatment outcome, y is the untreated control and 
##' z is the background noise.
##' @param nBoot Optional Number of time to bootstrap.
##' @return Store and print graphs/data of synergy/antagonism analyses.
##' @examples
##' fl <- system.file("/raw/testData.xls", package="COMBIA")
##' wellR <- list(c("l3:l10","m3:m10","b3:k10", "l13:l20","m13:m20","b13:k20", 
##'            "l23:l30","m23:m30","b23:k30", "l33:l40","m33:m40","b33:k40") )
##' sh <- 1
##' mdl <- "bliss"
##' xConc <- c(0.00,  0.20,  0.39,  0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50) 
##' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
##' xDrug <- "A"
##' yDrug <- "B"
##' cellLine <-"Cell"
##' calculateSynergy(filename = c(fl), sheet=1, model="bliss", inputFormates="macsynergy", 
##'                  wells=wellR, yConcentration= yConc, xConcentration=xConc,
##'                  xDrug=xDrug, yDrug=yDrug,cellLine=cellLine, nBoot=1000)
##'
##' @author Muhammad kashif
##' @export
calculateSynergy <- function( filename, sheet, model, inputFormates, platetype="384", keyposition =2,
                              selectionkey="65000", platekey=7051, minThersholdForCVCal=0.15, minThersholdForCV=0.3,
                              wells, yConcentration,xConcentration,
                              xDrug, yDrug,cellLine,     survivalFunc= function (x,y,z) {(x-z)/(y-z)}, nBoot=2500)
  {
    # Local variable declaration
    replicateValues <- 0;
    rawDataPreProcessed <- 0;  
    noOfRows <- 0   #  Variable store no of rows to be used to format data from matrix to vector and vice versa 
    noOfCols <- 0   #  Variable store no of cols to be used to format data from matrix to vector and vice versa
    inputFileFormates <- inputFormates
    # Extract number of rows and columns of an experiment
    noOfRows <- length(yConcentration)  # Calculate no of rows in a single data file/experiment
    noOfCols <- length(xConcentration)  # Calculate no of cols in a single data file/experiment
  
    # Apply appropriate functions based on formates and calculate values
    if (inputFileFormates=="fmca")
      {
        rawDataPreProcessed <- readFMCAValues( file=filename, platetype=platetype, keyposition =keyposition,      
                                                selectionkey=selectionkey, platekey=platekey, wells=wells,   
                                                minThersholdForCVCal=minThersholdForCVCal,
                                                minThersholdForCV=minThersholdForCV,  yConcentration, xConcentration
                                            )
      } else if(inputFileFormates=="macsynergy"){
          # wells in excel format
          wellRangesExcel <- wells
          rawDataPreProcessed <- readMacSynergyValues(filename, sheet, nrow=41, wellRangesExcel,
                                                      minThersholdForCVCal, minThersholdForCV)
      
      } else if(inputFileFormates=="others"){
          # wells in excel format
          wellRangesExcel <- wells
          rawDataPreProcessed <- readOtherValues( file=filename, sheet,  wellRangesExcel,platetype,
                                            minThersholdForCVCal, minThersholdForCV,  
                                            survivalFunc, yConcentration=yConcentration, xConcentration=xConcentration )
       
      }

   # Apply user selected synergy/antagonism detection model
   if (model=="bliss")
      { 
      
        # synergyBlissCalculation consists of three lists are resturned 
        # first row of raw Bliss synergy/antagonism
        # 2nd row lower bound of BIs
        # 3rd row upper bound of BIs
        # 2nd list is BI of the maximum significant synergistic concentration combination 
        # 3rd list is BI of the maximum significant antagonistic concentration combination 
        synergyBlissCalculationLists <-  applyBliss(noOfRows,noOfCols, rawDataPreProcessed, nBoot)
    
        # Process lists and save them in user readable formate
        processedData <- synergySignificant(synergyBlissCalculationLists,noOfRows,noOfCols, xDrug, yDrug, cellLine )
        # Synergy/Analysis plot 
        synAntPlot(processedData, xConcentration, yConcentration, xDrug, yDrug, cellLine)
    
      } else if(model=="loewe"){   
        # synergyBlissCalculation consists of three lists are resturned 
        # first row of raw Bliss synergy/antagonism
        # 2nd row lower bound of BIs
        # 3rd row upper bound of BIs
        # 2nd list is BI of the maximum significant synergistic concentration combination 
        # 3rd list is BI of the maximum significant antagonistic concentration combination 
        synergyLoeweCalculationLists <- applyLoewe(rawDataPreProcessed, xConcentration, yConcentration, nBoot)
        # Process lists and save them in user readale formate 
        processedData_Loewe <- synergySignificant(synergyLoeweCalculationLists,noOfRows,noOfCols, xDrug, yDrug, cellLine )
        # Synergy/Analysis plot 
        synAntPlot(processedData_Loewe, xConcentration, yConcentration, xDrug, yDrug, cellLine )
      }
  } # End of calculateSynergy Main function




################################Floustar data reading##############################################
# This function can read data from the file that may or may not be generated by automated 
# wet lab experimental plateforms, e-g fluostar, automated  robotics etc. 
# It calculates survival indices of the specified wells, one only needs to call "readFluostarPlates" 
# function with proper arguments to execute everything (Example are provided with sample values).
# It is also important to note that wells argument of readFluostarPlates should always be in triplet form that is 
# 1-control wells range, 2-empty wells range and 3-case wells range.
# It can also read the data from the file where a plate is read only one time, still it cope with variations if an experiment is
# repeated twice or many times in adjacent rows in the file. Another flexibility of is its ability to calculate Survival values of those
# experiments in which single plate and no repeated row is used. 
#
####################################################################################################

##' Read raw experimental data from a file.
##' This function reads the data from specified (excel,log, txt etc) file and store it in a data frame. 
##' @param filename Filename.ext.
##' @param separator Any character (, ; ' etc.) that is used as a separator in specified file.
##' @param sheet Need to use only when reading excel files. It is the number of the excel sheet to be read in a worksheet.
##' @param noofrows_skip Number of the rows in the file that should be skipped before starting the data reading.
##' @param readplates Number of the plates that you want to read from a set of plates in a file. This parameter can only 
##' be used with excel files. Otherwise it will be ignored.
##' @param numberofrowsperplate It is calculated on the basis of type of plates, that is, number of rows per plates are 17 for 
##' 384 well plates (16 lines from plates + 1 header lines) and 9 for 96 well plates (8 lines from plates + 1 header lines).
##' @param platetype Either 384 or 96. 
##' @return Data frame of data in the file.
##' @examples 
##' f <- system.file("/raw/optima.log", package="COMBIA")
##' fileDF <- readFile(filename = f, separator = ",",sheet=4, noofrows_skip=3,   
##' readplates=1, numberofrowsperplate=17, platetype="384")  
##' @author Muhammad Kashif
##' @export
readFile <- function(filename, separator, sheet, noofrows_skip, readplates, numberofrowsperplate, platetype)
  {
    # Change for Version3 Date 2011-05-27 
    # an argument platetype is added to specify the colulmns
    columnrange <- 0 
    if (platetype == "384"){
        columnrange <- 1:24
      } else if (platetype =="96"){
          columnrange <- 1:12
      }  
      
    # Change for Version2 Date 2011-02-17
    filenamechunks <- unlist( strsplit( filename, ".", fixed = TRUE))
    if ( filenamechunks[length(filenamechunks)] == "xls"){
        #if plate types are not specified
        if (platetype == ""){
            rawdata <- read.xls( filename, sheet = sheet, skip = noofrows_skip, nrows = readplates * numberofrowsperplate, header = FALSE, fill = TRUE)
        } else {
            rawdata <- (read.xls( filename, sheet = sheet, skip = noofrows_skip, nrows = readplates * numberofrowsperplate, header = FALSE, fill = TRUE)[,columnrange])
        }
        return(rawdata)   
    }else{
    
    # this part of code can read .log file and store in a data structure####
    
    # 1. Classical parser was not written because we were interested in barcode only that is always at second position of 
    # header and it can be extracted with methods provided by R. I guess we dont have any well designed and fairly
    # complex grammer for parser. Therefore regular expressions and regular expression like structures are the logical choice.
    
    # 2. Second choice of storing plate data was, one vector and one matrix (As per concept of database tables), 
    # one(vector) for storing headers and 
    # other one for storing the data. As we already know each plate has a line of header immediately followed by 16 lines of 
    # data(following the geometry of plates), therefore, it is not needed to have a link between these vector and matric through 
    # pointer like structure because 
    # 1st header * 1 and header > 1 * 17 give the same thing.
    
    
    # 3. Still another way of reading the data from .txt or .log is to use the built in functions of R that are related to file  
    # reading and string operations: they are
    # srcfile() for file name
    # getsrclines() read source file line by line
    # unlist(strsplit()) to split the line 
    
    # But it is clear from the above functions that they will still need datastructure to store their results 
    # after some  manipulations.
    
    # 4. Dataframe was the choice made because of there ability to store the numeric amd non numeric data. Here only one built in 
    # function (read.table) is sufficient to read and store data.
    
    #if plate types are not specified
    if (platetype == ""){
        rawdata <- data.frame(read.table( filename, sep = separator, fill = TRUE, skip = noofrows_skip))
      } else {  
        rawdata <- data.frame(read.table( filename, sep = separator, fill = TRUE, skip = noofrows_skip)[,columnrange])
      }
      return(rawdata)
    }
  }


##' Extracts the barcodes from  a dataset, each experimental plate has a barcode in a dataset.
##' Barcodes are extracted from the header of the each plate data from  position specified by keyposition argument.
##' @param keyposition Position of keyvalue in the header of plate.
##' @param rawdata An object(dataframe) of rawdata.
##' @param numberofrowsperplate  This argument is not needed when you call function "readFluostarPlates". The number of rows depend upon the
##' geometry of the plates. These are 16 in case of 384 well plates.
##' @param doubleplateexperiment This parameter has TRUE or FALSE value only. It is TRUE when an experiment is read
##' twice and one only want to choose only one of them. 
##' @return Barcodes.
##' @examples
##' f <- system.file("/raw/optima.log", package="COMBIA")
##' fileDF <- readFile(filename = f,  separator = ",", noofrows_skip=0,
##'                    platetype="384")  
##' Generatedbarcode <- extractKey(keyposition = 2, rawdata = fileDF, 
##'                               numberofrowsperplate = 17, 
##'                               doubleplateexperiment = FALSE) 
##' @author Muhammad Kashif
##' @export
extractKey <- function(keyposition, rawdata, numberofrowsperplate, doubleplateexperiment){
              # Generation of the missing barcode is not possible mainly there are limitless
              # ways plates can be arranged on FMCA day.               
  numberofrows  <- nrow(rawdata)  # total number of rows
  rawbarcode    <- seq(length = (numberofrows / numberofrowsperplate), from = 0, to = 0)   # initializing vector of raww barcodes
  rawbarcode[1] <- as.vector(rawdata [1, keyposition] )  # raw barcode of first plate      
  
  cnt=1     			
  while (cnt <= ((numberofrows / numberofrowsperplate)- 1) )
    {  
      # Loop to extract raw barcodes
      rawbarcode[cnt + 1] <- as.vector(rawdata[ (cnt * numberofrowsperplate) + 1, keyposition ])
      cnt=cnt+1 					
    }
  
    # Perform empty barcode check 
    if( any(rawbarcode == "NOREAD") )  
      {
        stop("Barcode is missing")
      }
  
  
#   #creating new barcode barcodevector<- c(2521,2521, 2531,2531, 2522,2522, 2532,2532,'NOREAD', 'NOREAD',2533, 2533,'NOREAD','NOREAD', 'NOREAD','NOREAD' )
#   barcodevector  <- rawbarcode
#   
#   #doubleplateexperiment<-FALSE
#   if (barcodevector[1] != "NOREAD") {
#         stepNext <- 2
#       if  (doubleplateexperiment == FALSE)
#         {
#         stepNext <- 1
#         }
#     
#         # Calculate  barcode of first element      
#         barcodeoffirstelement <- as.integer( barcodevector[ 1 ])        
#         cntr<-1
#         for(i in 1:(length(barcodevector)/stepNext) )
#           {# i=1
#           if (i%%2 ==1){
#             autoBarCode[cntr:(cntr+1)] <- barcodeoffirstelement + (i%/%2)
#             } else{ 
#               autoBarCode[cntr:(cntr+1)] <- barcodeoffirstelement + (10 + ((i-1)%/%2))
#             } 
#           cntr<-cntr+stepNext;
#           } #for 
#     }  else { 
#     }
  
  return(rawbarcode)
  }

# Note: Hash tables are used to store the FMCA data because it consists of a key (barcode) and a corresponding data.
# If COMBIA is used for large amount of FMCA data then the data in the hash table will be more transparent, logical and 
# easy to access than matrices.

##' Select one readout from the twice read data and store selected values.
##' In those cases where an experiment is read twice, only one readout from each pair is selected. 
##' The selection is based on presence of minimum number of "selection key" or read with 
##' highest mean is picked in case of tie.
##' @param rawdata An object (dataframe) of rawdata.
##' @param processedbarcode A vector containing all barcodes of data. In this case it is the output of function "extractKey".
##' @param numberofrowsperplate  This argument is not needed when you call function "readFluostarPlates". 
##' The number of rows depends upon the geometry of the plates. These are 16 in case of 384 well plates.
##' @param selectionkey key Value, that is used to select one of double read outs.
##' @param doubleplateexperiment This parameter value can either be TRUE or FALSE. It is TRUE when an experiment is read twice.
##' @return A hash table of selected readouts. Note: Hash tables are used to store the FMCA data because it consists of a key 
##' (barcode) and a corresponding data. If COMBIA is used for large amount of FMCA data then the data in the hash table will 
##' be more transparent, logical and easy to access than matrices.
##' @examples
##' f <- system.file("/raw/optima.log", package="COMBIA")
##' fileDF <- readFile(filename = f, separator = ",", noofrows_skip=0,
##'                    platetype="384") 
##' Generatedbarcode <-  extractKey(keyposition = 2,
##'     rawdata = fileDF, numberofrowsperplate = 17, doubleplateexperiment=TRUE) 
##' hashedplates <-  selectPlate(rawdata=fileDF,
##'     processedbarcode = Generatedbarcode, numberofrowsperplate=17,
##'     selectionkey="65000", doubleplateexperiment = FALSE )
##' @author Muhammad Kashif
##' @export

selectPlate <- function(rawdata, processedbarcode, numberofrowsperplate, selectionkey, doubleplateexperiment)
  {
    startindex <- 0 # variable starting from zero 
    dataplate1 <- 0 # variable storing data of plate1
    dataplate2 <- 0 # variable storing data of plate2
  
    # Variables for hashing function
    barcodekey <- 0
    hashedplates <- hash(keys = unique(processedbarcode), values = seq(length = length(unique(processedbarcode)), 0, 0))
    while(startindex < nrow(rawdata) ){
    # selecting one of the two consecutive plates ......
      dataplate1 <- as.integer(as.matrix(rawdata[(startindex + 2) : (startindex + numberofrowsperplate), ]))
      startindex <- startindex + numberofrowsperplate
      # Change for Version2 Date 2011-02-17 it is introduced here to control double and single plate experiments
      if (doubleplateexperiment == TRUE){
          dataplate2 <- as.integer(as.matrix(rawdata[(startindex + 2) : (startindex + numberofrowsperplate), ]))
          startindex <- startindex + numberofrowsperplate
            # selection on the basis of number of 65000 values
            if (length(dataplate1[dataplate1 == selectionkey] ) == length(dataplate2[dataplate2 == selectionkey])){
                if(mean(dataplate1) > mean(dataplate2)){
                  # plate1 selected
                  print("Plate 1 selected as per mean criteria")
                  barcodekey <- processedbarcode[(startindex / numberofrowsperplate)]
                  hashedplates[[as.character( barcodekey )]] <- dataplate1
                  } else{
                    # plate 2 selected
                    print("Plate 2 selected as per mean criteria")
                    barcodekey <- processedbarcode[(startindex / numberofrowsperplate)]
                    hashedplates[[as.character( barcodekey )]] <- dataplate2
                  }  
            } else {     
                    if( length(dataplate1[dataplate1 == selectionkey]) < length(dataplate2[dataplate2 == selectionkey])){
                        print("Plate 1 selected as per least 65000 values")
                        barcodekey<- processedbarcode[(startindex / numberofrowsperplate)]
                        hashedplates[[as.character(barcodekey)]] <- dataplate1
                      } else{
                          #plate2 selected 
                          print("Plate 2 selected as per least 65000 values")
                          barcodekey <- processedbarcode [(startindex / numberofrowsperplate)]
                          hashedplates[[as.character(barcodekey)]] <- dataplate2
                      }  
            }
        } else{ 
            # Change for Version2 Date 2011-02-17 #if single plate data is read
            # print("Plate 1 selected in single plate dataread")
            barcodekey <- processedbarcode[(startindex / numberofrowsperplate)]
            hashedplates[[as.character(barcodekey)]] <- dataplate1
        }
      }
    return(hashedplates)
  }



rangemean <- function(platebarcode, range, platetype, hashedplates, printcv )
  {
    # Start of function to extract range meanings, it works on the basis of plate labels and don't follow the excel style.
    plateforSI <- hashedplates[[as.character(platebarcode)]]
    if (platetype == "384"){
        dim(plateforSI) <- c(16,24)
      } else if (platetype == "96"){
      dim(plateforSI) <- c(8,12)
     }
  
    rangelist <- unlist(strsplit(range, ":"))
    rowstart <- substr(rangelist[1], 0, 1)
    columnstartindex <- as.integer(substr(rangelist[1], 2, 4))
    rowstartindex <- which( letters[1 : 26] == rowstart)
  
    rowend <- substr(rangelist[2 ], 0, 1)
    columnendindex <- as.integer(substr(rangelist[ 2 ], 2, 4 ))
    rowendindex <- which( letters[1:26] == rowend)
  
    welldata <- plateforSI[rowstartindex : rowendindex, columnstartindex : columnendindex]
    # FMCA Quality checks
    # Relation between empty and control   >5
    # Control well CV <30%
    #calculate Cv of wells and print it
  
    if (printcv==TRUE)
      {
        print(paste("CV% for control well is :", 100 * cVCal(welldata) ) )
      }
    return (mean( welldata ))
  }



# Calculates survival values for a range of input wells.
# Survival value is calculated as,  treated well - empty well/control well - empty well.    
# @param hashedplates A hash table of selec plates. It is the output of function "selectPlate".
# @param platekey It is the key of the plate whose S.I is needed to be calculated.
# @param platetype It is the type of plate (386 and 96). 
# @param rowsperexperiment It is the argument that specifies if the same experiment is reptead and how many times in a plate. If an experiment is 
# repeated twice in adjacent rows then average of its values will be used in the SI calculation. 
# @param wells This argument can take a list of arguments in the triplet form. Where first argument of triplet is the range of control wells,      
# second argument is the range of empty wells while third one is the range of case wells. It is made so that in labs plates layouts can differ  
# greatly. By using this  triplet scheme one can handel a number of palte layouts.
# @return A matrix with S.I showing values where they are actually exist on the plate.
# @examples
# f <- system.file("/raw/optima.log", package="COMBIA")
# fileDF <- readFile(filename = f, separator = ",", noofrows_skip=0,
#                     platetype="384") 
# Generatedbarcode <- extractKey(keyposition = 2,
#                             rawdata = fileDF, numberofrowsperplate = 17, 
#                             doubleplateexperiment=FALSE) 
# hashedplates <-  selectPlate(rawdata=fileDF,
#                             processedbarcode=Generatedbarcode, 
#                             numberofrowsperplate=17,
#                             selectionkey="65000", 
#                             doubleplateexperiment = FALSE  )
# survivalindeces <- calculatesi(hashedplates = hashedplates, 
#                                 platekey = "7051", platetype = "384",rowsperexperiment=1,
#                                 wells = c( "c8:h8","c1:n1","c3:c7",    "c8:h8","c1:n1","c9:c11", 
#                                 "c8:h8","c1:n1","e3:e7",     "c8:h8","c1:n1","e9:e11",
#                                 "c8:h8","c1:n1","g3:g7",     "c8:h8","c1:n1","g9:g11") 
#                               )    
# @author Muhammad Kashif
# @export 
calculatesi <- function(hashedplates, platekey, platetype, rowsperexperiment, wells )
  {
    # Change for Version2 Date 2011-02-17
    # three parameter cntrlrange, emptrange and caserange were removed and new "wells" is introduced to handel the 
    # variable number of these arguments.
    # rowsperexperiment argument stored the no of rows repeated per experiment 
  
    # Change for Version3 Date 2011-05-27 :: tolower function is added
    wellranges <- tolower(wells)
    if ( (length(wellranges)==0) | ( (length(wellranges)%%3)!=0 )){
        stop("Incorrect control, empty and case well ranges")
     }
  
    plateforSI <- 0 ## variable will store values of the plate underprocessing
    SI <- 0
    if (platetype == "384"){ ### this if statement will be helpful when we will be using same code for other plates may be 96 well
        plateforSI <- hashedplates[[as.character(platekey)]]  # extracting values of the plate based on barcode
        dim(plateforSI) <- c(16, 24)
        SI <- seq(length = 384, 0, 0)
        dim(SI) <- c(16, 24)
    } else if (platetype == "96" ){
          plateforSI <- hashedplates[[as.character(platekey)]]  # extracting values of the plate based on barcode
          dim(plateforSI) <- c(8, 12)
          SI <- seq(length = 96, 0, 0)
        dim(SI) <- c(8, 12)
    }  
    
      # Change for Version2 Date 2011-02-17
    cnt <- 0;
    cntrlMeanStore <- 0;
    emptyMeanStore <- 0;
    meanStoreCounter <- 0
    while (cnt < (length(wellranges)) ){
        #mean of the specified rangesControl and empty are calculated here
        cntrlmean <- rangemean(platekey, wellranges[cnt +1 ], platetype, hashedplates, printcv=TRUE)
        print(paste("Controlmean:", cntrlmean))
        # Change 280514, print  CV of CONTROLS
        # FMCA Quality checks
        # Ration between empty and control   >5
        # Control well CV <30%
        meanStoreCounter <- meanStoreCounter+1
        cntrlMeanStore[meanStoreCounter] <- cntrlmean
        emptmean  <- rangemean(platekey, wellranges[ cnt +2 ], platetype, hashedplates, printcv=FALSE)
        emptyMeanStore <- emptmean
        # These lines can be a part of the function but for clearity they are written seprate.
        range <- unlist(strsplit( wellranges[ cnt +3 ] ,":")) #split the range i-e from A2:D6 in to A2 and D6
        columnstartindex <- as.integer(substr(range[1], 2, 4)) #extract strating column
        rowstartindex <- which( letters[1 : 26] == substr(range[1], 0, 1)) #extract strating row
        columnendindex <- as.integer(substr(range[2], 2, 4))  #extract ending column
        rowendindex<- which( letters[1 : 26] == substr(range[2], 0, 1)) #extract ending row
          while(rowstartindex <= rowendindex)
            {
              if (rowsperexperiment > 1){
                SI[rowstartindex, columnstartindex : columnendindex] <- 
                    (( colMeans( matrix(plateforSI[rowstartindex: (rowstartindex + rowsperexperiment - 1),  columnstartindex : columnendindex],
                              nrow = length(rowstartindex: (rowstartindex + rowsperexperiment - 1)) , ncol = length(columnstartindex : columnendindex) )) 
                              - emptmean) / (cntrlmean - emptmean))
              }else{
                SI[rowstartindex, columnstartindex : columnendindex] <- (  plateforSI[rowstartindex,  
                                                                              columnstartindex : columnendindex]  -                      
                                                                     emptmean) / (cntrlmean - emptmean)
              }
          rowstartindex <- rowstartindex + rowsperexperiment
          } # end of the while
      cnt <- cnt+3 
    } # end of while loop of Change for Version2 Date 2011-02-17
  
    # FMCA Quality checks
    # Ratios between empty and control   >5
    # Control well CV <30%
    cntrlMeanStore;
    emptyMeanStore;
    print(paste("Ratio between empty and control", 100* (emptyMeanStore/ mean(cntrlMeanStore))  ) )
    return(SI)
  }


##' Read a file containing raw read outs and calculates survival values.
##' This function calculates survival values as 
##' (survival value= treated - empty/ untreated control- empty) of the wells specified
##' in the wells argument. Wells argument should  be in triplet form, they are,  untreated control well range, empty well range and treated well range.
##' The function also handles the double plate experiments in which one plate is read twice and only one of them is selected.
##' @param filename filename.ext.
##' @param separator It is the separation character within the file assigned to filename.  
##' @param noofrows_skip Number of the rows to skip in in data file.
##' @param sheet Sheet number conataining data in a workbook.
##' @param readplates Number of the plates to read from a set of plates from an excel file, This feature is only useful with spreadsheet files.
##' @param platetype Two types of plate formats (384 and 96) are supported. 
##' @param doubleplateexperiment This parameter value can either be TRUE or FALSE. It is TRUE when an experiment is read twice.
##' @param keyposition Position of barcode in the header. 
##' @param selectionkey Value, that is used to select one of double read outs.
##' @param platekey Barcode of data plate.    
##' @param rowsperexperiment It is the argument that specifies if the same experiment is repeated 
##' and how many times in a plate. If an experiment is repeated twice in adjacent rows then average 
##' of its values will be used in the survival value calculation. 
##' @param wells Well ranges to read data. 
##' @return Matrix containing survival values.
##' @examples
##' f <- system.file("/raw/optima.log", package="COMBIA")
##' platematrix <-readFluostarPlates(filename=f, platetype="384", keyposition=2,      
##' selectionkey="65000",platekey=7051,
##'  wells=c( "c8:h8","c1:n1","c3:c7",    "c8:h8","c1:n1","c9:c11", 
##'         "c8:h8","c1:n1","g3:g7",    "c8:h8","c1:n1","g9:g11",
##'
##'        "i8:n8","c1:n1","i3:i7",     "i8:n8","c1:n1","i9:i11",
##'        "i8:n8","c1:n1","k3:k7",     "i8:n8","c1:n1","k9:k11",
##'        "i8:n8","c1:n1","m3:m7",     "i8:n8","c1:n1","m9:m11",
##' 
##'        "c18:h18","c1:n1","c14:c17",     "c18:h18","c1:n1","c19:c22",
##'        "i18:n18","c1:n1","m14:m17",     "i18:n18","c1:n1","m19:m22",
##'        "i18:n18","c1:n1","m13:m13"
##'
##'       )                  
##'    )
##' @author Muhammad Kashif
##' @export  
readFluostarPlates <- function(filename, separator=",", noofrows_skip=0, sheet="1", readplates=1, platetype, doubleplateexperiment=FALSE, keyposition,    
                               selectionkey, platekey, rowsperexperiment=1, wells )
    {
      numberofrowsperplate <- 0  
      if (platetype == "384"){
          numberofrowsperplate <- 17
      } else if (platetype == "96") {
      numberofrowsperplate <- 9
    }
  
    rawdata <- readFile(filename, separator, sheet, noofrows_skip, readplates, numberofrowsperplate, platetype)
    Generatedbarcode <- extractKey(keyposition, rawdata, numberofrowsperplate, doubleplateexperiment)
    hashedplates <-  selectPlate(rawdata, processedbarcode = Generatedbarcode, numberofrowsperplate, selectionkey, doubleplateexperiment)
    survivalindeces <- calculatesi(hashedplates, platekey, platetype, rowsperexperiment, wells )    
  
  }









