#' Descretize Quantitative Variables.
#'
#' Takes in a data frame which contains only Quantitative varables in columns. Standadize the variables.
#' Discretize quantitative variables and returns discretized quantitative variables. Discretization was performed by
#' equal width bining algorithm.
#' @param myDataQuant A data frame which includes quantitative variables in columns.
#' @param noice Noice indicator. If noice = TRUE data standerdization is done by
#' deviding the difference between data point and median of the variable by the range of the variable.
#' If noice = FALSE data standerdization is done by deviding the difference between data point and
#' mean of the variable by the standard deviation of the variable.
#' @return A data frame consists of discretized quantitative variables.
#' @examples
#' QuantVars <- data.frame(Qnvar1 = c(1.5,3.2,4.9,5), Qnvar2 = c(4.8,2,1.1,5.8))
#' Discretized <- discretizeQuant(QuantVars)
#' @import stats
#' @export




discretizeQuant <- function(myDataQuant, noice = TRUE){

  i = 1
  noOfObj <-nrow(myDataQuant)
  noOfVar <- ncol(myDataQuant)
  NormData <- data.frame(data = NA)

  if(noice == TRUE){
    for(i in 1:ncol(myDataQuant)){
      NormData_i <- (myDataQuant[,i] - stats::median(myDataQuant[,i]))/(max(myDataQuant[,i]) - min(myDataQuant[,i]))
      NormData <- cbind(NormData, NormData_i)
    }
  }

  if(noice == FALSE){
    for(i in 1:ncol(myDataQuant)){
      NormData_i <- (myDataQuant[,i] - mean(myDataQuant[,i]))/ stats::sd(myDataQuant[,i])
      NormData <- cbind(NormData, NormData_i)
    }
  }

  NormData <- data.frame(NormData[,-1])
  colnames(NormData) <- colnames(myDataQuant)

  #decretize by equal no of Intervals
  dectrTzed <- data.frame()
  dectrTzedFac <- matrix(data = NA, ncol = noOfVar, nrow = noOfObj)
  noOfInt <- ceiling(sqrt(noOfObj))
  IntWidth <- (max(NormData[,i]) - min(NormData[,i])) / noOfInt

  for(i in 1:ncol(NormData)){
    for (j in 1:noOfObj){
      k = 1
      while(k <= noOfInt){
        if(NormData[j,i] <= min(NormData[,i]) + k*IntWidth){
          dectrTzed[j,i] <- k
          break
        }
        if(k == noOfInt)
          dectrTzed[j,i] <- k
        k = k + 1
      }
    }
  }
  i = 1

  #Find Significance
  for(i in 1:ncol(myDataQuant)){
    dectrTzedFac[,i] <- as.factor(dectrTzed[,i])
  }
  dectrTzedFac <- data.frame(dectrTzedFac)
  colnames(dectrTzedFac) <- colnames(myDataQuant)

  return(dectrTzedFac)

}


#' Calculate Conditional Probabilities.
#'
#' Takes in a data frame which contains only qualitative variables. Discretized quantitative variables
#' , a mixture of qualitative variables and discretized quantitative variables are also accepted.
#' Calculates conditional probabilities for each pair of attribute values in  the data frame.
#' Returns a data frame consists of J, A, B and C in columns where Pr(A|B) = C and J is the column number in the input
#' data frame corresponding to the values in A.
#' @param myDataAll A data frame which includes qualitative variables OR discretized quantitative variables
#' OR a mixture of qualitative variables and discretized quantitative variables in columns.
#' @return A data frame with four columns J, A, B and C in columns where Pr(A|B) = C and J is the column number in the input
#' data frame corresponding to the values in A.
#' @examples
#' QualiVars <- data.frame(Qlvar1 = c("A","B","A","C"), Qlvar2 = c("Q","Q","R","Q"))
#' CalcForQuali <- calcCondProb(QualiVars)
#' QuantVars <- data.frame(Qnvar1 = c(1.5,3.2,4.9,5), Qnvar2 = c(4.8,2,1.1,5.8))
#' Discretized <- discretizeQuant(QuantVars)
#' CalcForQuant <- calcCondProb(Discretized)
#' AllQualQuant <- data.frame(QualiVars, Discretized)
#' CalcForAll <- calcCondProb(AllQualQuant)
#' @export




calcCondProb <- function(myDataAll){

  ProbData <- data.frame()
  NoOfVars <- ncol(myDataAll)

  for(i in 1:NoOfVars){
    levelsOfDataFirst <- data.frame(levels(as.factor(myDataAll[,i])))
    NofFirst <- nrow(levelsOfDataFirst)

    j = 1
    while(j <= NoOfVars){
      if(j == i)
        j = j + 1
      if(j > NoOfVars)
        break

      levelsOfDataSec <- data.frame(levels(as.factor(myDataAll[,j])))
      NofSec <- nrow(levelsOfDataSec)

      for(k in 1:NofFirst){
        filteredfRST = myDataAll[myDataAll[,i] == levelsOfDataFirst[k,1], ]
        denomVal <- nrow(filteredfRST)
        for(a in 1:NofSec){
          filteredSec <- filteredfRST[filteredfRST[,j] == levelsOfDataSec[a,1], ]
          if (is.na(filteredSec[1,1]) == FALSE){
            numeVal <- nrow(filteredSec)
          }else
            numeVal <- 0
          jVal = j
          occur <- levelsOfDataSec[a,1]
          given <- levelsOfDataFirst[k,1]
          condProbVal <- numeVal/denomVal
          newRowToBind <- data.frame(jVal, occur, given, condProbVal)
          ProbData <- rbind(ProbData,newRowToBind)
        }
      }
      j = j + 1
    }
  }
  return(ProbData)
}




#' Calculate Distance  Between Attribute Values.
#'
#' Takes in a data frame which contains only qualitative variables. Discretized quantitative variables
#' , a mixture of qualitative variables and discretized quantitative variables are also accepted.
#' Calculates distance between each pair of attribute values for a given attribute. This calculation
#' is done according to the method proposed by Ahmad & Dey (2007).
#' @param myDataAll A data frame which includes qualitative variables OR discretized quantitative variables
#' OR a mixture of qualitative variables and discretized quantitative variables in columns.
#' @details distBetPairs is an implementtion of the method proposed by Ahmad & Dey (2007)
#' to find the distance between two catogorical values corresponding to a qualitative
#' variable. This distance measure considers distribution of values in the data set.
#' This function is also used to find the distance between discretized values
#' corresponding to quantitative variables which are used in calculating the significance
#' of quantitative attributes. See Ahmad & Dey (2007) for more datails.
#' @return A data frame with four columns J, A, B and C in columns where Distance(A, B) = C
#' and J is the column number in the input data frame corresponding to the values in A.
#' @references Ahmad, A., & Dey, L. (2007). A k-mean clustering algorithm for mixed numeric and categorical data. Data & Knowledge Engineering, 63(2), 503-527.
#' @examples
#' QualiVars <- data.frame(Qlvar1 = c("A","B","A","C"), Qlvar2 = c("Q","Q","R","Q"))
#' library(dplyr)
#' distForQuali <- distBetPairs(QualiVars)
#' QuantVars <- data.frame(Qnvar1 = c(1.5,3.2,4.9,5), Qnvar2 = c(4.8,2,1.1,5.8))
#' Discretized <- discretizeQuant(QuantVars)
#' distForQuant <- distBetPairs(Discretized)
#' AllQualQuant <- data.frame(QualiVars, Discretized)
#' distForAll <- distBetPairs(AllQualQuant)
#' @export



distBetPairs <- function(myDataAll){

  #library(dplyr)

  DistData <- data.frame()
  NoOfVars <- ncol(myDataAll)
  ConditioanalProbs <- calcCondProb(myDataAll)

  for(i in 1:NoOfVars){

    levelsOfDataFirst <- data.frame(levels(as.factor(myDataAll[,i])))
    NofFirst <- nrow(levelsOfDataFirst)

    for (j in 1:NofFirst){

      for (k in j+1:NofFirst){
        a = 1
        DistBetPair_ia <- data.frame()

        while(a <= NoOfVars){
          if(a == i)
            a = a + 1

          if(a > NoOfVars)
            break

          filterForX <-  dplyr::filter(ConditioanalProbs, ConditioanalProbs$jVal == a & ConditioanalProbs$condProbVal != 0 & ConditioanalProbs$given == as.character(levelsOfDataFirst[j,1]))
          wVals <- data.frame(filterForX$occur)
          filterForY <- dplyr::filter(ConditioanalProbs, ConditioanalProbs$jVal == a & ConditioanalProbs$given == as.character(levelsOfDataFirst[k,1]))
          b = 1
          while(b <= nrow(wVals)){
            filterForY <- dplyr::filter(filterForY, filterForY$occur != as.character(wVals[b,1]))
            b = b + 1
          }
          SumConPropX <- sum(filterForX$condProbVal) + sum(filterForY$condProbVal) - 1
          DistBetPair_ia <- rbind(DistBetPair_ia, SumConPropX)
          a = a + 1
        }

        VarIndex <- i
        xVal <- levelsOfDataFirst[j,1]
        yVal <- levelsOfDataFirst[k,1]
        DistBetPair <- mean(DistBetPair_ia[,1])
        newRowToBind <- data.frame(VarIndex, xVal, yVal, DistBetPair)
        DistData <- rbind(DistData, newRowToBind)
        DistData <- dplyr::filter(DistData, is.na(yVal) == FALSE)
      }
    }
  }
  return(DistData)
}





#' Calculate Significance of Quantitative Attributes.
#'
#' Takes in two data frames where first contains only qualitative attributes and the other
#' contains only quantitative attributes. Function calculates significance of quantitative
#' attributes based on the method proposed by Ahmad & Dey (2007).
#' @param myDataQuali A data frame which includes only qualitative variables in columns.
#' @param myDataQuant A data frame which includes only quantitative variables in columns.
#' @details signifOfQuantVars is an implementtion of the method proposed by Ahmad & Dey (2007)
#' to calculate the significance of quantitative attributes. Signinficance of an attribute is an
#' important fact to consider in the process of clustering. To calculate the significance
#' quantitative attributes are discreized first. These  significace values are used in calculating
#' distance between any two numeric values of aquantitative attribute. See Ahmad & Dey (2007) for
#' more datails.
#' @return A data frame with two columns A and B where A represents variable number and B
#' represents significane of corresponding variable.
#' @references Ahmad, A., & Dey, L. (2007). A k-mean clustering algorithm for mixed numeric and categorical data. Data & Knowledge Engineering, 63(2), 503-527.
#' @examples
#' QualiVars <- data.frame(Qlvar1 = c("A","B","A","C"), Qlvar2 = c("Q","Q","R","Q"))
#' QuantVars <- data.frame(Qnvar1 = c(1.5,3.2,4.9,5), Qnvar2 = c(4.8,2,1.1,5.8))
#' SigOfQuant <- signifOfQuantVars(QualiVars, QuantVars)
#' @export


signifOfQuantVars <- function(myDataQuali, myDataQuant){

  DiscretQuant <- discretizeQuant(myDataQuant)
  allData <- data.frame(myDataQuali, DiscretQuant)
  DistForAll <- distBetPairs(allData)

  QuantVarStart <- ncol(myDataQuali) + 1
  LoopStopVal <- ncol(allData)
  SignfData = data.frame()

  ct = 1
  for(i in QuantVarStart:LoopStopVal){
    VarIndex = DistForAll[i,1]
    SigDataForQuant <- dplyr::filter(DistForAll, VarIndex == i)
    SignfData[ct,1] <- ct
    SignfData[ct,2] <- mean(SigDataForQuant$DistBetPair)
    ct = ct + 1
  }
  return(SignfData)
}




#' Calculate Dissimilarity Matrix for Mixed Attributes.
#'
#' Takes in two data frames where first contains only qualitative attributes and the other
#' contains only quantitative attributes. Function calculates the dissimilarity matrix
#' based on the method proposed by Ahmad & Dey (2007).
#' @param myDataQuali A data frame which includes only qualitative variables in columns.
#' @param myDataQuant A data frame which includes only quantitative variables in columns.
#' @details calcDissimMat is an implementtion of the method proposed by Ahmad & Dey (2007)
#' to calculate the dissimilarity matrix at the presence of both qualitative and quantitative
#' attributes. This approach finds dissimilarity of qualitative and quantitative attributes seperately
#' and the final dissimilarity matrix is formed by combining both. See Ahmad & Dey (2007) for
#' more datails.
#' @return A dissimilarity matrix. This can be used as an input to pam, fanny, agnes and diana functions.
#' @references Ahmad, A., & Dey, L. (2007). A k-mean clustering algorithm for mixed numeric and categorical data. Data & Knowledge Engineering, 63(2), 503-527.
#' @examples
#' QualiVars <- data.frame(Qlvar1 = c("A","B","A","C","C","A"), Qlvar2 = c("Q","Q","R","Q","R","Q"))
#' QuantVars <- data.frame(Qnvar1 = c(1.5,3.2,4.9,5,2.8,3.1), Qnvar2 = c(4.8,2,1.1,5.8,3.1,2.2))
#' DisSimMatCalcd <- calcDissimMat(QualiVars, QuantVars)
#'
#' agnesClustering <- cluster::agnes(DisSimMatCalcd, diss = TRUE, method = "ward")
#' silWidths <- cluster::silhouette(cutree(agnesClustering, k = 2), DisSimMatCalcd)
#' mean(silWidths[,3])
#' plot(agnesClustering)
#'
#' PAMClustering <- cluster::pam(DisSimMatCalcd, k=2, diss = TRUE)
#' silWidths <- cluster::silhouette(PAMClustering, DisSimMatCalcd)
#' plot(silWidths)
#' @export




calcDissimMat <- function(myDataQuali, myDataQuant){

  #library(dplyr)

  dissimForQuant <- data.frame()
  dissimForQuali <- data.frame()
  dissimForAll <- data.frame()
  disSimForAttr <- numeric(0)
  NoOfObj <- nrow(myDataQuali)
  NoClolQuant <- ncol(myDataQuant)
  NoClolQuali <- ncol(myDataQuali)

  SignfDataFound <- signifOfQuantVars(myDataQuali, myDataQuant)

  for (a in 1:NoOfObj){
    for (b in 1:NoOfObj){
      for (i in 1:NoClolQuant){
        disSimForAttr[i] <- (SignfDataFound[i,2]*(myDataQuant[a,i] - myDataQuant[b,i]))^2
      }
      dissimForQuant[a,b] <- sum(disSimForAttr)
    }
  }

  disSimForAttr <- numeric(0)
  DistanceData <- distBetPairs(myDataQuali)

  for (a in 1:NoOfObj){
    for (b in 1:NoOfObj){
      for (i in 1:NoClolQuali){
        DistanceFound <- dplyr::filter(DistanceData, (DistanceData$xVal == as.character(myDataQuali[a,i]) & DistanceData$yVal == as.character(myDataQuali[b,i])) | (DistanceData$xVal == as.character(myDataQuali[b,i]) & DistanceData$yVal == as.character(myDataQuali[a,i])))

        if(nrow(DistanceFound) == 0){
          DistValUniq = 0
        }else
          DistValUniq <- DistanceFound$DistBetPair

        disSimForAttr[i] <- DistValUniq ^ 2
      }
      dissimForQuali[a,b] <- sum(disSimForAttr)
    }
  }
  dissimForAll <- dissimForQuant[,] + dissimForQuali[,]

  agnesClustering <- cluster::agnes(dissimForAll, diss = TRUE, method = "ward")

  return(dissimForAll)
}




