bws.dataset <- function(respondent.dataset, response.type = 1, choice.sets, design.type = 1, item.names = NULL, row.renames = TRUE)
{
# Name : bws.dataset
# Title: Creating a data set suitable for BWS analysis on the basis of counting and modeling approaches 
# Arguments:
#  respondent.dataset  a data frame containing a respondent data set
#  response.type       a value describing the format of response variables
#                        = 1 if the row number format is used
#                        = 2 if the item number format is used
#  choice.sets         a data frame or matrix containing choice sets
#  design.type         a value describing how to design the choice sets
#                        = 1 if a OMED is used
#                        = 2 if a BIBD is used
#  item.names          a vector containing the names of items
#  row.renames         a logical variable describing whether or not the row names of a data set created by this function are changed


# set variables, vectors, and matrices

  respondentDataset <- data.matrix(respondent.dataset)
  choicesets        <- data.matrix(choice.sets)
  numQuestions      <- nrow(choicesets)
  numRespondents    <- nrow(respondentDataset)

  if(design.type == 1)  # OMED
  {
    choicesets           <- choicesets - 1
    numItems             <- ncol(choicesets)
    frequencyItem        <- apply(choicesets, 2, table)[2, ]
    names(frequencyItem) <- c(1:length(frequencyItem))
  } 
  else                  # BIBD
  {
    numItems      <- length(table(choicesets))
    frequencyItem <- c(table(choicesets))
  }

  itemsInSet <- vector("list", numQuestions)
  if(design.type == 1)  # OMED
  {
    for(i in 1:numQuestions)
    {
      itemsInSet[[i]] <- which(choicesets[i, ] == 1)
    }
  }
  else                  # BIBD
  {
    for(i in 1:numQuestions)
    {
      itemsInSet[[i]] <- choicesets[i, ]
    }
  }

  numItemsInSet         <- sapply(itemsInSet, length)
  numPossiblePairsInSet <- numItemsInSet*(numItemsInSet - 1)


# expand respondent data set:
#  respondent data set is expanded to that in a format in which
#  each row shows a possible pair of the best and worst items

  # initial set
  expandRespondentDataset <- matrix(0, nrow = sum(numPossiblePairsInSet) * numRespondents,
                                       ncol = 5)
  m <- 0
  
  for(i in 1:numRespondents)
  {
    for(j in 1:numQuestions)
    {
      for(k in 1:numPossiblePairsInSet[j])
      {
        m <- m + 1
        # respondent's identification number
        expandRespondentDataset[m, 1]   <- respondentDataset[i, 1]
        # question number
        expandRespondentDataset[m, 2]   <- c(j)
        # number of possible pairs
        expandRespondentDataset[m, 3]   <- k
        # item numbers selected as best and worst
        expandRespondentDataset[m, 4:5] <- respondentDataset[i, c(j * 2, j * 2 + 1)]
        if(response.type == 1)
        {
          expandRespondentDataset[m, 4:5] <- itemsInSet[[j]][expandRespondentDataset[m, 4:5]]
        }
      }
    }
  }
  
  colnames(expandRespondentDataset) <- c("ID", "Q", "PAIR", "RES.B", "RES.W")

  
# make design matrix from choice sets

  # initial set
  designMatrix           <- matrix(0, nrow = sum(numPossiblePairsInSet), ncol = 4 + numItems)
  variableNames          <- paste("ITEM", 1:numItems, sep = "")
  colnames(designMatrix) <- c("Q", "PAIR", "BEST", "WORST", variableNames)
  lastRow                <- 0
  
  # create Q, PAIR, BEST, and WORST variables
  for(i in 1:numQuestions)
  {
    # create all combinations of items in choice set
    temp <- expand.grid(WORST = itemsInSet[[i]], BEST = itemsInSet[[i]])
    # exclude combinations of same item
    temp <- subset(temp, temp$BEST != temp$WORST)
    # combine Q and PAIR with possible pairs
    temp <- cbind(i, c(1:nrow(temp)), temp$BEST, temp$WORST)
    # store design matrix corresponding to i-th question in designMatrix
    designMatrix[(1 + lastRow):(lastRow + nrow(temp)), 1:4] <- temp
    lastRow <- lastRow + nrow(temp)
  }
  
  # assign values to ITEMj variables according to values of BEST and WORST:
  #  ITEMj = 1 if BEST = j; -1 if WORST = j; and 0 otherwise 
  for(i in 1:nrow(designMatrix))
  {
    designMatrix[i, c(designMatrix[i, 3] + 4, designMatrix[i, 4] + 4)] <- c(1, -1)
  }
  
  designMatrix <- as.data.frame(designMatrix)


# make and return data set for analysis

  # merge respondent data set with design matrix
  dataset <- merge(expandRespondentDataset, designMatrix, by = c("Q", "PAIR"))
  dataset <- dataset[order(dataset$ID, dataset$Q, dataset$PAIR), ]

  # create RES variable that is used as dependent (status) variable in clogit()
  TRUEorFALSE.B <- dataset$RES.B == dataset$BEST
  TRUEorFALSE.W <- dataset$RES.W == dataset$WORST
  dataset$RES   <- (TRUEorFALSE.B + TRUEorFALSE.W) == 2
  
  # create STR variable that is used as stratification variable in clogit()
  dataset$STR <- dataset$ID * 100 + dataset$Q
  
  # relabel item variables
  if(is.null(item.names) == FALSE)
  {
    colnames(dataset)[8:(7 + numItems)] <- item.names
  }

  # change row names
  if(row.renames == TRUE)
  {
    rownames(dataset) <- c(1:nrow(dataset))
  }

  # assign attributes to data set
  attributes(dataset)$nitems       <- numItems
  attributes(dataset)$nrespondents <- numRespondents
  attributes(dataset)$fitem        <- frequencyItem

  return(dataset)
}

