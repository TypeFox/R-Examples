
require("stringr")
library("foreign")
require("discretization")
require("functional")

#' mdl
#'
#' This function takes in a data frame to be discretized. The data type of the columns are important.
#' Only the numeric columns are discretized. 
#'
#' Here we use Fayyad's mdl discretization method.
#' Discretizing data by MDL method as implemented in the package 'discretization' 
#'
#' @param data data.frame. This data frame is discretized and returned.
#' @return  list of cut points and the discretized data frame 
#' @author saiteja ranuva
mdl <-  function(data){
    p <- length(data[1,])-1
    y <- data[,(p+1)]
    xd <- data
    for (i in 1:p){
      x <- data[,i]
      cuts1 <- cutPoints(x,y)
      cuts <- c(-Inf,cuts1,Inf)
      xd[,i] <- as.integer(cut(x,cuts,include.lowest = TRUE))
    }
    return (list(cutp=cuts,Disc.data=xd))
  }


#' cutX
#'
#' This function takes in a data frame to be discretized. The data type of the columns are important.
#' Only the numeric columns are discretized. 
#'
#' This uses the cut points generated while discretizing training data 
#' to discretize test data
#'
#' @param data data.frame. This data frame is discretized and returned.
#' @param cutp list - A list of cutp points obtained from training data
#' @return data data.frame. This the discretized data frame
#' @author saiteja ranuva


cutX <-
  function(data,cutp) {
    len <- length(data)
    y <- data[, len]
    discCol <- which(sapply(data, is.numeric))
    i=0
    cut <- list()
    for (col in discCol) {
      notNACol <- which(!is.na(data[, col]))
      if(length(notNACol)!=0){
        i=i+1        
        q <- as.integer(cut(data[, col][notNACol],cutp[[i]],include.lowest = TRUE))
        data[notNACol, col] <- q
      }
    }
    return(data)    
  }

 
#' discretizer
#'
#' This function takes in a data frame to be discretized. The data type of the columns are important.
#' Only the numeric columns are discretized. 
#'
#' Here we use Fayyad's mdl discretization method.
#' Discretizing data by MDL method as implemented in the package 'discretization' 
#'
#' @param data data.frame. This data frame is discretized and returned.
#' @return data data.frame. This the discretized data frame
#' @author saiteja ranuva
discretizer <- function(data) {
  df <- data
  len <- length(data)
  y <- data[, len]
  discCol <- which(sapply(data, is.numeric))
  i=0
  cut <- list()
  for (col in discCol) {
    notNACol <- which(!is.na(data[, col]))
    if(length(notNACol)!=0){
      i=i+1
      d <- data.frame(a = data[, col][notNACol], b = y[notNACol])
      q <- mdl(d)$Disc.data
      cut[[i]]<-mdl(d)$cutp
      data[notNACol, col] <- q[, 1]
    }
  }
  return(list(data=data,cutp=cut))
}


# Calcualtes all the required indexes for the AODE
#' fun1
#'
#' Calculates all the array and matrix indices required.
#'
#' Details - add later
#'
#' @param aode list. this is the list which has all the required variables in it.
#' @return a list 
#' @author sai teja ranuva
#' @seealso \code{\link{setVar}}
indexCalc <- function(aode) {
  featureAttrValues <- aode[["featureAttrValues"]]
  featureClasses <- aode[["featureClasses"]]
  indexI <- 0
  nOfDiscreteAttrValues <- aode[['nOfDiscreteAttrValues']]
  
  # def:attrIndex - gives each attribute an index(considering one memory location 
  #  for each attr value
  # Eg: attr1 - {1,2,3} attr2-{a,b} attr3-{1,2,3,4} attr4-{1,2}
  # then attrIndex for this will be 
  # [0,3,5,9]
  attrIndex <- vector(length = aode[["nOfAttr"]])
  
  # def:reverseAttrIndex - for each attr it stores how many more attr values exist
  # after it
  # Eg: attr1 - {1,2,3} attr2-{a,b} attr3-{1,2,3,4} attr4-{1,2}
  # reverseAttrIndex will be [8,6,2,0]
  reverseAttrIndex <- vector(length = aode[["nOfAttr"]])
   
  
  rIndexJ <- 0  
  
  #local variables for some use
  tempAttrSum <- 0   
  p <- 0
  
  # calculating reverseAttrIndex
  for (attrI in 1:aode[["nOfAttr"]]) {
    attrIndex[attrI] <- indexI
    p <- nOfDiscreteAttrValues[attrI]
    indexI <- indexI + p
    n <- length(nOfDiscreteAttrValues)
    if (attrI <= (aode[["nOfAttr"]] - 1)) {
      tempAttrSum <- sum(nOfDiscreteAttrValues[(attrI + 1):n])
    }
    reverseAttrIndex[attrI] <- tempAttrSum
  }
  
  # def:totalAttrValues - stores the total count of attr values
  totalAttrValues <- sum(nOfDiscreteAttrValues)
  
  
  temp <- 0
  
  # def:reverseCumulativeAttrIndex - takes a lot of words to explain exactly 
  # it is simply an index to get to the correct memory location of counts 
  # and prob matrices
  # Eg: attr1 - {1,2,3} attr2-{a,b} attr3-{1,2,3,4} attr4-{1,2}
  # [0,24,36,44]   
  reverseCumulativeAttrIndex <- vector(length = length(reverseAttrIndex))
  p <- 0
  q <- 0
  r <- 0
  for (index in 1:length(reverseAttrIndex)) {
    if (index == 1) {
      r <- 0
    } else {
      r <- reverseAttrIndex[index - 1]
    }
    reverseCumulativeAttrIndex[index] <- r * p + q
    p <- length(featureAttrValues[[index]])
    q <- reverseCumulativeAttrIndex[index]
  }
  
  # return as a list
  return(list(attrIndex = attrIndex, reverseAttrIndex = reverseAttrIndex, reverseCumulativeAttrIndex = reverseCumulativeAttrIndex))
}

#' setVar
#'
#' sets the required variables 
#'
#' calculates the required space for indices and allocates
#'
#' @param aode list. this is the list which has all the required variables in it.
#' @return aode list updated value of the list
#' @author sai teja ranuva
setVar <- function(aode) {
  
  #self explanatory
  data <- aode[["train"]]
  index <- indexCalc(aode)
  aode[["attrIndex"]] <- index$attrIndex
  aode[["reverseAttrIndex"]] <- index$reverseAttrIndex
  aode[["reverseCumulativeAttrIndex"]] <- index$reverseCumulativeAttrIndex
  
  
  # def:nOfNA - keeps the count of missing values for each attribute
  aode[["nOfNA"]] <- apply(data, 2,Compose(is.na, sum))
  #aode[['nOfNA']] <- rep(0,aode[['nOfAttr']])
  p <- matrix(0,aode[['totalClasses']],aode[["nOfAttr"]])
  y <- 1:aode[['totalClasses']]
  for(y in 1:aode[['totalClasses']]){
    for (xi in 1:aode[["nOfAttr"]]) {
      q <- !(data[,aode[["nOfAttr"]]+1]!=y)
      r <- is.na(data[, xi])
      p[y,xi] <- sum(q & r)
    }
  }    
  aode[["nOfNAyxi"]] <- p
  
  p <- vector(length = length(aode[["featureAttrValue"]]))
  for (xi in 1:aode[["nOfAttr"]]) {
    q <- is.na(data[, aode[["nOfAttr"]] + 1])
    r <- is.na(data[, xi])
    p[xi] <- sum(q | r)
  }
  aode[["nOfNAyxi"]] <- p
  
  m <- aode[['train']] 
  c <- aode[['nOfAttr']] +1
  aode[['cyxjpNA']] <- array(rep(0,aode[['nOfAttr']]*aode[['totalClasses']]*aode[['totalAttrValues']]),dim=c(aode[['nOfAttr']],aode[['totalClasses']],aode[['totalAttrValues']]))
  for(xi in 1:aode[["nOfAttr"]]){
    for(y in 1:aode[['totalClasses']]){
      for (xj in 1:aode[["nOfAttr"]]) {
        if(xj==xi) next
        for(vxj in 1:aode[['nOfDiscreteAttrValues']][xj]){
          yxj <- m[m[, c] == aode[["featureClasses"]][y] & m[, xj] == aode[['attrValues']][vxj] & is.na(m[, xi]),]
          # if yxi is a vector (happens if yxi has only one row), 
          # convert it into matrix
          if (is.vector(yxj))  
            yxi = matrix(yxi, nrow = 1)
          yxj <- yxj[rowSums(is.na(yxj)) != (aode[["nOfAttr"]] + 1), ]
          if (is.vector(yxj)) 
            yxj = matrix(yxj, nrow = 1)
          if(is.null(yxj))
            n=0
          else
            n <- nrow(yxj)
          aode[['cyxjpNA']][xi,y,aode[['atrIndex']][xj]+vxj] <- n
        }
      }
    }
  }
  #print(aode[['cyxjpNA']])
  # def:nOfNAyxi - keeps the count of missing of either (class or an attribute)
  
  
  
  
  # Allocate space for counts of count(y,xi)
  aode[["c_pyxi"]] <- matrix(0, nrow = aode[["totalClasses"]], ncol = aode[["totalAttrValues"]], 
                             byrow = TRUE)
  
  # Allocating space for count(y,xi,xj)
  aode[["c_pyxixj"]] <- matrix(0, nrow = aode[["totalClasses"]], ncol = (aode[["reverseCumulativeAttrIndex"]][length(aode[["reverseCumulativeAttrIndex"]])]), 
                               byrow = TRUE)
  
  # if subsumption flag is set 
  if (aode[["s"]]) {
    # allocate space for count(xi,xj)
    aode[["cxixj"]] <- rep(0, aode[["reverseCumulativeAttrIndex"]][length(aode[["reverseCumulativeAttrIndex"]])])
  }
  
  # allocate space for count(xi)
  aode[["cxi"]] <- rep(0, aode[["totalAttrValues"]])
  return(aode)
}



#' training
#'
#' Calculates the count matrices.
#'
#' class parent count and class parent child count are calculated
#'
#' @param aode list. this is the list which has all the required variables in it.
#' @return aode list. updated value
#' @author sai teja ranuva 
training <- function(aode) {
  c <- aode[["nOfAttr"]] + 1  # column number of class
  m <- aode[["train"]]        #matrix of train data   
  for (y in 1:length(aode[["featureClasses"]])) {        
    for (xi in 1:(aode[["nOfAttr"]])) {
      if (aode[["nOfDiscreteAttrValues"]][xi] == 0) 
        next
      attrValues <- aode[["featureAttrValues"]][[xi]]
      for (vxi in 1:length(attrValues)) {
        yxi <- m[m[, c] == aode[["featureClasses"]][y] & m[, xi] == attrValues[vxi],]
        
        # if yxi is a vector (happens if yxi has only one row), 
        # convert it into matrix
        if (is.vector(yxi))  
          yxi = matrix(yxi, nrow = 1)
        yxi <- yxi[rowSums(is.na(yxi)) != (aode[["nOfAttr"]] + 1), ]
        if (is.vector(yxi)) 
          yxi = matrix(yxi, nrow = 1)
        
        # temp variable to store a combo index
        ind <- aode[["reverseCumulativeAttrIndex"]][xi] + (vxi - 1) * aode[["reverseAttrIndex"]][xi]
        
        # attempt to reduce unnecessary looping , could be useless
        if (is.null(yxi)) {
          aode[["c_pyxi"]][y, aode[["attrIndex"]][xi] + vxi] = 0
          for (xj in (xi + 1):(aode[["nOfAttr"]]) - 1) {
            if (xj <= xi) 
              break
            p <- length(aode[["featureAttrValues"]][[xj]])
            for (vxj in 1:p) aode[["c_pyxixj"]][y, ind + aode[["attrIndex"]][xj] + 
                                                  vxj] <- 0
          }
          next
        }
        
        n <- as.integer(nrow(yxi))
        # update counts                
        aode[["cxi"]][aode[["attrIndex"]][xi] + vxi] <- aode[["cxi"]][aode[["attrIndex"]][xi] + 
                                                                        vxi] + n
        aode[["c_pyxi"]][y, aode[["attrIndex"]][xi] + vxi] <- n
        
        
        for (xj in (xi + 1):(aode[["nOfAttr"]])) {
          if (aode[["nOfDiscreteAttrValues"]][xi] == 0) 
            next
          if (xj > aode[["nOfAttr"]]) 
            break
          attrValuesJ <- aode[["featureAttrValues"]][[xj]]
          if (xj != xi) {
            p <- 0
            p <- length(aode[["featureAttrValues"]][[xj]])
            for (vxj in 1:p) {
              yxixj <- yxi[yxi[, xj] == attrValuesJ[vxj], ]
              if (is.vector(yxixj)) 
                yxixj <- matrix(yxixj, nrow = 1)
              yxixj <- yxixj[rowSums(is.na(yxixj)) != (aode[["nOfAttr"]] + 1), ]
              if (is.vector(yxixj)) 
                yxixj <- matrix(yxixj, nrow = 1)
              if(is.null(yxixj))
                n=0
              else
                n <- nrow(yxixj)
              
              # update counts
              aode[["cxixj"]][ind + aode[["attrIndex"]][xj] - aode[["attrIndex"]][xi + 
                                                                                    1] + vxj] = aode[["cxixj"]][ind + aode[["attrIndex"]][xj] - aode[["attrIndex"]][xi + 
                                                                                                                                                                      1] + vxj] + n
              aode[["c_pyxixj"]][y, ind + aode[["attrIndex"]][xj] - aode[["attrIndex"]][xi + 
                                                                                          1] + vxj] <- as.integer(n)
              
            }
            
          }
        }
      }
    }
  }
  return(aode)
}


#' calWeight
#'
#' This function calculates mutual information between the attr and class.
#'
#' This function is called when weighted flag is set
#'
#' @param aode list. this is the list which has all the required variables in it.
#' @return aode list. updated value
#' @author saiteja ranuva
calWeight <- function(aode) {
  featureAttrValues <- aode[["featureAttrValues"]]
  totalClasses <- aode[["totalClasses"]]
  nOfFeatureVectors <- dim(aode[["train"]])[1]
  mestimate <- aode[["m"]]
  train <- aode[["train"]]
  pxi <- vector(length = aode[["totalAttrValues"]])
  py <- vector(length = aode[["totalClasses"]])
  w <- vector(length = length(featureAttrValues))
  
  for (y in 1:totalClasses) {
    n<-0
    if(is.null(train[train[, dim(train)[2]] == y, ]))
      n <- 0 
    else {
      if(is.vector(train[train[, dim(train)[2]] == y, ]))
        n=1
      else
        n = nrow(train[train[, dim(train)[2]] == y, ])
    }
    py[y] <- (n + (mestimate/aode[["totalClasses"]])) /(dim(train)[1] - sum(is.na(train[, dim(train)[2]])) + mestimate)
    
  }
  for (xi in 1:aode[["nOfAttr"]]) {
    vxi <- seq(1:aode[['nOfDiscreteAttrValues']][xi])
    pxi[aode[["attrIndex"]][xi] + vxi] <- (aode[['cxi']][aode[['attrIndex']][xi]+vxi] + (mestimate/length(featureAttrValues[[xi]]))) * 
      1/((dim(train)[1]) + mestimate)
  }
  nonZFlag <- TRUE  # to check if not all weights are 0, if so we give equal weight to all attr
  for (xi in 1:aode[["nOfAttr"]]) {
    sum <- 0
    for (y in 1:totalClasses) {
      for (vxi in 1:aode[["nOfDiscreteAttrValues"]][xi]) {
        p_yxi <- aode[["pyxi"]][y, aode[["attrIndex"]][xi] + vxi]
        sum <- sum + p_yxi * log((p_yxi/(py[y] * pxi[aode[["attrIndex"]][xi] + vxi])), 2)
      }
    }
    if (sum != 0) 
      nonZFlag <- FALSE
    w[xi] <- sum
  }
  # returning the weights and the non zero flag
  return(list(w = w, nonZFlag = nonZFlag))
}

##' 
##'
##' 
##' @title aode
##' 
##' @details This is the training phase of the algorithm. Necessary count and probability tables are generated which will used for the prediction purpose.
##' 
##' @description This function builds the model using the AODE algorithm which can then be used classification.
##' 
##' @param train data.frame : training data. It should be a data frame. AODE works only discretized data. It would be better to discreetize the data frame before passing it to this function.However, aode discretizes the data if not done before hand. It uses an R package called discretization for the purpose. It uses the well known MDL discretization technique.(It might fail sometimes) 
##' @param mestimate optional numeric
##' @param weighted optional boolean
##' @param subsumption optional boolean
##' @param S optional numeric subsumption constant
##' @return An object of class AODE
##' @examples require("datasets")
##' aode(iris,mestimate=1) 
##' aode(iris)
##' aode(iris,weighted=TRUE)
##' @author saiteja ranuva
##' @export
aode <- function(train, mestimate = 1, weighted = FALSE, 
                 subsumption = FALSE, S = 100) {
  
  discCol <- which(sapply(train[, 1:(dim(train)[2] - 1)], is.numeric))
  featureNamesList <- names(train)
  ntrain <- dim(train)[1]
  x <- discretizer(train)
  train <- x$data
  cutp <- x$cutp
  nOfAttr <- dim(train)[2] - 1
  train <- data.matrix(train)
  
  class <- factor(train[,nOfAttr+1])
  totalClasses <- length(levels(class))
  featureClasses <- as.integer(levels(class))
  
  featureAttrValues <- apply(train,2,Compose(factor,levels,as.integer))
  nOfDiscreteAttrValues <- sapply(featureAttrValues, length)
  nOfFeatureVectors <- dim(train)[1]
  bayes <- list(train = train, cutp=cutp ,m = mestimate, w = weighted, s = subsumption, 
                S = S, nOfFeatureVectors = nOfFeatureVectors, featureNamesList = featureNamesList, 
                nOfAttr = nOfAttr, featureAttrValues = featureAttrValues, featureClasses = featureClasses, 
                attrIndex = NULL, reverseAttrIndex = NULL, reverseCumulativeAttrIndex = NULL, nOfDiscreteAttrValues = nOfDiscreteAttrValues, 
                totalClasses = totalClasses, totalAttrValues = sum(nOfDiscreteAttrValues), c_pyxi = NULL, c_pyxixj = NULL, 
                nOfNA = NULL, nOfNAyxi = NULL, cxi = NULL, ind = NULL, cpyxi = NULL, cxixj = NULL, 
                pyxi = NULL, pyxixj = NULL, pyxixjopp = NULL, weight = NULL, probInitAODE = 1, p = NULL, cyxjpNA=NULL)
  bayes[["p"]] <- vector()
  for (i in 1:bayes[["nOfAttr"]]) {
    if (length(bayes[["featureAttrValues"]][[i]]) == 0) {
      bayes[["p"]] <- c(bayes[["p"]], i)
    }
  }
  bayes <- setVar(bayes)
  t = system.time({
    bayes <- training(bayes)
  })
  
  # def:ind - combo index
  bayes[["ind"]] <- rep(0, bayes[["totalAttrValues"]])
  
  # def:pyxi - probablity of y=y,xi=xi
  bayes[["pyxi"]] <- matrix(0, nrow = bayes[["totalClasses"]], ncol = bayes[["totalAttrValues"]], 
                            byrow = TRUE)
  
  # def:pyxixj - probablity of y=y,xi=xi,xj=xj , xi<xj
  bayes[["pyxixj"]] <- matrix(0, nrow = bayes[["totalClasses"]], ncol = (bayes[["reverseCumulativeAttrIndex"]][length(bayes[["reverseCumulativeAttrIndex"]])]), 
                              byrow = TRUE)
  
  # def:pyxixjopp - probablity of y=y.xi=xi,xj=xj , xj<xi
  bayes[["pyxixjopp"]] <- matrix(0, nrow = bayes[["totalClasses"]], ncol = (bayes[["reverseCumulativeAttrIndex"]][length(bayes[["reverseCumulativeAttrIndex"]])]), 
                                 byrow = TRUE)
  
  #calculating prob from counts 
  for (y in 1:bayes[["totalClasses"]]) {
    for (xi in 1:bayes[["nOfAttr"]]) {
      if (bayes[["nOfDiscreteAttrValues"]][xi] == 0) 
        next
      for (vxi in 1:bayes[["nOfDiscreteAttrValues"]][xi]) {
        if (xi < bayes[["nOfAttr"]]) 
          bayes[["ind"]][bayes[["attrIndex"]][xi] + vxi] <- bayes[["reverseCumulativeAttrIndex"]][xi] + 
          (as.numeric(vxi) - 1) * bayes[["reverseAttrIndex"]][xi] - bayes[["attrIndex"]][xi + 
                                                                                           1]
        temp <- bayes[["c_pyxi"]][y, bayes[["attrIndex"]][xi] + vxi] + bayes[["m"]]/(bayes[["totalClasses"]] * 
                                                                                       bayes[["nOfDiscreteAttrValues"]][xi])
        bayes[["pyxi"]][y, bayes[["attrIndex"]][xi] + vxi] <- temp/(bayes[["nOfFeatureVectors"]] + 
                                                                      bayes[["m"]] - bayes[["nOfNAyxi"]][xi])
        for (xj in (xi + 1):(bayes[["nOfAttr"]])) {
          if (xj > bayes[["nOfAttr"]]) 
            break
          if (bayes[["nOfDiscreteAttrValues"]][xj] == 0) 
            next
          if (xj > bayes[["nOfAttr"]]) 
            break
          
          p <- 0
          p <- bayes[["nOfDiscreteAttrValues"]][xj]
          for (vxj in 1:p) {
            k <- bayes[["ind"]][bayes[["attrIndex"]][xi] + vxi] + bayes[["attrIndex"]][xj] + 
              vxj
            num <- bayes[["c_pyxixj"]][y, k] + (bayes[["m"]]/bayes[["nOfDiscreteAttrValues"]][xj])
            
            #(C(xp,y) - ?(xc)) + m
            #den <- (bayes[["c_pyxi"]][y, bayes[["attrIndex"]][xi] + vxi] - bayes[["nOfNA"]][xj]) + bayes[["m"]]
            den <- (bayes[["c_pyxi"]][y, bayes[["attrIndex"]][xi] + vxi] - bayes[["cyxjpNA"]][xi,y,bayes[['attrIndex']][xj]+vxj]) + bayes[["m"]]
            if(den==0 || is.infinite(num)){
              print('this is a bug, please report it - with the data set if possible')
            }
            bayes[["pyxixj"]][y, k] <- num/den
            num <- bayes[["c_pyxixj"]][y, k] + (bayes[["m"]]/bayes[["nOfDiscreteAttrValues"]][xi])
            den <- bayes[["c_pyxi"]][y, bayes[["attrIndex"]][xj] + vxj] + bayes[["m"]] - bayes[["cyxjpNA"]][xi,y,bayes[['attrIndex']][xj]+vxj]
            bayes[["pyxixjopp"]][y, k] <- num/den
          }
        }
      }
    }
  }
  
  # def:weight - vector containing mutual information of attr and class
  
  bayes[["weight"]] <- rep(1, bayes[["nOfAttr"]]) # initializing the weights with 1
  if (weighted) {
    x <- calWeight(bayes)
    if (x[["nonZFlag"]]) 
      bayes[["weight"]] <- rep(1/bayes[["totalClasses"]], bayes[["nOfAttr"]]) else bayes[["weight"]] <- x[["w"]]
  }
  
  
  obj <- list()
  attr(obj,"aode") <- bayes
  class(obj) <- "AODE"
  
  return(obj)
  
}   


#' @title predict
#' @description  This is a generic function. This function predicts the class of the test data and returns a vector of predicted values. 
#' @details Written in line with the E1071 package. 
#'
#' @method predict AODE 
#' @S3method predict AODE 
#' @param object object of class AODE
#' @param test test data frame. If the training data was discretized, then the same cut points shall be used to discretize the test data. So obviously if the training was not discretized, test data should also not be discretized.
#' @param ... extra arguments which might be needed in future
#' @return class vector containing the predicted class distribution of the test data
#' @examples data<-iris 
#' ode<-aode(data)
#' predict(ode,iris)
#' @author sai teja ranuva
#' @export
predict.AODE <- function(object,test,...){
  aode <- attr(object,"aode")
  test <- cutX(test,aode[['cutp']])
  test <- data.matrix(test)
  n <- aode[["nOfAttr"]] + 1
  y <- test[, n]
  x <- test[, 1:(n - 1)]
  class <- apply(x, 1, distributionForInstance, aode = aode)
  return(class)
}



#' distributionForInstance
#'
#' predicts class of a given instance based on the model
#'
#' details to be added
#'
#' @param x instance to be classified
#' @param aode list. this is the list which has all the required variables in it.
#' @return class integer predicted class of the instance
#' @author sai teja ranuva
distributionForInstance <- function(x, aode) {
  Aode <- rep(0, aode[["totalClasses"]])
  classVal <- aode[["nOfAttr"]] + 1
  ode <- matrix(0, nrow = aode[["totalClasses"]], ncol = aode[["nOfAttr"]], byrow = TRUE)
  
  #' for subsumption resolution
  specialGenArr <- rep(-1, aode[["nOfAttr"]])
  
  y <- x[classVal]
  x <- x[1:aode[["nOfAttr"]]]
  x[aode[["p"]]] <- NA
  #' if subsumption flag is set, do subsumntion resolution
  if (aode[["s"]]) {
    for (l in 1:aode[["nOfAttr"]]) {
      if (is.na(x[l])) 
        next
      for (m in 1:aode[["nOfAttr"]]) {
        if (is.na(x[m]) || l == m || specialGenArr[m] == l) 
          next
        if (l < m) 
          count <- aode[["cxixj"]][aode[["reverseCumulativeAttrIndex"]][l] + (x[l] - 
                                                                                1) * aode[["reverseAttrIndex"]][l] + aode[["attrIndex"]][m] - aode[["attrIndex"]][l + 
                                                                                                                                                                    1] + x[m]] else count <- aode[["cxixj"]][aode[["reverseCumulativeAttrIndex"]][m] + (x[m] - 
                                                                                                                                                                                                                                                          1) * aode[["reverseAttrIndex"]][m] + aode[["attrIndex"]][l] - aode[["attrIndex"]][m + 
                                                                                                                                                                                                                                                                                                                                              1] + x[l]]
        if (aode[["cxi"]][aode[["attrIndex"]][m] + x[m]] > aode[["S"]]) {
          if (aode[["cxi"]][aode[["attrIndex"]][m] + x[m]] == count) {
            if ((aode[["cxi"]][aode[["attrIndex"]][m] + x[m]] == aode[["cxi"]][aode[["attrIndex"]][l] + 
                                                                                 x[l]]) && (l < m)) {
              next
            } else {
              specialGenArr[l] = m
              break
            }
          }
        }
      }
    }
  }
  #' neglect xi's in the subsumption vector, initializing to NA => ignoring 
  x[which(specialGenArr != -1)] <- NA
  
  #' for convenience
  t <- aode[["attrIndex"]]
  y <- seq(1, aode[["totalClasses"]])
  n <- as.integer(aode[["nOfAttr"]])
  xi <- which(!is.na(x))
  ode[y, xi] <- aode[["pyxi"]][y, aode[["attrIndex"]][xi] + x[xi]] * aode[["weight"]][xi] * 
    aode[["probInitAODE"]]
  #' just for convenience
  ind <- aode[["ind"]]
  #' to track which of x attr values are NA(missing)
  r <- which(is.na(x))
  #' to vectorize the calculations
  y <- seq(1, aode[["totalClasses"]])
  
  #' using sapply instead of for loop here makes the code more efficient. but a bug is getting 
  #' introduced. have to look into it.
  for (xi in 1:(aode[["nOfAttr"]] - 1)) {
    if (is.na(x[xi]) || aode[["nOfDiscreteAttrValues"]][xi] == 0) 
      next
    p <- aode[["attrIndex"]][xi] + x[xi]
    
    #' sloppy code , need to replace with efficient one
    xj <- (xi + 1):aode[["nOfAttr"]]
    xj <- xj[!xj %in% aode[["p"]]]
    xj <- xj[!xj %in% r]
    
    vec <- ind[p] + aode[["attrIndex"]][xj] + x[xj]
    ode[y, xj] <- ode[y, xj] * aode[["pyxixjopp"]][y, ind[p] + aode[["attrIndex"]][xj] + 
                                                     x[xj]]
    
    if (length(vec) == 1) 
      mat <- matrix(aode[["pyxixj"]][y, vec]) else mat <- aode[["pyxixj"]][y, vec]
    ode[y, xi] <- ode[y, xi] * apply(mat, 1, prod)
  }
  
  for (y in 1:aode[["totalClasses"]]) {
    Aode[y] <- sum(ode[y, ])/sum(aode[["weight"]])
  }
  
  class <- which(Aode - max(Aode) == 0)[1]
  return(class)
}

