# 
# Copyright (c) 2010, 2014, IBM Corp. All rights reserved. 
#     
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
#
plot.idaLm <- function(x, names = TRUE, max_forw = 50, max_plot = 15, order = NULL,
                     lmgON = FALSE, backwardON = FALSE, ...){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot idaLm objects.")
  }
  
  idaLm <- x
  if(!("idaLm" %in% class(idaLm))){
    stop("Input is not an idaLm object.")
  }
  if(is.null(idaLm$CovMat)){
    stop("Function needs an idaLm object with a Covariance Matrix")
  }
  if(is.numeric(max_forw)){
    if(!(max_forw%%1 == 0 && max_forw > 1)){
      stop("max_forw must be an integer bigger than 1.")
    }    
  }else{
    stop("max_forw must be numeric.")
  }
  if(is.numeric(max_plot)){
    if(!(max_plot%%1 == 0 && max_plot > 0)){
      stop("max_plot must be an integer bigger than 0.")
    }    
  }else{
    stop("max_plot must be numeric.")
  }
  if(!is.logical(names)){
    stop("The parameter \"names\" must be logical.")
  }
  if(!is.logical(lmgON)){
    stop("The parameter \"lmgON\" must be logical.")
  }
  if(!is.logical(backwardON)){
    stop("The parameter \"backwardON\" must be logical.")
  }
  
  ######Get Covariance Matrix##########
  mat <- idaLm$CovMat
  nrowMat <- nrow(mat)
  if(rownames(mat)[nrowMat] != "Intercept"){
    stop("Plotting a model without intercept is currently not supported.")
  }else{  
   n = mat[nrowMat, nrowMat]
  }
  
  card <- idaLm$card
  colNames <- names(card)
  
  if(is.vector(order)){
    check <- order %in% colNames
    if(all(check)){
      #by adding up the cardinalities of the columns we get the index of the last column for every
      #predictor. +1 because of the target variable in the first row/column
      lastcol <- cumsum(card)[!(colNames %in% order)] + 1
      delcard <- card[!(colNames %in% order)]
      #we use the ":" function to create a vector for all indices that belong to one predictor
      deletecols <- unlist(mapply(FUN = function(x, y){do.call(":",list(x-y+1,x))}, lastcol, delcard))
        
      card <- card[colNames %in% order]
      mat <- mat[-deletecols, -deletecols]
     
      colNames <- names(card)
    }else{
      error <- "Some of the Attributes of the vector are not in the model:"
      error <- paste(error, paste(order[!check], collapse = ","))
      stop(error)
    }
  }
  
  YTY = mat[1, 1]
  XTX = mat[-1, -1]
  XTY = mat[1, -1]
  
  nomCols <- names(card)[card > 1]

#################Function that calculates relative Importance################
  getRelImp <- function(colNames, card, nomCols, numrow, XTX, XTY, YTY, order, lmgON,
                        backwardON){
    #This function gives us more Information about the impact of the attributes on the
    #linear model. It Calculates the values "Usefulness","First","Forward","Backward","Order"
    #See in each part of the programm for further explanation.
    nrowXTX <- nrow(XTX)          
    mean_y  <- XTY[nrowXTX] / numrow #mean value of y
    Var_y   <- YTY - mean_y * mean_y * numrow #sum (y-mean_y)^2, Variance of y
    
    #Create new Matrix from the information of the covariance matrix
    SXY <- matrix(nrow = 1, ncol = nrowXTX)
    
    SXY[1, ] <- XTY - XTX[nrowXTX, ] * mean_y #sum(x_i*y-mean_x_i*mean_y)=sum(x_i-mean_x_i)*(y-mean_y)
    
    SXY[1, nrowXTX] <- 0 #Intercept is constant so this is always zero
    SXY <- SXY / Var_y    #Here we divide by Var_y
        
    #R2 gives the value of the model and is between 0(bad) and 1(good).
    #Calculating R2 is [sum ((predicted yvalue) - mean_y)^2 ]/[sum (y-mean_y)^2]
    #Calculation formula with the Covariance matrix. (XTX)^-1 %*% XTY are the coefficients
    #of the linear model. R2=SXY%*%(XTX)^-1%*%XTY
    calc <- function(SXY, XTX, XTY){
      coeff <- NULL
      try({coeff = solve(XTX, XTY)}, silent=T)#trying to solve directly
      if(is.null(coeff)) {
        coeff <- ginv(XTX, tol = 1e-16) %*% XTY#when XTX not invertible. Calculating a pseudo-inverse
      }
      SXY %*% coeff #This is R2. We will also be calling this method for subsets of
    }             #the model.
    
    totalR2 <- calc(SXY, XTX, XTY)#Value of the Model with all Predictors given by the Model
    
    #trivial case of one predictor
    if(length(colNames) < 2){
      
      if(is.vector(order)){
          res <- data.frame(totalR2)
          names(res) <- c('Model_Values')
          rownames(res) <- order
      }else{     
        if(backwardON){
          res <- data.frame(0)
          names(res) <- c('Backward_Values')
          rownames(res) <- colNames
        }else{
          res <- data.frame(totalR2)
          names(res) <- c('Forward_Values')
          rownames(res) <- colNames
        }
      }
      
      if(lmgON){
        res2 <- data.frame(res, totalR2, totalR2, totalR2)
        names(res2) <- c(names(res), 'Usefulness', 'First', 'LMG')
      }else{
        res2 <- data.frame(res, totalR2, totalR2)
        names(res2) <- c(names(res), 'Usefulness', 'First')
      }      
      res2           
    }else{
      #Saves the position of the attributes in the covariance matrix
      posofAttr <- list()
      R2 <- c()
      First <- c()
      disposition <- 0

      for(i in 1:length(colNames)){         
            length <- card[colNames[i]]-1
            tempPos <- (i+disposition) : (i+disposition+length)
            posofAttr[[colNames[i]]] <- tempPos
            
            lastXTX <- XTX[-tempPos, -tempPos]
            lastXTY <- XTY[-tempPos]
            lastSXY <- SXY[1, -tempPos]
            R2 <- c(R2, calc(lastSXY, lastXTX, lastXTY))
            
            firstXTX <- XTX[c(tempPos, nrowXTX), c(tempPos, nrowXTX)]
            firstXTY <- XTY[c(tempPos, nrowXTX)]
            firstSXY <- SXY[1, c(tempPos, nrowXTX)]            
            First <- c(First, calc(firstSXY, firstXTX, firstXTY))
            
            disposition <- disposition + length
      }

      names(First) <- colNames
      Usefulness <- totalR2-R2 #Compare how much information/improvement each variable has
      names(R2) <- colNames
      #This function calculates forward or backward depending on direction.
      #It will find in each step the attribute that improves the model best(forward)
      #or the attribute that can be taken away with minimum value reduction(backward)
      stepbystep <- function(StartingVal, posofAttr, SXY, XTX, XTY, max_forw,
                             direction = "forward"){
        
        tempcolNames <- names(StartingVal)
        #switch the best attribute of the first step infront
        bestpos <- match(max(StartingVal), StartingVal)[1]
        swap <- tempcolNames[1]
        tempcolNames[1] <- names(StartingVal)[bestpos]
        tempcolNames[bestpos] <- swap    
                
        indices <- posofAttr[[tempcolNames[1]]]
        
        #Adds the Intercept line to the Subset of the covariance Matrix we are looking at
        if(direction == "forward"){
          indices <- c(indices,nrowXTX)
        }
        
        #Saves the value of the first attribute
        modelval <- StartingVal[bestpos]
        max_forw <- min(length(tempcolNames), max_forw)
        #starting at 2 cause StartingVal gave us the best attribute in the first step
        for(i in 2:max_forw){
          #tells us the value of the model with one additional/fewer attribute so we can find the best
          tempvalues <- c()
          for(j in i:length(tempcolNames)){
            #to create the model with one additional/fewer attribute we need to
            #add the position of this attribute to our Submodel
            tempindices <- indices
            tempindices <- c(tempindices,posofAttr[[tempcolNames[j]]])
            
            #Creates the covariance matrices of the model from the original cov.matrix
            if(direction == "forward"){
              tempXTX <- XTX[tempindices, tempindices]
              tempSXY <- SXY[tempindices]
              tempXTY <- XTY[tempindices]
            }else{
              tempXTX <- XTX[-tempindices, -tempindices]
              tempSXY <- SXY[-tempindices]
              tempXTY <- XTY[-tempindices]
            }
            #Calculation of the Submodel-Value
            tempvalues <- c(tempvalues, calc(tempSXY, tempXTX, tempXTY))
          }
          names(tempvalues) <- tempcolNames[i:length(tempcolNames)]
          
          #finds the best model with one additional/fewer attribute
          modelval[i] <- max(tempvalues)
          bestpos <- match(modelval[i], tempvalues)[1]
          
          #switch the best attribute of this step at the ith position. So the first i
          #are in the order of the heuristic
          swap <- tempcolNames[i]
          tempcolNames[i] <- tempcolNames[i+bestpos-1]
          tempcolNames[i+bestpos-1] <- swap    
          
          #change indices to add/take away the best attribute of that step to/from our model
          indices <- c(indices, posofAttr[[tempcolNames[i]]])
        }
        data.frame(tempcolNames[1:max_forw], modelval)
      } 
      
      
      ###LMG### This function averages above all the improvements you get when adding
      #the predictor to each possible subset of the Model. Note that this function
      #has EXPONENTIAL runningtime over the number of attributes!
      lmg <- function(posofAttr, First, SXY, XTX, XTY){
        colNames <- names(First)
        
        #Create covariance-matrix for the subset to calculate R^2 by using the calc method.
        calc.index <- function(indices,SXY,XTX,XTY){
          indices <- c(indices, nrow(XTX)) #Need to add the Intercept row
          tempXTX <- XTX[indices, indices]
          tempSXY <- SXY[indices]
          tempXTY <- XTY[indices]
          calc(tempSXY, tempXTX, tempXTY)
        }
        
        #Adds up the Improvements for adding it First
        NrAttr <- length(colNames)
        lmgRes <- factorial(NrAttr-1)*First
        
        #dataold is the matrix with all subsets with k-1 elements
        #We model a set as an orderd vector.
        
        dataold <- t(matrix(1:NrAttr))
        
        #valold[j] gives the Value of the Model with the Predictors given by dataold[, j]
        valold <- First
        
        #iterating over the cardinality of the subsets we are looking at
        for(cardsubset in 2:NrAttr){
          nrofSets <- choose(NrAttr, cardsubset)
          valnew <- vector(mode = "integer", length =  nrofSets)
          
          datanew <- matrix(, ncol = nrofSets, nrow = cardsubset)
          count <- 0 #counts the numbers of new subsets created for the next step. Indexing
          
          #iterating over the predictors
          for(NrPred in 1:NrAttr){
            
            #iterating over the subsets with cardsubset-1 elements and adding the Nrpred-th Predictor
            for(NrSubset in 1:ncol(dataold)){
              
              #if the predictor is in the subset there is no sense in adding it
              tmpSet <- dataold[, NrSubset]
              if(!(NrPred %in% tmpSet)){
                #Calculates the Value (R^2) of the Model with the Subset dataold[, NrSubset] added by the colNames[NrPred] Predictor
                newSet <- c(tmpSet, NrPred)
                temp <- calc.index(unlist(posofAttr[newSet]), SXY, XTX, XTY)
                
                #prevents us from creating the same set twice. Only if its orderd we will create it
                if(NrPred > max(tmpSet)){
                  #Adds a subset with cardsubset elements
                  count <- count + 1
                  datanew[, count] <- newSet
                  valnew[count] <- temp
                }
                #adds the weighted Improvement for this subset
                lmgRes[NrPred] = lmgRes[NrPred] + factorial(cardsubset-1) * factorial(NrAttr-cardsubset)*(temp-valold[NrSubset])
              }
            }
          }

          #preparing for the next step
          dataold <- datanew
          valold <- valnew
        }
        #By the weighted Calculation we actually averaged about every Permutation of Attributes
        #So we divide by the number of Permutations: factorial(NrAttr)
        lmgRes <- lmgRes/(factorial(NrAttr))
        lmgRes <- data.frame(lmgRes)

        lmgRes #Has the same order as the Attributes in the covariance matrix.
      }
      
      ###order
      ###If you want to know the improvement of the model if you add attributes in a specific
      ###order, you can give this vector as a parameter.
      orderfct <- function(order,posofAttr,XTX,XTY,SXY,nrowXTX,First){ 
        orderval <- First[order[1]]
        indices <- c(posofAttr[[order[1]]], nrowXTX)
        if(length(order) > 1){
          for(i in 2:length(order)){
            indices <- c(indices, posofAttr[[order[i]]])
            
            #Creates the covariance matrizes of the model from the original cov.matrix
            tempXTX <- XTX[indices, indices]
            tempSXY <- SXY[indices]
            tempXTY <- XTY[indices]
            orderval <- c(orderval, calc(tempSXY, tempXTX, tempXTY))
          }
        }
        names(orderval) <- NULL
        res <- data.frame(orderval)
        
        res
      }
      
      ###CALL#####
      #Sanity check if order is a vector with attributes that belong to the model.
      #Then it will calculate order,backward or forward depending on parameters
      if(is.vector(order)){
          res <- as.data.frame(orderfct(order, posofAttr, XTX, XTY, SXY, nrowXTX, First))
          names(res) <- c('Model_Values')
          rownames(res) <- order
      }else{
        if(backwardON){
          backward <- stepbystep(R2, posofAttr, SXY, XTX, XTY, max_forw = max_forw, direction = "backward")
          res <- data.frame(backward[, 2])
          names(res) <- c('Backward_Values')
          rownames(res) <- backward[, 1]
        }else{
          forward <- stepbystep(First, posofAttr, SXY, XTX, XTY, max_forw = max_forw)
          res <- data.frame(forward[, 2])
          names(res) <- c('Forward_Values')
          rownames(res) <- forward[, 1]
        }
      }
  
      if(lmgON){
        lmgval <- lmg(posofAttr, First, SXY, XTX, XTY)
      }
      
      #Order First and Usefulness into the order given by forward (or backward or order)
      neworder <- match(rownames(res), names(First))
    
      First <- First[neworder]
      Usefulness <- Usefulness[neworder]
      
      if(lmgON){
        lmgval <- lmgval[neworder,]
        res2 <- data.frame(res, Usefulness, First, lmgval)
        names(res2) <- c(names(res), 'Usefulness', 'First', 'LMG')
      }else{
        res2 <- data.frame(res, Usefulness, First)
        names(res2) <- c(names(res), 'Usefulness', 'First')
      }
      res2
    }
  }

  ###################Calling getRelImp############################
  if(lmgON){
    if(length(colNames) > 15){
      stop("The method recommends to not use LMG for more than 15 Attributes due to runtime.")
    }
  }
  tryCatch({
  relImp=getRelImp(colNames, card, nomCols, n, XTX, XTY, YTY, order, lmgON, backwardON)
  }, error = function(e){stop(e)})

  
###################Creating the Plot################################
    plot3 <- function(relImp, max_plot, names){
      reshaped <- relImp
      ModelVal <- colnames(reshaped)[1]
      
      reshaped$Attributes <- rownames(relImp)
      max_plot <- min(nrow(reshaped), max_plot)
      reshaped <- reshaped[1:max_plot, ]
      
      reshaped <- reshaped[, c('Attributes', 'Usefulness', 'First', ModelVal)]
      
      
      if(!names){
        reshaped$Attributes <- as.character(1:nrow(reshaped))
        reshaped$Attributes <- factor(reshaped$Attributes, levels = reshaped$Attributes)
      }else{
        reshaped$Attributes <- factor(reshaped$Attributes, levels = reshaped$Attributes)
        reshaped <- reshaped[1:nrow(reshaped), ]
      }
      
      #Using reshape function to switch table from wide to narrow format
      reshaped <- reshape(reshaped, varying=c( ModelVal, 'First', 'Usefulness'), 
                          v.names = c('Explanation'), direction="long")
      #A way to change the names in the Legend
      reshaped$time[reshaped$time == 1] <- 'Model_Value'
      reshaped$time[reshaped$time == 2] <- 'First'
      reshaped$time[reshaped$time == 3] <- 'Usefulness'
     
      #Change to factor so ggplot2 will not order it alphabetically
      reshaped$time <- factor(reshaped$time)
      
      
      #workaround some bug where another device was still open so the plot was not printed
      #close all devices.
      #try({graphics.off()}, silent=TRUE) #too destructive
      
      #Print plot into png
      png(filename = "RelImpPlot.png")
      plot1 <- ggplot(reshaped, aes_string(x = 'Attributes', y = 'Explanation', fill = "time",
                                           width=0.4)) +
              scale_x_discrete(expand=c(0.1,0))+
              scale_fill_manual(values=c("#B42600", "#033676", "#E7C500"),
                                breaks=c("Model_Value","First",  "Usefulness"),
                                c('Model_Value', 'First', 'Usefulness'))+
              geom_bar(stat = "identity", position = "identity")+
              guides(fill = guide_legend(reverse = FALSE, title = NULL))
      
      if(names == TRUE){
        plot1 <- plot1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      }

      print(plot1)
      dev.off()
      print(plot1)

      if(!names){  
        name_df <- as.data.frame(rownames(relImp[1:max_plot, ]))
        names(name_df) <- c('Attributes')
        cat("The order of the attributes in the plot is:\n")
        print(name_df)
      }
    }
    ############plot the lmg values if lmgON=TRUE############
    plotlmg <- function(relImp, max_plot, names){
      max_plot <- min(max_plot, nrow(relImp))
      data <- relImp
      data$Attributes <- rownames(relImp)
      data$Explanation <- relImp$LMG
      data <- data[1:max_plot, ]
      if(!names){
        data$Attributes <- as.character(1:nrow(data))
        data$Attributes <- factor(data$Attributes, data$Attributes)
      }else{
        data$Attributes <- factor(data$Attributes, data$Attributes)
        data <- data[1:nrow(data), ]
      }
      #try({graphics.off()}, silent=TRUE) #too destructive
      png(filename="RelImpPlot.png")
      plot1 <- ggplot(data, aes_string(x = 'Attributes', y = 'Explanation')) + 
        geom_bar(stat = "identity", fill="#FF6F00", width=0.4)+
        scale_x_discrete(expand=c(0.1,0))
       
      if(names == TRUE){
       plot1 <- plot1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      }
      print(plot1)
      dev.off()
      print(plot1)
      
      if(!names){  
        name_df <- as.data.frame(rownames(relImp[1:max_plot, ]))
        names(name_df) <- c('Attribute')
        cat("Die Attribute im Plot haben die folgende Reihenfolge:\n")
        print(name_df)
      }
    }
  
    relImp[relImp > 1] <- 1
    relImp[relImp < 0] <- 0
    if(lmgON){
      plotlmg(relImp, max_plot, names)
    }else{
      plot3(relImp, max_plot, names)
    }
  relImp#Returns the Results of the Importance Calculation
}