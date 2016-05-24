
##########################################################################
##function to fit a gdm object from a sitepair table
gdm <- function (data, geo=FALSE, splines=NULL, knots=NULL){
  #################
  ##lines used to quickly test function
  #data<-table
  #geo<-F
  #splines<-splines
  #knots<-knots
  #################
  options(warn.FPU = FALSE)
    
  ##adds error checking to gdm function
  ##checks to see if in site-pair format from formatsitepair function
  if(class(data)[1] != "gdmData"){
    warning("data class does not include type 'gdmData'. Make sure your data is in site-pair format or the gdm model will not fit.")
  }
  ##checks to makes sure data is a matrix or data frame
  if(!(class(data)[1]=="gdmData" | class(data)[1]=="matrix" | class(data)[1]=="data.frame")){
    stop("data argument needs to be a matrix or a data frame")
  }
    
  ##sanity check on the data table
  if(ncol(data) < 6){
    stop("Not enough columns in data. (Minimum need: Observed, weights, X0, Y0, X1, Y1)")
  }  
  if(nrow(data) < 1){
    stop("Not enough rows in data")
  }
    
  ##checks that geo has either TRUE or FALSE
  if(!(geo==TRUE | geo==FALSE)){
    stop("geo argument must be either TRUE or FALSE")
  }
  ##makes sure splines is a numeric vector
  if(is.null(splines)==FALSE & class(splines)!="numeric"){
    stop("argument splines needs to be a numeric data type")
  }
  ##checks knots inputs
  if(is.null(knots)==FALSE & class(knots)!="numeric"){
    stop("argument knots needs to be a numeric data type")
  }
    
  ##check that the response data is [0..1]
  rtmp <- data[,1]
  if(length(rtmp[rtmp<0]) > 0){
    stop("Response data has negative values. Must be between 0 - 1.")
  }
  if (length(rtmp[rtmp>1]) > 0){
    stop("Response data has values greater than 1. Must be between 0 - 1.")
  }

# gdm ---------------------------------------------------------------------


  ##current data format is response,weights,X0,Y0,X1,Y1 before any predictor data (thus 6 leading columns)
  LEADING_COLUMNS <- 6
  if(geo){
    nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2 + 1		
  }else{
    nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2
  }
  
  ##checks to make sure at least one predictor is available
  if(nPreds < 1){
    stop("Data has NO predictors")
  }

  ##setup the predictor name list, and removes the "s1." and "s2." to make resulting names more intuitive
  if(geo==TRUE){
    if(nPreds > 1){
      predlist <- c("Geographic", sapply(strsplit(names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds-1)], "s1."), "[[", 2))
    }else{
      predlist <- c("Geographic")
    }
  }else{
    predlist <- sapply(strsplit(names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds)], "s1."), "[[", 2)
  }

  ##deal with the splines and knots
  if(is.null(knots)){
    ##generate knots internally from the data
    if(is.null(splines)){
      nSplines <- 3
      quantvec <- rep(0, nPreds * nSplines)
      splinvec <- rep(nSplines, nPreds)
        
      if(geo==TRUE){
        ##get knots for the geographic distance
        v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
        quantvec[1] <- min(v)
        quantvec[2] <- median(v)
        quantvec[3] <- max(v)
          
        if(nPreds > 1){
          ##get knots for the environmental predictors
          for(i in seq(from = 1, to = nPreds-1, by = 1)){
            v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])                 
            index = i * nSplines
            quantvec[index+1] <- min(v)
            quantvec[index+2] <- median(v)
            quantvec[index+3] <- max(v)
          }
        }
      }else{
        ## get knots for the environmental predictors after skipping geographic preds
        for(i in seq(from = 1, to = nPreds, by = 1)){
          v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])                   
          index = (i-1) * nSplines
          quantvec[index+1] <- min(v)
          quantvec[index+2] <- median(v)
          quantvec[index+3] <- max(v)
        }
      }
    }else{
      ##otherwise check that the supplied splines vector has enough data and minumum spline values of 3
      if(length(splines) != nPreds){
        stop(paste("Number of splines does not equal the number of predictors. 
              Splines argument has", length(splines), "items but needs", nPreds, "items."))
      }
        
      ##count the total number of user defined splines to dimension the knots vector
      quantvec <- rep(0, sum(splines))
      splinvec <- splines
        
      if(geo==T){
        if(splines[1] < 3){
          stop("Must have at least 3 splines per predictor")
        }
          
        ## get knots for the geographic distance
        v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
        quantvec[1] <- min(v)		## 0% knot
        quantvec[splines[1]] <- max(v)	## 100% knot
          
        quant_increment <- 1.0 / (splines[1]-1)
        this_increment <- 1
        for (i in seq(from = 2, to = (splines[1]-1), by = 1)){  
          ## mid % knots
          quantvec[i] <- quantile(v,quant_increment*this_increment)
          this_increment <- this_increment + 1
        }
          
        if(nPreds > 1){
          ##get knots for the environmental predictors
          current_quant_index <- splines[1]
          for(i in seq(from = 1, to = nPreds-1, by = 1)){
            num_splines <- splines[i+1]
            if(num_splines < 3){
              stop("Must have at least 3 splines per predictor")
            }
                
            v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])                 
            quantvec[current_quant_index+1] <- min(v)	            ## 0% knot
            quantvec[current_quant_index+num_splines] <- max(v)	    ## 100% knot
              
            quant_increment <- 1.0 / (num_splines-1)
            this_increment <- 1
            for(i in seq(from = 2, to = (num_splines-1), by = 1)){  
              ##mid % knots
              quantvec[current_quant_index+i] <- quantile(v,quant_increment*this_increment)
              this_increment <- this_increment + 1
            }
              
            current_quant_index <- current_quant_index + num_splines
          }
        }
      }else{
        ##get knots for the environmental predictors
        current_quant_index <- 0
        for(i in seq(from = 1, to = nPreds, by = 1)){
          num_splines <- splines[i]
          if(num_splines < 3){
            stop("Must have at least 3 splines per predictor")
          }
              
          v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])          
          quantvec[current_quant_index+1] <- min(v)	        ## 0% knot
          quantvec[current_quant_index+num_splines] <- max(v)	## 100% knot
            
          quant_increment <- 1.0 / (num_splines-1)
          this_increment <- 1
          for(i in seq(from = 2, to = (num_splines-1), by = 1)){  
            ##mid % knots
            quantvec[current_quant_index+i] <- quantile(v,quant_increment*this_increment)
            this_increment <- this_increment + 1
          }
          current_quant_index <- current_quant_index + num_splines
        }
      }
    }
  }else{
    ##user defined knots supplied as an argument
    if(is.null(splines)){
      ##check that there are nPreds * 3 knots in the user defined vector
      if(length(knots) != (nPreds * 3)){
        stop(paste("When knots are supplied by the user, there should be", (nPreds * 3), "items in the knots argument, not", length(knots), "items."))
      }
        
      ## now check that each of the three knots for each predictor are in ascending order
      for(i in seq(from = 1, to = nPreds, by = 1)){
        index = i * 3
        if((knots[index-1] < knots[index-2]) || 
            (knots[index] < knots[index-2]) || 
            (knots[index] < knots[index-1])){
          stop(paste("Knots for ", predlist[i], "are not in ascending order."))
        }
      }
        
      nSplines <- 3
      quantvec <- knots
      splinvec <- rep(nSplines, nPreds)
    }else{
      ##check that there are sum(splines) knots in the user defined vector
      if(length(knots) != sum(splines)){
        stop(paste("When knots are supplied by the user, there should be", sum(splines), "items in the knots argument, not", length(knots), "items."))
      }
        
      ##now check that each of the knots for each predictor are in ascending order
      index = 0
      for(i in seq(from = 1, to = nPreds, by = 1)){
        for(j in seq(from = 2, to = splines[i], by = 1)){
          if(knots[index+j] < knots[index+j-1]){
            stop(paste("Knots for ", predlist[i], "are not in ascending order."))
          }
        }
        index <- index + splines[i]
      }
        
      quantvec <- knots
      splinvec <- splines
    }
  }
    
  p1 <- 0
  p2 <- 0
  p3 <- 0
  p4 <- 0
  p5 <- rep(0,times=length(quantvec))
  p6 <- rep(0,times=nrow(data))
  p7 <- rep(0,times=nrow(data))
  p8 <- rep(0,times=nrow(data))

  ##Call the dll function, fitting the gdm model
  z <- .C( "GDM_FitFromTable",
            paste(getwd()),
            as.matrix(data),
            as.integer(geo),
            as.integer(nPreds), 
            as.integer(nrow(data)), 
            as.integer(ncol(data)),
            as.integer(splinvec),
            as.double(quantvec),                 
            gdmdev = as.double(p1),
            nulldev = as.double(p2),
            expdev = as.double(p3),
            intercept = as.double(p4),         
            coeffs = as.double(p5),
            response = as.double(p6),
            preddata = as.double(p7),
            ecodist = as.double(p8), 
            PACKAGE = "gdm")
    
  m <- match.call(expand.dots = F)
  
  ##creates the gdm object, and fills its parts
  gdmModOb <- structure(list(dataname = m[[2]],
                              geo = geo,
                              sample = nrow(data),
                              gdmdeviance = z$gdmdev,
                              nulldeviance = z$nulldev,
                              explained = z$expdev,
                              intercept = z$intercept,
                              predictors = predlist,
                              coefficients = z$coeffs,
                              knots = quantvec,
                              splines = splinvec,
                              creationdate = date(),
                              observed = z$response,
                              predicted = z$preddata,
                              ecological = z$ecodist))
  ##sets gdm object class  
  class(gdmModOb) <- c("gdm", "list")
  
  ##reports a warning should the model "fit", yet the sum of coefficients = 0
  if(sum(gdmModOb$coefficients)==0){
    warning("Problem with model fitting, no solution obtained. Sum of spline coefficients = 0. Deviance explained = 0.")
    ##sets the deviance explained to 0, to reflect that the model didn't fit correctly
    gdmModOb$explained <- 0
  }
  
  ##returns gdm object
  return(gdmModOb)
}

##########################################################################


##########################################################################
##function to plot the splines of a gdm object
plot.gdm <- function (x, plot.layout = c(2,2), plot.color = "blue", 
                      plot.linewidth=2.0, include.rug=FALSE, rug.sitepair=NULL, ...){
  #################
  ##lines used to quickly test function
  #x <- gdm.1
  #plot.layout <- c(2,2)
  #plot.color <- "green"
  #plot.linewidth <- 2.0
  #include.rug <- T
  #rug.sitepair <- gdmTab
  #################
  ##checks to make sure that a site-pair table has been included
  
  
  options(warn.FPU = FALSE)
  PSAMPLE <- 200
  preddata <- rep(0,times=PSAMPLE)
  
  ##establish what plot layout to use
  thisplot <- 0
  one_page_per_plot <- FALSE
  if ((plot.layout[1]==1) && (plot.layout[2]==1)){
    one_page_per_plot <- TRUE
  }else{
    par(mfrow=plot.layout)
  }
  
  ##plots the compositional dissimilarity spline plot
  plot(x$ecological, x$observed, xlab="Predicted Ecological Distance", ylab="Observed Compositional Dissimilarity", type="n", ylim=c(0,1))
  points(x$ecological, x$observed, pch=20, cex=0.25, col=plot.color)
  overlayX <- seq( from=min(x$ecological), to=max(x$ecological), length=PSAMPLE )
  overlayY <- 1 - exp( - overlayX )
  lines( overlayX, overlayY, lwd=plot.linewidth )
  thisplot <- thisplot + 1
  
  ##determines rather or not to put mulitiple plots on one page or not
  if(one_page_per_plot){
    dev.new()
    dev.next()
  }
  ##plots the second compositional dissimilarity spline plot
  plot(x$predicted, x$observed, xlab="Predicted Compositional Dissimilarity", ylab="Observed Compositional Dissimilarity", type="n", ylim=c(0,1))
  points( x$predicted, x$observed, pch=20, cex=0.25, col=plot.color )
  overlayX <- overlayY <- seq( from=min(x$predicted), to=max(x$predicted), length=PSAMPLE )
  lines( overlayX, overlayY, lwd=plot.linewidth ) 
  thisplot <- thisplot + 1
  
  ##determine the max of all the predictor data, to be used in the plotting below
  preds <- length(x$predictors)
  predmax <- 0
  splineindex <- 1
  for(i in 1:preds){  
    ##only if the sum of the coefficients associated with this predictor is > 0.....
    numsplines <- x$splines[i]
    if(sum(x$coefficients[splineindex:(splineindex+numsplines-1)]) > 0){
      ## get predictor plot Y-data                            
      z <- .C( "GetPredictorPlotData", 
               pdata = as.double(preddata),
               as.integer(PSAMPLE),
               as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
               as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
               as.integer(numsplines), 
               PACKAGE = "gdm" )
      
      v <- max(z$pdata)
      if(v > predmax){
        predmax <- v
      } 
    }
    ##update the spline index
    splineindex <- splineindex + numsplines
  }
  
  ##plot the predictors with non-zero sum of coefficients      
  splineindex <- 1
  for(i in 1:preds){  
    #i <- 2
    ##only if the sum of the coefficients associated with this predictor is > 0.....
    numsplines <- x$splines[i]
    if(sum(x$coefficients[splineindex:(splineindex+numsplines-1)]) > 0){
      ##plots one graph per page, unless specified otherwise
      if (one_page_per_plot){
        dev.new()
        dev.next()
      }else{
        thisplot <- thisplot + 1
        if(thisplot > (plot.layout[1] * plot.layout[2])){	
          #dev.new()                              				
          #x11()
          #dev.next()					
          thisplot <- 1
          par(mfrow=plot.layout)	
        }
      }
      
      ##get predictor plot Y-data    
      z <- .C( "GetPredictorPlotData", 
               pdata = as.double(preddata),
               as.integer(PSAMPLE),
               as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
               as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
               as.integer(numsplines),
               PACKAGE = "gdm")
      
      if(x$geo & i==1){
        varNam <- "Geographic Distance"
        ##calculates rug plot data
        if(include.rug==TRUE){
          rugData <- unique(sqrt(((rug.sitepair$s1.xCoord-rug.sitepair$s2.xCoord)^2)+((rug.sitepair$s1.yCoord-rug.sitepair$s2.yCoord)^2)))
        }
      } else{
        varNam <- x$predictors[i]
        ##gets rug plot data
        if(include.rug==TRUE){
          varDat <- grep(varNam, colnames(rug.sitepair))
          rugData <- unique(c(rug.sitepair[,c(varDat[1])], rug.sitepair[,c(varDat[2])]))
        }
      }
      
      plot(seq(from=x$knots[[(i*3)-2]],to=x$knots[[(i*3)]], length=PSAMPLE), z$pdata, 
           xlab=varNam, ylab=paste("f(", varNam, ")", sep="" ), ylim=c(0,predmax), type="l")
      if(include.rug==TRUE){
        rug(rugData)
      }
    }
    ##update the spline index
    splineindex <- splineindex + numsplines
  }       
}
##########################################################################


##########################################################################
##function to either predict the biological dissimilarities between sites, 
##or to predict the dissimilarity of the same sites between two time periods,
##based on a gdm
predict.gdm <- function (object, data, time=FALSE, predRasts=NULL, ...){
  #################
  ##lines used to quickly test function
  ##object = gdm model
  ##data = a sitepair table
  #object <- gdm1
  #data <- climExtCurr
  #time <- T
  #predRasts <- climExtFuture
  #################
  options(warn.FPU = FALSE)
  
  ##if making a time prediction, makes sure all data is in the correct format,
  ##and then transforms the raster data into data tables in order to utalize
  ##the C predict utility
  if(time==TRUE){
    ##checks to make sure the inputs are correct
    if(is.null(predRasts)==TRUE){
      stop("Prediction rasters needed when time is TRUE")
    }
    if(class(data)!="RasterStack" & class(data)!="RasterLayer" & class(data)!="RasterBrick"){
      stop("Prediction data need to be a raster object when time is TRUE")
    }
    if(class(predRasts)!="RasterStack" & class(predRasts)!="RasterLayer" & class(predRasts)!="RasterBrick"){
      stop("predRasts need to be a raster object when time is TRUE")
    }
    if(nlayers(data)!=nlayers(predRasts)){
      stop("Current and future raster objects must have the same number of layers")
    }
    if(nlayers(data)!=length(object$predictors)-1 | nlayers(predRasts)!=length(object$predictors)-1){
      stop("Number of predictor variables does not equal the number used to fit the model")
    }
    for(i in 1:nlayers(data)){
      if(names(data)[i]!=names(predRasts)[i]){
        stop("Layer names do not match the variables used to fit the model")
      }
    }
    ##tests to see if raster data is stackable
    tryRasts <- try(stack(data[[1]], predRasts[[1]]), silent=TRUE)
    if(class(tryRasts)=="try-error"){
      stop("Current and Prediction rasters do not stack, make sure rasters are spatially congruant.")
    }
    
    ##sets up sitepair table with current and future data
    predLayer <- data[[1]]
    currXY <- as.data.frame(na.omit(rasterToPoints(data, progress='text')))
    predXY <- as.data.frame(na.omit(rasterToPoints(predRasts, progress='text')))
    cells <- cellFromXY(predLayer, cbind(currXY$x, currXY$y))
    dummData <- rep.int(0, nrow(currXY))
    data <- cbind(dummData, dummData, currXY[,1:2], currXY, predXY[,-c(1,2)])
    ##adds s1 or s2 to the variables name of the data
    t1var <- paste("s1.", colnames(currXY)[-c(1,2)], sep="")
    t2var <- paste("s2.", colnames(predXY)[-c(1,2)], sep="")
    ##sets the correct names to the data
    colnames(data) <- c("distance", "weights", "s1.xCoord", "s1.yCoord", 
                        "s2.xCoord", "s2.yCoord", t1var, t2var)
  }
  
  ##makes the prediction based on the data object
  predicted <- rep(0,times=nrow(data))
  z <- .C( "GDM_PredictFromTable",
           as.matrix(data),
           as.integer(object$geo),
           as.integer(length(object$predictors)), 
           as.integer(nrow(data)), 
           as.double(object$knots),
           as.integer(object$splines),
           as.double(c(object$intercept,object$coefficients)),
           preddata = as.double(predicted),
           PACKAGE = "gdm")
  
  ##if a time prediciton, maps the predicted values to a raster and returns
  ##the layer, otherwise returns a dataframe of the predicted values
  if(time==FALSE){
    return(z$preddata)
  }else{
    predLayer[cells] <- z$preddata
    return(predLayer)
  }
}
##########################################################################


##########################################################################
##function to transform input data into biological space based on a given gdm
gdm.transform <- function (model, data){
  #################
  ##lines used to quickly test function
  #model <- gdm.1 
  #data <- envTrans
  #data <- cropRasts[[3:nlayers(cropRasts)]]
  #################
  options(warn.FPU = FALSE)
  rastDat <- NULL
  dataCheck <- class(data)
  
  ##error checking of inputs
  ##checks to make sure a gdm model is given
  if(class(model)[1]!="gdm"){
    stop("model argument must be a gdm model object")
  }
  ##checks to make sure data is a correct format
  if(!(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick" | dataCheck=="data.frame")){
    stop("Data to be transformed must be either a raster object or data frame")
  }
  
  ##checks rather geo was T or F in the model object
  geo <- model$geo

  ##turns raster data into dataframe
  if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
    ##converts the raster object into a dataframe, for the gdm transformation
    rastDat <- data
    data <- rasterToPoints(rastDat)
    ##determines the cell number of the xy coordinates
    rastCells <- cellFromXY(rastDat, xy=data[,1:2]) 
    
    ##checks for NA in the 
    checkNAs <- as.data.frame(which(is.na(data), arr.ind=T))
    if(nrow(checkNAs)>0){
      warning("After extracting raster data, NAs found from one or more layers. Removing NAs from data object to be transformed.")
      data <- na.omit(data)
      rastCells <- rastCells[-c(checkNAs$row)]
    }
    
    ##if geo was not T in the model, removes the coordinates from the data frame
    if(geo==FALSE){
      data <- data[,3:ncol(data)]
    }
  }
  
  sizeVal <- 10000000
  ##sets up the data to be transformed into pieces to be transformed
  holdData <- data
  fullTrans <- matrix(0,nrow(holdData),ncol(holdData))
  rows <- nrow(holdData)
  istart <- 1
  iend <- min(sizeVal,rows)
  ##to prevent errors in the transformation of the x and y values when geo is a predictor,
  ##extracts the rows with the minimum and maximum x and y values, these rows will be added
  ##onto the "chuck" given to transform, and then immediately removed after the transformation,
  ##this makes sure that the c++ code will always have access to the minimum and maximum 
  ##x and y values
  if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
    xMaxRow <- holdData[which.max(holdData[,"x"]),]
    xMinRow <- holdData[which.min(holdData[,"x"]),]
    yMaxRow <- holdData[which.max(holdData[,"y"]),]
    yMinRow <- holdData[which.min(holdData[,"y"]),]
  }
  
  ##transform the data based on the gdm
  ##part of a loop to prevent memory errors 
  while(istart < rows){
    ##Call the dll function
    data <- holdData[istart:iend,]
    ##adds coordinate rows to data to be transformed
    if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
      data <- rbind(xMaxRow, xMinRow, yMaxRow, yMinRow, data)
    }
    transformed <- matrix(0,nrow(data),ncol(data))
    z <- .C( "GDM_TransformFromTable",
             as.integer(nrow(data)), 
             as.integer(ncol(data)),
             as.integer(model$geo),
             as.integer(length(model$predictors)), 
             as.integer(model$splines),             
             as.double(model$knots),             
             as.double(model$coefficients),
             as.matrix(data),
             trandata = as.double(transformed),
             PACKAGE = "gdm")
    
    ## Convert transformed from a vector into a dataframe before returning...
    nRows <- nrow(data)
    nCols <- ncol(data)
    
    ## z$trandata is the transformed data vector created
    myVec <- z$trandata
    pos <- 1
    ##fills out dataframe with transformed values
    for (i in seq(from = 1, to = nCols, by = 1)) {
      tmp <- myVec[seq(from=pos, to=pos+nRows-1)]
      transformed[,i] <- tmp
      pos <- pos + nRows
    }
    
    ##remove the coordinate rows before doing anything else
    if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
      transformed <- transformed[-c(1:4),]
    }
    
    ##places the transformed values into the readied data frame 
    fullTrans[istart:iend,] <- transformed
    istart <- iend + 1
    iend <- min(istart + (sizeVal-1), rows)
  }
  
  ##if wanted output data as raster, provides maps raster, or output table
  if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
    ##maps the transformed data back to the input rasters
    rastLay <- rastDat[[1]]
    rastLay[] <- NA
    outputRasts <- stack()
    for(nn in 1:ncol(fullTrans)){
      #print(nn)
      #nn=1
      holdLay <- rastLay
      holdLay[rastCells] <- fullTrans[,nn]
      #holdLay[rastCells] <- holdData[,nn]
      
      outputRasts <- stack(outputRasts, holdLay)
    }
    ##renames raster layers to be the same as the input
    if(geo){
      names(outputRasts) <- c("xCoord", "yCoord", names(rastDat))
    } else {
      names(outputRasts) <- names(rastDat)
    }
    
    ##get the predictors with non-zero sum of coefficients      
    splineindex <- 1
    predInd <- NULL
    for(i in 1:length(model$predictors)){  
      #i <- 1
      ##only if the sum of the coefficients associated with this predictor is > 0.....
      numsplines <- model$splines[i]
      if(sum(model$coefficients[splineindex:(splineindex+numsplines-1)])>0){
        predInd <- c(predInd, i)
      }
      splineindex <- splineindex + numsplines
    }
    if(geo){
      predInd <- c(1,2,predInd[-1]+1)
    }
    
    outputRasts <- outputRasts[[predInd]]
    
    ##returns rasters
    return(outputRasts)
  }else{
    if(is.null(rastDat)){
      ##if not raster data, sends back the transformed data
      colnames(fullTrans) <- colnames(data)
      return(fullTrans)
    }else{
      ##returns only the transformed variable data as a table, and the cells with which to map to
      colnames(fullTrans) <- colnames(data)
      return(list(fullTrans, rastCells))
    }
  }
}
##########################################################################


##########################################################################
##function to print a summary of a gdm object
summary.gdm <- function (object, ...){
  print( "", quote=F )    
  print( "", quote=F )    
  print( "GDM Modelling Summary", quote=F );
  print( paste( "Creation Date: ", object$creationdate ), quote=F );
  print( "", quote=F )    
  ##        call <- match.call()
  m <- match.call(expand.dots = F)
  print( paste( "Name: ", m[[2]] ), quote=F )
  print( "", quote=F )    
  print( paste( "Data: ", object$dataname ), quote=F )
  print( "", quote=F )    
  print( paste( "Samples: ", object$sample ), quote=F )
  print( "", quote=F )    
  print( paste( "Geographical distance used in model fitting? ", object$geo ), quote=F )
  print( "", quote=F )    
  print( paste( "NULL Deviance: ", object$nulldeviance ), quote=F )
  print( paste( "GDM Deviance: ", object$gdmdeviance ), quote=F )  
  print( paste( "Deviance Explained: ", object$explained ), quote=F )
  print( "", quote=F )    
  print( paste( "Intercept: ", object$intercept ), quote=F )
  print( "", quote=F )    
  thiscoeff <- 1
  thisquant <- 1
  for(i in 1:length(object$predictors)){
    print( paste( "Predictor ",i,": ",object$predictors[[i]], sep="" ), quote=F )            
    print( paste( "Splines: ",object$splines[[i]], sep="" ), quote=F )
    numsplines <- object$splines[[i]]
    for(j in 1:numsplines){
      if ( j == 1 ) print( paste( "Min Knot: ",object$knots[[thisquant]], sep="" ), quote=F )          
      else if ( j == numsplines ) print( paste( "Max Knot: ",object$knots[[thisquant]], sep="" ), quote=F )
      else print( paste( round(100/(numsplines-1),digits=2),"% Knot: ",object$knots[[thisquant]], sep="" ), quote=F )
      thisquant <- thisquant + 1
    }
    for(j in 1:numsplines){
      print( paste( "Coefficient[",j,"]: ",object$coefficients[[thiscoeff]], sep="" ), quote=F )
      thiscoeff <- thiscoeff + 1
    }
    print( "", quote=F )                
  }   
}
##########################################################################


##########################################################################
##Takes species data from a variety of commonly used formats and transforms it into a
##site-pair table for GDM
formatsitepair <- function(bioData, bioFormat, dist="bray", abundance=FALSE, 
                           siteColumn=NULL, XColumn, YColumn, sppColumn=NULL, 
                           abundColumn=NULL, sppFilter=0, predData, distPreds=NULL, 
                           weightType="equal", custWeights=NULL, samples=NULL){
  ###########################
  ##lines used to quickly test function
  #bioData <- sppTab
  #bioFormat <- 2
  #dist <- "bray"
  #abundance <- F
  #siteColumn <- "site"
  #XColumn <- "Long"
  #YColumn <- "Lat"
  #sppColumn <- "species"
  #sppFilter <- 10
  #abundColumn <- NULL
  #predData <- envTab
  #distPreds <- list(newGDMdissim)
  #weightType <- "equal"
  #custWeights <- NULL
  #samples <- 80
  #################
  #bioData <- exFormat2a
  #bioFormat <- 4
  #dist <- "bray"
  #abundance <- F
  #siteColumn <- "site"
  #XColumn <- NULL
  #YColumn <- NULL
  #sppColumn <- NULL
  #sppFilter <- 0
  #abundColumn <- NULL
  #predData <- envTab
  #distPreds <- list(as.matrix(gdmDissim))
  #weightType <- "equal"
  #custWeights <- NULL
  #samples <- NULL
  ###########################
  ##input error checking
  ##makes sure bioData is in an acceptable format
  if(!(class(bioData)[1]=="data.frame" | class(bioData)[1]=="matrix" | class(bioData)[1]=="gdmData")){
    "bioData should be a data frame or matrix in one of the acceptable formats"
  }
  ##makes sure predData is in an acceptable format
  if(!(class(predData)=="data.frame" | class(predData)=="matrix" | class(predData)=="RasterStack" | class(predData)=="RasterLayer" | class(predData)=="RasterBrick")){
    "predData should be a data frame or matrix in one of the acceptable formats"
  }
  ##if bioFormat is not an acceptable number, exit function
  if(bioFormat %in% c(1:4)){} else{
    stop("Acceptable values for the bioFormat argument are: 1, 2, 3, or 4")
  }
  ##checks that geo has either TRUE or FALSE
  if(!(abundance==TRUE | abundance==FALSE)){
    stop("abundance argument must be either TRUE or FALSE")
  }
  ##if samples is not a number, then exit function
  if(is.null(samples)==FALSE){
    if(is.numeric(samples)==FALSE | samples<0){
      stop("samples argument must be a positive integer")
    }
  }
  
  ##makes sure that sppFilter is a number, if not exit function
  if((is.numeric(sppFilter)==FALSE & is.null(sppFilter)==FALSE) | sppFilter<0){
    stop("sppFilter argument must be a positive integer")
  }
  ##makes sure a proper weightType is used
  if(weightType %in% c("equal", "richness", "custom")){} else{
    stop("Acceptable values for the weightType argument are: equal, richness, or custom")
  }
  ##if weightType == custom, makes sure a custWeights is attached
  if(weightType=="custom" & is.null(custWeights)==T){
    stop("custom weightType provided with no custWeights")
  }
  ##with bioFormat 2, makes sure spp data has been included
  if(bioFormat==2 & is.null(sppColumn)==TRUE){
    stop("Need to define sppColumn argument when bioFormat==2") 
  }
  ##makes sure that the sppColumn name can be found in the bioData with bioFormat 2
  #if(bioFormat==2 & (sppColumn %in% names(bioData))){} else{
  #  stop("Cannot find sppColumn in bioData - check name?")
  #}
  ##makes sure that a site column is provided when using table type 2 and raster environmental data
  if(bioFormat==2 & is.null(siteColumn)==TRUE){
    if(!(class(predData)=="RasterStack" | class(predData)=="RasterLayer" | class(predData)=="RasterBrick")){
      stop("A siteColumn needs to be provided in either the bioData or predData inputs")
    }
  }
  ##when a site column is provided
  if(is.null(siteColumn)==FALSE){
    ##makes sure the site column name is of type character
    if(class(siteColumn)!="character"){
      stop("siteColumn argument needs to be as type character")
      ##checks to see if siteColumn exists in the bioData for bioFormats 1 and 2    
    }else if(!(siteColumn %in% colnames(bioData)) & (bioFormat==1 | bioFormat==2)){
      stop("Cannot find a match for siteColumn in the columns of bioData.")
    }
    ##if the siteColumn is provided with input type 3, remove it
    #if(bioFormat==3 & siteColumn %in% colnames(bioData)){
    #  wSite <- which(colnames(bioData)==siteColumn)
    #  bioData <- bioData[-wSite]
    #}
  }
  ##checks to make sure that the coordinate columns are characters and can be found in either the biological or 
  ##environmental data
  if(bioFormat!=4){
    if(class(XColumn)!="character"){
      stop("XColumn argument needs to be as type character")
    }else if(class(YColumn)!="character"){
      stop("YColumn argument needs to be as type character")
    }else if(!(XColumn %in% colnames(bioData) | XColumn %in% colnames(predData))){
      stop("XColumn not found in either the bioData or predData arguments")
    }else if(!(YColumn %in% colnames(bioData) | YColumn %in% colnames(predData))){
      stop("YColumn not found in either the bioData or predData arguments")
    }
  }
  ##checks table type 3 specific requirements
  if(bioFormat==3){
    if(weightType=="richness"){
      stop("Cannot weight by site richness when supplying the biological data as a distance matrix.")
    }else if(nrow(bioData)!=(ncol(bioData)-1)){
      stop("Biological dissimilarity has differing number of rows to number of columns. Did you forget to add a column for site ID's?")
    }
  }
  
  ##warns if distPreds are not matrices
  for(mat in distPreds){
    if(class(mat)!="matrix" & class(mat)!="data.frame"){
      warning("One or more of the provided distance predictor matrices are not of class 'matrix'.")
    }
  }
  ##if a custom weight vector is provided, makes sure it is a vector
  if(is.null(custWeights)==FALSE & (class(custWeights)!="data.frame" & class(custWeights)!="matrix")){
    stop("argument custWeights needs to be a data frame or matrix data type")
  }
  
  ##sets up variables to be used later
  toRemove <- NULL       ##removed from sppfilter
  removeRand <- NULL     ##remove from random sample out sites
  distData <- NULL
  ##checks input data format
  ##species data as site-species matrix
  if(bioFormat==1 | bioFormat==2){
    ##first, if bioFormat=2, then transforms it into a site-spp matrix (bioFormat 1)
    if(bioFormat==2){
      ##makes sure that the sppColumn name can be found in the bioData with bioFormat 2
      if((sppColumn %in% names(bioData))){} else{
        stop("Cannot find sppColumn in bioData - check name?")
      }
      
      ##insert a site column if one was not given
      if(is.null(siteColumn)){
        colnames(bioData)[which(colnames(bioData)==XColumn)] <- "myXness"
        colnames(bioData)[which(colnames(bioData)==YColumn)] <- "myYness"
        bioData <- transform(bioData, siteUltimateCoolness=as.numeric(interaction(bioData$myXness, bioData$myYness, drop=TRUE)))
        siteColumn <- "siteUltimateCoolness"
        colnames(bioData)[which(colnames(bioData)=="myXness")] <- XColumn
        colnames(bioData)[which(colnames(bioData)=="myYness")] <- YColumn 
      }
      
      ##insert presence if abundance was not given
      if(is.null(abundColumn)){
        warning("No abundance column specified, assuming species data are presence")
        bioData["reallysupercoolawesomedata"] <- 1
        abundColumn <- "reallysupercoolawesomedata"
      }
      
      ##rename the siteColumn and sppColumn in order to cast the data into a siteXspp matrix
      preCastBio <- bioData
      colnames(preCastBio)[which(colnames(preCastBio)==siteColumn)] <- "siteUltimateCoolness"
      colnames(preCastBio)[which(colnames(preCastBio)==sppColumn)] <- "spcodeUltimateCoolness"
      castData <- reshape2::dcast(preCastBio, siteUltimateCoolness~spcodeUltimateCoolness, value.var=abundColumn)
      ##adds coordinates to the cast data
      uniqueCoords <- unique(preCastBio[which(colnames(preCastBio) %in% c("siteUltimateCoolness", XColumn, YColumn))])
      bioData <- merge(castData, uniqueCoords, by="siteUltimateCoolness")
      colnames(bioData)[which(colnames(bioData)=="siteUltimateCoolness")] <- siteColumn
    }
    
    ##checks to see if the coordinates can be found in bioData, if not, checks to see if 
    ##they can be found in envData
    if((XColumn %in% colnames(bioData))==FALSE | (YColumn %in% colnames(bioData)==FALSE)){
      xCol <- which(colnames(predData)==XColumn)
      yCol <- which(colnames(predData)==YColumn)
      locs <- predData[c(xCol,yCol)]
    }else{
      xCol <- which(colnames(bioData)==XColumn)
      yCol <- which(colnames(bioData)==YColumn)
      locs <- bioData[c(xCol,yCol)]
    }
    
    ##checks unique sites against rasters
    if(class(predData)=="RasterStack" | class(predData)=="RasterLayer" | class(predData)=="RasterBrick"){
      ##when using rasters, uses the cell as the site 
      warning("With raster prediction data, defaulting sites to cells.")
      ##gets the cell location of the given coordinates
      cellID <- as.data.frame(cellFromXY(predData, locs))
      colnames(cellID)[which(colnames(cellID)=="cellFromXY(predData, locs)")] <- "cellName"
      ##if none of the points intersected with the prediction raster
      if(nrow(cellID)==sum(is.na(cellID$cellName))){
        stop("None of the given points intersect with the given raster data. Double check that your geography is correct and that the given XColumn and YColumn values are correct.")
      }
      cellLocs <- as.data.frame(xyFromCell(predData, cellID$cellName))
      ##temporarily keeps old site in to identify what to remove from other objects
      rastBioData <- cbind(cellID, cellLocs, bioData[-c(which(colnames(bioData) %in% c(XColumn, YColumn)))])
      
      ##if custom weights selected, gives weight table new cell site names
      if(weightType=="custom" & !is.null(custWeights)){
        nameTab <- unique(rastBioData[c("cellName", siteColumn)])
        tempWeightTab <- merge(x=nameTab, y=custWeights, by=siteColumn)
        siteNum <- which(colnames(tempWeightTab)=="cellName")
        custWeights <- tempWeightTab[-siteNum]
        colnames(custWeights)[1] <- siteColumn
      }
      ##removes original site column from the rastBioData table
      siteNum <- which(colnames(rastBioData)==siteColumn)
      rastBioData <- rastBioData[-siteNum] 
      
      ##aggregates species data by cell
      cellNum <- which(colnames(rastBioData)=="cellName")
      bioData <- aggregate(rastBioData, rastBioData[cellNum], FUN=mean)
      bioData <- bioData[-cellNum]
      
      ##extracts raster data into environmental prediction data table
      rastEx <- as.data.frame(extract(predData, bioData$cellName))
      
      ##renames bioData columns which have been updated from rasters
      colnames(bioData)[which(colnames(bioData)=="cellName")] <- siteColumn
      colnames(bioData)[which(colnames(bioData)=="x")] <- XColumn
      colnames(bioData)[which(colnames(bioData)=="y")] <- YColumn
      
      ##updates locs object based on raster coordinates
      xCol <- which(colnames(bioData)==XColumn)
      yCol <- which(colnames(bioData)==YColumn)
      locs <- bioData[c(xCol,yCol)] 
      
      ##recreates the predData
      siteCol <- which(colnames(bioData)==siteColumn)
      predData <- cbind(bioData[siteCol], locs, rastEx)
    }
    
    ##filters out sites with low species count 
    ##first isolates the species data
    siteCol <- which(colnames(bioData)==siteColumn)
    xCol <- which(colnames(bioData)==XColumn)
    yCol <- which(colnames(bioData)==YColumn)
    sppDat <- bioData[-c(siteCol, xCol, yCol)]
    ##totals the number of species per site
    sppDat[sppDat>=1] <- 1
    sppDat[sppDat==0] <- 0
    sppDat[is.na(sppDat)] <- 0
    sppTotals <- cbind(as.data.frame(bioData[siteCol]), apply(sppDat, 1, function(m){sum(as.numeric(m))}))
    ##filters out sites with less species than filter
    filterBioDat <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] >= sppFilter)
    toRemove <- bioData[,siteCol][which(!(bioData[,siteCol] %in% filterBioDat[,siteColumn]))]
    ##reassembles bioData after filtering
    spSiteCol <- filterBioDat[1]
    bioData <- unique(merge(spSiteCol, bioData, by=siteColumn))
    
    ##removes random sampling of sites
    if(is.null(samples)==FALSE){
      if(nrow(bioData)<samples){
        warning("After species filter, fewer records remaining than specified in samples, continuing without subsampling")
      }else{
        fullSites <- bioData[,siteCol]
        randRows <- sample(1:nrow(bioData), samples)
        ##actual selection of the random rows to keep
        bioData <- bioData[c(randRows),]
        #removeRand <- fullLength[-(randRows)]
        ##records the sites that have been removed, for distPreds later in function
        removeRand <- fullSites[which(! (fullSites %in% bioData[,siteCol]))]
      } 
    }
    
    ##identifies and removes filtered out sites and sampled sites from predData
    ##renames siteColumn in order to access objects correctly
    colnames(bioData)[colnames(bioData)==siteColumn] <- "gettingCoolSiteColumn"
    colnames(predData)[colnames(predData)==siteColumn] <- "gettingCoolSiteColumn"
    predData <- unique(predData)
    predData <- predData[which(predData$gettingCoolSiteColumn %in% bioData$gettingCoolSiteColumn),]
    
    ##remove custom weights from any sites removed by species filtering and sampling
    if(weightType=="custom" & !is.null(custWeights)){
      colnames(custWeights)[colnames(custWeights)==siteColumn] <- "gettingCoolSiteColumn"
      custWeights <- custWeights[which(predData$gettingCoolSiteColumn %in% custWeights[,"gettingCoolSiteColumn"]),]
      custWeights <- custWeights[order(custWeights[,"gettingCoolSiteColumn"]),]
      colnames(custWeights)[colnames(custWeights)=="gettingCoolSiteColumn"] <- siteColumn
    }
    
    ##rename site columns for results
    colnames(bioData)[colnames(bioData)=="gettingCoolSiteColumn"] <- siteColumn
    colnames(predData)[colnames(predData)=="gettingCoolSiteColumn"] <- siteColumn
    
    ##as a final check, makes sure bioData and predData sites are in same order
    predSite <- which(names(predData) == siteColumn)
    bioSite <- which(names(bioData)==siteColumn)
    predData <- predData[order(predData[,predSite]),]
    bioData <- bioData[order(bioData[,bioSite]),]
    
    ##sets up species data for calculating dissimilarity
    bx <- which(names(bioData)==XColumn)
    by <- which(names(bioData)==YColumn)
    sppData <- bioData[-c(bioSite, bx, by)]
    
    ##creates distance matrix
    if(abundance==F){
      sppData[sppData>=1] <- 1
      sppData[sppData==0] <- 0
      sppData[is.na(sppData)] <- 0
      distData <- vegan::vegdist(sppData, dist, binary=T)
    }else{
      sppData[is.na(sppData)] <- 0
      distData <- vegan::vegdist(sppData, dist, binary=F)
    }
    
    ########################################################################
    ##species data as site-site distance matrix
  }else if(bioFormat==3){
    ##orders bioData to match ordering of predData below
    holdSite <- bioData[,which(siteColumn %in% colnames(bioData))]
    bioData <- bioData[,-which(siteColumn %in% colnames(bioData))]
    orderedData <- as.matrix(as.dist(bioData[order(holdSite),order(holdSite)]))
    
    ##site-site distance already calculated
    distData <- lower.tri(as.matrix(orderedData), diag=FALSE)
    distData <- as.vector(orderedData[distData])
    predData <- unique(predData)
    ##orders the prediction data by site
    predData <- predData[order(predData[siteColumn]),]
    ########################################################################
    ##site pair table, already preped 
  }else if(bioFormat==4){
    ##site-pair distance value
    outTable <- bioData
    ########################################################################
    
  }else{
    ##return error, bioFormat argument out of bounds
    stop(paste("bioFormat argument of '", as.character(bioFormat), "' is not an accepted input value", sep=""))
  }
  ########################################################################
  
  ##With the dissim distance calculated, creates and fills the table in gdm format
  if(bioFormat!=4){
    ##creates base site-pair table
    outTable <- as.data.frame(createsitepair(dist=distData, spdata=bioData, envInfo=predData, dXCol=XColumn, dYCol=YColumn, siteCol=siteColumn, weightsType=weightType, custWeights=custWeights))
  }else{
    outTable <- bioData
  }
  
  ##first checks that the size of the dissimilarity matrices, if any were provided
  ##then pastes any dissimilarity matrices onto the created site-pair table
  if(length(distPreds)>0){
    
    baseMat <- distPreds[[1]]
    ##checks to size of each dissimilarity matrix, to make sure they are all the same
    lapply(distPreds, function(mat, mat1){
      if((dim(mat1)[1]!=dim(mat)[1]) & (dim(mat1)[2]!=dim(mat)[2])){
        stop("The dimensions of your predictor matrices are not the same.")
      }
    }, mat1=baseMat)
    #print(dim(baseMat))
    
    ##hold site columns
    holdSiteCols <- lapply(distPreds, function(dP){dP[,which(siteColumn %in% colnames(dP))]})
    #remove site column from matrices
    distPreds <- lapply(distPreds, function(dP){dP[,-which(siteColumn %in% colnames(dP))]})
    #print(dim(distPreds[[1]]))
    ##orders the distance matices of distPreds
    distPreds <- mapply(function(dP, hSC){as.matrix(as.dist(dP[order(hSC),order(hSC)]))}, dP=distPreds, hSC=holdSiteCols, SIMPLIFY=FALSE)
    #print(dim(distPreds[[1]]))
    ##orders the site columns to match the distance matrices
    orderSiteCols <- lapply(holdSiteCols, function(hSC){hSC[order(hSC)]})

    ##compares sites with sites removed, for one reason or another
    rmSites <- c(toRemove, removeRand)
    ##removes the sites removed above when creating the site-pair table
    if(length(rmSites)>0){
      rmIndex <- lapply(orderSiteCols, function(hSC, tR){which((hSC %in% tR))}, tR=rmSites)
      distPreds <- mapply(function(mat, tR){mat[-c(tR), -c(tR)]}, mat=distPreds, tR=rmIndex, SIMPLIFY=FALSE)
    }
    ##set new baseMat
    baseMat <- distPreds[[1]]
    #print(dim(baseMat))
    
    ##checks the size of the dissimilarity matrices against the size of distData
    baseMatDat <- lower.tri(as.matrix(baseMat),diag=FALSE)
    baseMatDat <- as.vector(baseMat[baseMatDat])  
    #print(nrow(outTable))
    #print(length(baseMatDat))
    if(nrow(outTable)!=length(baseMatDat)){
      stop("The dimensions of the distance predictor matrices do not match the biological data.")
    }
    
    ##addes any associated distance predictors to the sitepair table
    for(num in 1:length(distPreds)){
      #num <- 2
      ##isolate matrix
      matrixDat <- lower.tri(as.matrix(distPreds[[num]], diag=FALSE))
      sweetFreakenPredMatrix <- as.vector(distPreds[[num]][matrixDat])
      ##add matrix to table
      if(ncol(outTable)>6){
        ##break table up to insert matrix data columns
        baseSitePair <- outTable[,1:6]
        otherSitePair <- outTable[,7:ncol(outTable)]
        otherNames <- colnames(otherSitePair)
        s1SitePair <- as.data.frame(otherSitePair[,1:(ncol(otherSitePair)/2)])
        colnames(s1SitePair) <- otherNames[1:(ncol(otherSitePair)/2)]
        s2SitePair <- as.data.frame(otherSitePair[,(ncol(otherSitePair)/2+1):ncol(otherSitePair)])
        colnames(s2SitePair) <- otherNames[(ncol(otherSitePair)/2+1):ncol(otherSitePair)]
        ##formats data from dissimilarity matrices 
        s1Zeros <- as.data.frame(rep(0,length(sweetFreakenPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, sep="")
        s2Mat <- as.data.frame(sweetFreakenPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep="")
        ##restructures the output table
        outTable <- cbind(baseSitePair, s1SitePair, s1Zeros, s2SitePair, s2Mat)
      }else{
        ##formats data from dissimilarity matrices 
        s1Zeros <- as.data.frame(rep(0,length(sweetFreakenPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, sep="")
        s2Mat <- as.data.frame(sweetFreakenPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep="")
        ##restructures the output table
        outTable <- cbind(outTable, s1Zeros, s2Mat)
      }
    }
  }
  
  ##return output table
  class(outTable) <- c("gdmData", "data.frame")
  return(outTable)
}
##########################################################################


##########################################################################
##Used in formatGDMData to transform data from a site-site distance matrix into
##a site pair format
createsitepair <- function(dist, spdata, envInfo, dXCol, dYCol, siteCol, 
                           weightsType, custWeights){
  ###########################
  ##lines used to quickly test function
  #dist = distData
  #spdata = bioData
  #envInfo = predData
  #dXCol = XColumn
  #dYCol = YColumn
  #siteCol = siteColumn
  #weightsType = weightType
  #custWeights = custWeights
  ###########################
  ##required libraries
  #require(raster)
  
  ##Create gdm ready table
  weightsType <- as.character(weightsType)
  distance <- as.vector(dist)
  ##calculates richness total, the sums of the two most populus sites
  if(weightsType[1]=="richness"){
    sppOnly <- spdata[-c(1,2,3)]
    sppSums <- rowSums(sppOnly)
    sppSiteSums <- cbind(spdata[1], sppSums)
    orderedSums <- sppSiteSums[order(-sppSiteSums[,2]),]
    richTotal <- orderedSums[1,2]+orderedSums[2,2]
  }
  
  ##Builds index needed for output gdm table format
  s1.xCoord <- s1.yCoord <- s2.xCoord <- s2.yCoord <- NULL
  s1 <- s2 <- NULL
  
  if((siteCol %in% colnames(envInfo))==T){
    count <- seq(length(unique(envInfo[,siteCol]))-1,1,-1)
  }else{
    count <- seq(length(unique(envInfo[,"siteUltimateCoolness"]))-1,1,-1)
  }
  s1 <- unlist(sapply(seq(length(count),1), function(y){c(s1, rep((max(count)-y)+1, times=y))}))
  s2 <- unlist(sapply(seq(length(count),1), function(y){c(s2, (max(count)-y+2):(max(count)+1))}))
  
  if(length(s1)!=length(distance)){
    stop("The length of distance values are not the same as the expected number of rows of the site-pair table, unable to proceed.")
  }
  
  if(weightsType[1]=="equal"){
    weights <- rep(1, times=length(distance))
  }else if(weightsType[1]=="custom"){
    weights <- (custWeights[s1, "weights"] + custWeights[s2, "weights"]) / 2
  }else{
    weights <- (sppSiteSums[s1, "sppSums"] + sppSiteSums[s2, "sppSums"]) / richTotal
  }
  gdmTable <- cbind(distance, weights)
  
  ##from environmental or species table, copy coordinates for site-pair table
  if((dXCol %in% colnames(envInfo))==T){
    if((siteCol %in% colnames(envInfo))==T){
      checkTab <- table(envInfo[siteCol])
    }else{
      checkTab <- table(envInfo["siteUltimateCoolness"])
    }
    
    if(sum(checkTab>1)>0){
      stop("A site has two or more unique entries of data associated with it. Check data for issues.")
    }
    s1.xCoord <- envInfo[s1, dXCol]
    s2.xCoord <- envInfo[s2, dXCol]
    s1.yCoord <- envInfo[s1, dYCol]
    s2.yCoord <- envInfo[s2, dYCol]  
  }else if((dXCol %in% colnames(spdata))==T){
    s1.xCoord <- spdata[s1, dXCol]
    s2.xCoord <- spdata[s2, dXCol]
    s1.yCoord <- spdata[s1, dYCol]
    s2.yCoord <- spdata[s2, dYCol]
  }else{
    stop("X,Y Coordinates not found with unique sites, unable to build site-pair table")
  }
  
  ##sets up output table
  gdmForm <- cbind(gdmTable, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord)
  xhold <- which(names(envInfo)==dXCol)
  yhold <- which(names(envInfo)==dYCol)
  sitehold <- which(names(envInfo)==siteCol)
  sitehold2 <- which(names(envInfo)=="siteUltimateCoolness")
  envInfo <- envInfo[-c(xhold, yhold, sitehold, sitehold2)]
  
  ##fills output table
  if(ncol(envInfo)>0){
    gdmTableFill <- cbind(gdmForm, envInfo[s1,1:ncol(envInfo)], envInfo[s2,1:ncol(envInfo)])
    names.s1 <- paste("s1.",names(envInfo[1:ncol(envInfo)]), sep="")
    names.s2 <- paste("s2.",names(envInfo[1:ncol(envInfo)]), sep="")
    colnames(gdmTableFill) <- c(colnames(gdmTableFill)[1:6], names.s1, names.s2)
  }else{
    gdmTableFill <- gdmForm
  }
  
  ##returns results
  return(gdmTableFill)
}
##########################################################################


##########################################################################
##Extracts Ispline data from a gdm model
isplineExtract <- function (model){
  ###########################
  #model = gdmOb
  ###########################
  ##error checking
  ##checks to make sure a gdm model is given
  if(class(model)[1]!="gdm"){
    stop("model argument must be a gdm model object")
  }
  
  ##Collects or sets simple data
  options(warn.FPU = FALSE)
  PSAMPLE <- 200
  preddata <- rep(0, times = PSAMPLE)
  pn <- model$predictors
  nPreds <- length(pn)
  yDat <- xDat <- matrix(0,PSAMPLE,nPreds)
  colnames(yDat) <- colnames(xDat) <- pn
  pmin <- 1
  pmax <- PSAMPLE
  
  ##cycles through each prodictor and fills the spline matrices
  splineindex <- 1
  for (i in 1:nPreds){ 
    #i<-1
    numsplines <- model$splines[i]
    z <- .C("GetPredictorPlotData", 
            pdata = as.double(preddata), 
            as.integer(PSAMPLE), 
            as.double(model$coefficients[splineindex:(splineindex + numsplines - 1)]), 
            as.double(model$knots[splineindex:(splineindex + numsplines - 1)]), 
            as.integer(numsplines),
            PACKAGE = "gdm")
    yDat[,i] <- z$pdata
    pmin <- pmin + PSAMPLE
    pmax <- pmax + PSAMPLE
    xDat[,i] <-  seq(from=model$knots[[(i*3)-2]],to=model$knots[[(i*3)]], length=PSAMPLE)
    splineindex <- splineindex + numsplines
  }
  
  ##lists and returns matrices
  outData <- list(x=xDat,y=yDat)
  return(outData)
}
##########################################################################

##########################################################################
plotUncertainty = function(spTable, leaveOut, bsIters, geo=FALSE, splines=NULL, 
                           knots=NULL, splineCol="blue", errCol="grey80", 
                           plot.linewidth=2.0, plot.layout=c(2,2), parallel=FALSE,
                           cores=2){
  ##A function to plot the uncertantiy of each variable from a GDM model
  #################
  #spTable <- sitePairTab          ##the input site-pair table to subsample from
  #leaveOut <- 0.1     ##percent of sites that should be randomly removed from site pair table 
  #bsIters <- 50       ##the number of time the site-pair table should be sampled
  #geo <- F              ##rather or not the gdm model takes geography into account, see gdm
  #splines <- NULL       ##splines gdm setting, see gdm
  #knots <- NULL         ##knots gdm setting, see gdm
  #splineCol <- "blue"    ##color of the center line
  #errCol <- "grey80"        ##color of the uncertainty polygon
  #plot.linewidth <- 2.0    ##line width of the center line
  #plot.layout <- c(3,3)    ##number of plots per page
  #parallel <- F            ##rather or not the sampling should happen in parallel processing, to speed it up
  #cores <- 6               ##number of cores to if parallel processing
  #################
  ##function breaks and warnings
  ##makes sure that table is a properly formated site-pair table
  if(class(spTable)[1] != "gdmData"){
    warning("table class does not include type 'gdmData'. Make sure your data is in site-pair format. See the formatsitepair fucntion for help.")
  }
  ##checks to makes sure table is a matrix or table frame
  if(!(class(spTable)[1]=="gdmData" | class(spTable)[1]=="matrix" | class(spTable)[1]=="data.frame")){
    stop("table argument needs to be a matrix or a table frame")
  }
  
  ##sanity check on the data table
  if(ncol(spTable) < 6){
    stop("Not enough columns in table. (Minimum need: Observed, weights, X0, Y0, X1, Y1)")
  }  
  if(nrow(spTable) < 1){
    stop("Not enough rows in table")
  }
  
  ##checks that geo has either TRUE or FALSE
  if(!(geo==TRUE | geo==FALSE)){
    stop("geo argument must be either TRUE or FALSE")
  }
  ##makes sure splines is a numeric vector
  if(is.null(splines)==FALSE & class(splines)!="numeric"){
    stop("argument splines needs to be a numeric data type")
  }
  ##checks knots inputs
  if(is.null(knots)==FALSE & class(knots)!="numeric"){
    stop("argument knots needs to be a numeric data type")
  }
  
  ##checks that parallel has either TRUE or FALSE
  if(!(parallel==TRUE | parallel==FALSE)){
    stop("parallel argument must be either TRUE or FALSE")
  }
  ##makes sure that cores has a value when parallel is true
  if(parallel==TRUE & is.null(cores)==TRUE){
    stop("If parallel==TRUE, the number of cores must be specified")
  }
  ##makes sure that cores is a positive integer 
  if((is.null(cores)==FALSE & is.numeric(cores)==FALSE) | cores<1){
    stop("argument cores needs to be a positive integer")
  }
  
  ##makes sure that bsIters is a positive integer 
  if((is.null(bsIters)==FALSE & is.numeric(bsIters)==FALSE) | bsIters<1){
    stop("argument bsIters needs to be a positive integer")
  }
  ##makes sure leaveOut is a number
  if(is.numeric(leaveOut)==FALSE){
    stop("leaveOut must be a number between 0 and 1")
  }
  ##makes sure that leaveOut is between 0 and 1
  if(leaveOut < 0){
    stop("leaveOut must be a number between 0 and 1")
  }
  if(leaveOut > 1){
    stop("leaveOut must be a number between 0 and 1")
  }
  
  ##assign k to prevent issues to cran checking
  k <- NULL
  
  ##runs parallel if desired by the users
  if(parallel==TRUE){
    ##loads libraries
    #require(foreach)
    #require(doParallel)
    #requireNamespace("foreach")
    #requireNamespace("parallel")
    #requireNamespace("doParallel")
    
    ##sets cores
    cl <- makeCluster(cores, outfile="")
    registerDoParallel(cl)
    ##first removes a number of sites according to input
    subSamps <- foreach(k=1:bsIters, .verbose=F, .packages=c("gdm")) %dopar%
      removeSitesFromSitePair(k, tab=spTable, rmFrac=leaveOut)
    ##models the subsamples
    gdmMods <- foreach(k=1:length(subSamps), .verbose=F, .packages=c("gdm")) %dopar%
                     #gdmMods <- try(foreach(k=1, .verbose=F, .packages=c("gdm")) %dopar%
                     gdm(subSamps[[k]], geo=geo, splines=splines, knots=knots)
    stopCluster(cl)
  }else{
    ##first removes a number of sites according to input
    subSamps <- lapply(1:bsIters, removeSitesFromSitePair, tab=spTable, rmFrac=leaveOut)
    ##models the subsamples
    gdmMods <- lapply(subSamps, gdm, geo=geo, splines=splines, knots=knots)
  }
  
  ##models the full gdm
  fullGDMmodel <- gdm(spTable, geo=geo, splines=splines, knots=knots)
    
  ##Extracts the splines for each model
  exUncertSplines <- lapply(gdmMods, isplineExtract)
  fullGDMsplines <- isplineExtract(fullGDMmodel)
  
  ##get the names of the predictor variables
  predVars <- colnames(exUncertSplines[[1]][[1]])
  
  ##establish what plot layout to use
  thisplot <- 0
  one_page_per_plot <- FALSE
  if ((plot.layout[1]==1) && (plot.layout[2]==1)){
    one_page_per_plot <- TRUE
  }else{
    par(mfrow=plot.layout)
  }
  
  ##sets the plotting minimum and maximum y-values 
  totalYmin <- Inf
  totalYmax <- -Inf
  
  ##determines the bounds of the plots
  for(p in 1:length(predVars)){
    #p=1
    predV <- predVars[p]
    
    ##gets the minimum and maximum, to set the ploting extent
    for(nm in 1:length(exUncertSplines)){
      #nm=1
      selPlot <- exUncertSplines[[nm]]
      
      spYmax <- max(selPlot[[2]][,predV])
      spYmin <- min(selPlot[[2]][,predV])
      
      totalYmax <- max(c(totalYmax, spYmax))
      totalYmin <- min(c(totalYmin, spYmin))
    }
  }
  
  ##plots by variable
  for(p in 1:length(predVars)){
    #p <- 1
    predV <- predVars[p]
    
    ##sets the plotting minimum and maximum x-values
    totalXmin <- Inf
    totalXmax <- -Inf
    ##gets the minimum and maximum, to set the ploting extent
    for(nm in 1:length(exUncertSplines)){
      #nm=1
      selPlot <- exUncertSplines[[nm]]
      
      spXmax <- max(selPlot[[1]][,predV])
      spXmin <- min(selPlot[[1]][,predV])
      
      if(spXmax > totalXmax){totalXmax = spXmax}
      if(spXmin < totalXmin){totalXmin = spXmin}
    }
    
    ##checks to make sure that there is some significance, if there is the the data is plotted
    if(totalYmax!=0){
      ##add mean of subsets to plot
      plotX <- NULL
      plotY <- NULL
      byVarMatX <- NULL
      byVarMatY <- NULL
      ##create matrices based on the variable and its x and y location for each model iteration
      for(nn in 1:length(exUncertSplines)){
        #nn=1
        plotX[[nn]] <- exUncertSplines[[nn]][[1]]
        plotY[[nn]] <- exUncertSplines[[nn]][[2]]
        byVarMatY <- cbind(byVarMatY, plotY[[nn]][,predV])
        byVarMatX <- cbind(byVarMatX, plotX[[nn]][,predV])
      }
      ##gets spline data from the full gdm model
      fullPlotX <- fullGDMsplines[[1]]
      fullPlotX <- fullPlotX[,predV]
      fullPlotY <- fullGDMsplines[[2]]
      fullPlotY <- fullPlotY[,predV]
      
      ##calculates the confidence intervals for plotting
      sdX <- apply(as.matrix(byVarMatX), 1, sd)
      sdY <- apply(as.matrix(byVarMatY), 1, sd)
      highBoundX <- fullPlotX + sdX
      lowBoundX <- fullPlotX - sdX
      highBoundY <- fullPlotY + sdY
      lowBoundY <- fullPlotY - sdY
      
      ##collects the data to be used in the rug plot
      if(predV=="Geographic"){
        ##calculates unique eucildian distance between sites
        rugData <- unique(sqrt(((spTable$s1.xCoord-spTable$s2.xCoord)^2)+((spTable$s1.yCoord-spTable$s2.yCoord)^2)))
      }else{
        ##gets unique values of variable data
        varDat <- grep(predV, colnames(spTable))
        rugData <- unique(c(spTable[,c(varDat[1])], spTable[,c(varDat[2])]))
      }
      
      ##plots one graph per page, unless specified otherwise
      if (one_page_per_plot){
        dev.new()
        dev.next()
      }else{
        thisplot <- thisplot + 1
        if(thisplot > (plot.layout[1] * plot.layout[2])){    			
          thisplot <- 1
          par(mfrow=plot.layout)	
        }
      }
      #settings <- par(pars)
      ##plots mean data line and polygon of uncertanty 
      plot(NULL, xlim=c(totalXmin, totalXmax), ylim=c(totalYmin, totalYmax), xlab=predV, ylab="Partial Ecological Distance")
      polygon(c(lowBoundX, rev(highBoundX)), c(lowBoundY, rev(highBoundY)), col=errCol, border=NA)
      lines(fullPlotX, fullPlotY, col=splineCol, lwd=plot.linewidth)
      rug(rugData)
    }
  }
}
##########################################################################

########################################################################## 
removeSitesFromSitePair <- function(n, tab, rmFrac){
  ##a function to remove a random number of sites from a sitepair table
  ##involves assigning an index to each site, picking the indicies to be 
  ##removed, then identifying which site pairs those indices are a part of
  ##and remove those site pairs from the table
  #################
  #n <- 1                    ##a number, not actually used, but needed to preform lapply
  #tab <- gdmTab             ##sitepair table
  #rmFrac <- leaveOut    ##fraction of sites to remove from sitepair table
  #################
  ##First create an index of each site in the JTable
  sortMat<-matrix(NA,(nrow(tab)*2),5)
  ##Fill the sorting table
  i_site<-1
  for(i_row in 1:nrow(tab)){
    sortMat[i_site,1] <- tab[i_row,3] 
    sortMat[i_site,2] <- tab[i_row,4]
    sortMat[i_site,3] <- 1
    sortMat[i_site,4] <- i_row
    i_site <- i_site+1
    sortMat[i_site,1] <- tab[i_row,5] 
    sortMat[i_site,2] <- tab[i_row,6]
    sortMat[i_site,3] <- 2
    sortMat[i_site,4] <- i_row
    i_site <- i_site+1  
  }
  
  ##Now sort the sorting table by Long then Lat
  sortMat <- sortMat[order(sortMat[,1], sortMat[,2]),]
  ##loop through and give each site an index based on whether it matches the 
  ##Lat & Long of the site above it
  sortMat[1,5]<-1
  for(i_row in 2:nrow(sortMat)){
    if(sortMat[i_row,1] == sortMat[(i_row-1),1]){
      if(sortMat[i_row,2] == sortMat[(i_row-1),2]){
        sortMat[i_row,5] <- sortMat[(i_row-1),5]
      }else{
        sortMat[i_row,5] <- sortMat[(i_row-1),5] + 1
      } ## end else
    }else{
      sortMat[i_row,5] <- sortMat[(i_row-1),5] + 1
    }
  }
  ##And finally create the indexTab, table to indexes of sites
  indexTab <- matrix(NA,nrow(tab),2)
  for(i_row in 1:nrow(sortMat)){
    indexTab[sortMat[i_row,4],sortMat[i_row,3]] <- sortMat[i_row,5]
  }
  
  ##determines the number of sites to remove
  numToRemove <- round(max(indexTab)*rmFrac)
  ##randomly determines the index of sites to remove
  rmSites <- sample(1:max(indexTab), numToRemove)
  rmIndexCol1 <- which(indexTab[,1] %in% rmSites)
  rmIndexCol2 <- which(indexTab[,2] %in% rmSites)
  ##creates sampled table
  sampTable <- tab[-c(unique(c(rmIndexCol1, rmIndexCol2))),]
  
  ##returns the sampled table
  return(sampTable)
}
##########################################################################


##########################################################################
##for single species site pair, when ready
##########################################################################


##########################################################################
gdm.varImp <- function(spTable, geo, splines=NULL, knots=NULL, fullModelOnly=FALSE, 
                                       nPerm=100, parallel=FALSE, cores=2){
  ##function to test the signifiance of the various variables of a gdm site-pair table
  #################
  #spTable <- sitePairTab          ##the input site-pair table to subsample from
  #load("M:/UAE/kavyaWorking/Code/GDM/GDMSitepairTable.RData")
  #spTable <- gdmSitea
  #geo <- T              ##rather or not the gdm model takes geography into account, see gdm
  #splines <- NULL       ##splines gdm setting, see gdm
  #knots <- NULL         ##knots gdm setting, see gdm
  #fullModelOnly <- F         ##not sure about this argument, acceptable values are TRUE and FALSE
  #nPerm <- 2
  #parallel <- T
  #cores <- 2
  #################
  ##things to be examined later
  ##spline and knot inputs while variable testing
  ##change values of wieghts when changing site-pair table site pairs
  #################
  ##assign k to prevent issues to cran checking
  k <- NULL
  
  ##number of variables in the site-pair table, adds 1 if geo is to be TRUE
  nVars <- (ncol(spTable)-6)/2
  ##collects variable names
  varNames <- colnames(spTable[,c(7:(6+nVars))])
  varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
  if(geo==TRUE){
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  
  ##First create a spTable to determine the index of each site in the site-pair table
  sortMatX <- sapply(1:nrow(spTable), function(i, spTab){c(spTab[i,3], spTab[i,5])}, spTab=spTable)
  sortMatY <- sapply(1:nrow(spTable), function(i, spTab){c(spTab[i,4], spTab[i,6])}, spTab=spTable)
  sortMatNum <- sapply(1:nrow(spTable), function(i){c(1,2)})
  sortMatRow <- sapply(1:nrow(spTable), function(i){c(i,i)})
  ##adds a column of NA for index to be added to
  fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY), as.vector(sortMatNum), as.vector(sortMatRow), rep(NA, length(sortMatX)))
  ##assigns sites by unique coordinates
  siteByCoords <- unique(fullSortMat[,1:2])
  ##number of sites to expect by coordinates
  numSites <- nrow(siteByCoords)
  ##assigns site index based on coordinates
  for(i in 1:numSites){
    fullSortMat[which(fullSortMat[,1]==siteByCoords[i,1] & fullSortMat[,2]==siteByCoords[i,2]),5] <- i
  }
  
  ##create index table to know where each site is in input site-pair table
  indexTab <- matrix(NA,nrow(spTable),2)
  for(iRow in 1:nrow(fullSortMat)){
    indexTab[fullSortMat[iRow,4],fullSortMat[iRow,3]] <- fullSortMat[iRow,5]
  }
  ## And remove the sorting table and supporting objects to free up memory
  rm(fullSortMat)
  rm(sortMatX)
  rm(sortMatY)
  rm(sortMatNum)
  rm(sortMatRow)
  rm(siteByCoords)
  
  ##create siteXvar table, to be able to rebuild site-pair table later in function
  exBySite <- lapply(1:numSites, function(i, index, tab){rowSites <- which(index[,1] %in% i)
  if(length(rowSites)<1){
    rowSites <- which(index[,2] %in% i)
  }
  exSiteData <- tab[rowSites[1],]
  return(exSiteData)
  }, index=indexTab, tab=spTable)
  ##extracts the data from the site-pair table by site
  siteData <- lapply(exBySite, function(row){row[grep("s1.", colnames(row))]})
  ##identifies the one site not in the first column of the index table                                                 
  outSite <- which(!(1:numSites %in% indexTab[,1]))
  ##removes the site that did not appear in the first column
  siteData[[outSite]] <- NULL
  ##transforms data from list to table
  siteData <- do.call("rbind", siteData)
  ##removes 's1.' from column names
  colnames(siteData) <- sapply(strsplit(colnames(siteData), "s1."), "[[", 2)
  
  ##adding the site removed above, with correct values
  outSiteData <- exBySite[[outSite]][grep("s2.", colnames(exBySite[[outSite]]))]
  colnames(outSiteData) <- sapply(strsplit(colnames(outSiteData), "s2."), "[[", 2)
  siteData <- rbind(siteData, outSiteData)
  
  ##sets up objects to be returned by the function
  modelTestValues <- matrix(NA,3,nVars,dimnames = list(c("Model deviance", "Percent deviance explained", "Model p-value"),c("fullModel", paste("fullModel-", seq(1,nVars-1), sep=""))))
  ##deviance reduction variable table, not yet sure why imporant
  devReductVars <- matrix(NA, nVars, nVars-1)
  rownames(devReductVars) <- varNames
  colnames(devReductVars) <- c("fullModel", paste("fullModel-", seq(1,nVars-2), sep=""))
  ##p value variable table, not yet sure why imporant
  pValues <- devReductVars
  
  ##assigns given site-pair table to new variable, to prevent changing the original input
  currSitePair <- spTable
  
  for(v in 1:nVars){
    #v <- 2
    print(varNames[v])

    ##runs gdm, first time on full site-pair table
    ##however as variables are removed the "full" site-pair table will have less varialbes in it
    fullGDM <- gdm(currSitePair, geo=geo, splines=splines, knots=knots)
    
    ##create a series of permutated site-pair tables, randomized site comparisons
    if(parallel == TRUE){
      #require(foreach)
      #require(doParallel)
      
      ##sets cores
      cl <- makeCluster(cores, outfile="")
      registerDoParallel(cl)
      
      permSitePairs <- foreach(k=1:nPerm, .verbose=F, .packages=c("gdm"), .export=c("permutateSitePair")) %dopar%
        permutateSitePair(currSitePair, siteData, indexTab, varNames)
      
      permGDM <- try(foreach(k=1:length(permSitePairs), .verbose=F, .packages=c("gdm")) %dopar%
                       gdm(permSitePairs[[k]], geo=geo, splines=NULL, knots=NULL))
      ##closes cores
      stopCluster(cl)
    }else{
      ##non-parallel version
      permSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn){permutateSitePair(csp,sd,it,vn)}, 
                              csp=currSitePair, sd=siteData, it=indexTab, vn=varNames)
      permGDM <- lapply(permSitePairs, gdm, geo=geo, splines=NULL, knots=NULL)
    }
    
    ##runs gdm on the permuted tables
    permModelDev <- sapply(permGDM, function(mod){mod$gdmdeviance})
    
    ##begins to fill in the output table with data from fully fitted model
    modelTestValues[1,v] <- fullGDM$gdmdeviance
    modelTestValues[2,v] <- fullGDM$explained
    #p-value
    modelTestValues[3,v] <- sum(permModelDev<=fullGDM$gdmdeviance)/nPerm
    
    ##ends the loop if only 1 variable was used in the model
    if(length(varNames)<2){
      #if(length(varNames)<1){
      break
    }
    
    ##begins running tests on variations
    ##runs model without geo if geo was part of the model
    if(geo==TRUE){
      noGeoGDM <- gdm(currSitePair, geo=FALSE, splines=NULL, knots=NULL)
      
      ##create a series of permutated site-pair tables, randomized site comparisons
      if(parallel == TRUE){
        ##sets cores
        cl <- makeCluster(cores, outfile="")
        registerDoParallel(cl)
        
        permSitePairs <- foreach(k=1:nPerm, .verbose=F, .packages=c("gdm"), .export=c("permutateSitePair")) %dopar%
          permutateSitePair(currSitePair, siteData, indexTab, varNames)
        ##closes cores
        stopCluster(cl)
      }else{
        ##non-parallel version
        permSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn){permutateSitePair(csp,sd,it,vn)}, 
                                csp=currSitePair, sd=siteData, it=indexTab, vn=varNames)
      }
      
      ##runs gdm on the permuted tables
      permGDM <- lapply(permSitePairs, gdm, geo=geo, splines=NULL, knots=NULL)
      ##permutations for geographic, and adds them to output objects
      permDevReduct <- sapply(permGDM, function(mod, ngeo){ngeo$gdmdeviance - mod$gdmdeviance}, ngeo=noGeoGDM)
      devReductVars[1,v] <- noGeoGDM$gdmdeviance - fullGDM$gdmdeviance
      pValues[1,v] <- sum(permDevReduct>=devReductVars[1,v])/nPerm
    }
    
    ##now tests all other variables 
    for(varChar in varNames){
      #varChar <- varNames[2]
      if(varChar!="Geographic"){
        ##select variable columns to be removed from original site-pair table
        testVarCols1 <- grep(paste("^s1.", varChar, "$", sep=""), colnames(currSitePair))
        testVarCols2 <- grep(paste("^s2.", varChar, "$", sep=""), colnames(currSitePair))
        testSitePair <- currSitePair[,-c(testVarCols1, testVarCols2)]
        ##run gdm for the missing variable
        noVarGDM <- gdm(testSitePair, geo=geo, splines=NULL, knots=NULL)
        
        ##create a series of permutated site-pair tables, randomized site comparisons
        if(parallel == TRUE){
          ##sets cores
          cl <- makeCluster(cores, outfile="")
          registerDoParallel(cl)
          
          noVarSitePairs <- foreach(k=1:nPerm, .verbose=F, .packages=c("gdm"), .export=c("permutateVarSitePair")) %dopar%
            permutateVarSitePair(currSitePair, siteData, indexTab, varChar)
          ##closes cores
          stopCluster(cl)
        }else{
          ##non-parallel version
          noVarSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn){permutateVarSitePair(csp,sd,it,vn)}, 
                                   csp=currSitePair, sd=siteData, it=indexTab, vn=varChar)
        }
        
        ##runs gdm on the permuted tables
        permGDM <- lapply(noVarSitePairs, gdm, geo=geo, splines=NULL, knots=NULL)
        ##permutations for geographic, and adds them to output objects
        permDevReduct <- sapply(permGDM, function(mod, nvar){nvar$gdmdeviance - mod$gdmdeviance}, nvar=noVarGDM)
        devReductVars[which(rownames(devReductVars) %in% varChar),v] <- noVarGDM$gdmdeviance - fullGDM$gdmdeviance
        #print(varChar)
        pValues[which(rownames(pValues) %in% varChar),v] <- sum(permDevReduct>=devReductVars[which(rownames(devReductVars) %in% varChar),v])/nPerm
      }
    }
    
    ##if fullModelOnly == TRUE, breaks script
    if(fullModelOnly==TRUE){
      break
    }
    
    ##based on the P-value, and then deviance reduction, select the variable to be omitted
    ##from future iterations of the testing
    tempPVals <- as.numeric(pValues[c(1:nVars),v])
    tempDevs <- as.numeric(devReductVars[c(1:nVars),v])
    tempPVals <- tempPVals[!is.na(tempPVals)]
    tempDevs <- tempDevs[!is.na(tempDevs)]
    #print(length(tempPVals))
    #print(length(tempDevs))
    varToOmit <- which.max(tempPVals)
    
    for(iCheck in 1:length(varNames)){
      if(tempPVals[iCheck] == tempPVals[varToOmit]){
        if(tempDevs[iCheck] < tempDevs[varToOmit]){
          varToOmit <- iCheck
        }
      }
    }
    
    ##removes variables 
    ##if selected, removes geo
    if(varToOmit==1 & geo==TRUE){
      geo <- FALSE
      varNames <- varNames[-1]
    }else{
      ##removes any non-geo variable
      nameToRemove <- varNames[varToOmit]
      #nameToRemove <- "bio1"
      ##remove from variables
      varNames <- varNames[-varToOmit]
      removeFromSitePs1 <- grep(paste("^s1.", nameToRemove, "$", sep=""), colnames(currSitePair))
      removeFromSitePs2 <- grep(paste("^s2.", nameToRemove, "$", sep=""), colnames(currSitePair))
      ##removes variable from important objects
      currSitePair <- currSitePair[,-c(removeFromSitePs1,removeFromSitePs2)]
    }
  }
  ##lists tables into one object
  outObject <- list(modelTestValues, devReductVars, pValues)
  return(outObject)
}
##########################################################################

##########################################################################
permutateSitePair <- function(spTab, siteVarTab, indexTab, vNames){
  ##a function to randomize the rows of a site-pair table
  #################
  #spTab <- currSitePair    ##site-pair table
  #siteVarTab <- siteData   ##siteXvar table
  #indexTab <- indexTab     ##table of the index of sites
  #vNames <- varNames       ##variables names
  #vNames <- c("awcA", "phTotal", "shcA", "solumDepth", "bio5", "bio19")
  #################
  ##randomizes the row order of the given siteXvar table
  randVarTab <- siteVarTab[sample(nrow(siteVarTab), nrow(siteVarTab)), ]
  
  ##sets up the coordinate values for the randomized site-pair table
  s1xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],1]})
  s1yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],2]})
  s2xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],1]})
  s2yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],2]})
  
  #print(vNames)
  ##extracts values of other variables 
  varLists <- lapply(vNames, function(vn, rvTab, spt, inT){if(vn!="Geographic"){
    ###################
    #vn <- vNames[[2]]
    #rvTab=randVarTab
    #spt=spTab
    #inT=indexTab
    ###################
    ##identifies variable columns in randVarTab
    randCols <- grep(paste("^", vn, "$", sep=""), colnames(rvTab))
    #print(randCols)
    ##identifies variable columns in site-pair table
    spCols <- grep(vn, colnames(spt))
    
    s1var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,1],randCols]})
    s2var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,2],randCols]})
    
    return(list(s1var, s2var))
  }
  }, rvTab=randVarTab, spt=spTab, inT=indexTab)
  ##'unravels' the varList into a data.frame of the variable portion of a site-pair table
  bySite <- lapply(1:2, function(i,vlist){sapply(vlist, function(vl,k){vl[[k]]},k=i)}, vlist=varLists)
  
  if(class(bySite[[1]])=="list"){
    site1Vars <- do.call("cbind", bySite[[1]])
    site2Vars <- do.call("cbind", bySite[[2]])
  }else{
    site1Vars <- bySite[[1]]
    site2Vars <- bySite[[2]]
  }
  ##sets up new site-pair table
  newSP <- as.data.frame(cbind(spTab$distance, spTab$weights, s1xCoord, s1yCoord, s2xCoord, s2yCoord, site1Vars, site2Vars))
  colnames(newSP) <- colnames(spTab)
  class(newSP) <- c(class(spTab))
  return(newSP)
}
##########################################################################
permutateVarSitePair <- function(spTab, siteVarTab, indexTab, vName){
  ##only randomizes the values for a particular variable 
  #################
  #spTab <- currSitePair    ##site-pair table
  #siteVarTab <- siteData   ##siteXvar table
  #indexTab <- indexTab     ##table of the index of sites
  #vName <- varChar         ##variables names
  #################
  ##randomizes the row order of the given siteXvar table
  randVarTab <- siteVarTab[sample(nrow(siteVarTab), nrow(siteVarTab)), ]
  
  ##identifies variable columns in randVarTab
  randCols <- grep(paste("^", vName, "$", sep=""), colnames(randVarTab))
  ##identifies variable columns in site-pair table
  spCols1 <- grep(paste("^s1.", vName, "$", sep=""), colnames(spTab))
  spCols2 <- grep(paste("^s2.", vName, "$", sep=""), colnames(spTab))
  
  ##extracts values based on new index position
  s1var <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],randCols]})
  s2var <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],randCols]})
  ##places values back into site-pair table
  spTab[,spCols1] <- s1var
  spTab[,spCols2] <- s2var
  
  return(spTab)
}
##########################################################################


