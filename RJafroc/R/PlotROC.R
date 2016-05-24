#' @import ggplot2
PlotROC <- function(dataset, plottingModalities, plottingReaders, legendPosition) {
  NL <- dataset$NL
  LL <- dataset$LL
  lesionNum <- dataset$lesionNum
  maxNL <- dim(NL)[4]
  dataType <- dataset$dataType
  modalityID <- dataset$modalityID
  readerID <- dataset$readerID
  
  if (length(dim(NL)) != 4 || length(dim(LL)) != 4) 
    stop("The dimension of NL or LL is not correct. ")
  
  if (!is.list(plottingModalities) && !is.list(plottingReaders)) {
    ROCPoints <- CalculateROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders)
  } else {
    if (is.list(plottingModalities) && is.list(plottingReaders) && length(plottingModalities) == length(plottingReaders)) {
      ROCPoints <- data.frame(FPF = NULL, TPF = NULL)
      for (i in 1:length(plottingReaders)) {
        if (length(plottingModalities[[i]]) == 1 && (length(plottingReaders[[i]]) == 1)) {
          tempROCPoints <- CalculateROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders[[i]])
          ROCPoints <- rbind(ROCPoints, tempROCPoints)
        } else {
          tempROCPoints <- CalculateAvgROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders[[i]])
          ROCPoints <- rbind(ROCPoints, tempROCPoints)
        }
      }
    } else if (is.list(plottingModalities) && length(plottingReaders) == 1) {
      plottingReaders <- plottingReaders[[1]]
      ROCPoints <- data.frame(FPF = NULL, TPF = NULL)
      for (i in 1:length(plottingModalities)) {
        if (length(plottingModalities[[i]]) == 1) {
          tempROCPoints <- CalculateROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders)
          ROCPoints <- rbind(ROCPoints, tempROCPoints)
        } else {
          tempROCPoints <- CalculateAvgROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders)
          ROCPoints <- rbind(ROCPoints, tempROCPoints)
        }
      }
    } else if (is.list(plottingReaders) && length(plottingModalities) == 1) {
      if (is.list(plottingModalities)) 
        plottingModalities <- plottingModalities[[1]]
      ROCPoints <- data.frame(FPF = NULL, TPF = NULL)
      for (i in 1:length(plottingReaders)) {
        if (length(plottingReaders[[i]]) == 1) {
          tempROCPoints <- CalculateROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders[[i]])
          ROCPoints <- rbind(ROCPoints, tempROCPoints)
        } else {
          tempROCPoints <- CalculateAvgROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders[[i]])
          ROCPoints <- rbind(ROCPoints, tempROCPoints)
        }
      }
    } else {
      stop("Lengths of trts and rdrs do not match.")
    }
  }
  
  classes <- unique(ROCPoints$class)
  if (missing(legendPosition)){
    tooSmall <- FALSE
    for (i in 1:length(classes)){
      indices <- ROCPoints$class == classes[i]
      maxTPF <- max(ROCPoints$TPF[indices])
      if (maxTPF < 0.25){
        tooSmall <- TRUE
        legendPosition <- "bottom"        
      }
      if (sum(indices) > 20 && all(ROCPoints$type[indices] == "individual")){
        typeTmp <- as.character(ROCPoints$type)
        typeTmp[indices] <- "continuous"
        ROCPoints$type <- typeTmp
      }
    }
    if (!tooSmall){
      legendPosition <- c(1, 0)
    }
  }
  
  if (legendPosition == "right" || legendPosition == "left" ){
    legDir <- "vertical"
  }else{
    legDir <- "horizontal"
  }
  
  if (!is.list(plottingModalities) && !is.list(plottingReaders)) {
    mr <- unlist(strsplit(as.character(ROCPoints$class), split = "\n"))
    dim(mr) <- c(2, length(mr)/2)
    ROCPoints <- cbind(ROCPoints, data.frame(Modality = mr[1, ], Reader = mr[2, ]))
    opratingPoints <- ROCPoints[ROCPoints$type == "individual" & !((ROCPoints$FPF == 0 & ROCPoints$TPF == 0) | (ROCPoints$FPF == 1 & ROCPoints$TPF == 1)), ]
    legendLength <- length(plottingReaders)
    shapeVector <- rep(NA, length(plottingReaders))
    for (n in 1:legendLength) {
      index <- which(ROCPoints$Reader == levels(ROCPoints$Reader)[n])[1]
      if (ROCPoints$type[index] == "individual") 
        shapeVector[n] <- 16
    }
    if (length(plottingModalities) == 1){
      ROCPlot <- with(ROCPoints, {
        ROCPlotTemp <- ggplot()
        mStrings <- unique(as.character(ROCPoints$Modality))
        for (i in 1:length(plottingModalities)) {
          ROCPlotTemp <- ROCPlotTemp + geom_line(data = ROCPoints, aes(x = FPF, y = TPF, color = class), size = 1)
        }
        ROCPlotTemp <- ROCPlotTemp + geom_point(data = opratingPoints, size = 4, aes(x = FPF, y = TPF, color = class)) +
          theme(legend.title = element_blank(), legend.position = legendPosition, legend.direction = legDir, legend.justification = c(1, 0)) + 
          guides(color = guide_legend(override.aes = list(shape = shapeVector))) + 
          scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
      })
    } else {      
      ROCPlot <- with(ROCPoints, {
        ROCPlotTemp <- ggplot()
        mStrings <- unique(as.character(ROCPoints$Modality))
        for (i in 1:length(plottingModalities)) {
          ROCPlotTemp <- ROCPlotTemp + geom_line(data = ROCPoints[ROCPoints$Modality == mStrings[i], ], aes(x = FPF, y = TPF, color = Reader, linetype = Modality), size = 1)
        }
        ROCPlotTemp <- ROCPlotTemp + geom_point(data = opratingPoints, size = 4, aes(x = FPF, y = TPF, color = Reader)) +
          theme(legend.title = element_blank(), legend.position = legendPosition, legend.direction = legDir, legend.justification = c(1, 0)) + 
          guides(color = guide_legend(override.aes = list(shape = shapeVector))) + 
          scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
      })
    }
    ROCPoints <- data.frame(FPF = ROCPoints$FPF, TPF = ROCPoints$TPF, class = ROCPoints$class, type = ROCPoints$type)
  } else {
    opratingPoints <- ROCPoints[ROCPoints$type == "individual" & !((ROCPoints$FPF == 0 & ROCPoints$TPF == 0) | (ROCPoints$FPF == 1 & ROCPoints$TPF == 1)), ]
    
    legendLength <- length(levels(ROCPoints$class))
    shapeVector <- rep(NA, length(levels(ROCPoints$class)))
    for (n in 1:legendLength) {
      index <- which(ROCPoints$class == levels(ROCPoints$class)[n])[1]
      if (ROCPoints$type[index] == "individual") 
        shapeVector[n] <- 16
    }
    ROCPlot <- with(ROCPoints, {
      ggplot(data = ROCPoints, aes(x = FPF, y = TPF, color = class)) + geom_line(size = 1) + geom_point(data = opratingPoints, size = 4) + 
        theme(legend.title = element_blank(), legend.position = legendPosition, legend.direction = legDir, legend.justification = c(1, 0)) + 
        guides(color = guide_legend(override.aes = list(shape = shapeVector))) + 
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    })
  }
  
  # print(ROCPlot)
  return(list(ROCPlot = ROCPlot, ROCPoints = ROCPoints))
} 
