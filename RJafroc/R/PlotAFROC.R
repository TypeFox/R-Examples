#' @import ggplot2
PlotAFROC <- function(dataset, plottingModalities, plottingReaders, legendPosition) {
  NL <- dataset$NL
  LL <- dataset$LL
  lesionNum <- dataset$lesionNum
  maxNL <- dim(NL)[4]
  dataType <- dataset$dataType
  modalityID <- dataset$modalityID
  readerID <- dataset$readerID
  
  if (length(dim(NL)) != 4 || length(dim(LL)) != 4) 
    stop("The dimension of NL or LL is not corretc. ")
  
  if (!is.list(plottingModalities) && !is.list(plottingReaders)) {
    AFROCPoints <- CalculateAFROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders)
  } else {
    if (is.list(plottingModalities) && is.list(plottingReaders) && length(plottingModalities) == length(plottingReaders)) {
      AFROCPoints <- data.frame(FPF = NULL, LLF = NULL)
      for (i in 1:length(plottingReaders)) {
        if (length(plottingModalities[[i]]) == 1 && (length(plottingReaders[[i]]) == 1)) {
          tempAFROCPoints <- CalculateAFROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders[[i]])
          AFROCPoints <- rbind(AFROCPoints, tempAFROCPoints)
        } else {
          tempAFROCPoints <- CalculateAvgAFROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders[[i]])
          AFROCPoints <- rbind(AFROCPoints, tempAFROCPoints)
        }
      }
    } else if (is.list(plottingModalities) && length(plottingReaders) == 1) {
      plottingReaders <- plottingReaders[[1]]
      AFROCPoints <- data.frame(FPF = NULL, LLF = NULL)
      for (i in 1:length(plottingModalities)) {
        if (length(plottingModalities[[i]]) == 1) {
          tempAFROCPoints <- CalculateAFROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders)
          AFROCPoints <- rbind(AFROCPoints, tempAFROCPoints)
        } else {
          tempAFROCPoints <- CalculateAvgAFROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities[[i]], plottingReaders)
          AFROCPoints <- rbind(AFROCPoints, tempAFROCPoints)
        }
      }
    } else if (is.list(plottingReaders) && length(plottingModalities) == 1) {
      if (is.list(plottingModalities)) 
        plottingModalities <- plottingModalities[[1]]
      AFROCPoints <- data.frame(FPF = NULL, LLF = NULL)
      for (i in 1:length(plottingReaders)) {
        if (length(plottingReaders[[i]]) == 1) {
          tempAFROCPoints <- CalculateAFROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders[[i]])
          AFROCPoints <- rbind(AFROCPoints, tempAFROCPoints)
        } else {
          tempAFROCPoints <- CalculateAvgAFROCPoints(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders[[i]])
          AFROCPoints <- rbind(AFROCPoints, tempAFROCPoints)
        }
      }
    } else {
      stop("Lengths of trts and rdrs do not match.")
    }
  }
  
  classes <- unique(AFROCPoints$class)
  if (missing(legendPosition)){
    tooSmall <- FALSE
    for (i in 1:length(classes)){
      indices <- AFROCPoints$class == classes[i]
      maxLLF <- max(AFROCPoints$LLF[indices])
      if (maxLLF < 0.25){
        tooSmall <- TRUE
        legendPosition <- "bottom"        
      }
      if (sum(indices) > 20 && all(AFROCPoints$type[indices] == "individual")){
        typeTmp <- as.character(AFROCPoints$type)
        typeTmp[indices] <- "continuous"
        AFROCPoints$type <- typeTmp
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
    mr <- unlist(strsplit(as.character(AFROCPoints$class), split = "\n"))
    dim(mr) <- c(2, length(mr)/2)
    AFROCPoints <- cbind(AFROCPoints, data.frame(Modality = mr[1, ], Reader = mr[2, ]))
    opratingPoints <- AFROCPoints[AFROCPoints$type == "individual" & !((AFROCPoints$FPF == 0 & AFROCPoints$LLF == 0) | (AFROCPoints$FPF == 1 & AFROCPoints$LLF == 1)), ]
    legendLength <- length(plottingReaders)
    shapeVector <- rep(NA, length(plottingReaders))
    for (n in 1:legendLength) {
      index <- which(AFROCPoints$Reader == levels(AFROCPoints$Reader)[n])[1]
      if (AFROCPoints$type[index] == "individual") 
        shapeVector[n] <- 16
    }
    if (length(plottingModalities) == 1){
      AFROCPlot <- with(AFROCPoints, {
        AFROCPlotTemp <- ggplot()
        mStrings <- unique(as.character(AFROCPoints$Modality))
        for (i in 1:length(plottingModalities)) {
          AFROCPlotTemp <- AFROCPlotTemp + geom_line(data = AFROCPoints, aes(x = FPF, y = LLF, color = class), size = 1)
        }
        AFROCPlotTemp <- AFROCPlotTemp + geom_point(data = opratingPoints, size = 4, aes(x = FPF, y = LLF, color = class)) +
          theme(legend.title = element_blank(), legend.position = legendPosition, legend.direction = legDir, legend.justification = c(1, 0)) + 
          guides(color = guide_legend(override.aes = list(shape = shapeVector))) + 
          scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
      })
    } else {      
      AFROCPlot <- with(AFROCPoints, {
        AFROCPlotTemp <- ggplot()
        mStrings <- unique(as.character(AFROCPoints$Modality))
        for (i in 1:length(plottingModalities)) {
          AFROCPlotTemp <- AFROCPlotTemp + geom_line(data = AFROCPoints[AFROCPoints$Modality == mStrings[i], ], aes(x = FPF, y = LLF, color = Reader, linetype = Modality), size = 1)
        }
        AFROCPlotTemp <- AFROCPlotTemp + geom_point(data = opratingPoints, size = 4, aes(x = FPF, y = LLF, color = Reader)) +
          theme(legend.title = element_blank(), legend.position = legendPosition, legend.direction = legDir, legend.justification = c(1, 0)) + 
          guides(color = guide_legend(override.aes = list(shape = shapeVector))) + 
          scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
      })
    }
    AFROCPoints <- data.frame(FPF = AFROCPoints$FPF, LLF = AFROCPoints$LLF, class = AFROCPoints$class, type = AFROCPoints$type)
  } else {
    opratingPoints <- AFROCPoints[AFROCPoints$type == "individual" & !((AFROCPoints$FPF == 0 & AFROCPoints$LLF == 0) | (AFROCPoints$FPF == 1 & AFROCPoints$LLF == 1)), ]
    
    legendLength <- length(levels(AFROCPoints$class))
    shapeVector <- rep(NA, length(levels(AFROCPoints$class)))
    for (n in 1:legendLength) {
      index <- which(AFROCPoints$class == levels(AFROCPoints$class)[n])[1]
      if (AFROCPoints$type[index] == "individual") 
        shapeVector[n] <- 16
    }
    
    AFROCPlot <- with(AFROCPoints, {
      ggplot(data = AFROCPoints, aes(x = FPF, y = LLF, color = class)) + geom_line(size = 1) + geom_point(data = opratingPoints, size = 4) + 
        theme(legend.title = element_blank(), legend.position = legendPosition, legend.direction = legDir, legend.justification = c(1, 0)) + 
        guides(color = guide_legend(override.aes = list(shape = shapeVector))) + 
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    })
  }
  
  # print(AFROCPlot)
  return(list(AFROCPlot = AFROCPlot, AFROCPoints = AFROCPoints))
} 
