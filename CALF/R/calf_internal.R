calf_internal <- function(data,
                          nMarkers,
                          randomize  = FALSE,
                          proportion = NULL,
                          times){
  # getting rid of global variable warning -------------------------- #
  x = NULL
  y = NULL
  refx = NULL
  refy = NULL

  # setting up some initial values -----------------------------------#
  if (any(apply(data, 2, is.numeric) == FALSE)) {
    stop("CALF ERROR: Data are not numeric. Please check that data were read in correctly.")
  }

  nVars <- ncol(data) - 1
  dNeg  <- data[ ,2:ncol(data)]
  dNeg  <- dNeg * - 1
  data  <- data.frame(data, dNeg)

  if (randomize == TRUE) data[ ,1] <- sample(data[ ,1])

  if (!is.null(proportion)){
    ctrlRows  <- which(data[ ,1] == 0)
    caseRows  <- which(data[ ,1] == 1)
    # calculate number of case and control to keep
    nCtrlKeep <- round(length(ctrlRows)*proportion, digits = 0)
    nCaseKeep <- round(length(caseRows)*proportion, digits = 0)
    # sample randomly rows of case and control to keep, record rows to keep
    keepRows  <- c(sample(ctrlRows)[1:nCtrlKeep], sample(caseRows)[1:nCaseKeep])
    # subset original data to keep these rows
    data      <- data[keepRows, ]
  }

  ctrl  <- data[data[ ,1] == 0, 2:ncol(data)]
  case  <- data[data[ ,1] == 1, 2:ncol(data)]
  indexNegPos <- rep(0, (nVars*2))
  # end of setting up some initial values ----------------------------#

  # initial loop to establish first optimal marker -------------------#
  allP <- numeric()
  for (i in 1:(nVars*2)){
    caseVar    <- case[ ,i]
    ctrlVar    <- ctrl[ ,i]
    p       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
    allP[i] <- p
  }
  # end of initial loop ----------------------------------------------#

  keepMarkers  <- names(case)[which.min(allP)]
  bestP        <- min(allP, na.rm = TRUE)
  keepIndex    <- which.min(allP)
  # second loop to add another marker --------------------------------#

  allP <- numeric()
  casePrev <- case[ ,keepIndex]
  ctrlPrev <- ctrl[ ,keepIndex]
  for (i in 1:(nVars*2)){
    if (i != keepIndex){
      caseVar <- casePrev + case[ ,i]
      ctrlVar <- ctrlPrev + ctrl[ ,i]
      p       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
    } else {
      p <- 1
    }
    allP[i] <- p
  }
  # end of second loop ----------------------------------------------#

  keepMarkers  <- append(keepMarkers, names(case)[which.min(allP)])
  bestP        <- append(bestP, min(allP, na.rm = TRUE))
  keepIndex    <- append(keepIndex, which.min(allP))

  # check if the latest p is lower than the previous p               #
  continue <- bestP[length(bestP)] < bestP[length(bestP)-1]

  # loop for third through nMarker ----------------------------------#

  while (continue == TRUE){
    allP     <- numeric()
    casePrev <- rowSums(case[ ,keepIndex], na.rm = TRUE)
    ctrlPrev <- rowSums(ctrl[ ,keepIndex], na.rm = TRUE)
    for (i in 1:(nVars*2)){
      if (!(i %in% keepIndex)){
        caseVar <- casePrev + case[ ,i]
        ctrlVar <- ctrlPrev + ctrl[ ,i]
        p       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
      } else {
        p <- 1
      }
      allP[i] <- p
    }
    keepMarkers  <- append(keepMarkers, names(case)[which.min(allP)])
    bestP        <- append(bestP, min(allP, na.rm = TRUE))
    keepIndex    <- append(keepIndex, which.min(allP))
    continue     <- bestP[length(bestP)] < bestP[length(bestP)-1]
    # stop the search when it hits the max number of markers
    if (length(keepMarkers) == nMarkers) continue <- FALSE
  }

  indexNegPos[keepIndex] <- ifelse(keepIndex >= nVars, -1, 1)
  finalIndex   <- ifelse(keepIndex <= nVars, keepIndex, keepIndex - nVars)
  finalMarkers <- data.frame(names(case)[finalIndex], indexNegPos[keepIndex])
  names(finalMarkers) <- c("Marker","Weight")

  ## AUC -------------------------------------------------------------#
  # create function value for each individual
  funcValue <- c(rowSums(case[,c(keepIndex)]), rowSums(ctrl[,c(keepIndex)]))
  # rank individual function values
  ranks       <- rank(funcValue)
  seqCaseCtrl <- c(rep(1, nrow(case)), rep(0, nrow(ctrl)))

  # set up plot -----------------------------------------------------#
  all <- data.frame(funcValue,
                    seqCaseCtrl,
                    ranks)
  all <- all[order(all$ranks),]
  all$refx <- seq(0,1,1/(nrow(all)-1))
  all$refy <- seq(0,1,1/(nrow(all)-1))
  initVal  <- all$seqCaseCtrl[1]
  moveRight <- ifelse(initVal == 0, nrow(case), nrow(ctrl))
  moveUp    <- ifelse(initVal == 0, nrow(ctrl), nrow(case))
  # moveLeft
  for (i in 2:nrow(all)){
    all$x[1] <- 0
    all$y[1] <- 0
    if (all$seqCaseCtrl[i] == initVal){
      all$x[i] = all$x[i-1]
      all$y[i] = all$y[i-1] + 1/(moveUp-1)
    } else {
      all$x[i] = all$x[i-1] + 1/(moveRight)
      all$y[i] = all$y[i-1]
    }
  }

  rocPlot <- ggplot(all, aes(x = x, y = y)) +
    geom_line(size = 1) +
    geom_line(aes(x = refx, y = refy, colour = "red"), size = 1.5) +
    scale_x_continuous(limits = c(0,1)) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("True Positive Rate (Sensitivity)") +
    xlab("False Positive Rate (1 - Specificity)")
  # set up plot -----------------------------------------------------#

  # compute arguments for AUC
  caseFunc  <- sum(ranks[1:nrow(case)]) - nrow(case)*(nrow(case)+1)/2
  ctrlFunc  <- sum(ranks[(nrow(case)+1):length(ranks)]) - nrow(ctrl)*(nrow(ctrl)+1)/2
  # compute AUC
  auc       <- round(max(ctrlFunc, caseFunc)/(caseFunc + ctrlFunc), digits = 4)
  est       <- list(selection  = finalMarkers,
                    auc        = auc,
                    randomize  = randomize,
                    proportion = proportion,
                    rocPlot    = rocPlot)
  class(est) <- "calf"
  return(est)
}
