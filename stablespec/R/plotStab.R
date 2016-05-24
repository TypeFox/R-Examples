#' Plot each of the stability of causal path and edge including the threshold
#' of stability and model complexity.
#' @title Plot of edge and causal path stability.
#' @param listOfFronts \code{\link{list}} of models including their fitness
#' and subset index, which is one of the outputs of \code{\link{stableSpec}}.
#' @param threshold threshold of stability selection. The default is 0.6.
#' @param stableCausal \code{\link{list}} of the stability of causal path
#' for the whole range of model complexities, which is one of the outputs
#' of \code{\link{stableSpec}}.
#' @param stableCausal_l1 \code{\link{list}} of the stability of
#' causal path of length 1
#' for the whole range of model complexities, which is one of the outputs
#' of \code{\link{stableSpec}}.
#' @param stableEdge \code{\link{list}} of the stability of edge
#' for the whole range of model complexities, which is one of the outputs
#' of \code{\link{stableSpec}}.
#' @param longitudinal \code{TRUE} for longitudinal data,
#' and \code{FALSE} cross-sectional data.
#' @return Plot of causal path and edge stability for
#' every pair of variables, and a plot of altogether combined
#' including the corresponding stability and complexity thresholds.
#' @examples
#' \donttest{
#' the_data <- adhd
#' numSubset <- 1
#' num_iteration <- 5
#' num_pop <- 10
#' mut_rate <- 0.075
#' cross_rate <- 0.85
#' longi <- FALSE
#' num_time <- 1
#' the_co <- "covariance"
#' # assummed that nothing causing variable Gender
#' cons_matrix <- matrix(c(2, 1, 3, 1, 4, 1, 5, 1, 6, 1), 5, 2, byrow=TRUE)
#' th <- 0.1
#' to_plot <- FALSE
#'
#' result_adhd <- stableSpec(theData=the_data, nSubset=numSubset,
#' iteration=num_iteration,
#' nPop=num_pop, mutRate=mut_rate, crossRate=cross_rate,
#' longitudinal=longi, numTime=num_time,
#' co=the_co, consMatrix=cons_matrix, threshold=th, toPlot=to_plot)
#'
#' plotStability(listOfFronts=result_adhd$listOfFronts, threshold=th,
#' stableCausal=result_adhd$causalStab,
#' stableCausal_l1=result_adhd$causalStab_l1,
#' stableEdge=result_adhd$edgeStab,
#' longitudinal=longi)
#' }
#' @author Ridho Rahmadi \email{r.rahmadi@cs.ru.nl}
#' @export
plotStability <- function(listOfFronts = NULL, threshold = NULL,
                          stableCausal = NULL, stableCausal_l1 = NULL,
                          stableEdge = NULL, longitudinal = NULL) {


  if (!is.null(listOfFronts)) {
    if (!is.list(listOfFronts)) {
      stop("Argument listOfFronts should be a list.")
    }
  } else {
    stop("Argument listOfFronts cannot be missing.")
  }

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || is.matrix(threshold)) {
      stop("Argument threshold should be positive numeric, e.g., 0.6.")
    }
  } else {
    threshold <- 0.6
  }

  if (!is.null(stableCausal)) {
    if (!is.list(stableCausal)) {
      stop("Argument stableCausal should be a list.")
    }
  } else {
    stop("Argument stableCausal cannot be missing.")
  }

  if (!is.null(stableCausal_l1)) {
    if (!is.list(stableCausal_l1)) {
      stop("Argument stableCausal_l1 should be a list.")
    }
  } else {
    stop("Argument stableCausal_l1 cannot be missing.")
  }

  if (!is.null(stableEdge)) {
    if (!is.list(stableEdge)) {
      stop("Argument stableEdge should be a list.")
    }
  } else {
    stop("Argument stableEdge cannot be missing.")
  }

  if (!is.null(longitudinal)) {
    if (!is.logical(longitudinal)) {
      stop("Argument longitudinal should be either logical TRUE or FALSE.")
    }
  } else {
    stop("Argument longitudinal cannot be missing.")
  }



  # get the number of variables
  numVar <- nrow(stableCausal[[1]])

  # get the total number of model complexities
  if (longitudinal) {
    numComp <- ((numVar / 2) * (numVar / 2)) +
      ((numVar / 2) * (numVar / 2 - 1) / 2)

    stringSize <- (numVar * numVar + (numVar * (numVar - 1)))

  } else {
    numComp <- numVar * (numVar-1) / 2
    stringSize <- (numVar*(numVar-1))
  }

  #get at which model complexity the minimum average of BIC is
  minBicAt <- getMinBic(listOfFronts, stringSize)

  # create a matrix of M by N, where M is
  # the number of model complexities
  # and N is the number of edges
  numCol <- numVar * numVar - numVar
  numColEdge <- numVar * (numVar - 1) / 2
  stableEdgeB <- stableEdge

  mat4PlotEdgeB <-  matrix(0, length(stableEdge), numCol)


  mat4PlotCausal <-  mat4PlotCausal_l1 <- matrix(0, length(stableCausal),
                                                 numCol)

  mat4PlotEdge <-  matrix(0, length(stableEdge), numCol)


  #get the title, later will be changed by the rownames of the data
  theTitle <- NULL
  for(i in 1:numVar) {
    for(j in 1:numVar) {
      if(i != j) {
        theTitle <- c(theTitle, paste('Var',j,'-->Var', i, sep=""))
      }
    }
  }

  for(i in 1:length(stableCausal)) {

    #make NA for diag and remove NA
    diag(stableCausal[[i]]) <-
      diag(stableCausal_l1[[i]]) <- diag(stableEdge[[i]]) <-
      stableEdgeB[[i]][upper.tri(stableEdgeB[[i]], diag=TRUE)] <- NA

    mat4PlotCausal[i, ] <- as.vector(stableCausal[[i]]
                                    [!is.na(stableCausal[[i]])])

    mat4PlotCausal_l1[i, ] <- as.vector(stableCausal_l1[[i]]
                                       [!is.na(stableCausal_l1[[i]])])

    mat4PlotEdge[i,] <- as.vector(stableEdge[[i]][!is.na(stableEdge[[i]])])

    #these are for plot stable edges all in one
    stableEdgeB[[i]][upper.tri(stableEdgeB[[i]], diag=TRUE)] <- NA

    mat4PlotEdgeB[i, ] <- as.vector(stableEdgeB[[i]][!is.na(stableEdgeB[[i]])])
  }

  #plot each causal path and edge in 2 by 4 tables
  par( mfrow = c(2, 4) )

  for(j in 1:numCol) {
    #initialize a plot
    plot(NULL, NULL, xlab = "Model Complexity", ylab="Selection Probability",
         xlim=c(length(stableCausal), 1), ylim=c(0, 1), xaxt = "n")

    #set up the axis
    axis(1, at=(numComp + 1):1, labels=0:numComp)

    #causal stability, color is red
    lines(mat4PlotCausal[, j], title(theTitle[j]),
          xlim=rev(range(1:length(stableCausal))), col=2, lwd=3, lty=3)
    #causal stability length 1, color is green
    lines(mat4PlotCausal_l1[, j],title(theTitle[j]),
          xlim=rev(range(1:length(stableCausal_l1))), col=3, lwd=3, lty=3)
    #edge stability, color is blue
    lines(mat4PlotEdge[, j],title(theTitle[j]),
          xlim=rev(range(1:length(stableEdge))), col=4, lwd=3, lty=3)
    #bic threshold
    lines(rep(minBicAt, 7), seq(-0.1, 1.1, 0.2))
    #stability threshold
    lines(c(numComp + 5, numComp / 2, numComp / 3, -1),
          rep(threshold, 4))
  }

  #plot all in one edge stability
  par( mfrow = c(1, 1))

  plot(NULL,NULL, main="Edge Stability", xlab = "Model Complexity",
       ylab="Selection Probability", xlim=c(length(stableEdge), 1),
       ylim=c(0, 1), xaxt = "n")

  axis(1, at=(numComp + 1):1, labels=0:numComp)

  for(k in 1:numColEdge) {
    lines(mat4PlotEdgeB[, k], xlim=rev(range(1:length(stableEdge))),
          col=k, lwd=3, lty=3)
  }
  #bic threshold
  lines(rep(minBicAt, 7), seq(-0.1, 1.1, 0.2))
  #stability threshold
  lines(c(numComp + 5, numComp / 2, numComp / 3, -1),
        rep(threshold, 4))

  #plot all in one causal path stability
  par(mfrow = c(1, 1))

  plot(NULL,NULL, main="Causal Path Stability", xlab = "Model Complexity",
       ylab="Selection Probability", xlim=c(length(stableCausal), 1),
       ylim=c(0, 1), xaxt = "n")

  axis(1, at=(numComp + 1):1, labels=0:numComp)

  for(k in 1:numCol) {
    lines(mat4PlotCausal[,k], xlim=rev(range(1:length(stableCausal))),
          col=k, lwd=3, lty=3)
  }
  #bic threshold
  lines(rep(minBicAt, 7), seq(-0.1, 1.1, 0.2))
  #stability threshold
  lines(c(numComp + 5, numComp / 2, numComp / 3, -1),
        rep(threshold, 4))
}
