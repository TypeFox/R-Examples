################################################################################
#' Recurrence Quantification Analysis (RQA)
#' @description
#' The Recurrence Quantification Analysis (RQA) is an advanced technique for the nonlinear
#' analysis that allows to quantify the number and duration of the recurrences in the 
#' phase space. 
#' @param time.series The original time series from which the phase-space reconstruction is performed.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the \emph{time.series}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param takens Instead of specifying the \emph{time.series}, the \emph{embedding.dim} and the \emph{time.lag}, the user
#' may specify directly the Takens' vectors. 
#' @param radius Maximum distance between two phase-space points to be considered a recurrence.
#' @param lmin Minimal length of a diagonal line to be considered in the RQA. Default \emph{lmin} = 2.
#' @param vmin Minimal length of a vertical line to be considered in the RQA. Default \emph{vmin} = 2.
#' @param save.RM Logical value. If TRUE, the recurrence matrix is stored as a sparse matrix. Note that
#' computing the recurrences in matrix form can be computationally expensive.
#' @param do.plot Logical. If TRUE, the recurrence plot is shown. However, plotting the recurrence matrix is computationally 
#'  expensive. Use with caution.
#' @param ... Additional plotting parameters.
#' @param distanceToBorder In order to avoid border effects, the \emph{distanceToBorder} points near the 
#' border of the recurrence matrix are ignored when computing the RQA parameters. Default, \emph{distanceToBorder} = 2.
#' @return A \emph{rqa}  object that consist of a list with the most important RQA parameters:
#' \itemize{
#'  \item \emph{recurrence.matrix}: A sparse symmetric matrix containing the recurrences of the phase space.
#'  \item \emph{REC}: Recurrence. Percentage of recurrence points in a Recurrence Plot.
#'  \item \emph{DET}: Determinism. Percentage of recurrence points that form diagonal lines.
#'  \item \emph{LAM}: Percentage of recurrent points that form vertical lines.
#'  \item \emph{RATIO}: Ratio between \emph{DET} and \emph{RR}.
#'  \item \emph{Lmax}: Length of the longest diagonal line.
#'  \item \emph{Lmean}: Mean length of the diagonal lines. The main diagonal is not taken into account.
#'  \item \emph{DIV}: Inverse of \emph{Lmax}.
#'  \item \emph{Vmax}: Longest vertical line.
#'  \item \emph{Vmean}: Average length of the vertical lines. This parameter is also referred to as the Trapping time.
#'  \item \emph{ENTR}: Shannon entropy of the diagonal line lengths distribution
#'  \item \emph{TREND}: Trend of the number of recurrent points depending on the distance to the main diagonal
#'  \item \emph{diagonalHistogram}: Histogram of the length of the diagonals.
#'  \item \emph{recurrenceRate}: Number of recurrent points depending on the distance to the main diagonal.
#' }
#' 
#' @references Zbilut, J. P. and C. L. Webber. Recurrence quantification analysis. Wiley Encyclopedia of Biomedical Engineering  (2006).
#' @examples
#' \dontrun{
#' rossler.ts =  rossler(time=seq(0, 10, by = 0.01),do.plot=FALSE)$x
#' rqa.analysis=rqa(time.series = rossler.ts, embedding.dim=2, time.lag=1,
#'                radius=1.2,lmin=2,do.plot=FALSE,distanceToBorder=2)
#' plot(rqa.analysis)
#' }
#' @author Constantino A. Garcia and Gunther Sawitzki
#' @rdname rqa
#' @export rqa
rqa=function(takens = NULL, time.series=NULL, embedding.dim=2, time.lag = 1,
             radius,lmin = 2,vmin = 2,distanceToBorder=2,
             save.RM = TRUE, do.plot=FALSE,...){
  if(is.null(takens)){
    takens = buildTakens( time.series, embedding.dim = embedding.dim, time.lag = time.lag)  
  } 
  ntakens = nrow(takens)
  # distance to the border of the matrix to use in the linear regression that estimates
  #the trend
  maxDistanceMD=ntakens-distanceToBorder
  if (maxDistanceMD <=1) maxDistanceMD=2 # this should not happen
  
  neighs=findAllNeighbours(takens,radius)
  if (save.RM || do.plot){
    neighs.matrix = neighbourList2SparseMatrix(neighs)
  }
  if (do.plot) {
    rec.plot = recurrencePlotFromMatrix(neighs.matrix,...)
  }
  hist=getHistograms(neighs,ntakens,lmin,vmin)
  # calculate the number of recurrence points from the recurrence rate. The recurrence
  # rate counts the number of points at every distance in a concrete side of the main diagonal.
  # Thus, sum all points for all distances, multiply by 2 (count both sides) and add the main
  # diagonal
  numberRecurrencePoints=sum(hist$recurrenceHist)+ntakens
  # calculate the recurrence rate dividing the number of recurrent points at a given
  # distance by all points that could be at that distance
  recurrence_rate_vector=hist$recurrenceHist[1:(ntakens-1)]/((ntakens-1):1)
  #percentage of recurrent points
  REC=(numberRecurrencePoints)/ntakens^2
  diagP=calculateDiagonalParameters(ntakens,numberRecurrencePoints,lmin,hist$diagonalHist,recurrence_rate_vector,maxDistanceMD)  
  #paramenters dealing with vertical lines
  vertP=calculateVerticalParameters(ntakens,numberRecurrencePoints,vmin,hist$verticalHist)
  #join all computations
  rqa.parameters = c(REC=REC,RATIO=diagP$DET/REC,
                     diagP,vertP,
                     list(diagonalHistogram=hist$diagonalHist,
                          recurrenceRate=recurrence_rate_vector))
  
  if (!save.RM){
    neighs.matrix = NULL
  }
  rqa.analysis = c(list(recurrence.matrix=neighs.matrix),
                   rqa.parameters)
  
  rqa.analysis = propagateTakensAttr(rqa.analysis, takens)
  attr(rqa.analysis, "radius") = radius
  class(rqa.analysis) = "rqa"
 
  rqa.analysis
}

#' @export
plot.rqa = function(x,...){
  if (!is.null(x$recurrence.matrix)){
    recurrencePlotFromMatrix(x$recurrence.matrix,
                             ...)
  }else{
    stop("The recurrence matrix has not been stored... Impossible to plot!")
  }
  
}
################################################################################
#' Recurrence Plot 
#' @description
#' Plot the recurrence matrix.
#' @details
#' WARNING: This function is computationally very expensive. Use with caution.
#' @param time.series The original time series from which the phase-space reconstruction is performed.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the \emph{time.series}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param takens Instead of specifying the \emph{time.series}, the \emph{embedding.dim} and the \emph{time.lag}, the user
#' may specify directly the Takens' vectors. 
#' @param radius Maximum distance between two phase-space points to be considered a recurrence.
#' @param ... Additional plotting parameters.
#' @references Zbilut, J. P. and C. L. Webber. Recurrence quantification analysis. Wiley Encyclopedia of Biomedical Engineering  (2006).
#' @author Constantino A. Garcia
#' @export recurrencePlot
#' @import Matrix
#' @useDynLib nonlinearTseries
recurrencePlot=function(takens = NULL, time.series, 
                        embedding.dim, time.lag,radius,
                        ...){
  if(is.null(takens)){
    takens = buildTakens( time.series,
                          embedding.dim = embedding.dim, 
                          time.lag = time.lag)  
  } 
  neighs.matrix = neighbourList2SparseMatrix(findAllNeighbours(takens,radius))
  recurrencePlotFromMatrix(neighs.matrix,...)
}

#private 
recurrencePlotFromMatrix=function(neighs.matrix,
                                  main="Recurrence plot",
                                  xlab="Takens vector's index",
                                  ylab="Takens vector's index",...){
  # need a print because it is a trellis object!!
  rec.plot = image(neighs.matrix,
                   main = main, xlab = xlab, ylab = ylab, 
                   ...)
  print(rec.plot)
  rec.plot
}

neighs2numericType = function(neighs){
  lapply(neighs,
         FUN = function(x){
           if(length(x)==0){
             numeric();
           }else{
             x
           }
         })
}

neighbourList2SparseMatrix = function(neighs){
  ntakens = length(neighs)
  neighs = neighs2numericType(neighs)
  neigh.len = sum(sapply(neighs, FUN=length)) + ntakens
  neighs.matrix = matrix(0,nrow= neigh.len ,ncol=2)
  .Call("nonlinearTseries_neighsList2SparseRCreator",neighs=as.list(neighs),ntakens=as.integer(ntakens),
        neighs_matrix=as.matrix(neighs.matrix),PACKAGE="nonlinearTseries")
  neighs.matrix
  sparseMatrix(neighs.matrix[,1],neighs.matrix[,2],dims = c(ntakens,ntakens),
               symmetric = T)
}

neighbourListToCsparseNeighbourMatrix = function(neighs){
  # sum 1 to columns to include the diagonal (i,i) elements
  neighs.len = sapply(neighs,length)
  max.neighs = 1 + max(neighs.len)
  neighs.matrix= matrix(-1,nrow=length(neighs),
                        ncol = max.neighs)
  neighs = neighs2numericType(neighs)
  .Call("nonlinearTseries_neighsList2Sparse",neighs=as.list(neighs),
        neighs_matrix = as.matrix(neighs.matrix),
        PACKAGE = "nonlinearTseries")
    
  list(neighs = neighs.matrix, nneighs = (neighs.len + 1) )
  
}


calculateVerticalParameters=function(ntakens,numberRecurrencePoints,vmin=2,verticalLinesHistogram){
  #begin parameter computations
  num=sum((vmin:ntakens)*verticalLinesHistogram[vmin:ntakens])
  LAM=num/numberRecurrencePoints
  Vmean=num/sum(verticalLinesHistogram[vmin:ntakens])
  if (is.nan(Vmean)) Vmean=0
  #pick the penultimate
  histogramWithoutZeros=which(verticalLinesHistogram>0)
  if (length(histogramWithoutZeros)>0) Vmax=tail(histogramWithoutZeros,1) else Vmax=0
  
  #results
  params=list(LAM=LAM,Vmax=Vmax,Vmean=Vmean)
  return(params)
}

calculateDiagonalParameters=function(ntakens,numberRecurrencePoints,lmin=2,lDiagonalHistogram,recurrence_rate_vector,maxDistanceMD){
  #begin parameter computations
  num=sum((lmin:ntakens)*lDiagonalHistogram[lmin:ntakens]);
  DET=num/numberRecurrencePoints
  Lmean=num/sum(lDiagonalHistogram[lmin:ntakens])
  aux.index=lmin:(ntakens-1)
  LmeanWithoutMain=(sum((aux.index)*lDiagonalHistogram[aux.index]))/(sum(lDiagonalHistogram[aux.index]))
  #pick the penultimate
  Lmax=tail(which(lDiagonalHistogram>0),2)[1]
  if (Lmax==ntakens) Lmax=0
  DIV=1/Lmax
  pl=lDiagonalHistogram/sum(lDiagonalHistogram)
  diff_0=which(pl>0)
  ENTR=-sum(pl[diff_0]*log(pl[diff_0]));
  
  # use only recurrent points with a distance to the main diagonal < maxDistance
  recurrence_rate_vector=recurrence_rate_vector[1:maxDistanceMD]
  mrrv=mean(recurrence_rate_vector)
  #auxiliar vector for the linear regresion: It is related to the general regression
  #formula xi-mean(x)
  auxiliarVector=(1:maxDistanceMD-(maxDistanceMD+1)/2);auxiliarVector2=auxiliarVector*auxiliarVector
  num=sum(auxiliarVector*((recurrence_rate_vector-mrrv)/2) ) # divide by two because we are having into account just one side of the main diag
  den=sum(auxiliarVector2)
  TREND=num/den
  #results
  params=list(DET=DET,DIV=DIV,Lmax=Lmax,Lmean=Lmean,LmeanWithoutMain=LmeanWithoutMain,ENTR=ENTR,TREND=TREND)
  return(params)
}

getHistograms=function(neighs,ntakens,lmin,vmin){
  
  # the neighbours are labeled from 0 to ntakens-1
  c.matrix = neighbourListToCsparseNeighbourMatrix(neighs)
  verticalHistogram = rep(0,ntakens)
  diagonalHistogram = rep(0,ntakens)
  recurrenceHistogram = rep(0,ntakens)
  # auxiliar variables
  hist = .C("getHistograms", neighs = as.integer(c.matrix$neighs),
            nneighs = as.integer(c.matrix$nneighs), ntakens = as.integer(ntakens), 
            vmin = as.integer(vmin), lmin = as.integer(lmin),
            verticalHistogram = as.integer(verticalHistogram),
            diagonalHistogram = as.integer(diagonalHistogram),
            recurrenceHistogram = as.integer(recurrenceHistogram),
            PACKAGE="nonlinearTseries" )
  
  
  return(list(diagonalHist=hist$diagonalHistogram,recurrenceHist=hist$recurrenceHistogram,
       verticalHist=hist$verticalHistogram))
  
}  





