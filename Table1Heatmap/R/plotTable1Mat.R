#   Copyright 2012, 2013, 2014, Philip C. Schouten
#
#   This file is part of Table1Heatmap.
#
#   Table1Heatmap is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Table1Heatmap is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with Table1Heatmap.  If not, see <http://www.gnu.org/licenses/>.


#' Get the number of patients by crosstable-ing all elements in a list. The output matrix
#' can be used as input for getValueMatPt, maskImageDiag and subsequently plotTable1Matrix

#' @title getTextMatPt
#' @param lstM named list with elements (clinical characteristic) of equal length containing a
#' value for clinical characteristic for every patient. NA's allowed.
#' @return textMat matrix containing crosstable of all clinical characteristics supplied in lstM
.getTextMatPt <- function(lstM) {
  # repeat for every variable
  for (i in 1:length(lstM)) {
    # table the rest of the list versus the variable, using NA's,
    # temporary list
    aggList <- lapply(lstM, function(x) table(x, lstM[[i]], useNA='always'))
    # create a temporary matrix
    aggMat <- do.call(rbind, aggList)

    if (! exists('textMat')) {
      textMat <- aggMat
      } else {
      textMat <- cbind(textMat, aggMat)
      }
    }


  reps <- unlist(lapply(lstM, function(x)  length(levels(x))+1))

  rownames(textMat) <- paste(rep(names(lstM),reps), unlist(lapply(lstM, function(x) c(levels(x),NA))),sep=":")
  colnames(textMat) <- paste(rep(names(lstM),reps), unlist(lapply(lstM, function(x) c(levels(x),NA))),sep=":")

  attr(textMat,'reps') <- reps
  attr(textMat, 'nPt') <- unique(unlist(lapply(lstM, length)))

  return(textMat)
}


#' Get the percentage of patients by crosstabled by getTextMatPt

#' @title getValueMatPt
#' @param textMat textMat obtained by \code{\link{getTextMatPt}}
#' @param nPt number of patients in the analysis, default is extract from the object
#' returned by  \code{\link{getTextMatPt}}
#' @return valueMat matrix containing textMat divided by the total number of patients
#' in the analysis.
.getValueMatPt <- function(textMat, nPt=attr(textMat,'nPt')) {
  # just calculate percentage
  valueMat <- apply(textMat, 2, function(x) x/nPt)

  return(valueMat)
  }


#' Get a matrix of p-values by applying a Fisher Exact (2x2 tables) or Chi-square
#' (> 2x2 tables) to crosstables of the list elements. The output can be used for
#' maskImageDiag and plotTable1Matrix

#' @title getValueMatP
#' @param lstM named list with elements (clinical characteristic) of equal length containing a
#' value for clinical characteristic for every patient. NA's allowed.
#' @return valueMat matrix with p-values obtained by Fisher Exact or Chi square tests
#' on clinical characteristics supplied in lstM
.getValueMatP <- function(lstM) {

    for (i in 1:length(lstM)) {
      aggList <- lapply(lstM, function(x) .checkDimAndTest(table(x, lstM[[i]])))
      aggMat <- matrix(unlist(lapply(aggList, function(x) x$p.value)),ncol=1)

      if (! exists('valueMat')) {
        valueMat <- aggMat
        } else {
        valueMat <- cbind(valueMat, aggMat)
        }
    }
    return(valueMat)

    }


#' Get a matrix of directionality of the crosstables of p-values by applying
#' a Fisher Exact (2x2 tables) to crosstables of the list elements. The output can be used for
#' maskImageDiag and plotTable1Matrix

#' @title getTextMatP
#' @param lstM named list with elements (clinical characteristic) of equal length containing a
#' value for clinical characteristic for every patient. NA's allowed.
#' @return valueMat matrix with directionality of the point estimate of obtained by
#' Fisher Exacton clinical characteristics supplied in lstM. ChiSquare not supported.
.getTextMatP <- function(lstM) {

    for (i in 1:length(lstM)) {
      aggList <- lapply(lstM, function(x) .checkDimAndTest(table(x, lstM[[i]])))

      estMat <- matrix(unlist(lapply(aggList, function(x) ifelse(!is.null(x$estimate), x$estimate,NA))), ncol=1)


      if (! exists('textMat')) {
      textMat <- estMat
      } else {
      textMat <- cbind(textMat, estMat)
      }
    }

    textMat <- textMat > 1

    rownames(textMat) <- names(lstM)
    colnames(textMat) <- names(lstM)



    return(textMat)
    }


#' Check the dimension of a crosstable and then apply a Fisher Exact (2x2 table) or
#' Chi-square test (>2x2 table)

#' @title checkDimAndTest
#' @param crosstab crosstable
#' @return output of euther chisq.test or fisher.test
.checkDimAndTest <- function(crosstab) {
  if (any(!( dim(crosstab) == c(2,2)))) {
    set.seed(1)
    chisq.test(crosstab, simulate.p.value=T)
    } else {
    fisher.test(crosstab)
    }
  }

#' Mask half matrix (lowerright corner) for plotting

#' @title maskImageDiag
#' @param mat matrix for which the right lower corner will be masked
#' @return mat masked matrix with NA's
.maskImageDiag <- function(mat) {
  # create a matrix with nrow(mat) rows filled TTTTTTTTF, TTTTTTTTFF - FFFFFFFFF
  mask <- ! sapply(1:nrow(mat), function(x) c(rep(T, x), rep(F,dim(mat)[1]-x)))

  mat[mask] <- NA

  return(mat)
}

#' Plot a heatmap of table 1, either by plotting p values, n of patients colored by
#' p-values or n of patients colored by percentage of total patients

#' @title plotTable1Heatmap
#' @param factorList named list with clinical variables coded as factors
#' @param method 'AssociationByP' : draw a heatmap with p-values of the association
#' between parameters and a direction of effect if a Fisher Exact test was performed.
#' 'CrosstableByP' : draw a heatmap with p-values of the association between parameters
#' and show the crosstables
#' 'CrosstableByN' : draw a heatmap with percentage of patients in the cell of a crosstable
#'  and show the crosstables
#' @param drawRaster draw horizontal and vertical lines marking characteristics
#' @param ... plotting parameters for image
#' @examples
#' 
#' 
#'  lst <- list(a=sample(c(TRUE,FALSE), 10, replace=TRUE), b=sample(c(TRUE,FALSE), 10, replace=TRUE),
#' c=sample(c(TRUE,FALSE), 10, replace=TRUE))  
#' lst <- lapply(lst, as.factor)
#' 
#' dev.new(height=10, width=10)
#'
#' plotTable1Heatmap(factorList=lst, method='AssociationByP', drawRaster=TRUE)
#' 
#' plotTable1Heatmap(factorList=lst, method='CrosstableByP', drawRaster=TRUE)
#' 
#' plotTable1Heatmap(factorList=lst, method='CrosstableByN', drawRaster=TRUE)
#' 
#' @export
plotTable1Heatmap <- function(factorList, method=c('AssociationByP', 'CrosstableByP', 'CrosstableByN')[1], drawRaster=NULL, ...) {



    textMat <- switch(method,
      AssociationByP =  .getTextMatP(factorList),
      CrosstableByP = .getTextMatPt(factorList),
      CrosstableByN = .getTextMatPt(factorList)
      )
      
    valueMat <- switch(method,
      AssociationByP =  .getValueMatP(factorList),
      CrosstableByP = .getValueMatP(factorList),
      CrosstableByN = .getValueMatPt(textMat)
      )

    par(mar=c(15.1,15.1,2.1,5.1))
    #par(xpd=NA)


    vm <- valueMat
    vt <- .maskImageDiag(textMat)


    if (method %in% c('AssociationByP', 'CrosstableByP' )) {
      if (is.null(drawRaster)) {
        drawRaster <- F
        }



      if (all(dim(vm) == dim(vt))) {
        vt[vm > 0.10] <- NA
        txt <-  c('-','+')[factor(vt, levels=c(F,T))]
        axlab   <- rownames(textMat)
        axtick <- 1:nrow(valueMat)-0.5
        setcex=3

        verticallines <- 1:nrow(vm)
        horizontallines <- 1:nrow(vm)
        }

      if (! all(dim(vm) == dim(vt))) {

         axlab <- rownames(textMat)
         axtick <- 1:nrow(textMat)-0.5
         setcex=1

         rowP <- list()

         for (i in 1:nrow(vm)) {
           rowP[[i]] <- vm[,i]
           }

        rowP <- lapply(rowP, function(x) rep(x, times=attr(vt,'reps')))
        rowP <- lapply(1:length(rowP), function(x) matrix(rep(rowP[[x]], times=attr(vt,'reps')[x]), ncol=attr(vt,'reps')[x]))
        setcex <- 1

        vm <- do.call(cbind, rowP)

        txt <- vt

        verticallines <- cumsum(attr(vt, 'reps'))
        horizontallines <- cumsum(attr(vt, 'reps'))
        }

      colFills <- c('orange',rev(blue2yellow(19)))
      textCols <-  c('black','white')[factor(vm>0.6)]

      }

      if  (method=='CrosstableByN') {
         if (! all(dim(vm) == dim(vt))) {
          stop('Pt invalid option')
          }

      if (is.null(drawRaster)) {
        drawRaster <- T
        }

      txt <- vt
      colFills=c(blue2yellow(19),'orange')
      textCols <- c('black','white')[factor(vm<0.4)]
      setcex=1


      axlab <- rownames(textMat)
      axtick <- 1:nrow(textMat)-0.5


      verticallines <- cumsum(attr(vt, 'reps'))
      horizontallines <- cumsum(attr(vt, 'reps'))

      }


    vm <- .maskImageDiag(vm)

    image(z=vm, x=(0 : nrow(vm)), y=c(0 : ncol(vm)),
      col=colFills, zlim=c(0,1),
      xaxt='n', yaxt='n', xlab='',ylab='', ...)

    if (drawRaster) {
      abline(v=verticallines, h=horizontallines)
    }


    axis(1, at=axtick, labels=axlab, las=2)
    axis(2, at=axtick, labels=axlab, las=2)

    ys <- rep(1:nrow(vm), each=ncol(vm))-0.5
    xs <- rep(1:nrow(vm), times=ncol(vm))-0.5

    #print(setcex)

    text(x=xs, y=ys, labels=txt, cex=setcex,
      col=textCols)

    legend('bottomright', fill=colFills,
      legend=levels(cut(0:1, breaks=c(0, seq(0.05,1,by=0.05)), include.lowest=T)), cex=0.90,
      bg='white')

    box()
}
