#' Performs a cluster analysis on the basis of the word-occurrence data.frame.
#'
#' Performs a cluster analysis on the basis of the word-occurrence data.frame
#'   using \code{\link[stats]{hclust}} function.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param myMethod The method to compute distances, see \code{\link[stats]{dist}}
#'   function.
#' @param gp A logical to specify if groups should be made.
#' @param nbGp An intger to specify the number of groups. Ignored if \code{gp=FALSE}.
#' @param getPlot If \code{TRUE}, save the cluster plot in the RESULTS directory.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @param ... Additional arguments from the \code{\link[stats]{hclust}} function.
#' @return An object of class \code{\link[stats]{hclust}}.
#' @examples
#' data("loremIpsum")
#' loremIpsum01 <- loremIpsum[1:100]
#' loremIpsum02 <- loremIpsum[101:200]
#' loremIpsum03 <- loremIpsum[201:300]
#' loremIpsum04 <- loremIpsum[301:400]
#' loremIpsum05 <- loremIpsum[401:500]
#' subDir <- "RESULTS"
#' dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
#' write(x = loremIpsum01, file = "RESULTS/loremIpsum01.txt")
#' write(x = loremIpsum02, file = "RESULTS/loremIpsum02.txt")
#' write(x = loremIpsum03, file = "RESULTS/loremIpsum03.txt")
#' write(x = loremIpsum04, file = "RESULTS/loremIpsum04.txt")
#' write(x = loremIpsum05, file = "RESULTS/loremIpsum05.txt")
#' wordOccuDF <- getwordOccuDF(mywd = paste0(getwd(), "/RESULTS"), excludeSW = FALSE)
#' file.remove(list.files(pattern = "loremIpsum"))
#' doCluster(wordF = wordOccuDF, myMethod = "ward.D2")
#' @export
doCluster <- function(wordF, myMethod = "ward.D2", gp = FALSE, nbGp = 5, getPlot = TRUE, mwidth = 800,
  mheight = 800, formatType = "png", ...){
  ## create RESULTS folder
  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
  ## make cluster analysis
  if(ncol(wordF) > 3){
    fitClust <-  stats::hclust(stats::dist(t(as.matrix(wordF[,2:length(wordF[1,])]))),
      method = myMethod, ...)

    if(getPlot == TRUE){
      if (gp == TRUE){
        ## strat graph
        R.devices::devEval(type = formatType, name = "HCLUST",
          aspectRatio = mheight / mwidth,
          scale = do.call(function(){if((mheight / mwidth) <= 1) {
            x <- max(mheight / 480, mwidth / 480)} else {
              x <- min(mheight / 480, mwidth / 480)}
            return(x)}, list())
          , path = file.path(getwd(), subDir), {
            try(graphics::plot(fitClust, hang = -1), silent = TRUE)
            groups <- try(stats::cutree(fitClust, k = nbGp), silent = TRUE)
            try(stats::rect.hclust(fitClust, k = nbGp, border = "red"), silent = TRUE)
          }
        )
        # end graph
      }else{
        R.devices::devEval(type = formatType, name = "HCLUST",
          aspectRatio = mheight / mwidth,
          scale = do.call(function(){if((mheight / mwidth) <= 1) {
            x <- max(mheight / 480, mwidth / 480)} else {
              x <- min(mheight / 480, mwidth / 480)}
            return(x)}, list())
          , path = file.path(getwd(), subDir), {
          try(graphics::plot(fitClust), silent = TRUE)
          }
        )
      }
    }
    return(fitClust)
  }
}

#' Performs a k-means cluster analysis on the basis of the word-occurrence data.frame.
#'
#' Performs a k-means cluster analysis on the basis of the word-occurrence data.frame
#'   using \code{\link[stats]{kmeans}} function.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param nbClust The number of clusters.
#' @param nbIter The number of iterations allowed.
#' @param algo The algoritm used (see \code{\link[stats]{kmeans}}).
#' @param getPlot If \code{TRUE}, save the k-means cluster plot in the RESULTS directory.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @param ... Additional arguments from the \code{\link[stats]{kmeans}} function.
#' @return An object of class kmeans (see \code{\link[stats]{kmeans}}).
#' @examples
#' data("loremIpsum")
#' loremIpsum01 <- loremIpsum[1:100]
#' loremIpsum02 <- loremIpsum[101:200]
#' loremIpsum03 <- loremIpsum[201:300]
#' loremIpsum04 <- loremIpsum[301:400]
#' loremIpsum05 <- loremIpsum[401:500]
#' subDir <- "RESULTS"
#' dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
#' write(x = loremIpsum01, file = "RESULTS/loremIpsum01.txt")
#' write(x = loremIpsum02, file = "RESULTS/loremIpsum02.txt")
#' write(x = loremIpsum03, file = "RESULTS/loremIpsum03.txt")
#' write(x = loremIpsum04, file = "RESULTS/loremIpsum04.txt")
#' write(x = loremIpsum05, file = "RESULTS/loremIpsum05.txt")
#' wordOccuDF <- getwordOccuDF(mywd = paste0(getwd(), "/RESULTS"), excludeSW = FALSE)
#' file.remove(list.files(pattern = "loremIpsum"))
#' doKmeansClust(wordF = wordOccuDF, nbClust = 2)
#' @export
doKmeansClust <- function(wordF, nbClust = 4, nbIter = 10, algo = "Hartigan-Wong", getPlot = TRUE,
  mwidth = 800, mheight = 800, formatType = "png", ...){
  ## create RESULTS folder
  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
  ## make kmeans-cluster analysis
  if(ncol(wordF)>3){
    dd <-(stats::dist(t(as.matrix(wordF[,2:length(wordF[1,])])), method="euclidian"))# ,colnames=fileNames
    kfit <- stats::kmeans(x=dd, centers=nbClust, iter.max=nbIter, algorithm=algo,...)
    if(getPlot==TRUE){
      R.devices::devEval(type = formatType, name = "KMEANCLUST",
        aspectRatio = mheight / mwidth,
        scale = do.call(function(){if((mheight / mwidth) <= 1) {
          x <- max(mheight / 480, mwidth / 480)} else {
            x <- min(mheight / 480, mwidth / 480)}
          return(x)}, list())
        , path = file.path(getwd(), subDir), {
        try(cluster::clusplot(as.matrix(dd), kfit$cluster, color = TRUE, 
          shade = TRUE, labels = 2, lines = 0), silent = TRUE)
        }
      )
    }
    return(kfit)
  }
}
