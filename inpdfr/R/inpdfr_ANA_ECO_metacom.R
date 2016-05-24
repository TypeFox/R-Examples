#' Performs a metacomunity analysis.
#'
#' Use the package \code{\link[metacom]{Metacommunity}} to analyse the word-occurrence data.frame,
#'    considering words as species and documents as communities.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param numSim Number of simulated null matrices, see \code{\link[metacom]{Metacommunity}}.
#' @param limit An integer to limit the number of words to use in the analysis.
#' @param getPlot If \code{TRUE}, save the plot in the RESULTS directory.
#' @param getTextSink If \code{TRUE}, save the console output in the RESULTS directory.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @return An object of class \code{\link[metacom]{Metacommunity}}.
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
#' wordOccuDF <- getwordOccuDF(mywd = paste0(getwd(), "/RESULTS"), 
#'   excludeSW = FALSE)
#' file.remove(list.files(pattern = "loremIpsum"))
#' doMetacomMetacom(wordF = wordOccuDF)
#' @export
doMetacomMetacom <- function(wordF, numSim = 10, limit = "Inf", getPlot = TRUE, 
  getTextSink = TRUE, mwidth = 800, mheight = 800, formatType = "png"){

  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)

  if(is.numeric(limit) && limit != Inf && limit != "Inf"){
    wordF <- wordF[1:limit,]
  }

  metacomDB <- wordF[,2:length(wordF[1,])]
  metacomDB[metacomDB >= 1] <- 1
  rownames(metacomDB) <- wordF[,1]
  metacomDB <- t(metacomDB)
  metaCom <- metacom::Metacommunity(comm = metacomDB, sims = numSim, allowEmpty = TRUE)

  if(getPlot == TRUE){
    R.devices::devEval(type = formatType, name = "metacom_Metacommunity",
     aspectRatio = mheight / mwidth,
     scale = do.call(function(){if((mheight / mwidth) <= 1) {
       x <- max(mheight / 480, mwidth / 480)} else {
         x <- min(mheight / 480, mwidth / 480)}
       return(x)}, list())
     , path = file.path(getwd(), subDir), {
      try(metacom::Imagine(comm = metacomDB, col = c(0, grDevices::grey(0.6))), silent = TRUE)
     }
    )
  }

  if(getTextSink == TRUE){
    sink('RESULTS/metacom_Metacommunity.txt')
    cat('\n#######################\n### STRUCTURE       ###\n#######################\n')
    try(print(metacom::IdentifyStructure(metaCom)), silent = TRUE)
    cat('\n#######################\n### SUMMARY         ###\n#######################\n')
    try(print(summary(metaCom)), silent = TRUE)
    cat('\n#######################\n### RESULTS         ###\n#######################\n')
    try(print(metaCom), silent = TRUE)
    sink()
  }

  try(print(paste0("Identified community structure: ", metacom::IdentifyStructure(metaCom))), silent = TRUE)

  metaComPkg <- metaCom
  return(metaComPkg)
}
