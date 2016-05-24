useRtopWithIntamap <- function() {
  if (!"intamap" %in% loadedNamespaces()) stop("Please load intamap with library(intamap) before rerunning this function")
  
    packageStartupMessage("Loading optional package: intamap \n")
    info = matrix(c("estimateParameters","spatialPredict","methodParameters",
             rep("rtop",3),rep(NA,3)),ncol = 3)
#    info[3,2] = "rtopVariogramModel"
    registerS3methods(info,package = "intamap",env = environment(rtopVariogram))
}

