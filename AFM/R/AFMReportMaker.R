require("fftwtools")
require("pracma")

require("data.table")

require("gstat")
require(sp)


require("stringr")

require(gridExtra)
require(ggplot2)

require(reshape2)

# for reporting
require(png)
require(grid)


if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c("r", "roughness", "dir.hor", "id","name","press","ang1","ang2","ang3","anis1","anis2","id"))


#' Generate a pdf report for all AFM images in a directory
#' 
#' A function to generate a pdf report for each \code{\link{AFMImage}} in a directory. Images should be in export Nanoscope format as the \code{\link{importFromNanoscope}} function will be used.
#'
#' @param imageDirectory a directory where are located image as Nanoscope export format
#' @param imageNumber (optional) an image number in the directory. If it is set only the selected image will be processed.
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' # A report will be generated for all the images in imageDirectory directory
#' # imageDirectory="c:/images"
#' imageDirectory=tempdir()
#' exit<-generateReportFromNanoscopeImageDirectory(imageDirectory)
#' 
#' # A report will be generated for the fifth image in the imageDirectory directory
#' exit<-generateReportFromNanoscopeImageDirectory(imageDirectory,5)
generateReportFromNanoscopeImageDirectory<-function(imageDirectory, imageNumber) {
  filesToProcess<-list.files(imageDirectory, include.dirs = FALSE, recursive = FALSE,full.names = TRUE,pattern = "\\.txt$")
  
  if (!missing(imageNumber)) {
    if (imageNumber<=length(filesToProcess)){
      filesToProcess<-filesToProcess[imageNumber]
    }else{
      print(paste("Selected image number",imageNumber,paste("exceeds the number of image in directory (",length(filesToProcess),")", sep="")))
      return(FALSE)
    }
  }
  
  for(fullfilename in filesToProcess){
    AFMImage<-importFromNanoscope(fullfilename)
    generateReport(AFMImage)
  }
  return(TRUE)
}

#' Generate an analysis report for one AFMImage
#' 
#' A function to analyse an \code{\link{AFMImage}} and save on disk the analysis. The analysis are saved in outputs directory located in the image directory.
#' All the rdata and image files in the reportDirectory directory are loaded to generate one report for one \code{\link{AFMImage}}.
#'
#' @param AFMImage an \code{\link{AFMImage}} to be analysed
#' @author M.Beauvais
#' @export
#' @examples
#' \dontrun{
#' library(AFM)
#' 
#' # Analyse the AFMImageOfRegularPeaks AFMImage sample from this package
#'   data("AFMImageOfRegularPeaks")
#'   AFMImage<-AFMImageOfRegularPeaks
#' 
#' # exportDirectory="C:/Users/my_windows_login" or exportDirectory="/home/ubuntu"
#'   exportDirectory=tempdir()
#'   AFMImage@@fullfilename<-paste(exportDirectory,"AFMImageOfRegularPeaks.txt",sep="/")
#'   
#' # Start to check if your sample is normaly distributed and isotropic.
#'   generateCheckReport(AFMImage)
#' # If the sample is normaly distributed and isotropic, generate a full report
#'   generateReport(AFMImage)
#' 
#' 
#' # Analyse your own AFM image from nanoscope analysis (TM) software tool
#'    anotherAFMImage<-importFromNanoscope("c:/users/my_windows_login/myimage.txt")
#' # Start to check if your sample is normaly distributed and isotropic.
#'    generateCheckReport(anotherAFMImage)
#' # If your sample is normaly distributed and isotropic, generate a full report
#'    generateReport(anotherAFMImage)
#' }
generateReport <- function(AFMImage) {
  
  sampleName<-basename(AFMImage@fullfilename)
  sampleDirectory<-dirname(AFMImage@fullfilename)
  print(paste("generating a full Report for", sampleName, "in", sampleDirectory))
  
  reportDirectory<-paste(sampleDirectory, "outputs", sep="/")
  createReportDirectory(reportDirectory)
  
  AFMImageAnalyser<-new("AFMImageAnalyser", AFMImage=AFMImage, fullfilename= AFMImage@fullfilename)
  AFMImageAnalyser<-analyse(AFMImageAnalyser=AFMImageAnalyser)
  putAnalysisOnDisk(AFMImageAnalyser=AFMImageAnalyser, AFMImage=AFMImage)
  
  #   # find rdata file for the AFMImage
  #   rdata_directoryfiles<-list.files(reportDirectory, 
  #                                    include.dirs = FALSE, recursive = FALSE, full.names = TRUE,
  #                                    pattern = paste(sampleName,"AFMImageAnalyser.rda$",sep="-"))
  
  #   if (length(rdata_directoryfiles)>0) {
  reportFullfilename<-paste(reportDirectory, paste(sampleName,"fullreport.pdf",sep="-"),sep="/")
  generateAFMImageReport(AFMImageAnalyser, reportFullfilename, isCheckReport = FALSE)
  #   }else{
  #     print("analysis not found...")
  #     print(paste(sampleName,"AFMImageAnalyser.rda",sep="-"))
  #   }
  
  print("done")
}

#' Generate a check report for one AFMImage
#' 
#' Generate a check report in pdf format in order to analyse the distribution and the isotropy of heights of the \code{\link{AFMImage}}.
#'
#' @param AFMImage an \code{\link{AFMImage}} imported from Nanoscope Analysis(TM) with \code{\link{importFromNanoscope}} or created manually \code{\link{AFMImage}}
#' @author M.Beauvais
#' @export
#' @examples
#' \dontrun{
#' library(AFM)
#' 
#' # Analyse the AFMImageOfRegularPeaks AFMImage sample from this package
#'   data("AFMImageOfRegularPeaks")
#'   AFMImage<-AFMImageOfRegularPeaks
#' # exportDirectory="C:/Users/my_windows_login" or exportDirectory="/home/ubuntu"
#'   exportDirectory=tempdir()
#'   AFMImage@@fullfilename<-paste(exportDirectory,"AFMImageOfRegularPeaks.txt",sep="/")
#'   
#' # Start to check if your sample is normaly distributed and isotropic.
#'   generateCheckReport(AFMImage)

#' # If the sample is normaly distributed and isotropic, generate a full report
#'   generateReport(AFMImage)
#'
#' # Analyse your own AFM image from nanoscope analysis (TM) software tool
#'    anotherAFMImage<-importFromNanoscope("c:/users/me/myimage.txt")
#' # Start to check if your sample is normaly distributed and isotropic.
#'    generateCheckReport(anotherAFMImage)
#' # If your sample is normaly distributed and isotropic, generate a full report
#'    generateReport(anotherAFMImage)
#' }
generateCheckReport <- function(AFMImage) {
  sampleName<-basename(AFMImage@fullfilename)
  sampleDirectory<-dirname(AFMImage@fullfilename)
  print(paste("Generating a check report for", sampleName, "in", sampleDirectory))
  
  reportDirectory<-paste(sampleDirectory, "outputs", sep="/")
  createReportDirectory(reportDirectory)
  
  AFMImageAnalyser<-checkIsotropy(AFMImage)
  putAnalysisOnDisk(AFMImageAnalyser, AFMImage)
  
  #   sampleName<-basename(AFMImage@fullfilename)
  #   rdata_directoryfiles<-list.files(reportDirectory, 
  #                                    include.dirs = FALSE, recursive = FALSE, full.names = TRUE,
  #                                    pattern = paste(sampleName,"AFMImageAnalyser.rda$",sep="-"))
  
  # if (length(rdata_directoryfiles)>0) {
  
  reportFullfilename<-paste(reportDirectory, paste(sampleName,"checkreport.pdf",sep="-"),sep="/")
  
  generateAFMImageReport(AFMImageAnalyser, reportFullfilename, isCheckReport = TRUE)
  #   }else{
  #     print("analysis not found...")
  #     print(paste(sampleName,"AFMImageAnalyser.rda",sep="-"))
  #   }
  print("done")
}

#' @title Generate an analysis report from an AFMImageAnalyser object
#' 
#' @description \code{generateAFMImageReport} generates a report from an AFMImageAnalyser object
#'
#' @param AFMImageAnalyser an \code{\link{AFMImageAnalyser}} to be used to produce report
#' @param reportFullfilename location on disk where to save the generated report
#' @param isCheckReport TRUE to generate a check report must be generated, FALSE to generate a full report
#' @author M.Beauvais
#' @export
generateAFMImageReport<-function(AFMImageAnalyser, reportFullfilename, isCheckReport){
  numberOfModelsPerPage=3
  
  if (missing(isCheckReport)) {
    isCheckReport = TRUE
  }
  
  AFMImage<-AFMImageAnalyser@AFMImage
  
  #   # load AFMImageAnalyser in rda format
  #   print(paste("loading", basename(oneAFMImageAnalyserFile)))
  #   
  #   x = load(file = oneAFMImageAnalyserFile)
  #   AFMImageAnalyser= get(x)
  #   rm(x)
  #   
  fullfilename <- AFMImageAnalyser@AFMImage@fullfilename
  sampleName<-basename(AFMImageAnalyser@AFMImage@fullfilename)
  reportDirectory=dirname(AFMImageAnalyser@fullfilename)
  createReportDirectory(reportDirectory)
  #   
  #   # load image in rda format
  #   afmImageFullfilename<-paste(dirname(oneAFMImageAnalyserFile) ,paste(sampleName, "AFMImage.rda", sep="-"),sep="/")
  #   print(paste("loading", basename(afmImageFullfilename)))
  #   x = load(file = afmImageFullfilename)
  #   AFMImage= get(x)
  #   rm(x)
  #   
  #   
  #   
  #   print(paste("processing image", sampleName))
  #   
  
  
  
  # save all images necessary for the report on disk
  putImagesFromAnalysisOnDisk(AFMImageAnalyser, AFMImage, reportDirectory)

  print(paste("creating", basename(reportFullfilename), "..."))
  pdf(reportFullfilename, width=8.27, height=11.69)
  
  # first page
  rglImagefullfilename<-get3DImageFullfilename(reportDirectory, sampleName)
  print(paste("reading", basename(rglImagefullfilename), "..."))
  img <- readPNG(rglImagefullfilename)
  
  roughnesses<-getRoughnessParameters(AFMImageAnalyser@AFMImage)  
  basicImageInfo<-data.table(name=c("Scan size",
                                    "Samples per line",
                                    "Lines",
                                    "Total Rrms",
                                    "Ra (mean roughness)"),
                             values=c(paste(AFMImageAnalyser@AFMImage@scansize,"nm"),
                                      paste(as.character(AFMImageAnalyser@AFMImage@samplesperline)),
                                      paste(as.character(AFMImageAnalyser@AFMImage@lines)),
                                      paste(round(roughnesses$totalRMSRoughness_TotalRrms, digits=4),"nm"),
                                      paste(round(roughnesses$MeanRoughness_Ra, digits=4),"nm")))
  imageInformationDTPlot<-getGgplotFromDataTable(basicImageInfo,
                                                 removeRowNames= TRUE,
                                                 removeColNames=TRUE)
  
  grid.newpage() # Open a new page on grid device
  pushViewport(viewport(layout = grid.layout(5, 4)))
  
  vp1<-viewport(layout.pos.row = 2:3, layout.pos.col = 1:4)
  grid.raster(img,vp=vp1)
  
  vp0<-viewport(layout.pos.row = 1, layout.pos.col = 2:3)
  grid.text(sampleName,  vp=vp0, gp=gpar(fontsize=20, col="black"))
  
  vp2<-viewport(layout.pos.row = 4:5, layout.pos.col = 1:4)
  print(imageInformationDTPlot,vp=vp2)
  
  
  # page for checking
  # normality / omni direction of samples
  if (!length(AFMImageAnalyser@variogramAnalysis@directionalVariograms)==0) {
    
    exportpng2FullFilename<-getDirectionalVarioPngFullfilename(reportDirectory, sampleName)
    print(paste("reading",basename(exportpng2FullFilename)))
    directionalVariograms<-readPNG(exportpng2FullFilename)
    
    grid.newpage() # Open a new page on grid device
    pushViewport(viewport(layout = grid.layout(4, 4)))
    
    qq <- checkNormalityQQ(AFMImage)
    m <- checkNormalityDensity(AFMImage)
    
    vp2<-  viewport(layout.pos.row = 1:2, layout.pos.col = 1:2)
    print(qq, vp = vp2)
    
    vp3<-viewport(layout.pos.row = 1:2, layout.pos.col = 3:4)
    print(m, vp = vp3) 
    
    vp4<-viewport(layout.pos.row = 3:4, layout.pos.col = 1:4)
    grid.raster(directionalVariograms,vp=vp4)
  }
  
  if (!isCheckReport) {
    # get variogram model evaluation
    if (!length(AFMImageAnalyser@variogramAnalysis@variogramModels)==0) {
      
      mergedDT<-getDTModelEvaluation(AFMImageAnalyser@variogramAnalysis)
      print(mergedDT)
      
      sillrangeDT<-getDTModelSillRange(AFMImageAnalyser@variogramAnalysis)
      setkey(sillrangeDT, "model")
      name<-press<-NULL
      sampleDT <- mergedDT[name==basename(AFMImageAnalyser@AFMImage@fullfilename)]
      setkey(sampleDT, "model")
      #sampleDT <- sampleDT[cor>0.98]
      sampleDT<-merge(sampleDT, sillrangeDT, by="model")
      sampleDT<-sampleDT[,name:=NULL]
      sampleDT <- unique(sampleDT)
      sampleDT <- sampleDT[order(-rank(cor), rank(press))]
      print(basename(AFMImageAnalyser@AFMImage@fullfilename))
      print(basename(AFMImageAnalyser@fullfilename))
      print(sampleDT)
      
      
      summarySampleDT<-copy(sampleDT)
      summarySampleDT$press<-round(sampleDT$press)
      summarySampleDT$sill<-round(sampleDT$sill)
      summarySampleDT$range<-round(sampleDT$range)
      
      print("plotting variogram table...")
      existsVariogramModel<-TRUE
      if (nrow(sampleDT)!=0) {
        plotBestVariogramModelsTable<-getGgplotFromDataTable(summarySampleDT, 
                                                             removeRowNames=TRUE, 
                                                             removeColNames=FALSE)
      }else{
        print("no good variogram table...")  
        existsVariogramModel<-FALSE
        sampleDT <- mergedDT[name==basename(fullfilename)]
        sampleDT <- unique(sampleDT)
        sampleDT <- sampleDT[order(-rank(cor), rank(press))]
        plotBestVariogramModelsTable<-getGgplotFromDataTable(summarySampleDT, 
                                                             removeRowNames=TRUE, 
                                                             removeColNames=FALSE)
      }
      #print(plotBestVariogramModelsTable)
    }

    # best variogram models page
    if (!length(AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram)==0) {
      # chosen sample
      grid.newpage() # Open a new page on grid device
      pushViewport(viewport(layout = grid.layout(7, 2)))
      
      sampleSpplotfullfilename<-getSpplotImagefullfilename(reportDirectory, sampleName)
      print(paste("reading",basename(sampleSpplotfullfilename)))
      sampleImg <- readPNG(sampleSpplotfullfilename)
      
      sampleSpplotfullfilename<-getVarioPngchosenFitSample(reportDirectory, sampleName)
      print(paste("reading",basename(sampleSpplotfullfilename)))
      chosenFitSampleImg <- readPNG(sampleSpplotfullfilename)
      
      
      vp0<-  viewport(layout.pos.row = 1, layout.pos.col = 1:2)
      grid.text("Variogram analysis",  vp=vp0, gp=gpar(fontsize=20, col="black"))
      
      vp1<-  viewport(layout.pos.row = 2:3, layout.pos.col = 1)
      grid.raster(sampleImg,vp=vp1)
      #vp3<-viewport(layout.pos.row = 9, layout.pos.col = 1)
      #grid.text("Original",  vp=vp3, gp=gpar(fontsize=10, col="black"))
      
      vp2<-  viewport(layout.pos.row = 2:3, layout.pos.col = 2)
      grid.raster(chosenFitSampleImg,vp=vp2)
      #vp4<-viewport(layout.pos.row = 9, layout.pos.col = 2)
      #grid.text("Sample",  vp=vp4, gp=gpar(fontsize=10, col="black"))    

      vp5<-viewport(layout.pos.row = 4:7, layout.pos.col = 1:2)
      print(plotBestVariogramModelsTable,vp=vp5)
      
      printVariogramModelEvaluations(AFMImageAnalyser, sampleDT, numberOfModelsPerPage)
    }
    
    
    
    # Roughness against length scale
    if (!length(AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)==0) {
      
      grid.newpage() # Open a new page on grid device
      pushViewport(viewport(layout = grid.layout(7, 2)))
      
      vp0<-viewport(layout.pos.row = 1, layout.pos.col = 1:2)
      grid.text("Roughness vs. lengthscale",  vp=vp0, gp=gpar(fontsize=20, col="black"))
      
      exportCsvFullFilename<-getRoughnessAgainstLengthscale(reportDirectory, sampleName)
      
      print(paste("reading",basename(exportCsvFullFilename)))
      samplePredictedImg <- readPNG(exportCsvFullFilename)
      vp1<-viewport(layout.pos.row = 2:4, layout.pos.col = 1)
      grid.raster(samplePredictedImg,vp=vp1)
      
      exportCsvFullFilename<-getRoughnessAgainstLengthscale10nm(reportDirectory, sampleName)
      print(paste("reading",basename(exportCsvFullFilename)))
      samplePredictedImg <- readPNG(exportCsvFullFilename)
      vp1<-viewport(layout.pos.row = 2:4, layout.pos.col = 2)
      grid.raster(samplePredictedImg,vp=vp1)
      
      
      for(i in c(0,1)) {
        exportpng2FullFilename=getRoughnessAgainstLengthscaleIntersection(reportDirectory, paste(sampleName,i*2,sep="-"))
        if (file.exists(exportpng2FullFilename)) {
          print("intersection inserted...")  
          
          img<-readPNG(exportpng2FullFilename)
          vp2<-viewport(layout.pos.row = 5:7, layout.pos.col = i+1)
          
          grid.raster(img,vp=vp2)
        }
      }
    }
    
    # export fractal dimension
    if (!length(AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)==0) {
      grid.newpage() # Open a new page on grid device
      pushViewport(viewport(layout = grid.layout(7, 4)))
      
      vp0<-viewport(layout.pos.row = 1, layout.pos.col = 1:4)
      grid.text("Fractal dimension analysis",  vp=vp0, gp=gpar(fontsize=20, col="black"))
      
      n=length(AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)
      print(n)
      if (n>0) {
        sampleDT <- data.table( 
          fd_method= c(sapply(1:n, function(i) AFMImageAnalyser@fdAnalysis@fractalDimensionMethods[[i]]@fd_method)),
          fd= c(sapply(1:n, function(i) AFMImageAnalyser@fdAnalysis@fractalDimensionMethods[[i]]@fd)),
          fd_scale= c(sapply(1:n, function(i) AFMImageAnalyser@fdAnalysis@fractalDimensionMethods[[i]]@fd_scale)))
        
        #       sampleDT <- data.table( AFMImageAnalyser@fdAnalysis@fractalDimensionMethods)
        #       setnames(sampleDT,c("method","scale"),c("fd_method","fd_scale"))
        
        print(sampleDT)
        #sampleDT <- sampleDT[,c(2,13,14,15), with = FALSE]
        setkey(sampleDT, "fd_method")
        sampleDT <- unique(sampleDT)
        
        name<-NULL
        plotFractalDimensionTable<-getGgplotFromDataTable(sampleDT[, name:=NULL])
        vp3<-viewport(layout.pos.row = 2:3, layout.pos.col = 1:4)
        print(plotFractalDimensionTable,vp=vp3)
        
        exportpng2FullFilename=getFractalDimensionsPngFullfilename(reportDirectory, sampleName, "isotropic")
        if (file.exists(exportpng2FullFilename)) {
          img<-readPNG(exportpng2FullFilename)
          vp4<-viewport(layout.pos.row = 4:5, layout.pos.col = 1:2)
          grid.raster(img,vp=vp4)
        }
        exportpng2FullFilename=getFractalDimensionsPngFullfilename(reportDirectory, sampleName, "squareincr")
        if (file.exists(exportpng2FullFilename)) {
          img<-readPNG(exportpng2FullFilename)
          vp5<-viewport(layout.pos.row = 4:5, layout.pos.col = 3:4)
          grid.raster(img,vp=vp5)
        }
        exportpng2FullFilename=getFractalDimensionsPngFullfilename(reportDirectory, sampleName, "filter1")
        if (file.exists(exportpng2FullFilename)) {
          img<-readPNG(exportpng2FullFilename)
          vp6<-viewport(layout.pos.row = 6:7, layout.pos.col = 1:2)
          grid.raster(img,vp=vp6)
        }
      }
    }
  }
  dev.off()
  rm(AFMImageAnalyser)
  
}


#' @title printVariogramModelEvaluations
#' 
#' @description \code{printVariogramModelEvaluations} generates a graphic element containing the evaluation of all variogram models
#'
#' @param AFMImageAnalyser an \code{\link{AFMImageAnalyser}} to be used to produce report
#' @param numberOfModelsPerPage numeric to specify the number of model evaluations per pages
#' @param sampleDT a data.table containg the evaluation information
#' @author M.Beauvais
#' @export
printVariogramModelEvaluations<-function(AFMImageAnalyser, sampleDT, numberOfModelsPerPage){
  
  #####################
  # new page for experimental variogram and models
  experimentalVariogramDT<-AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram
  experimentalVariogramDT$name<-"Experimental"
  #drops <- c("dir.hor","dir.ver","id","np")
  experimentalVariogramDT<-experimentalVariogramDT[,c("dir.hor","dir.ver","id","np"):=NULL]
  #names(experimentalVariogramDT)
  
  sampleName<-basename(AFMImageAnalyser@AFMImage@fullfilename)
  
  sampleSpplotfullfilename<-getSpplotImagefullfilename(tempdir(), sampleName)
  saveSpplotFromAFMImage(AFMImageAnalyser@AFMImage, sampleSpplotfullfilename, withoutLegend=TRUE)
  #print(paste("reading",basename(sampleSpplotfullfilename)))
  #sampleImg<-getSpplotFromAFMImage(AFMImage=AFMImageAnalyser@AFMImage, expectedWidth=80, expectHeight=60, withoutLegend=TRUE)
  sampleImg <- readPNG(sampleSpplotfullfilename)
  
  
  
  allVarioModels<-str_sub(sampleDT$model,-3)
  for (i in seq(1:length(allVarioModels))) {
    indexInPage<-i%%numberOfModelsPerPage
    
    if (indexInPage==1) {
      # Open a new page
      grid.newpage() 
      pushViewport(viewport(layout = grid.layout(numberOfModelsPerPage*2, 3)))
    }
    if (indexInPage==0)indexInPage=numberOfModelsPerPage
    
    #print(indexInPage)
    
    #plot experimental variogram and model variogram
    vp1<-viewport(layout.pos.row = (indexInPage-1)*2+1, layout.pos.col = 2)
    grid.raster(sampleImg,vp=vp1)
    
    #     pushViewport(vp1)
    #     print(sampleImg,newpage=FALSE)
    #     popViewport(1)
    
    
    
    #print(i)
    
    allVariogramModelEvaluation<-AFMImageAnalyser@variogramAnalysis@variogramModels
    for (j in seq(1:length(allVariogramModelEvaluation))) {
      if (allVariogramModelEvaluation[j][[1]]@fit.v[2]$model==allVarioModels[i]) break;
    }
    #print(j)
    #print(allVariogramModelEvaluation[j][[1]]@fit.v[2]$model)
    #predictedfullfilename<-getSpplotPredictedImageFullfilename(reportDirectory, sampleName, allVarioModels[i])
    
    
    modelName<-allVariogramModelEvaluation[j][[1]]@fit.v[2]$model
    part_valid_pr<-AFMImageAnalyser@variogramAnalysis@variogramModels[[i]]@mykrige
    cuts<-AFMImageAnalyser@variogramAnalysis@cuts
    withoutLegend<-TRUE
    colLimit<-length(cuts)+3
    cols <- getSpplotColors(colLimit) 
    
    predictedspplotfullfilename<-getSpplotPredictedImageFullfilename(tempdir(), sampleName, modelName)
    saveSpplotFromKrige(predictedspplotfullfilename, modelName, part_valid_pr,cuts, withoutLegend = TRUE)  
    # TODO save on disk as png and read
    
    # read image on disk
    #print(paste("reading", basename(predictedspplotfullfilename)))
    samplePredictedImg <- readPNG(predictedspplotfullfilename)
    #samplePredictedImg<-spplot(part_valid_pr["var1.pred"], cuts=cuts, col.regions=cols,key=list(lines=FALSE, col="transparent"))
    
    vp2<-viewport(layout.pos.row = (indexInPage-1)*2+1, layout.pos.col = 3)
    grid.raster(samplePredictedImg,vp=vp2)
    #     pushViewport(vp2)
    #     print(samplePredictedImg,newpage=FALSE)
    #     popViewport(1)
    
    ang1<-ang2<-ang3<-anis1<-anis2<-name<-NULL
    fit.v<-allVariogramModelEvaluation[j][[1]]@fit.v
    vgm1<-vgm(fit.v[2]$psill, fit.v[2]$model, fit.v[2]$range, kappa = fit.v[2]$kappa, anis = c(fit.v[2]$anis1, fit.v[2]$anis2), add.to = vgm(fit.v[1]$psill, fit.v[1]$model, fit.v[1]$range,kappa = fit.v[1]$kappa,anis = c(fit.v[1]$anis1, fit.v[1]$anis2)))
    newModelDT<-data.table(vgm1)
    setnames(newModelDT, "psill", "sill" )
    newModelDT<-rbind(newModelDT, sampleDT[i], fill=TRUE)
    newModelDT<- newModelDT[, ang1:=NULL]
    newModelDT<- newModelDT[, ang2:=NULL]
    newModelDT<- newModelDT[, ang3:=NULL]
    newModelDT<- newModelDT[, anis1:=NULL]
    newModelDT<- newModelDT[, anis2:=NULL]
    
    
    plotVariogramModelTable<-getGgplotFromDataTable(newModelDT[,name:=NULL])
    vp4<-viewport(layout.pos.row = (indexInPage-1)*2+1+1, layout.pos.col = 2:3)
    print(vp=vp4, plotVariogramModelTable, row.names= FALSE, include.rownames=FALSE)
    
    # variogram from model
    myvgm<-experimentalVariogramDT
    experimentalVariogramDTnrow=nrow(myvgm)
    class(myvgm) = c("gstatVariogram", "data.frame")
    myvgm$np=rep(1,experimentalVariogramDTnrow)
    myvgm$dir.hor=rep(0,experimentalVariogramDTnrow)
    myvgm$dir.ver=rep(0,experimentalVariogramDTnrow)
    myvgm$id=rep(factor("var1"),experimentalVariogramDTnrow)
    
    begin<-(indexInPage-1)*2+1
    vp3<-viewport(layout.pos.row = begin:(begin+1), layout.pos.col = 1, width=200, height=200)
    vgLine <- rbind(
      cbind(variogramLine(vgm1, maxdist = max(myvgm$dist)), id = "Raw")
    )
    p1<-ggplot(myvgm, aes(x = dist, y = gamma, colour = id)) +  geom_line(data = vgLine) + geom_point()   
    p1 <- p1 + ylab("semivariance")
    p1 <- p1 + xlab("distance (nm)")
    p1 <- p1 + ggtitle("Semivariogram")
    p1 <- p1 + guides(colour=FALSE)
    p1 <- p1 + expand_limits(y = 0)
    print(p1,vp=vp3)
    
  }
  
}

getGgplotFromDataTable<-function(DT, removeRowNames, removeColNames) {
  if (missing(removeRowNames)) removeRowNames<-TRUE
  if (missing(removeColNames)) removeColNames<-FALSE
  
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.8)),
    colhead = list(fg_params=list(cex = 0.9)),
    rowhead = list(fg_params=list(cex = 0.9)))
  
  qplotFromDataTable<- qplot(1:10, 1:10, geom = "blank") + 
    theme_bw() +
    theme(line = element_blank(), text = element_blank())
  if ((removeRowNames)&&(removeColNames)) {
    qplotFromDataTable<- qplotFromDataTable + annotation_custom(grob = tableGrob(DT,  theme = mytheme, rows = NULL, cols=NULL))
  }else{
    if (removeRowNames) {
      qplotFromDataTable<- qplotFromDataTable + annotation_custom(grob = tableGrob(DT,  theme = mytheme, rows = NULL))
    }else{
      if (removeColNames) {
        qplotFromDataTable<- qplotFromDataTable + annotation_custom(grob = tableGrob(DT,  theme = mytheme, cols = NULL))
      }else{
        qplotFromDataTable<- qplotFromDataTable + annotation_custom(grob = tableGrob(DT,  theme = mytheme))
      }
    }
  }
  
  return(qplotFromDataTable)
}

#' Put the images from all analysis on disk
#' 
#' A function to put on disk all the images from variogram, PSD Analysis of an \code{\link{AFMImage}}
#' An AFM Image 3D representation is saved on disk thanks to the \code{\link{rgl}} package.
#' On Unix system, it is necessary to have a X server connection to be able to use the \code{\link{rgl}} package.
#'
#' @param AFMImageAnalyser an \code{\link{AFMImageAnalyser}}
#' @param AFMImage an \code{\link{AFMImage}}
#' @param exportDirectory where the images will be stored
#' @export
#' @author M.Beauvais
putImagesFromAnalysisOnDisk<-function(AFMImageAnalyser, AFMImage, exportDirectory) {
  exportFractalDimImagesForReport(AFMImage, exportDirectory)
  exportPSDImagesForReport(AFMImageAnalyser, AFMImage, exportDirectory)          
  exportVariogramImagesForReport(AFMImageAnalyser, AFMImage, exportDirectory)
  export3DImageForReport(AFMImage, exportDirectory)
}

exportVariogramImagesForReport<- function(AFMImageAnalyser, AFMImage, exportDirectory) {
  
  class(AFMImageAnalyser)="AFMImageAnalyser"
  class(AFMImageAnalyser@variogramAnalysis)="AFMImageVariogramAnalysis"
  
  sampleName<-basename(AFMImage@fullfilename)
  
  # ssplot of real sample for comparison with predicted sample from each variogram model
  spplotImagefullfilename<-getSpplotImagefullfilename(exportDirectory, sampleName)
  saveSpplotFromAFMImage(AFMImage, spplotImagefullfilename, withoutLegend=TRUE)
  
  # directional variograms files
  if (!length(AFMImageAnalyser@variogramAnalysis@directionalVariograms)==0) {
    
    exportCsvFullFilename<-getDirectionalVarioCsvFullfilename(exportDirectory, sampleName)
    print(paste("saving", basename(exportCsvFullFilename)))
    tryCatch({
      write.table(AFMImageAnalyser@variogramAnalysis@directionalVariograms, exportCsvFullFilename, sep=",") 
    }, error = function(e){
      print("error",e)
    }) 
    
    dvarios<-AFMImageAnalyser@variogramAnalysis@directionalVariograms
    dist<-gamma<-dir.hor<-NULL
    p1 <- ggplot(dvarios, aes(x=dist, y=gamma, color=as.factor(dir.hor) , shape=as.factor(dir.hor)))
    p1 <- p1 + geom_point()
    p1 <- p1 + ylab("semivariance")
    p1 <- p1 + xlab("distance (nm)")
    p1 <- p1 + ggtitle("Semivariogram")
    p1 <- p1 + expand_limits(y = 0)
    p1 <- p1 + guides(colour=FALSE)
    #print(p1)
    exportpng2FullFilename<-getDirectionalVarioPngFullfilename(exportDirectory, sampleName)
    print(paste("saving", basename(exportpng2FullFilename)))
    png(filename=exportpng2FullFilename, units = "px", width=800, height=800)
    print(p1)
    dev.off()
  }
  # omnidirectional variogram files 
  if (!length(AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram)==0) {
    
    exportCsvFullFilename<-getOmnidirectionalVarioCsvFullfilename(exportDirectory, sampleName)
    print(paste("saving", basename(exportCsvFullFilename)))
    AFMImageVariogram<-AFMImageAnalyser@variogramAnalysis@omnidirectionalVariogram
    class(AFMImageVariogram)=c("gstatVariogram","data.frame")
    tryCatch({
      write.table(AFMImageVariogram, exportCsvFullFilename, sep=",") 
    }, error = function(e){
      print("error",e)
    }) 
    
    
    myvgm<-AFMImageVariogram
    experimentalVariogramDTnrow=nrow(myvgm)
    class(myvgm) = c("gstatVariogram", "data.frame")
    myvgm$np=rep(1,experimentalVariogramDTnrow)
    myvgm$dir.hor=rep(0,experimentalVariogramDTnrow)
    myvgm$dir.ver=rep(0,experimentalVariogramDTnrow)
    myvgm$id=rep(factor("var1"),experimentalVariogramDTnrow)
    
    dist<-gamma<-id<-NULL
    p1<-ggplot(myvgm, aes(x = dist, y = gamma, colour = id)) + geom_point()   
    p1 <- p1 + ylab("semivariance")
    p1 <- p1 + xlab("distance (nm)")
    p1 <- p1 + ggtitle("Semivariogram")
    p1 <- p1 + expand_limits(y = 0)
    p1 <- p1 + guides(colour=FALSE)
    
    exportpng2FullFilename<-getOmnidirectionalVarioPngFullfilename(exportDirectory, sampleName)
    print(paste("saving", basename(exportpng2FullFilename)))
    png(filename=exportpng2FullFilename, units = "px", width=800, height=800)
    print(p1)
    dev.off()
    
    # chosen sample plot
    TheData<-as.data.frame(AFMImage@data)
    TheData=na.omit(TheData)
    part_model <- TheData[AFMImageAnalyser@variogramAnalysis@chosenFitSample, ]
    coordinates(part_model) = ~x+y
    proj4string(part_model)=CRS("+init")
    is.projected(part_model)
    
    pchosenFitSample<-spplot(part_model, col.regions="black",contour=TRUE,key=list(lines=FALSE, col="transparent"))
    expectedWidth = 400
    expectHeight = 300
    
    exportpngFullFilename<-getVarioPngchosenFitSample(exportDirectory, sampleName)
    print(paste("saving", basename(exportpngFullFilename)))
    png(filename=exportpngFullFilename, units = "px", width=expectedWidth, height=expectHeight)
    print(pchosenFitSample)
    dev.off()
    
    # save images from variogram modeling
    totalVariogramModels=length(AFMImageAnalyser@variogramAnalysis@variogramModels)
    #print(totalVariogramModels)
    fullfilename<-AFMImage@fullfilename
    cuts<-AFMImageAnalyser@variogramAnalysis@cuts
    for (i in seq(1,totalVariogramModels)) {
      #print(AFMImageAnalyser@variogramAnalysis@variogramModels[[i]]@res)
      testedModel<-AFMImageAnalyser@variogramAnalysis@variogramModels[[i]]@model
      print(testedModel)
      if (testedModel=="Wav2") { 
        vgm<-vgm( 5, "Exp", 1, add.to = vgm(5, "Wav", 1, nugget = 2.5))
      }else{
        vgm<-vgm(5,testedModel,1,0)
      }
      mykrige<-AFMImageAnalyser@variogramAnalysis@variogramModels[[i]]@mykrige
      
      predictedspplotfullfilename<-getSpplotPredictedImageFullfilename(exportDirectory, sampleName, testedModel)
      saveSpplotFromKrige(predictedspplotfullfilename, vgm,mykrige,cuts, withoutLegend = TRUE)  
      predictedAFMImage<-getAFMImageFromKrige(AFMImage, vgm, mykrige)
      class(predictedAFMImage) = c("AFMImage")
      #displayIn3D(predictedAFMImage,1024, full2Dfilename)
      export3DImageForReport(predictedAFMImage, exportDirectory)
      
    }
  }
  
}

exportPSDImagesForReport<-function(AFMImageAnalyser, AFMImage, exportDirectory) {
  #class(AFMImageAnalyser)="AFMImageAnalyser"
  #class(AFMImageAnalyser@psdAnalysis)="AFMImagePSDAnalysis"
  filename<-basename(AFMImage@fullfilename)
  
  # export Roughness against lengthscale graph
  if (!length(AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale)==0) {
    
    data<-AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale
    r<-roughness<-NULL
    p1 <- ggplot(data, aes(x=r, y=roughness, colour= basename(filename)))
    p1 <- p1 + geom_point()
    p1 <- p1 + geom_line()
    p1 <- p1 + ylab("roughness (nm)")
    p1 <- p1 + xlab("lengthscale (nm)")
    p1 <- p1 + guides(colour=FALSE)  
    pngFilename=paste(filename,"roughness-against-lengthscale.png",sep="-")
    exportpngFullFilename<-paste(exportDirectory, pngFilename, sep="/")
    print(paste("saving", basename(exportpngFullFilename)))
    png(filename=exportpngFullFilename, units = "px", width=800, height=800)
    print(p1)
    dev.off()
    
    # focus on the first 10nm
    newdata<-data[r<10,]
    r<-roughness<-NULL
    p1 <- ggplot(newdata, aes(x=r, y=roughness, colour= basename(filename)))
    p1 <- p1 + geom_point()
    p1 <- p1 + geom_line()
    p1 <- p1 + ylab("roughness (nm)")
    p1 <- p1 + xlab("lengthscale (nm)")
    p1 <- p1 + guides(colour=FALSE)
    pngFilename=paste(filename,"roughness-against-lengthscale-10nm.png",sep="-")
    exportpngFullFilename<-paste(exportDirectory, pngFilename, sep="/")
    print(paste("saving", basename(exportpngFullFilename)))
    png(filename=exportpngFullFilename, units = "px", width=800, height=800)
    print(p1)
    dev.off()
    
    # save intersections images
    if (!length(AFMImageAnalyser@psdAnalysis@intersections)==0) {
      saveOnDiskIntersectionForRoughnessAgainstLengthscale(AFMImageAnalyser, exportDirectory)
    }
  }
}

export3DImageForReport<-function(AFMImage, exportDirectory) {
  sampleName<-basename(AFMImage@fullfilename)
  rglImagefullfilename<-get3DImageFullfilename(exportDirectory, sampleName)
  if (displayIn3D(AFMImage, 1024, rglImagefullfilename)) {
    rgl.close()
  }
}

createReportDirectory<-function(reportDirectory) {
  if (!file.exists(reportDirectory)){
    print(paste("creating report directory",reportDirectory))
    dir.create(file.path(reportDirectory), showWarnings = FALSE)
  }
  if (!isReportDirectoryWritePermissionCorrect(reportDirectory)) {
    stop(paste("Error: can't write to output directory", reportDirectory))
  }
  print(paste("report directory is", reportDirectory))  
}

isReportDirectoryWritePermissionCorrect<-function(reportDirectory) {
  tryCatch({
    fullfilename=paste(reportDirectory, "permCheck.txt", sep="/")
    fileConn<-file(fullfilename)
    writeLines(c("Hello","World"), fileConn)
    close(fileConn)
    file.remove(fullfilename)
    return(TRUE)  
  }, error = function(e){
    close(fileConn)
    return(FALSE)
  }) 
}
