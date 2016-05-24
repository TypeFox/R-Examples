#
# function to create Splot 'sendplt' object 
#

# The Splot object holds all information needed to create static images



initSplot <- function(mat,
                      plot.calls,
                      Iflag=NA,
                      figTypes=NA,

                      mai.mat=NA,
                      mai.prc=FALSE,
                      plot.extras = NA,
                      
                      source.plot=NA,
                      image.size="800x1100",

                      pointsize=12,
                      res=NA,
                      ps.paper="letter",
                      ps.width=8,
                      ps.height=11,
                      
                      returnVl=TRUE,
                      saveFlag=FALSE,
                      saveName="Splot.RData"
                      
                      ){




  #set up splot object
  Splot = list()
  Splot$mat = mat
  
  # track number of figures in layout
  nfig = max(mat, na.rm=TRUE)
  Splot$nfig = nfig

  #
  # check plot.calls
  #
  # if not of right length fix 
  if( length(plot.calls) != nfig) {
    if(length(plot.calls) >  nfig){
      # if there are more plot calls than specified figures in layout
      # only use plot calls up to the number of figures
      cat(paste("Note:  the length of plot.calls,", length(plot.calls), ", is more than the number of figures,", nfig,".\n     Only using first ", nfig, "plot.calls \n"))
      Splot$plot.calls = plot.calls[1:nfig]      
    }
    if(length(plot.calls) <  nfig){
      # if there are less plot calls than specified figures in layout
      # add empty figures
      cat(paste("Note:  the length of plot.calls,", length(plot.calls), ", is less than the number of figures,", nfig,".\n     Adding empty plots to additional fields \n"))
      plot.calls = c(plot.calls, rep("plot(0,0, pch='', axes=FALSE, xlab='', ylab='')", (nfig-length(plot.calls))))
      Splot$plot.calls = plot.calls     
    }       
  }else{
    Splot$plot.calls = plot.calls
  }
 
  #
  # track which of layout figures are interactive
  #
  # if NA, assume not interactive 
  if( (length(Iflag) == 1) & is.na(Iflag[1]) ){
    cat("Note: Iflag not specified. Assuming figure is not interactive \n")
    Splot$Iflag = FALSE
  }
  # if length of Iflag does not equal the number of figures continue with all as FALSE
  # remember to update this in the makeImap function
  if( length(Iflag) != nfig) {
    cat(paste("Note: length of Iflag,",length(Iflag),", does not equal the number of figures, ", nfig,". \n      Continuing with all Iflags as FALSE \n"))
    Splot$Iflag = rep(FALSE, nfig)
  }
  if( length(Iflag) == nfig) {
    Splot$Iflag = Iflag
  }



  #
  # check figTypes
  #

  # if length of figTypes does not equal the number of figures continue with all as FALSE
  # remember to update this in the makeImap function
  if( length(figTypes) != nfig) {
    cat(paste("Note: length of figTypes,",length(figTypes),", does not equal the number of figures, ", nfig,". \n      Continuing with all figTypes as NA \n"))
    figTypes = rep(NA, nfig)
  } 
  Splot$figTypes = figTypes


  
  #
  # check dimensions of mai.mat if not NA
  #
  if(length(mai.mat) > 1){
    if(!is.null(dim(mai.mat))){
      # columns represent bottom, left, top, right 
      if(dim(mai.mat)[2] != 4){
        cat(paste("Warning:  mai.mat column dimension not correct. \n          mai.mat should be a nx4 matrix, where rows indicate a figure in the layout \n          and columns represent bottom, left, top, right figure margins \n          There should be ",nfig,"rows\n          Continuing with mai.mat as NA \n"))
  
        mai.mat = NA
      }
      # row should be for each figure 
      if(dim(mai.mat)[1] != nfig){
        cat(paste("Warning:  mai.mat row dimension not correct. \n          mai.mat should be a nx4 matrix, where rows indicate a figure in the layout \n          and columns represent bottom, left, top, right figure margins \n          There should be ",nfig,"rows\n          Continuing with mai.mat as NA \n"))
        mai.mat = NA
      }
    }
    # should be a matrix or data.frame 
    if(is.null(dim(mai.mat))){
        cat(paste("Warning:  mai.mat should be a nx4 matrix, where rows indicate a figure in the layout \n          and columns represent bottom, left, top, right figure margins \n          There should be ",nfig,"rows\n          Continuing with mai.mat as NA \n"))
        mai.mat=NA
      }
  }
  Splot$mai.mat = mai.mat
  
  # mai.prc if mai.mat precentage of original size (TRUE) or set values (FALSE)
  Splot$mai.prc = mai.prc

  
  #
  # check length of plt.extras
  #
  if(length(plot.extras)==1){
    # if plt.extras is NA no plotting for all plots
    if(is.na(plot.extras)) plot.extras = as.list(rep(NA, length(plot.calls)))
  }   
  # if there are more plots than plot extras
  # specify no plotting for remaining plots
  if(nfig > length(plot.extras)){
    cat("Note: There are more figures than extra plot calls. \n      Assuming extra plot calls are in order and \n      adding NAs for remaining figures \n" )
    dif = nfig - length(plot.extras)
    for(i in 1:dif){
      plot.extras[(length(plot.extras)+1)] = NA
    }    
  }
  # if there are more plot extras than figures
  # subset  
  if(nfig < length(plot.extras)){
    cat(paste("Note: There are more plot extras than figures \n      plt.extras should be a list of length",nfig," but is of lenght",length(plot.extras),"\n      Using the first",nfig,"plt.extra values \n" ))
    plot.extras = plot.extras[1:nfig]
  }
  Splot$plot.extras = plot.extras


  #
  # if source.plot is NA use appropriate 
  #
  platform = .Platform$OS.type


  # if single entry check for acceptable type
  if(length(source.plot)==1){

    if(!is.na(source.plot) & ((source.plot != "ps") & (source.plot != "png") & (source.plot != "jpeg") & (source.plot != "tiff") & (source.plot != "pdf"))){
      source.plot = NA
    }
    if(!is.na(source.plot) & ((source.plot != "png") & (source.plot != "jpeg"))){
      source.plot = c("png", source.plot)
    }                                     
    if(is.na(source.plot[1]))  source.plot = "png"

    pngF = length(which(source.plot == "png")) > 0
    psF = length(which(source.plot == "ps")) > 0
    tiffF = length(which(source.plot == "tiff")) > 0
    jpegF = length(which(source.plot == "jpeg")) > 0
    pdfF = length(which(source.plot == "pdf")) > 0

    source.plot = NA
    if(pngF) source.plot = c(source.plot, "png")
    if(jpegF) source.plot = c(source.plot, "jpeg")
    if(tiffF) source.plot = c(source.plot, "tiff")
    if(psF) source.plot = c(source.plot, "ps")
    if(pdfF) source.plot = c(source.plot, "pdf")
    
    source.plot = source.plot[-1]
    
  }else{
    # if vector check acceptable types 
    pngF = length(which(source.plot == "png")) > 0
    psF = length(which(source.plot == "ps")) > 0
    tiffF = length(which(source.plot == "tiff")) > 0
    jpegF = length(which(source.plot == "jpeg")) > 0
    pdfF = length(which(source.plot == "pdf")) > 0
  
    if(!pngF & !jpegF) pngF = TRUE
    
    source.plot = NA
    if(pngF) source.plot = c(source.plot, "png")
    if(jpegF) source.plot = c(source.plot, "jpeg")
    if(tiffF) source.plot = c(source.plot, "tiff")
    if(psF) source.plot = c(source.plot, "ps")
    if(pdfF) source.plot = c(source.plot, "pdf")
    source.plot = source.plot[-1] 
    
  }


  
  Splot$platform = platform
  Splot$source.plot = source.plot
  Splot$image.size = image.size
  Splot$pointsize = pointsize
  Splot$res = res
  # if a postscript is used store information for device
  if(source.plot[1] == "ps"){
    Splot$ps.paper = ps.paper
    Splot$ps.width = ps.width
    Splot$ps.height = ps.height
  }

  # set up object to store interactive mappings 
  Splot$iMaps = list()
  for(ifig in 1:Splot$nfig){
    eval.js(paste("Splot$iMaps$Figure", ifig," = list()", sep=""))
  }
  
  # set up object to store info for mappings 
  Splot$iTypes = list()
  for(ifig in 1:Splot$nfig){
    eval.js(paste("Splot$iTypes$Figure", ifig," = list()", sep=""))
  }

  
  
  # specify class 
  class(Splot) <- "Splot"

  # save and return
  if(saveFlag) save(Splot, file=saveName, compress=TRUE)
  if(returnVl) return(Splot)

}# end initSplot 




