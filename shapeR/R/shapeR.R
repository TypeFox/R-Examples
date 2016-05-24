
#' @import pixmap
#' @import gplots
#' @import wavethresh
#' @import jpeg
#' @import methods
#' @import vegan
#' @import MASS


#' @title shapeR
#' @description Collection and analysis of otolith shape data
#' @name shapeR
#' @docType package
#' @author Lisa Anne Libungan & Snaebjorn Palsson
NULL

#' @title A constructor for the shapeR class
#' @description a shapeR class
#' @slot project.path Path to the project where the images are stored
#' @slot info.file Info file containing fish and otolith information
#' @slot master.list.org The contents of the \code{info.file}
#' @slot master.list The contents of the \code{info.file} with added shape parameters and descriptors
#' @slot outline.list.org A list of all the original otolith outlines
#' @slot outline.list A list of all the otolith outlines. It returns a list of smoothed if contour smoothing (using\code{smoothout}) has been conducted.
#' @slot filter A logical vector selecting the otoliths used for analysis
#' @slot wavelet.coef.raw The wavelet coefficients for all the otolith outlines
#' @slot wavelet.coef The wavelet coefficients after aligning with the \code{info.file}. The data is generated when \link{enrich.master.list} is run
#' @slot wavelet.coef.std The standardized wavelet coefficients. The data is generated when \link{stdCoefs} is run
#' @slot wavelet.coef.std.removed The index of the removed wavelet coefficients after standardization. The data is generated when \link{stdCoefs} is run
#' @slot fourier.coef.raw The Fourier coefficients for all the otolith outlines
#' @slot fourier.coef The Fourier coefficients for after aligning with the info file. The data is generated when \link{enrich.master.list} is run
#' @slot fourier.coef.std The standardized Fourier coefficients. The data is generated when \link{stdCoefs} is run
#' @slot fourier.coef.std.removed The index of the removed Fourier coefficents after standardization. The data is generated when \link{stdCoefs} is run
#' @slot shape.coef.raw The uncalibrated shape measurements for all the otoliths. The shape parameters are: otolith.area, otolith.length, otolith.width, otolith.perimeter
#' @slot shape.coef The shape measurements for after aligning with the info file. The shape parameters have been calibrated using the calibration parameter as registered in the datafile as the column 'cal'.
#' @slot shape.std The standardized shape measurements. The data is generated when \link{stdCoefs} is run
#' @slot shape.std.removed The index of the removed shape measurements after standardization. The data is generated when \link{stdCoefs} is run
#' 
#' @return a \code{\linkS4class{shapeR}} object
#' @references Libungan LA and Palsson S (2015) ShapeR: An R Package to Study Otolith Shape Variation among Fish Populations. PLoS ONE 10(3): e0121102. \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121102}
#' @seealso \url{https://github.com/lisalibungan/shapeR}
#' @seealso \code{\link{plotWavelet}}
#' @seealso \code{\link{plotFourier}}
#' @seealso \code{\link{plotWaveletShape}}
#' @seealso \code{\link{plotFourierShape}}
#' @seealso \code{\link[vegan]{capscale}}
#' @seealso \code{\link{cluster.plot}}
#' @seealso \code{\link{setFilter}}
#' @seealso \code{\link[MASS]{lda}}
#' @seealso \code{\link{detect.outline}}
#' @seealso \code{\link{generateShapeCoefficients}}
#' @seealso \code{\link{enrich.master.list}}
#' @examples
#'  
#'\dontrun{
#'
#'# This example has two sections: (1) Demonstration of how a shapeR object 
#'# is analyzed and (2) How to create a shapeR object from an archive of 
#'# image files. 
#'
#'#-----------------------------------------
#'# Section 1: Analyzing a shapeR object
#'
#'data(shape)
#'
#'#Standardize coefficients
#'shape = stdCoefs(shape,"pop","length_cm")
#'
#'#Visualize Wavelet and Fourier coefficients
#'plotWavelet(shape,level=5,class.name= "pop",useStdcoef=TRUE)
#'plotFourier(shape,class.name= "pop",useStdcoef=TRUE)
#'
#'#Examine the mean shapes
#'plotWaveletShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1) 
#'plotFourierShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1)
#' 
#'#Canonical analysis
#'library(vegan)
#'cap.res = capscale(getStdWavelet(shape) ~ getMasterlist(shape)$pop)
#'anova(cap.res)
#'
#'#Visualize the canonical scores
#'eig=eigenvals(cap.res,constrained=TRUE)
#'eig.ratio = eig/sum(eig)
#'
#'cluster.plot(scores(cap.res)$sites[,1:2],getMasterlist(shape)$pop
#',plotCI=TRUE
#',xlab=paste("CAP1 (",round(eig.ratio[1]*100,1),"%)",sep="")
#',ylab=paste("CAP2 (",round(eig.ratio[2]*100,1),"%)",sep="")
#',main="Canonical clustering"
#')
#' 
#'#Only analyze Icelandic and Norwegian samples
#'shape = setFilter(shape, getMasterlist(shape, useFilter = FALSE)$pop %in% c("NO","IC"))
#'
#'#Classifier on standardized wavelet
#'lda.res.w = lda(getStdWavelet(shape),getMasterlist(shape)$pop,CV=TRUE)
#'ct.w = table(getMasterlist(shape)$pop,lda.res.w$class)
#'
#'diag(prop.table(ct.w, 1))
#'
#'# Total percent correct
#'sum(diag(prop.table(ct.w)))
#'
#'cap.res = capscale(getStdWavelet(shape) ~ getMasterlist(shape)$pop)
#'anova(cap.res)
#'
#'#Classifier on canoncial values
#'lda.res.w = lda(scores(cap.res)$sites,getMasterlist(shape)$pop,CV=TRUE)
#'ct.w = table(getMasterlist(shape)$pop,lda.res.w$class)
#'
#'diag(prop.table(ct.w, 1))
#'
#'# Total percent correct
#'sum(diag(prop.table(ct.w)))
#'
#'#-----------------------------------------
#'# Section 2: Creating a shapeR object from image files
#'
#'# The following example requires the user to download an archive of JPEG 
#'# files from https://github.com/lisalibungan/shapeR/
#'# place the ShapeAnalysis directory inside the working directory.
#'
#'shape = shapeR("~/ShapeAnalysis/","FISH.csv")
#'shape = detect.outline(shape,write.outline.w.org = TRUE)
#'
#'shape = generateShapeCoefficients(shape)
#'shape = enrich.master.list(shape)
#'}
#' @exportClass shapeR
#' @rdname shapeR
setClass("shapeR", 
         representation(project.path="character",
                        info.file="character",outline.list.org="list",outline.list="list",filter="vector",
                        fourier.coef="matrix",wavelet.coef="matrix",
                        fourier.coef.std="matrix",wavelet.coef.std="matrix",
                        fourier.coef.std.removed="vector",wavelet.coef.std.removed="vector",
                        fourier.coef.raw="matrix",#fourier.coef.raw.info="matrix",
                        wavelet.coef.raw="matrix",
                        shape.coef.raw="matrix",  
                        shape.coef="matrix", shape.std="matrix",shape.std.removed="vector",
                        master.list.org = "data.frame",
                        master.list="data.frame")
         ,prototype(outline.list.org=list(list()),outline.list=list(list()))
         )

#' @export shapeR
#' @param project.path The base project path where the images are stored
#' @param info.file The information file which store the information on the fish and otoliths. This is the base for the master.list
#' @param ... Additional parameters to be passed to 'read.csv' for reading the info.file
#' @rdname shapeR

shapeR <- function(project.path,info.file,...)
{
  shape = new("shapeR",project.path=project.path,info.file=info.file)
  shape@master.list.org = read.csv(paste(project.path,info.file,sep="/"),header=T,...)
  return(shape)
}



#' @exportMethod show
#' @title Show a shapeR object
#' @param object a shapeR oject
#' @description Show the project.path and info.file, the number of outlines that have been read and which fundamental methods have been run.
#' @rdname show
setMethod("show", "shapeR", 
          function(object){
            cat(class(object), "instance with", sum(sapply(object@outline.list,length)), 
                "outlines in project path:\n",object@project.path,"\n")
            if( length(object@shape.coef.raw) == 0 )
              cat("generateShapeCoefficients has not been run\n")
            else
              cat("generateShapeCoefficients has been run\n")
            
            if(length(object@master.list) == 0)
              cat("enrich.master.list has not been run\n")
            else
              cat("enrich.master.list has been run\n")
              
          }
)

setValidity("shapeR",
            function(object) {
              #Check if the project.path is valid
              if(!file.exists(object@project.path))
                FALSE
              #Check if the info file exists
              if(!file.exists(object@info.file))
                FALSE
              
              TRUE
            }
)


#' @title Get wavelet/Fourier coefficients and basic shape variables
#' @description Generates shape variables based on Fourier/wavelet reconstruction. Wavelet coefficients for wavelet. Basic shape parameters are also collected (area, length, width, perimeter).
#' @usage generateShapeCoefficients(object,...)
#' @param object \code{\linkS4class{shapeR}} object 
#' @param ... Additional parameters to be passed to the \code{\link{wd}} function of the \code{\link{wavethresh}} package for the wavelet decomposition of the otolith outlines 
#' @return A \code{\linkS4class{shapeR}} object with values in slots:
#' \itemize{
#'   \item wavelet.coef.raw
#'   \item fourier.coef.raw
#'   \item shape.coef.raw
#' }
#' @examples
#'\dontrun{
#'data(shape)
#'shape = generateShapeCoefficients(shape)}
#' @rdname generateShapeCoefficients
#' @export generateShapeCoefficients
#' @author Lisa Anne Libungan & Snaebjorn Palsson
#' @seealso \code{\link[wavethresh]{wavethresh}}
#' @references Nason, G. (2012). \code{\link{wavethresh}}: Wavelets statistics and transforms. R package, version 4.5.
#' @references Claude, J. (2008). Morphometrics with R. Springer. 316 p.
generateShapeCoefficients <- function(object,...) {
  f.coef = matrix(NA,nrow=0,ncol= 32*4) # Maximum number of harmonics
  #f.coef.info = matrix(NA,nrow=0,ncol= 5)
  f.coef.names = c()
  #Set column names for 32 Fourier harmonics
  colnames(f.coef) = c(paste("FD",1:32,'a',sep=""),paste("FD",1:32,'b',sep=""),paste("FD",1:32,'c',sep=""),paste("FD",1:32,'d',sep=""))
            
  w.coef = matrix(NA,nrow=0,ncol= sum(2^(0:10)) +1 )
  w.coef.names = c()
  #Set column names for 10 wavelet levels
  colnames(w.coef) = c("mean.radii","Ws0c1", paste("Ws1c",1:2^1,sep=""), paste("Ws2c",1:2^2,sep=""),
                       paste("Ws3c",1:2^3,sep=""),paste("Ws4c",1:2^4,sep=""),paste("Ws5c",1:2^5,sep=""),
                       paste("Ws6c",1:2^6,sep=""),paste("Ws7c",1:2^7,sep=""), paste("Ws8c",1:2^8,sep=""),
                       paste("Ws9c",1:2^9,sep=""),paste("Ws10c",1:2^10,sep="") 
                       )
            
            
  stocks = names(object@outline.list)
            
  sp.coef = matrix(NA,nrow=0,ncol=4)
  colnames(sp.coef) = c("otolith.area","otolith.length","otolith.width","otolith.perimeter")
            
  total=0
  for( dname in stocks)
    total=total+length(object@outline.list[[dname]])
              
  .shapeR.check.outline.list(object) #Stops execution if no outlines
  
  pb <- txtProgressBar(min = 0, max = total, initial = 0, char = "=",style=3,
                       width = NA, title="generateShapeCoefficients")          
  ptotal=0
  tryCatch( {
    for( dname in stocks){
      num.fname = length(object@outline.list[[dname]])
      i=0
                
      for(fname in names(object@outline.list[[dname]])){
        i=i+1
        #Outline. X and Y coordinates
        utlina = object@outline.list[[dname]][[fname]]
                  
        #Normalized elliptic Fourier
        f.coef.names = c(f.coef.names, paste(dname,fname,sep=";") ) 
                  
        #Fourier coeffs
        N=.shapeR.NEF(cbind(utlina$X,utlina$Y),n=32) # n is a parameter in the NEF function. Gives information about how many harmonics is needed. Possible to set as 32, 64, or skip it.
        values = matrix(c(N$A,N$B,N$C,N$D),nrow=1)
        colnames(values) = colnames(f.coef)
        f.coef = .shapeR.rbind.fill.matrix(f.coef,values)
        f.coef[is.na(f.coef)] = 0
        
        #f.coef.info = .shapeR.rbind.fill.matrix(f.coef.info,matrix(c(N$size,N$theta,N$psi,N$ao,N$co),nrow=1))
                  
        #Wavelet
        w.coef.names = c(w.coef.names, paste(dname,fname,sep=";") ) 
                  
        #Outline. X and Y coordinates
        utlina = object@outline.list[[dname]][[fname]]
                  
        #Wavelet coeffs
        #First normalize the otolith for area then find centroid
        mynd = utlina 
        mynd.fl = .shapeR.otolith.image.parameters(cbind(mynd$X,mynd$Y))$area
        mynd = .shapeR.rotate.svd(mynd)
        mynd = .shapeR.centroid(list(X=mynd[,1],Y=mynd[,2]))
        mynd$X = mynd$X/sqrt(mynd.fl)
        mynd$Y = mynd$Y/sqrt(mynd.fl)
                  
        #Polar coordinates
        ppolar= .shapeR.regularradius(mynd$X,mynd$Y,1024)
                  
        # wd transforms the polar coordinates to wavelet object with the package wavethresh (object is a container with various information)
        # wd calculates also the wavelet coefficients
                  
        wwd = wd(ppolar$radii,...)
                  
        # accessD is a function in wavethresh that returns the wavelet coefficients from the wavelet object 
        # The wavelet coeffs increase by the power of two for each level that we go
        # Level 1: 2^1, level 2: 2^2, level 3: 2^3, level 4: 2^4, level 5: 2^5 = 32 coeffs. We use coeffs to level 5 in our 
        # calculation because it is a sufficient number to reconstruct the otolith outline. Level 6 would be a overshoot in 
        # our calculations 
                  
        w.coef = rbind(w.coef, c(mean(ppolar$radii), accessD(wwd, level=0), accessD(wwd, level=1),
                                 accessD(wwd, level=2),
                                 accessD(wwd, level=3),
                                 accessD(wwd, level=4),
                                 accessD(wwd, level=5),
                                 accessD(wwd, level=6),
                                 accessD(wwd, level=7),
                                 accessD(wwd, level=8),
                                 accessD(wwd, level=9)
                                 )
                       )
        
        otolith.sp = .shapeR.otolith.image.parameters(cbind(utlina$X,utlina$Y))
        sp.coef = rbind(sp.coef,c(otolith.sp$area,otolith.sp$width,otolith.sp$height,otolith.sp$perimeter))
                  
                  
        ptotal=ptotal+1
        setTxtProgressBar(pb, ptotal, label=paste(dname,fname,sep="/"))        
      }
    }
    rownames(f.coef) = f.coef.names
    rownames(w.coef) = w.coef.names
    rownames(sp.coef) = f.coef.names
    #names(f.coef.info) = c('size','theta','psi','ao','co')
                
    object@fourier.coef.raw = f.coef
    #object@fourier.coef.raw.info = f.coef.info
    object@wavelet.coef.raw = w.coef
    object@shape.coef.raw = sp.coef
    return(object)
    },error = function(e) 
    {message(paste("Error. directory:",dname,", filename:",fname,sep=" "));
     print(e)
     return(object)
    },
    finally={
      close(pb)}
    )

}


#' @title Detect otolith outline
#' @description Determine the outline of otolith images in jpeg format which have been stored in the \code{Fixed} folder.
#' @usage detect.outline(object, threshold=0.2, mouse.click=FALSE, 
#'            display.images=FALSE, write.outline.w.org=FALSE)
#' @param object \code{\linkS4class{shapeR}} object
#' @param threshold Grayscale threshold. Value between 0 and 1.
#' @param mouse.click If TRUE, the user clicks where the starting point for the otolith contour extraction algorithm should start. Default is the center of the image. Could be good to set as TRUE if the otolith detection produces an error.
#' @param display.images If TRUE, each image is displayed and the user can visualize how the outline is captured
#' @param write.outline.w.org If TRUE, the outline is written on top of the original image using the function \code{\link{write.image.with.outline}}, and can be seen in the \code{Original_with_outline} folder
#' @return A \code{\linkS4class{shapeR}} object with otolith outlines in the slot outline.list
#' @examples
#'\dontrun{
#'#Use test data from Libungan and Palsson (2015):
#'shape = shapeR("ShapeAnalysis/","FISH.csv")
#'shape = detect.outline(shape, threshold=0.2,write.outline.w.org = TRUE)}
#' @rdname detect.outline
#' @export detect.outline
#' @details Based on the Conte function (Claude 2008)
#' @author Lisa Anne Libungan & Snaebjorn Palsson 
#' @references Claude, J. (2008). Morphometrics with R. Springer. 316 p.
#' @references Urbanek, S. (2014). \code{\link{jpeg}}: Read and write JPEG images. R package  version 0.1-8.
#' @references Bivand, R., Leisch, F. & Maechler, M. (2011) \code{\link{pixmap}}: Bitmap Images (''Pixel Maps''). R package version 0.4-11.
#' @references Libungan LA and Palsson S (2015) ShapeR: An R Package to Study Otolith Shape Variation among Fish Populations. PLoS ONE 10(3): e0121102. \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121102}
detect.outline <- 
  function(object,threshold=0.2,mouse.click=FALSE,display.images=FALSE,write.outline.w.org = FALSE)
{
    fix.image.path = paste(object@project.path,"Fixed",sep="/")
    w.dir = getwd()
            
    setwd(fix.image.path)
    dnames = list.dirs('.',full.names=TRUE)
    dnames = gsub("[.|/]","", dnames[-1]) 
            
    total=0
    total.images=0
    for (dname in dnames)
    {
      fnames = list.files(dname,pattern="\\.(jpg|jpeg|png)",ignore.case=T)
      for (fname in fnames)
      {
        strip.fname = gsub(".(jpg|jpeg|png)","",fname,ignore.case=T) 
        total.images = total.images+1                
        if((!(dname %in% names(object@outline.list)) || !(strip.fname %in% names(object@outline.list[[dname]])) ))
        {
          total=total+1 
        }
      }
    }
    
    if(total==0)
    {
      message(paste("All otolith",total.images,"outlines have been detected."))
      return(object)
    }
    
    pb <- txtProgressBar(min = 0, max = total, initial = 0, char = "=",style=3,
                         width = NA, title="detect.outline")
            
    ptotal=0
    tryCatch(
      {
        for (dname in dnames)
        {
          fnames = list.files(dname,pattern="\\.(jpg|jpeg)",ignore.case=T)
          for (fname in fnames)
          {
            strip.fname = gsub(".(jpg|jpeg)","",fname,ignore.case=T) 
            if((!(dname %in% names(object@outline.list)) || !(strip.fname %in% names(object@outline.list[[dname]])) ))
            {
              Rc = .shapeR.find.outline(paste(dname,fname,sep="/"),threshold=threshold,main=paste(dname,fname),mouse.click=mouse.click,display.images=display.images)
              M = cbind(Rc$X,Rc$Y) # Bind together x and y coordinates
              object@outline.list.org[[dname]][[strip.fname]] <- Rc
              object@outline.list[[dname]][[strip.fname]] <- Rc
              if(write.outline.w.org){
                write.image.with.outline(object,folder=dname,fname=strip.fname,doProgress=F)
              }
              
              ptotal=ptotal+1
              setTxtProgressBar(pb, ptotal, label=paste(dname,fname,sep="/"))
            }
          }
                  
        }
      },error = function(e) 
        {message(paste("Error. directory:",dname,", filename:",fname,sep=" "));
         message("The otolith image is possibly not in the center of the image. Try running detect.outline again for this one otolith with mouse.click = TRUE. Click with the mouse on the otolith. If that does not work, you may need to sharpen the image in the Fixed folder using an image manipulation program.")
         print(e)
         e
         },
      finally={
        close(pb);
        setwd(w.dir);
        return(object)
      }
      )
}


#' @title Write outlines on top of the original images for quality checking
#' @description A function which writes the outlines which were extracted from the images in the folder "Fixed" on top of the corresponding images in the "Original" folder. Viewing the resulted images in the folder "Original_with_outlines" is a good quality check to ensure the correct outline has been extracted. If the outline is not correct, then the image can be fixed in an image software, such as GIMP (www.gimp.org), placed in the "Fixed" folder and then the \code{detect.outline} step is repeated.
#' The function \code{\link{detect.outline}} calls this function if the parameter \code{write.outline.w.org} is set to \code{TRUE}.
#' @usage write.image.with.outline(object, folder = NA, fname = NA, doProgress = T)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param folder The folder name where the image is stored
#' @param fname Image file name. Not including the extension ".jpg"
#' @param doProgress If TRUE, a progressbar is shown
#' @examples
#'\dontrun{
#'#Use test data from Libungan and Palsson (2015) and run the following lines:
#'shape = shapeR("ShapeAnalysis/","FISH.csv")
#'shape = detect.outline(shape,write.outline.w.org = FALSE)
#'write.image.with.outline(shape)}
#' @rdname write.image.with.outline
#' @export write.image.with.outline
#' @author Lisa Anne Libungan
#' @references Libungan LA and Palsson S (2015) ShapeR: An R Package to Study Otolith Shape Variation among Fish Populations. PLoS ONE 10(3): e0121102. \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121102}
write.image.with.outline <- 
  function(object,folder=NA,fname=NA,doProgress=T)
{
    image.path = paste(object@project.path,"Original",sep="/")
    write.image.path = paste(object@project.path,"Original_with_outline",sep="/")
    if(!file.exists(write.image.path)){
      dir.create(write.image.path)
    }
            
    w.dir = getwd()
    setwd(image.path)
            
    stocks = c(folder)
    if(is.na(folder))
      stocks = names(object@outline.list)
            
    total=0
    for (dname in stocks)
    {
      total = total + length(object@outline.list[[dname]])
    }
    
    if(doProgress)
      pb <- txtProgressBar(min = 0, max = total, initial = 0, char = "=",style=3,
                         width = NA, title="write.image.with.outline")
    
    ptotal=0
            
            
    par(mfrow = c(1,1))
            
    tryCatch( {
      for( dname in stocks){
        fnames = c(fname)
        if(is.na(fname))
          fnames = names(object@outline.list[[dname]])
                
        dname_woutline = paste(write.image.path,"/",dname,sep="")
        if(!file.exists(dname_woutline)){
          dir.create(dname_woutline)
        }
                
        for(f.name in fnames){
          skra = paste(dname,"/",f.name,".jpg",sep="")
                  
          if(file.exists(skra) ){
            #M=read.pnm(skra)
            M<- readJPEG(skra)
            if(length(dim(M))>2){
              M<- suppressWarnings(pixmapRGB(M[,,1:3]))
              M<- as(M, "pixmapGrey")
            }
            else{
              M<- suppressWarnings(pixmapGrey(M))
            }
            png(paste(dname_woutline,"/",f.name,".png",sep=""),res=200,width=1000,height=1000)
            plot(M,main=paste(dname,f.name,sep="   "))
            with(object@outline.list[[dname]][[f.name]],lines(X,Y,col="red",lwd=2))
            dev.off()
            ptotal=ptotal+1

            if(doProgress)
              setTxtProgressBar(pb, ptotal, label=paste(dname,f.name,sep="/"))
          }
        }
      }
      },
      finally={
        setwd(w.dir);
        if(doProgress)
          close(pb);
      }
      )
}

#' @title Show the extracted outline on top of the original image
#' @description A function which displayes the outlines which were extracted from the image in the "Fixed" folder on top of the corresponding image in the "Original" folder.
#' @usage show.original.with.outline(object, folder, fname)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param folder The folder name where the image is stored
#' @param fname Image file name. Not including the extension ".jpg"
#' @references Libungan LA and Palsson S (2015) ShapeR: An R Package to Study Otolith Shape Variation among Fish Populations. PLoS ONE 10(3): e0121102. \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121102}
#' @examples
#'\dontrun{
#'#Follow the example in Libungan and Palsson (2015) and run the following lines:
#'show.original.with.outline(shape,"IC","403_2")}
#' @rdname show.original.with.outline
#' @export show.original.with.outline
#' @author Lisa Anne Libungan
show.original.with.outline <- 
  function(object,folder,fname)
  {
    image.path = paste(object@project.path,"Original",sep="/")
    
    w.dir = getwd()
    setwd(image.path)
          
    par(mfrow = c(1,1))
    
    tryCatch( {

      skra = paste(folder,"/",fname,".jpg",sep="")
        
      if(file.exists(skra) ){
        #M=read.pnm(skra)
        M<- readJPEG(skra)
        if(length(dim(M))>2){
          M<- suppressWarnings(pixmapRGB(M[,,1:3]))
          M<- as(M, "pixmapGrey")
        }
        else{
          M<- suppressWarnings(pixmapGrey(M))
        }
        plot(M,main=paste(folder,fname,sep="   "))
        with(object@outline.list[[folder]][[fname]],lines(X,Y,col="red",lwd=2))
        
      }
    },
    finally={
      setwd(w.dir);
    }
    )
  }

#' @title Contour smoothing
#' @description Remove high frequency pixel noise around the otolith outline
#' @usage smoothout(object, n)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param n The number of iterations. The default value is 100.
#' @examples
#'
#'\dontrun{ 
#'data(shape)
#'shape = smoothout(shape,n=100)
#'
#'# Plot smoothed outline on top of original outline for comparison
#'outline.org=shape@@outline.list.org[["IC"]][["403_2"]]
#'outline=shape@@outline.list[["IC"]][["403_2"]]
#'plot(outline.org$X,outline.org$Y,type='l',xlab="",ylab="",lwd=2,axes=F)
#'lines(outline$X,outline$Y,col="red",lwd=2)
#'legend("bottomleft",c('Original','Smoothed'),lty=1,col=c('black','red'),lwd=2)}
#'
#' @rdname smoothout
#' @export smoothout
#' @author Lisa Anne Libungan
#' @references Haines, A.J., Crampton, J.S. (2000). Improvements to the method of Fourier shape analysis as applied in morphometric studies. Palaeontology 43: 765-783.
#' @references Claude, J. (2008) Morphometrics with R. Springer. 316 p.
smoothout <- function(object,n=100)
  {
    stocks = names(object@outline.list)
    for( dname in stocks){
      num.fname = length(object@outline.list[[dname]])
      
      for(fname in names(object@outline.list[[dname]])){
        #Outline. X and Y coordinates
        outline = object@outline.list.org[[dname]][[fname]]
        
        M = .shapeR.smoothout(outline,n)
        
        object@outline.list[[dname]][[fname]] <- list(X=M$X,Y=M$Y) # Bind together x and y coordinates
        
      }
    }
    
    return(object)
  }



#' @title Link information in the info.file to the coefficients obtained from the otolith images
#' @description Link the original info file to the otolith coefficients
#' @usage enrich.master.list(object, folder_name = "folder", pic_name = "picname", 
#'                    calibration = "cal", include.wavelet = TRUE, include.fourier = TRUE, 
#'                    n.wavelet.levels = 5, n.fourier.freq = 12,...)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param folder_name Should contain the first letters of the area and the serie or station number of the sample, for example: "IC"
#' @param pic_name Should contain the serie number of a given sample and fish number, for example "403_2" (not including the .jpg extension)
#' @param calibration The name of the column where the pixels to measurement calibration is located
#' @param include.wavelet If TRUE, the wavelet coefficient are included in the master.list
#' @param include.fourier If TRUE then the Normalized Elliptic Fourier coefficients are included in the master.list
#' @param n.wavelet.levels Integer saying how many levels of wavelet levels should be included
#' @param n.fourier.freq Integer saying how many Fourier frequency levels should be included
#' @param ... Additional parameter for \code{read.csv} for reading the info.file
#' @return A \code{\linkS4class{shapeR}} object with values in slots:
#' \itemize{
#'   \item wavelet.coef
#'   \item fourier.coef
#'   \item shape
#'   \item filter
#'   \item master.list
#' }
#' @examples
#'\dontrun{
#'data(shape)
#'shape = generateShapeCoefficients(shape)
#'
#'shape = enrich.master.list(shape)}
#'
#' @rdname enrich.master.list
#' @export enrich.master.list
#' @author Lisa Anne Libungan
enrich.master.list <- 
          function(object,folder_name = "folder",pic_name="picname",calibration="cal",include.wavelet=TRUE,include.fourier=TRUE,n.wavelet.levels=5,n.fourier.freq=12,...)
          { 
            master.list = object@master.list.org
            master.list$key = paste(unlist(master.list[folder_name]),unlist(master.list[pic_name]),sep=";")
            
            .shapeR.check.outline.list(object) #Stops execution if no outlines
            
            .shapeR.check.shape.coefs(object) #Stops execution of shape coefficients
            
            shape.coef = object@shape.coef.raw
            org.dim = dim(master.list)
            master.list = merge(master.list,data.frame(key=rownames(shape.coef),shape.coef),by.x="key",by.y="key",suffixes = c("",".a"),all.x=TRUE)
            
            cal.n = master.list[,calibration]
            master.list$otolith.area =  master.list$otolith.area * 1/(cal.n^2)
            master.list$otolith.length =   master.list$otolith.length * 1/cal.n
            master.list$otolith.width =  master.list$otolith.width * 1/cal.n
            master.list$otolith.perimeter =  master.list$otolith.perimeter * 1/cal.n

            object@shape.coef = as.matrix( master.list[,(org.dim[2]+1):dim(master.list)[2]])
            
            if(include.wavelet){
              wavelet.coef = object@wavelet.coef.raw[,0:(sum(2^(1:n.wavelet.levels))+1)+1]
              org.dim = dim(master.list)
              master.list = merge(master.list,data.frame(key=rownames(wavelet.coef),wavelet.coef),by.x="key",by.y="key",suffixes = c("",".a"),all.x=TRUE)
              
              object@wavelet.coef = as.matrix( master.list[,(org.dim[2]+2):dim(master.list)[2]])
            }
            
            if(include.fourier)
            {
              f.seq = seq(1,32*4,by=32)
              f.seq=c(sapply(0:(n.fourier.freq-1),function(x) f.seq+x))[-c(1,2,3)]
              org.dim = dim(master.list)
              master.list = merge(master.list,data.frame(key=rownames(object@fourier.coef.raw[,f.seq]),object@fourier.coef.raw[,f.seq]),by.x="key",by.y="key",suffixes = c("",".a"),all.x=TRUE)
              object@fourier.coef = as.matrix( master.list[,(org.dim[2]+1):dim(master.list)[2]])
            }
                        
            object@master.list = master.list
            object = setFilter(object)
            
            return(object)
          }



#' @title Mean otolith shape based on wavelet reconstruction
#' @description A function for showing the mean otolith shape based on wavelet reconstruction
#' @usage plotWaveletShape(object, class.name,show.angle=FALSE,lty=1:5,col=1:6,...)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param class.name A string as the column name in the master list
#' @param show.angle If TRUE angles are shown on the plot
#' @param lty,col Vector of line types and colors. Values are used cyclically.
#' @param ... Additional parameters to be passed to 'plot'
#' @examples
#'data(shape)
#'plotWaveletShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1)
#' @rdname plotWaveletShape
#' @export plotWaveletShape
#' @author Lisa Anne Libungan
#' @references Nason, G. (2012) \code{\link{wavethresh}}: Wavelets statistics and transforms, version 4.5. R package.
plotWaveletShape <- 
  function(object,class.name,show.angle=FALSE,lty=1:5,col=1:6,...)
  {
    first=T
    
    classes = object@master.list[class.name][,1]
    lev.classes = levels(factor(classes[object@filter]))

    k = length(lev.classes)
    if (length(lty) < k) 
      lty <- rep_len(lty, k)
    if (length(col) < k) 
      col <- rep_len(col, k)
    
    for(ctr in lev.classes)
    {
      i = which(ctr == lev.classes)
              
      ind = classes== ctr & object@filter
      
      m.w.coef = apply(object@wavelet.coef[ind,,drop=F],2,mean)

      iw.utl = .shapeR.inverse.wavelet(m.w.coef,mean(object@master.list$mean.radii[ind]))
              
      if(first==T){
        plot(iw.utl$X,iw.utl$Y,type="l",col=col[i],xlab="",ylab="",xlim=range(iw.utl$X)*1.3,ylim=range(iw.utl$Y)*1.2, lty=lty[i],axes=F,frame.plot=F,...)
        first=F
      }
      else
        lines(iw.utl$X,iw.utl$Y,col=col[i],lty=lty[i],...)
    }
    legend("bottomleft",lev.classes,cex=0.9,lty=lty[1:k],col=col[1:k],...)
    
    if(show.angle)
    {
      centr.mr.x = mean(iw.utl$X)
      centr.mr.y = mean(iw.utl$Y)
      
      abline(h=centr.mr.y,lty=2,col="grey")
      abline(v=centr.mr.x,lty=2,col="grey")
      y.r = range(iw.utl$Y)
      x.r = range(iw.utl$X)
      y.dev = (y.r[2] - y.r[1])*.05
      x.dev = (x.r[2] - x.r[1])*.07
      text(x.dev,max(iw.utl$Y)*1.25,"90",cex=1.2)
      text(x.dev,min(iw.utl$Y)*1.2,"270",cex=1.2)
      text(max(iw.utl$X)*1.2,y.dev,"0",cex=1.2)
      text(min(iw.utl$X)*1.2,y.dev,"180",cex=1.2)
    }
    
  }



#' @title Mean otolith shape based on Fourier reconstruction
#' @description A function for showing the mean otolith shape based on Fourier reconstruction
#' @usage plotFourierShape(object, class.name, show.angle = FALSE,lty=1:5,col=1:6, ...)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param class.name A string as the column name in the master list
#' @param show.angle If TRUE angles are shown on the plot
#' @param lty,col Vector of line types and colors. Values are used cyclically.
#' @param ... Additional parameters to be passed to 'plot'
#' @examples
#'data(shape)
#'plotFourierShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1)
#' @rdname plotFourierShape
#' @export plotFourierShape
#' @author Lisa Anne Libungan
#' @references Libungan LA and Palsson S (2015) ShapeR: An R Package to Study Otolith Shape Variation among Fish Populations. PLoS ONE 10(3): e0121102. \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121102}
plotFourierShape <- 
  function(object,class.name,show.angle=FALSE,lty=1:5,col=1:6,...)
  {
    coef.f = cbind(-1,0,0,object@fourier.coef)
    first=T

    classes = object@master.list[class.name][,1]
    lev.classes = levels(factor(classes[object@filter]))
    
    k = length(lev.classes)
    if (length(lty) < k) 
      lty <- rep_len(lty, k)
    if (length(col) < k) 
      col <- rep_len(col, k)
    
    for(ctr in lev.classes)
    {
      i = which(ctr == lev.classes)
         
      ind = classes== ctr & object@filter
      tmp.Mx1 = apply(coef.f[ind,,drop=F],2,mean)
      i.seq = seq(1,12*4,by=4)
      iw.utl = .shapeR.iefourier(tmp.Mx1[i.seq],tmp.Mx1[i.seq+1],tmp.Mx1[i.seq+2],tmp.Mx1[i.seq+3],12,64*4)
              
      if(first==T){
                
        plot(iw.utl$x,iw.utl$y,type="l",col=col[i],xlab="",ylab="",xlim=range(iw.utl$x)*1.3,ylim=range(iw.utl$y)*1.3, lty=lty[i],axes=F,frame.plot=F,...)
        first=F
      }
      else
        lines(iw.utl$x,iw.utl$y,col=col[i],lty=lty[i],...) 
    }
    legend("bottomleft",lev.classes,cex=0.9,lty=lty,col=col,...)

    if(show.angle == T){
      centr.mr.x = mean(iw.utl$x)
      centr.mr.y = mean(iw.utl$y)
      
      abline(h=centr.mr.y,lty=2,col="grey")
      abline(v=centr.mr.x,lty=2,col="grey")
      y.r = range(iw.utl$y)
      x.r = range(iw.utl$x)
      y.dev = (y.r[2] - y.r[1])*.05
      x.dev = (x.r[2] - x.r[1])*.05
      text(x.dev,max(iw.utl$y)*1.25,"90",cex=1.2)
      text(x.dev,min(iw.utl$y)*1.2,"270",cex=1.2)
      text(max(iw.utl$x)*1.2,y.dev,"0",cex=1.2)
      text(min(iw.utl$x)*1.2,y.dev,"180",cex=1.2)
      
    
    }
            
  }
          

#' @title Standardize coefficients
#' @description Function to standardized the wavelet and Fourier coefficients for a specific parameter such as the fish length.
#' For each country/population a regression coefficient is calculated as a function of fish length.
#' If the slope is significantly different from zero, a correction is made according to Lleonart et al 2000.
#' First ANCOVA is performed: variable ~ pop*length_cm, following a method by Longmore et al 2010.
#' If there is a significant interaction between population and length_cm, then the coefficients are not used and automatically discarded. If there is no interaction, the 
#' coefficients are kept and standardized with regards to fish length.
#' @usage stdCoefs(object, classes=NA, std.by, std.type = "mean", p.crit = 0.05,bonferroni= FALSE)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param classes The classes to be grouped for standardization. Should be the same as used for the statistical tests
#' @param std.by The parameter to be used for standardization. Typically the length of the fish from the master.list.
#' @param std.type The tuning of the standardization. The standardization can be sensitive to what value all the fishes are standardized to. Possible values are:
#' \itemize{
#' \item min Standardized as the minimum value of std.by
#' \item mean Standardized as the mean value of std.by
#' \item max Standardized as the maximum value of std.by
#' }
#' @param p.crit An argument used to select the threshold critera for omitting coefficients which show interaction with fish length. If p.crit = 0.05, all coefficients which have p<0.05 are omitted. If p.crit = 0.01, only coefficients with p<0.01 are omitted. 
#' @param bonferroni A logical parameter for performing Bonferroni for multiple testing 
#' @examples
#'data(shape)
#'shape = stdCoefs(shape,classes="pop","length_cm")
#' @rdname stdCoefs 
#' @export stdCoefs
#' @author Lisa Anne Libungan
#' @references Lleonart, J., Salat, J. & Torres, G.J. (2000) Removing allometric effects of body size in morphological analysis. Journal of Theoretical Biology, 205, 85-93.
#' @references Longmore, C., Fogarty, K., Neat, F., Brophy, D., Trueman, C., Milton, A. & Mariani, S. (2010) A comparison of otolith microchemistry and otolith shape analysis for the study of spatial variation in a deep-sea teleost, \emph{Coryphaenoides rupestris}. Environmental Biology of Fishes, 89, 591-605.
#' @references Reist, J.D. (1985) An Empirical-Evaluation of Several Univariate Methods That Adjust for Size Variation in Morphometric Data. Canadian Journal of Zoology-Revue Canadienne De Zoologie, 63, 1429-1439.
stdCoefs <- 
  function(object,classes=NA,std.by,std.type="mean",p.crit=0.05,bonferroni= FALSE)
  {
    filter = object@filter
    in.classes = rep(NA,sum(filter))
    if(!is.na(classes))
      in.classes = object@master.list[classes][filter,1]
    
    ind = which(filter)
    message("Wavelet standardization. ", appendLF=F)
    tmp.wavelet.coef.info = .shapeR.coef.standardize.f(object@wavelet.coef[filter,],in.classes,object@master.list[std.by][filter,1],std.type,is.na.classes = is.na(classes),p.crit=p.crit,bonferroni=bonferroni)
    tmp.wavelet.coef.std = tmp.wavelet.coef.info$coef.standard
    tmp.wavelet.r.y = tmp.wavelet.coef.info$r.y
    message("Fourier standardization. ", appendLF=F)
    tmp.fourier.coef.info = .shapeR.coef.standardize.f(object@fourier.coef[filter,],in.classes,object@master.list[std.by][filter,1],std.type,is.na.classes = is.na(classes),p.crit=p.crit,bonferroni=bonferroni)
    tmp.fourier.coef.std = tmp.fourier.coef.info$coef.standard
    tmp.fourier.r.y = tmp.fourier.coef.info$r.y
    message("Measurement standardization. ",appendLF=F)
    tmp.shape.info = .shapeR.coef.standardize.f(object@shape.coef[filter,],in.classes,object@master.list[std.by][filter,1],std.type,is.na.classes = is.na(classes),p.crit=p.crit,bonferroni=bonferroni)
    tmp.shape.std = tmp.shape.info$coef.standard
    tmp.shape.r.y = tmp.shape.info$r.y
    
    
    object@wavelet.coef.std = matrix(nrow=length(filter),ncol=ncol(tmp.wavelet.coef.std))
    object@fourier.coef.std = matrix(nrow=length(filter),ncol=ncol(tmp.fourier.coef.std))
    object@shape.std = matrix(nrow=length(filter),ncol=ncol(as.matrix(tmp.shape.std)))
    
    object@wavelet.coef.std[ind,] = tmp.wavelet.coef.std
    object@fourier.coef.std[ind,] = tmp.fourier.coef.std
    object@shape.std[ind,] = tmp.shape.std
    
    object@wavelet.coef.std.removed = if(!is.null(tmp.wavelet.r.y)) tmp.wavelet.r.y else vector()
    object@fourier.coef.std.removed = if(!is.null(tmp.fourier.r.y)) tmp.fourier.r.y else vector()
    object@shape.std.removed = if(!is.null(tmp.shape.r.y)) tmp.shape.r.y else vector()

    if(length(object@shape.std.removed) == 0 ||(length(object@shape.std.removed) == 1 && is.na(object@shape.std.removed)))
      colnames(object@shape.std) = colnames(object@shape.coef)
    else
      colnames(object@shape.std) = colnames(object@shape.coef)[-object@shape.std.removed]
    
    return(object)
  }

#' @title Read updated master list
#' @description Reads an updated master list. This is important to run if you want to ensure that a updated master list is used in the analysis.
#' @usage read.master.list(object, ...)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param ... Additional parameter for \code{read.csv} for reading the info.file
#' @return \code{\linkS4class{shapeR}} object with values in slots:
#' \itemize{
#'   \item master.list.org
#' }
#' @rdname read.master.list
#' @export read.master.list
#' @author Lisa Anne Libungan
read.master.list <- function(object,...) {
  object@master.list.org = read.csv(paste(object@project.path,object@info.file,sep="/"),header=T,...)
  return(object)
}

#' @title Get standardized wavelet coefficients, filtered according to filter
#' @description Returns the standardized wavelet coefficients determined in  \code{stdCoefs}.
#' Returns only values as set in the slot \code{filter}
#' @usage getStdWavelet(object)
#' @param object \code{\linkS4class{shapeR}} object 
#' @return The standardized wavelet coefficients for all valid fish as determined by the slot \code{filter}
#' @rdname getStdWavelet
#' @export getStdWavelet
#' @author Lisa Anne Libungan
getStdWavelet <- function(object) {
  return(object@wavelet.coef.std[object@filter,])
}

#' @title Get wavelet coefficients, filtered according to filter
#' @description Returns the wavelet coefficients determined in \code{\link{generateShapeCoefficients}}.
#' Returns only values as set in the slot \code{filter}
#' @usage getWavelet(object)
#' @param object \code{\linkS4class{shapeR}} object 
#' @return The wavelet coefficients for all valid fish as determined by the slot \code{filter}
#' @rdname getWavelet
#' @export getWavelet
#' @author Lisa Anne Libungan
getWavelet <- function(object) {
  return(object@wavelet.coef[object@filter,])
}

#' @title Get simple shape variables, filtered according to filter
#' @description Returns shape variables length, width, perimeter and area determined in \code{\link{generateShapeCoefficients}}.
#' Returns only values as set in the slot \code{filter}.
#' These variables can only be obtained if the calibration measurements in pixels have been registered in the csv data file in a column labelled 'cal' (see example data file). To get the calibration measurements, use a image manipulation program and measure 1mm on the calibration measurement stick (that was taken for that particular dataset) and register how many pixels 1mm is into the column 'cal'.
#' @usage getMeasurements(object)
#' @param object \code{\linkS4class{shapeR}} object 
#' @return A data frame with all valid fish as determined by the slot \code{filter} and with columns:
#' \itemize{
#'   \item otolith.area
#'   \item otolith.length
#'   \item otolith.width
#'   \item otolith.perimeter
#' }
#' @examples
#'data(shape)
#'# Calculate the mean otolith area for each fish population
#'# The results are in square mm since the calibration ('cal') column 
#'# in the data file is in pixels (1 mm/pixel).
#'tapply(getMeasurements(shape)$otolith.area, getMasterlist(shape)$pop,mean)
#' @rdname getMeasurements
#' @export getMeasurements
#' @author Lisa Anne Libungan
getMeasurements <- function(object) {
  return(as.data.frame(object@shape.coef[object@filter,]))
}

#' @title Get simple shape variables after standardization, filtered according to filter
#' @description Returns the simple shape variables determined in \code{stdCoefs}.
#' Returns only values as set in the slot \code{filter}
#' @usage getStdMeasurements(object)
#' @param object \code{\linkS4class{shapeR}} object 
#' @return A data frame with all valid fish as determined by the slot \code{filter}. Returns only variables that have not been removed after standardization.
#' @examples
#'data(shape)
#'#Calculate the mean standardized otolith length for each fish population
#'tapply(getStdMeasurements(shape)$otolith.length,
#'getMasterlist(shape)$pop,mean)
#' @rdname getStdMeasurements
#' @export getStdMeasurements
#' @author Lisa Anne Libungan
getStdMeasurements <- function(object) {
  stdM = as.data.frame(object@shape.std[object@filter,])
  names(stdM) = colnames(object@shape.std[object@filter,,drop=FALSE])
  return(stdM)
}

#' @title Get standardized Fourier coefficients, filtered according to filter
#' @description Returns the standardized Fourier coefficients determined in  \code{\link{stdCoefs}}.
#' Returns only values as set in the slot \code{filter}
#' @usage getStdFourier(object)
#' @param object \code{\linkS4class{shapeR}} object 
#' @return The standardized Fourier coefficients for all valid fish as determined by the slot \code{filter}
#' @rdname getStdFourier
#' @export getStdFourier
#' @author Lisa Anne Libungan
getStdFourier <- function(object) {
  return(object@fourier.coef.std[object@filter,])
}

#' @title Get Fourier coefficients, filtered according to filter
#' @description Returns the Fourier coefficients determined in \code{\link{stdCoefs}}.
#' Returns only values as set in \code{\link{setFilter}}
#' @usage getFourier(object)
#' @param object \code{\linkS4class{shapeR}} object 
#' @return The Fourier coefficients for all fish as determined by \code{\link{setFilter}}
#' @rdname getFourier
#' @export getFourier
#' @author Lisa Anne Libungan
getFourier <- function(object) {
  return(object@fourier.coef[object@filter,])
}

#' @title Get filtered master.list values
#' @description Returns selected values from \code{master.list}
#' @usage getMasterlist(object, useFilter = TRUE)
#' @param object \code{\linkS4class{shapeR}} object
#' @param useFilter If TRUE, the master.list values are filtered by the slot \code{filter}. FALSE = no filtering.
#' @return The \code{master.list} is filtered by the slot \code{filter} if the \code{useFilter} is TRUE, else no filtering is done.
#' @rdname getMasterlist
#' @export getMasterlist
#' @author Lisa Anne Libungan
getMasterlist <- function(object, useFilter = TRUE) {
  if(useFilter)
    return(object@master.list[object@filter,])
  else
    return(object@master.list)
}

#' @title Set a filter to analyze the shape data 
#' @description Sets a filter on \code{master.list}. Here it is possible to filter the \code{master.list} by specific ages, maturity stages, areas, etc.
#' If no value is set, all data with shape parameters are used
#' @usage setFilter(object, filter)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param filter A vector restricting the new filter value. Only otoliths having shape parameters are selected.
#' @return A \code{\linkS4class{shapeR}} object with the slot \code{filter} set.
#' @examples
#'data(shape)
#'#Filter only Icelandic and Norwegian samples
#'shape = setFilter(shape,
#'getMasterlist(shape, useFilter = FALSE)$pop %in% c("NO","IC")) 
#'table(getMasterlist(shape)$pop)
#'#Reset filter
#'shape = setFilter(shape)
#'table(getMasterlist(shape)$pop)
#' @rdname setFilter
#' @export setFilter
#' @author Lisa Anne Libungan
setFilter <- function(object,filter=NULL) {
  if(is.null(filter)){
    object@filter = !is.na(object@master.list$Ws0c1) | !is.na(object@master.list$FD1d)
  }
  else{
    object@filter = object@filter & filter
  }
  
  return(object)
}



#' @title Mean and standard deviation of the wavelet coefficients
#' @description The mean and standard deviation of the wavelet coefficients
#' @usage plotWavelet(object, level, start.angle = 0, class.name=NULL,useStdcoef=FALSE,...)
#' @param object A \code{\linkS4class{shapeR}} object
#' @param level The wavelet level to be shown
#' @param start.angle The angle in degrees, the x-axis should start on
#' @param class.name Column name in master list for partitioning the data into groups and showing the ratio of variation among groups.
#' @param useStdcoef Choose "TRUE" or "FALSE" if coefficients should be standardized or not
#' @param ... Additional parameters to be passed to 'plot'
#' @examples 
#'data(shape)
#'shape = stdCoefs(shape,classes="pop","length_cm")
#'plotWavelet(shape,level=5,class.name= "pop",useStdcoef=TRUE)
#' @rdname plotWavelet
#' @export plotWavelet
#' @author Lisa Anne Libungan
  plotWavelet <-
  function(object,level,start.angle=0,class.name=NULL,useStdcoef=FALSE,...){
    ind = rep(T,2^level)
    level.ind = (sum(2^(c(1:(level-1)))):sum(2^(c(1:level))))[-1]+1
    cnt.less = 0

    if(useStdcoef){
      .shapeR.check.shape.coefs.std(object)
      ind.o = unlist(sapply(object@wavelet.coef.std.removed , function(x) which(x == level.ind)))
      if(length(ind.o)>0)
      {
        ind[ind.o] = F
      }
      cnt.less = sum(object@wavelet.coef.std.removed < min(level.ind))
    }

    sw.coef = object@wavelet.coef[object@filter,level.ind]
    if(useStdcoef)
      sw.coef = object@wavelet.coef.std[object@filter,(1:sum(ind)) +(level.ind[1]-cnt.less-1)]
    p.angle = ((seq(-pi,pi-0.0002,by=2*pi/2^level)+2*pi)%%(2*pi))[ind]*180/pi
    a.ind = order(p.angle)
    temp.angle = p.angle[a.ind]-start.angle
    ind.pos = which(temp.angle>=0)
    a.ind = c(a.ind[ind.pos],a.ind[-ind.pos])
    xaxis = p.angle[a.ind]
    if(!is.null(class.name))
      par(mar=c(5,4,4,5)+.1)
    else
      par(mar=c(5, 4, 4, 2)+.1)
    suppressWarnings(plotCI(1:length(xaxis),apply(sw.coef,2,mean)[a.ind],type="p",uiw=apply(sw.coef,2,sd)[a.ind],xaxt="n",xlab="Angle",ylab="Mean +/-sd",...))
    xaxis2 = c(seq(start.angle,360,by=10),seq(0,start.angle,by=10))
    xaxis2 = xaxis2[-length(xaxis2)]
    axis(side = 1, at = which(ceiling(xaxis2)%%60==0)*length(xaxis)/length(xaxis2), labels=xaxis2[which(ceiling(xaxis2)%%60==0)])
    #axis(side = 1, at = which(ceiling(xaxis)%%5==0), labels = xaxis[which(ceiling(xaxis)%%5==0)])

    tryCatch(
      {
        if(!is.null(class.name)){
          
          .shapeR.check.shape.coefs.std(object)
          classes = getMasterlist(object)[class.name][,1]
          sVar = .shapeR.variation.among(classes,sw.coef)
          
          par(new=TRUE)
          plot(1:length(xaxis),sVar[a.ind],type="l", col="black",xaxt="n",yaxt="n",xlab="",ylab="")
          axis(4)
          mtext("ICC",side=4,line=3)
        }
      }
    ,finally=par(new=FALSE))
          
}


#' @title Mean and standard deviation of the Fourier coefficients
#' @description The mean and standard deviation of the Fourier coefficients
#' @usage plotFourier(object, coef.index=NULL,class.name=NULL,useStdcoef=FALSE, ...)
#' @param object \code{\linkS4class{shapeR}} object
#' @param coef.index An index vector for which fourier coefficents to be shown. Default is \code{NULL} and all coefficients are shown.
#' @param class.name Column name in master list for partitioning the data into groups and showing the ratio of variance among to the sum of variance among and variance within.
#' @param useStdcoef Boolean saying if to use the standardized coefficients or not
#' @param ... Additional parameters to be passed to 'plot'
#' @examples
#'data(shape)
#'shape = stdCoefs(shape,classes="pop","length_cm")
#'plotFourier(shape,class.name= "pop",useStdcoef=TRUE)
#' @rdname plotFourier
#' @export plotFourier
#' @author Lisa Anne Libungan
plotFourier <- 
  function(object,coef.index=NULL,class.name=NULL,useStdcoef=FALSE,...){ 
    if(is.null(coef.index) && useStdcoef==FALSE)
      coef.index = 1:dim(object@fourier.coef)[2]
    if(is.null(coef.index) && useStdcoef==TRUE)
      coef.index = 1:dim(object@fourier.coef.std)[2]
    
    
    coef = object@fourier.coef[object@filter,coef.index]
    if(useStdcoef)
      coef = object@fourier.coef.std[object@filter,coef.index]

    if(!is.null(class.name))
      par(mar=c(5,4,4,5)+.1)
    else
      par(mar=c(5, 4, 4, 2)+.1)
    
    suppressWarnings(plotCI(apply(coef,2,mean),type="p",uiw=apply(coef,2,sd),xlab="coefficient",ylab="Mean +/-sd",...))
    lines(apply(coef,2,mean),lty=2)

    tryCatch(
    {
      if(!is.null(class.name)){
        .shapeR.check.shape.coefs.std(object)
        classes = getMasterlist(object)[class.name][,1]
        sVar = .shapeR.variation.among(classes,coef)
        
        par(new=TRUE)
        plot(coef.index,sVar,type="l", col="black",xaxt="n",yaxt="n",xlab="",ylab="")
        axis(4)
        mtext("ICC",side=4,line=3)
      }
    }
    ,finally=par(new=FALSE))

  
  }


#' @title Estimate the outline reconstruction based on Fourier/wavelet compared to the outlines that have not been transformed
#' @description Estimate outline reconstruction using a different number of coefficients of wavelet and Fourier compared to the original otolith
#' @usage estimate.outline.reconstruction(object, ...)
#' @param object \code{\linkS4class{shapeR}} object
#' @param ... Additional parameters to be passed to 'plot' and 'points'
#' @return A list containing values
#' \itemize{
#'  \item w.dev.m a list for number of coefficients for mean error of wavelet reconstruction
#'  \item w.dev.sd a list for number of coefficients for standard deviation of wavelet reconstruction
#'  \item f.power.total Fourier power for number of Fourier harmonics
#' }
#' @examples
#'\dontrun{
#'data(shape)
#'estimate.outline.reconstruction(shape)}
#'
#' @rdname estimate.outline.reconstruction
#' @export estimate.outline.reconstruction
#' @author Lisa Anne Libungan
#' @references Claude, J. (2008) Morphometrics with R. Springer. 316 p.
estimate.outline.reconstruction <- 
  function(object,... ){
    
    .shapeR.check.outline.list(object) #Stops execution if no outlines
    
    .shapeR.check.shape.coefs(object) #Stops execution of shape coefficients
    
    total=0
    for( dname in names(object@outline.list))
      total=total+length(object@outline.list[[dname]])
    
    
    pb <- txtProgressBar(min = 0, max = total, initial = 0, char = "=",style=3,
                         width = NA, title="estimate.outline.reconstruction")          
    ptotal=0
    #Go over coefficients
    r.names = rownames(object@wavelet.coef.raw)
    t.deviation = data.frame(level=c(),w.deviation=c(),dname=c(),fname=c())
    
    t.deviation.f = data.frame(harmonic=c(),f.deviation=c(),dname=c(),fname=c())
    
    for(i in 1:length(r.names)){
      rname.split = unlist(strsplit(r.names[i],";"))
      dname = rname.split[1]
      fname = rname.split[2]
              
      outline = object@outline.list[[dname]][[fname]]
      if(is.null(outline))
        next
      im.area = .shapeR.otolith.image.parameters(cbind(outline$X,outline$Y))$area
      mynd = .shapeR.rotate.svd(outline)
      mynd= .shapeR.centroid(list(X=mynd[,1],Y=mynd[,2]))
      mynd$X = mynd$X/sqrt(im.area)
      mynd$Y = mynd$Y/sqrt(im.area)
      ppolar= .shapeR.regularradius(mynd$X,mynd$Y,1024)
      # Done like in Morph book, pages 227-228
      for(n.levels in 2:9){
        iw.utl = .shapeR.inverse.wavelet(object@wavelet.coef.raw[i, 1:sum(2^(0:n.levels)) +1],object@wavelet.coef.raw[i,1],n.levels=n.levels,return.polar=TRUE)
        w.deviation = sqrt( sum((ppolar$radii-iw.utl$radii)^2)/sum(ppolar$radii^2)) #Normalized error
        t.deviation = rbind(t.deviation,data.frame(level=n.levels,w.deviation=w.deviation,dname=dname,fname=fname))
      }
       
      N<-.shapeR.efourier(cbind(outline$X,outline$Y),n=32)

      for(n.harm in 1:32)
      {
        iw.utl = .shapeR.iefourier(N$an[1:n.harm],N$bn[1:n.harm],N$cn[1:n.harm],N$dn[1:n.harm],k=n.harm,n=length(outline$X),ao=N$ao,co=N$co)
        f.deviation = sqrt( sum((outline$X-iw.utl$x)^2+(outline$Y-iw.utl$y)^2)/sum(outline$X^2+outline$Y^2)) #Normalized error
        t.deviation.f = rbind(t.deviation.f,data.frame(harmonic=n.harm,f.deviation=f.deviation,dname=dname,fname=fname))
      }
      
      ptotal=ptotal+1
      setTxtProgressBar(pb, ptotal)        
    }
    close(pb)
    
    w.dev.m =aggregate(t.deviation$w.deviation*100,by=list(t.deviation$level),mean)
    w.dev.sd =aggregate(t.deviation$w.deviation*100,by=list(t.deviation$level),sd)
    
    f.dev.m =aggregate(t.deviation.f$f.deviation*100,by=list(t.deviation.f$harmonic),mean)
    f.dev.sd =aggregate(t.deviation.f$f.deviation*100,by=list(t.deviation.f$harmonic),sd)
            
    return(list(w.dev.m=w.dev.m,w.dev.sd=w.dev.sd,f.dev.m=f.dev.m,f.dev.sd=f.dev.sd))
  }

#' @title Plot outline reconstruction
#' @description Show graphs of the reconstruction using different number of levels of wavelet reconstruction and Fourier power using different number of Fourier harmonics. 
#' Uses the output from \link{estimate.outline.reconstruction} 
#' @usage outline.reconstruction.plot(outline.rec.list,ref.w.level=5,
#'                             ref.f.harmonics=12,max.num.harmonics=32,...)
#' @param outline.rec.list The output from \link{estimate.outline.reconstruction}
#' @param ref.w.level Reference level for graphical purposes. The default is 5 as is the default of \linkS4class{shapeR}.
#' @param ref.f.harmonics Reference Fourier harmonize. The default is 12 as is the default in \linkS4class{shapeR}.
#' @param max.num.harmonics Maxinum number of Fourier harmonics to be shown
#' @param ... Additional parameters to be passed to 'plot'
#' @examples
#'\dontrun{data(shape)
#'est.list = estimate.outline.reconstruction(shape)
#'outline.reconstruction.plot(est.list,panel.first = grid())}
#' @rdname outline.reconstruction.plot
#' @export outline.reconstruction.plot
#' @author Lisa Anne Libungan
outline.reconstruction.plot <- function(outline.rec.list,ref.w.level=5, ref.f.harmonics=12,max.num.harmonics=32,...)
{  
  with(outline.rec.list, 
    {par(mfrow=c(1,2))
     ylim.v = c(0,max(w.dev.m[,2]+w.dev.sd$x*qnorm(0.95),f.dev.m[,2]+f.dev.sd$x*qnorm(0.95)))
    suppressWarnings(plotCI(w.dev.m[,1],w.dev.m[,2],uiw=w.dev.sd$x*qnorm(0.95),type="o",err='y',ylim=ylim.v,xlab="Wavelet levels",ylab="Deviation from outline %",...))
    abline(v=ref.w.level,col="red")
    
    suppressWarnings(plotCI(f.dev.m[1:max.num.harmonics,1],f.dev.m[1:max.num.harmonics,2],uiw=f.dev.sd$x[1:max.num.harmonics]*qnorm(0.95),type="o",err='y',ylim=ylim.v,xlab="Fourier harmonics",ylab="Deviation from outline %",...))
    abline(v=ref.f.harmonics,col="red") 
    }
  )
}
#' @title Plot data clusters
#' @description Plots data clusters
#' @usage cluster.plot( ddata, classes, main="", col.stock=NULL,
#'               plotCI = FALSE, conf.level = 0.68, ...)
#' @param ddata Matrix of points
#' @param classes A factor including the cluster values
#' @param main Title for the plot
#' @param col.stock Colors for the plotted classes
#' @param plotCI Plot means with confidence intervals
#' @param conf.level The confidence interval for the standard error of the mean
#' @param ... Additional parameters to be passed to 'plot' or 'ldahist' if one dimension
#' @examples
#'data(shape)
#'library(vegan)
#'cap.res = capscale(getStdWavelet(shape) ~ getMasterlist(shape)$pop)
#'
#'eig=eigenvals(cap.res,constrained=TRUE)
#'eig.ratio = eig/sum(eig)
#'
#'cluster.plot(scores(cap.res)$sites[,1:2],getMasterlist(shape)$pop
#',plotCI=TRUE
#',xlab=paste("CAP1 (",round(eig.ratio[1]*100,1),"%)",sep="")
#',ylab=paste("CAP2 (",round(eig.ratio[2]*100,1),"%)",sep="")
#',main="Canonical clustering"
#')
#'
#' @rdname cluster.plot
#' @export cluster.plot
#' @author Lisa Anne Libungan
#' @references Oksanen, J., Blanchet, F.G., Kindt, R., Legendre, P., Minchin, P.R., O'Hara, R.B., Simpson, G.L., Solymos, P., Stevens, M.H.H. and Wagner, H. (2013). \code{\link{vegan}}: Community Ecology Package. R package version 2.0-10.
cluster.plot <- function( ddata,classes,main="",col.stock = NULL,plotCI=FALSE,conf.level=0.68,...)
{  
  dim.dd = dim(ddata)
  
  if(is.null(dim.dd) || dim.dd[2]==1){
    ld1.mean = tapply(ddata,classes,mean)
  }
  else
  {
    ld1.mean = tapply(ddata[,1],classes,mean)
    ld2.mean = tapply(ddata[,2],classes,mean)
  }
    
  if(is.null(dim.dd) || dim.dd[2]==1){
    ld1.se = tapply(ddata,classes,function(x) sd(x)/sqrt(length(x)))
  }
  else
  {    
    ld1.se = tapply(ddata[,1],classes,function(x) sd(x)/sqrt(length(x)))    
    ld2.se = tapply(ddata[,2],classes,function(x) sd(x)/sqrt(length(x)))    
  }
        

  if(is.null(col.stock))
  {
    col.stock=1:length(levels(classes))
    names(col.stock) = levels(classes)
  }
  
  if(is.null(dim.dd) || dim.dd[2]==1)
    MASS::ldahist(ddata, classes, type="histogram",...)
  else
  {
    i=1
    for( lev in levels(classes)){
      ind = which(lev==classes)
      data = ddata[ind,]
      if(i==1)
        plot(data,cex=0.6,col=col.stock[[lev]],pch=lev,main=main,...)
      else
        points(data,cex=0.6,col=col.stock[[lev]],pch=lev,...)
        
      i=i+1
        
    }
          
    text(ld1.mean,ld2.mean,labels = levels(classes), cex=1.5,font=2)
      
    if(plotCI==TRUE){
      r = length(unique(classes)) # Number of means
      nis =dim(ddata)[1] # Number of samples
      df = nis -r
        
      suppressWarnings(plotCI(ld1.mean,ld2.mean,uiw=ld1.se*qtukey(conf.level,r,df),err='x',pch=NA,add=TRUE))
      suppressWarnings(plotCI(ld1.mean,ld2.mean,uiw=ld2.se*qtukey(conf.level,r,df),err='y',pch=NA,add=TRUE))
        
    }
  }  
}




#' @title Remove otolith outline
#' @description A function for removing an otolith outline from the file 'outline.list'. Typically done if the image is of bad quality and needs to be enhanced in a image processing software
#' @usage remove.outline(object, folder = "", fname = "")
#' @param object A \code{\linkS4class{shapeR}} object
#' @param folder The folder name where the outline that needs to be removed is stored
#' @param fname The file name of the outline to be removed
#' @return \code{\linkS4class{shapeR}} object
#' @examples
#'\dontrun{
#'#Use test data from example in Libungan and Palsson (2015):
#'shape = shapeR("ShapeAnalysis/","FISH.csv")
#'shape = detect.outline(shape)
#'#If otolith outline in folder IC named 403_1 needs to be removed
#'shape = remove.outline(shape, "IC", "403_1")}
#' @rdname remove.outline
#' @export remove.outline
#' @author Lisa Anne Libungan
#' @references Libungan LA and Palsson S (2015) ShapeR: An R Package to Study Otolith Shape Variation among Fish Populations. PLoS ONE 10(3): e0121102. \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121102}
remove.outline <- 
  
  function(object,folder="",fname=""){
    
    if(fname != ""){
      ind = which(names(object@outline.list[[folder]])==fname)
      if(length(ind)>0)
      {
        object@outline.list[[folder]] = object@outline.list[[folder]][-ind]
        message("The outline was removed")
      }
      else
      {
        message("Did not find the outline. It was not removed.")
      }
    }
    return(object)
  }


#' An example shapeR instance including 160 images.
#' 
#' The shape coefficients have been generated.
#' The wavelet coefficients have been standardized using pop and length_cm.   
#'
#' The class slot's are as follows:
#' \itemize{
#'   \item project.path. A path as "ShapeAnalysis/"
#'   \item info.file. A file as FISH.csv. The information is stored in the data frame master.list
#'   \item outline.list. A list with three elements (IC, NO, SC) which give a list of the otolith outlines
#'   \item filter. A logical vector showing which elements of the master list have valid otoliths
#'   \item fourier.coef. A matrix of the Normalized Elliptic Fourier coefficients
#'   \item wavelet.coef. A matrix of the wavelet coefficients
#'   \item shape. A matrix of shape variables after scaling according to calibration otolith.area, otolith.length, otolith.width, otolith.perimeter.
#'   \item fourier.coef.std. A matrix which will contain standardized Fourier coefficients
#'   \item wavelet.coef.std. A matrix which will contain standardized wavelet coefficients
#'   \item shape.coef.raw. A matrix of shape variables before scaling according to calibration otolith.area, otolith.length, otolith.width, otolith.perimeter.
#'   \item master.list. The contents of the info.file
#' }
#'
#' @docType data
#' @keywords datasets
#' @source \url{https://github.com/lisalibungan/shapeR}
#' @name shape
#' @usage data(shape)
#' @format A shapeR class including 160 images
NULL

#' An example data file
#' 
#' The file's columns are:
#' \itemize{
#'   \item country
#'   \item station
#'   \item pop
#'   \item stockID
#'   \item day
#'   \item month
#'   \item year
#'   \item lat
#'   \item lon
#'   \item fishno
#'   \item length_cm
#'   \item weight_g
#'   \item age
#'   \item sex
#'   \item maturity
#'   \item folder
#'   \item picname
#'   \item cal
#' }
#'
#' @docType data
#' @keywords datasets
#' @usage data(FISH)
#' @source \url{https://github.com/lisalibungan/shapeR}
#' @name FISH
#' @format An example data file
NULL