reportHTML.repeatedTrain <- function(object,
                                     filename,
                                     extension="html",
                                     directory=getwd(),
                                     Title,
                                     binSize=0.025,
                                     method=c("gsspsth0","gsspsth","gampsth"),
                                     stimTimeCourse=NULL,
                                     colCI=2,
                                     doTimeTransformation=TRUE,
                                     k=100,
                                     bs="tp",
                                     doGamCheck=FALSE,
                                     ...)
### reportHTML method for repeatedTrain objects.
### A raster plot is built first with a superposed
### gampsth.
### A summary of the gam object obtained by fitting
### an inhomogenous Poisson model to the data is printed.
### A summary of the relative change of the estimated
### intensity over each of the bins used is printed.
### The  smoothed psth is plotted with 95% bands.
### If doGamCheck is TRUE a gam.check plot is built.
### if doTimeTransformation is TRUE the time transformation
### is performed using the inhomogenous Poisson model and
### an Ogata's test plot is generated.
{
  
  objectN <- deparse(substitute(object))
  nbTrials <- length(object)
  
  ## check is object is a repeatedTrain object
  if (!is.repeatedTrain(object)) object <- as.repeatedTrain(object)

  if (missing(filename))
    filename <- paste(objectN,"analysis")
  
  if (missing(Title))
    Title <- filename

  
  HTMLInitFile(outdir=directory,
               filename=filename,
               extension=extension,
               Title=Title)

  fullName <- paste(directory,"/",
                    filename,".",
                    extension,sep="")

  saveName <- paste(directory,"/",
                    filename,".rda",
                    sep="")
  
  HTML.title(filename,
             file=fullName,
             HR=2)

  ## get the gampsth of the train
  PoissonF <- switch(method[1],
                     gampsth = gampsth(object,
                       binSize=binSize,
                       k=k,
                       bs=bs,
                       plot=FALSE,
                       ...),
                     gsspsth = gsspsth(object,
                       binSize=binSize,
                       plot=FALSE,
                       ...),
                     gsspsth0 = gsspsth0(object,
                       binSize=binSize,
                       plot=FALSE,
                       ...)
                     )

  PoissonFS <- summary(PoissonF)

  ## add a spike train plot after time transformation
  HTML.title(paste("Raster plot of ",
                   filename,
                   " with superposed PSTH",sep=""),
             file=fullName,HR=3)

  
  rpFigName <- paste(filename,"_rp.png",sep="")
  figFname <- paste(directory,"/",rpFigName,sep="")
  png(figFname,width=800,height=800)
  plot(object,
       stimTimeCourse=stimTimeCourse,
       main="")
  lines(PoissonF$mids,
        PoissonF$freq*nbTrials/max(PoissonF$freq),
        col=2,
        lwd=2)
  dev.off()
  HTMLInsertGraph(rpFigName,
                  file=fullName,
                  WidthHTML=800,
                  HeightHTML=800)


  HTMLbr(1,file=fullName)
  HTMLhr(file=fullName)

  HTML.title(paste("Fit summary of ",filename,sep=""),
             file=fullName,HR=3)
  txtA <- paste("<p>Fit performed after binning the \"collapsed train\" with",
                " a bin size of: ",
                binSize,
                ", assuming ",
                "no trial effect and an inhomogenous Poisson model. ",
                sep="")
  txtB <- switch(method[1],
                 gampsth=paste("Fit with gampsth and:k=",
                   k," and bs=\"",bs,"\".</p>",sep=""),
                 gsspsth=paste("Fit with gsspsth\".</p>",sep=""),
                 gsspsth0=paste("Fit with gsspsth0\".</p>",sep="")
                 )
  cat(paste(txtA,txtB,sep=""),
      file=fullName,
      append=TRUE
      )
  HTML(PoissonFS,file=fullName)

  HTMLbr(1,file=fullName)
  rcS <- summary(100*abs(diff(PoissonF$freq)/PoissonF$freq[-1]))
  cat(paste("<p>The statistics of the relative change <b>x 100</b> of the predicted rate",
            " over a bin of size: ",
            binSize,
            " are:</p>",
            sep=""),
      file=fullName,
      append=TRUE
      )
  HTML(rcS,file=fullName)
  
  HTMLbr(1,file=fullName)
  HTML.title(paste("Predicted firing rate under the inhomogenous ",
                   "Poisson hypothesis for ",
                   filename,sep=""),
             file=fullName,HR=3)
  
  frFigName <- paste(filename,"_fr.png",sep="")
  figFname <- paste(directory,"/",frFigName,sep="")
  png(figFname,width=800,height=800)
  plot(PoissonF,
       ylab="Instantaneous firing rate (Hz)",
       colCI=colCI,
       stimTimeCourse=stimTimeCourse,
       main="")
  dev.off()
  HTMLInsertGraph(frFigName,
                  file=fullName,
                  WidthHTML=800,
                  HeightHTML=800)

  if (method[1] == "gampsth" && doGamCheck) {
    ## Do gam check
    HTMLbr(1,file=fullName)
    HTMLhr(file=fullName)
    
    HTML.title(paste("GAM goodness of fit diagnostics of ",
                     filename,sep=""),
               file=fullName,HR=3)

    gcFigName <- paste(filename,"_gc.png",sep="")
    figFname <- paste(directory,"/",gcFigName,sep="")
    png(figFname,width=800,height=800)
    gam.check(gamObj(PoissonF))
    dev.off()
    HTMLInsertGraph(gcFigName,
                    file=fullName,
                    WidthHTML=800,
                    HeightHTML=800)
  } ## End of conditional on method[1] == "gampsth" && doGamCheck

  
  if (doTimeTransformation) {
    ## Do time transformation under the Poisson model
    HTMLbr(1,file=fullName)
    HTMLhr(file=fullName)
    
    HTML.title(paste("Ogata's tests on the time transformed spike trains of ",
                     filename,sep=""),
               file=fullName,HR=3)
  
    Lambda <- lapply(object,
                     function(l) {
                       class(l) <- NULL
                       L <- PoissonF$LambdaFct(l)
                       diff(L)
                     }
                     )
  
    Lambda <- cumsum(unlist(Lambda))
    
    class(Lambda) <- c("transformedTrain","spikeTrain")
    otFigName <- paste(filename,"_TTot.png",sep="")
    figFname <- paste(directory,"/",otFigName,sep="")
    png(figFname,width=800,height=800)
    plot(Lambda,
         which=c(1,2,4,5),
         ask=FALSE)
    dev.off()
    HTMLInsertGraph(otFigName,
                    file=fullName,
                    WidthHTML=800,
                    HeightHTML=800)
  } ## End of conditional on doTimeTransformation
  
  HTMLEndFile()

  fctCall <- match.call()
  
  if (doTimeTransformation) {
    save(PoissonF,Lambda,fctCall,file=saveName)
  } else {
    save(PoissonF,fctCall,file=saveName)
  }
  
}
