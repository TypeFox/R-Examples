reportHTML.gam <- function(object,
                           filename,
                           extension="html",
                           directory=getwd(),
                           Title,
                           neuron,
                           neuronEvts,
                           ...) {
  
  objectN <- deparse(substitute(object))
  
  ## if neuronEvts given check that it is a named list
  if (!missing(neuronEvts)) {
    if (!inherits(neuronEvts,"list")) neuronEvts <- list(neuronEvts)
    if (is.null(names(neuronEvts)))
      stop("neuronEvts should be a named list.")
  } ## End of conditional on !missing(neuronEvts)

  ## if neuron is missing try to figure out its value from the
  ## neuron variable of the "data" component of object
  if (missing(neuron)) {
    if (!("neuron" %in% names(object$data)))
      stop("neuron is not specified and cannot be obtained from object.")
    neuron <- with(object$data,unique(neuron)[1])
  } ## End of conditional on missing(neuron)
  
  if (missing(filename))
    filename <- paste(objectN,"GAM analysis")
  
  if (missing(Title))
    Title <- filename

  
  HTMLInitFile(outdir=directory,
               filename=filename,
               extension=extension,
               Title=Title)

  fullName <- paste(directory,"/",
                    filename,".",
                    extension,sep="")

  ## Write object's summary
  HTML.title("GAM fit summary:",file=fullName,HR=3)
  objectS <- summary(object) 
  HTML(objectS,file=fullName)

  
  HTMLhr(file=fullName)

  ## add a spike train plot after time transformation
  HTML.title("Spike train after time transformation:",
             file=fullName,HR=3)

  Lambda <- transformedTrain(object)
  st <- Lambda
  class(st) <- "spikeTrain"

  stFigName <- paste(filename,"_TTst.png",sep="")
  figFname <- paste(directory,"/",stFigName,sep="")
  png(figFname,width=500,height=500)
  plot(st,
       xlab="Transformed time",
       main=""
       )
  dev.off()
  HTMLInsertGraph(stFigName,
                  file=fullName,
                  WidthHTML=500,
                  HeightHTML=500)

  ## add a renewal test plot after time transformation
  HTMLhr(file=fullName)
  HTML.title("Renewal test after time transformation:",
             file=fullName,HR=3)
  
  rtFigName <- paste(filename,"_TTrt.png",sep="")
  figFname <- paste(directory,"/",rtFigName,sep="")
  png(figFname,width=800,height=800)
  renewalTestPlot(Lambda)
  dev.off()
  HTMLInsertGraph(rtFigName,
                  file=fullName,
                  WidthHTML=800,
                  HeightHTML=800)

  ## add a Ogata tests after time transformation
  HTMLhr(file=fullName)
  HTML.title("Ogata's tests after time transformation:",
             file=fullName,HR=3)
  
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

  ## Check if the Ogata tests are passed and if yes and if
  ## neurons is given add interaction plots
  LambdaS <- summary(Lambda)
  otTRUE <- LambdaS[[1]][2] &&
  LambdaS[[2]][2] &&
  (0.005 <= LambdaS[[3]][2]) &&
  (0.995 >= LambdaS[[3]][2])
  
  if (otTRUE && !missing(neuronEvts)) {
    HTMLhr(file=fullName)
    HTML.title("Interaction tests after time transformation:",
               file=fullName,HR=3)

    ## load other neurons data
    pre.Lambda <- lapply(neuronEvts,
                         function(l) transformedTrain(object,l)
                         )
    names(pre.Lambda) <- names(neuronEvts)
    
    itFigName <- paste(filename,"_TTit.png",sep="")
    figFname <- paste(directory,"/",itFigName,sep="")
    png(figFname,width=300*length(neuronEvts),height=1000)
    layout(matrix(1:(2*length(neuronEvts)),
                  nrow=2)
           )

    sapply(seq(neuronEvts),
           function(nIdx) {
             theFRT <- pre.Lambda[[nIdx]] %frt% Lambda
             plot(theFRT,
                  main=paste(names(neuronEvts)[nIdx],"/ N",neuron),
                  ask=FALSE,which=1)
             plot(theFRT,
                  main=paste(names(neuronEvts)[nIdx],"/ N",neuron),
                  ask=FALSE,which=2)
           }
           )
    dev.off()
    HTMLInsertGraph(itFigName,
                    file=fullName,
                    WidthHTML=300*length(neuronEvts),
                    HeightHTML=1000)
  } ## end of conditional on otTRUE && !missing(neuronEvts)

  HTMLbr(1,file=fullName)
  HTMLhr(file=fullName)
  
  ## End with the call used to generate the object
  HTML.title(paste("Call used to generate object ",
                   objectN, ":",sep=""),
             file=fullName,HR=3)
  thisC <- object$call
  HTML(thisC,file=fullName)
  HTMLbr(1,file=fullName)
  
  HTMLEndFile()

}
