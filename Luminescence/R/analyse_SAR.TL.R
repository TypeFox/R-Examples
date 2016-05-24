#' Analyse SAR TL measurements
#'
#' The function performs a SAR TL analysis on a
#' \code{\linkS4class{RLum.Analysis}} object including growth curve fitting.
#'
#' This function performs a SAR TL analysis on a set of curves. The SAR
#' procedure in general is given by Murray and Wintle (2000). For the
#' calculation of the Lx/Tx value the function \link{calc_TLLxTxRatio} is
#' used.\cr\cr \bold{Provided rejection criteria}\cr\cr
#' \sQuote{recyling.ratio}: calculated for every repeated regeneration dose
#' point.\cr \sQuote{recuperation.rate}: recuperation rate calculated by
#' comparing the Lx/Tx values of the zero regeneration point with the Ln/Tn
#' value (the Lx/Tx ratio of the natural signal).  For methodological
#' background see Aitken and Smith (1988)\cr
#'
#' @param object \code{\linkS4class{RLum.Analysis}}(\bold{required}): input
#' object containing data for analysis
#'
#' @param object.background currently not used
#'
#' @param signal.integral.min \link{integer} (\bold{required}): requires the
#' channel number for the lower signal integral bound (e.g.
#' \code{signal.integral.min = 100})
#'
#' @param signal.integral.max \link{integer} (\bold{required}): requires the
#' channel number for the upper signal integral bound (e.g.
#' \code{signal.integral.max = 200})
#'
#' @param sequence.structure \link{vector} \link{character} (with default):
#' specifies the general sequence structure. Three steps are allowed (
#' \code{"PREHEAT"}, \code{"SIGNAL"}, \code{"BACKGROUND"}), in addition a
#' parameter \code{"EXCLUDE"}. This allows excluding TL curves which are not
#' relevant for the protocol analysis.  (Note: None TL are removed by default)
#'
#' @param rejection.criteria \link{list} (with default): list containing
#' rejection criteria in percentage for the calculation.
#'
#' @param dose.points \code{\link{numeric}} (optional): option set dose points manually
#'
#' @param log \link{character} (with default): a character string which
#' contains "x" if the x axis is to be logarithmic, "y" if the y axis is to be
#' logarithmic and "xy" or "yx" if both axes are to be logarithmic. See
#' \link{plot.default}).
#'
#' @param \dots further arguments that will be passed to the function
#' \code{\link{plot_GrowthCurve}}
#'
#' @return A plot (optional) and an \code{\linkS4class{RLum.Results}} object is
#' returned containing the following elements:
#' \item{De.values}{\link{data.frame} containing De-values and further
#' parameters} \item{LnLxTnTx.values}{\link{data.frame} of all calculated Lx/Tx
#' values including signal, background counts and the dose points.}
#' \item{rejection.criteria}{\link{data.frame} with values that might by used
#' as rejection criteria. NA is produced if no R0 dose point exists.}\cr\cr
#' \bold{note:} the output should be accessed using the function
#' \code{\link{get_RLum}}
#' @note \bold{THIS IS A BETA VERSION}\cr\cr None TL curves will be removed
#' from the input object without further warning.
#' @section Function version: 0.1.4
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @seealso \code{\link{calc_TLLxTxRatio}}, \code{\link{plot_GrowthCurve}},
#' \code{\linkS4class{RLum.Analysis}}, \code{\linkS4class{RLum.Results}}
#' \code{\link{get_RLum}}
#'
#' @references Aitken, M.J. and Smith, B.W., 1988. Optical dating: recuperation
#' after bleaching.  Quaternary Science Reviews 7, 387-393.
#'
#' Murray, A.S. and Wintle, A.G., 2000. Luminescence dating of quartz using an
#' improved single-aliquot regenerative-dose protocol. Radiation Measurements
#' 32, 57-73.
#' @keywords datagen plot
#' @examples
#'
#'
#' ##load data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##transform the values from the first position in a RLum.Analysis object
#' object <- Risoe.BINfileData2RLum.Analysis(TL.SAR.Data, pos=3)
#'
#' ##perform analysis
#' analyse_SAR.TL(object,
#'                signal.integral.min = 210,
#'                signal.integral.max = 220,
#'                log = "y",
#'                fit.method = "EXP OR LIN",
#'                sequence.structure = c("SIGNAL", "BACKGROUND"))
#'
#' @export
analyse_SAR.TL <- function(
  object,
  object.background,
  signal.integral.min,
  signal.integral.max,
  sequence.structure = c("PREHEAT", "SIGNAL", "BACKGROUND"),
  rejection.criteria = list(recycling.ratio = 10, recuperation.rate = 10),
  dose.points,
  log = "",
  ...
){

  # CONFIG  -----------------------------------------------------------------

  ##set allowed curve types
  type.curves <- c("TL")

  ##=============================================================================#
  # General Integrity Checks ---------------------------------------------------

  ##GENERAL

  ##MISSING INPUT
  if(missing("object")==TRUE){
    stop("[analyse_SAR.TL] No value set for 'object'!")
  }

  if(missing("signal.integral.min") == TRUE){
    stop("[analyse_SAR.TL] No value set for 'signal.integral.min'!")
  }

  if(missing("signal.integral.max") == TRUE){
    stop("[analyse_SAR.TL] No value set for 'signal.integral.max'!")
  }

  ##INPUT OBJECTS
  if(is(object, "RLum.Analysis") == FALSE){
    stop("[analyse_SAR.TL] Input object is not of type 'RLum.Analyis'!")
  }


  # Protocol Integrity Checks --------------------------------------------------

  ##Remove non TL-curves from object by selecting TL curves
  object@records <- get_RLum(object, recordType = type.curves)

  ##ANALYSE SEQUENCE OBJECT STRUCTURE

  ##set vector for sequence structure
  temp.protocol.step <- rep(sequence.structure,length(object@records))[1:length(object@records)]

  ##grep object strucute
  temp.sequence.structure <- structure_RLum(object)

  ##set values for step
  temp.sequence.structure[,"protocol.step"] <- temp.protocol.step

  ##remove TL curves which are excluded
  temp.sequence.structure <- temp.sequence.structure[which(
    temp.sequence.structure[,"protocol.step"]!="EXCLUDE"),]

  ##check integrity; signal and bg range should be equal
  if(length(
    unique(
      temp.sequence.structure[temp.sequence.structure[,"protocol.step"]=="SIGNAL","x.max"]))>1){

    stop(paste(
      "[analyse_SAR.TL()] Signal range differs. Check sequence structure.\n",
      temp.sequence.structure
    ))
  }

  ##check if the wanted curves are a multiple of the structure
  if(length(temp.sequence.structure[,"id"])%%length(sequence.structure)!=0){

    stop("[analyse_SAR.TL()] Input TL curves are not a multiple of the sequence structure.")

  }



  # # Calculate LnLxTnTx values  --------------------------------------------------

  ##grep IDs for signal and background curves
  TL.preheat.ID <- temp.sequence.structure[
    temp.sequence.structure[,"protocol.step"] == "PREHEAT","id"]

  TL.signal.ID <- temp.sequence.structure[
    temp.sequence.structure[,"protocol.step"] == "SIGNAL","id"]

  TL.background.ID <- temp.sequence.structure[
    temp.sequence.structure[,"protocol.step"] == "BACKGROUND","id"]


  ##calculate LxTx values using external function

  for(i in seq(1,length(TL.signal.ID),by=2)){

    temp.LnLxTnTx <- get_RLum(
      calc_TLLxTxRatio(
        Lx.data.signal = get_RLum(object, record.id=TL.signal.ID[i]),
        Lx.data.background = get_RLum(object, record.id=TL.background.ID[i]),
        Tx.data.signal = get_RLum(object, record.id=TL.signal.ID[i+1]),
        Tx.data.background = get_RLum(object, record.id = TL.background.ID[i+1]),
        signal.integral.min,
        signal.integral.max))

    ##grep dose
    temp.Dose <- object@records[[TL.signal.ID[i]]]@info$IRR_TIME


    temp.LnLxTnTx <- cbind(Dose=temp.Dose, temp.LnLxTnTx)

    if(exists("LnLxTnTx")==FALSE){

      LnLxTnTx <- data.frame(temp.LnLxTnTx)

    }else{

      LnLxTnTx <- rbind(LnLxTnTx,temp.LnLxTnTx)

    }
  }

  ##set dose.points manual if argument was set
  if(!missing(dose.points)){
    temp.Dose <- dose.points
    LnLxTnTx$Dose <- dose.points

  }

  # Set regeneration points -------------------------------------------------

  #generate unique dose id - this are also the # for the generated points
  temp.DoseID <- c(0:(length(temp.Dose)-1))
  temp.DoseName <- paste("R",temp.DoseID,sep="")
  temp.DoseName <- cbind(Name=temp.DoseName,Dose=temp.Dose)

  ##set natural
  temp.DoseName[temp.DoseName[,"Name"]=="R0","Name"]<-"Natural"

  ##set R0
  temp.DoseName[temp.DoseName[,"Name"]!="Natural" & temp.DoseName[,"Dose"]==0,"Name"]<-"R0"

  ##find duplicated doses (including 0 dose - which means the Natural)
  temp.DoseDuplicated<-duplicated(temp.DoseName[,"Dose"])

  ##combine temp.DoseName
  temp.DoseName<-cbind(temp.DoseName,Repeated=temp.DoseDuplicated)

  ##correct value for R0 (it is not really repeated)
  temp.DoseName[temp.DoseName[,"Dose"]==0,"Repeated"]<-FALSE

  ##combine in the data frame
  temp.LnLxTnTx<-data.frame(Name=temp.DoseName[,"Name"],
                            Repeated=as.logical(temp.DoseName[,"Repeated"]))


  LnLxTnTx<-cbind(temp.LnLxTnTx,LnLxTnTx)
  LnLxTnTx[,"Name"]<-as.character(LnLxTnTx[,"Name"])


  # Calculate Recycling Ratio -----------------------------------------------

  ##Calculate Recycling Ratio

  if(length(LnLxTnTx[LnLxTnTx[,"Repeated"]==TRUE,"Repeated"])>0){

    ##identify repeated doses
    temp.Repeated<-LnLxTnTx[LnLxTnTx[,"Repeated"]==TRUE,c("Name","Dose","LxTx")]

    ##find concering previous dose for the repeated dose
    temp.Previous<-t(sapply(1:length(temp.Repeated[,1]),function(x){
      LnLxTnTx[LnLxTnTx[,"Dose"]==temp.Repeated[x,"Dose"] &
                 LnLxTnTx[,"Repeated"]==FALSE,c("Name","Dose","LxTx")]
    }))

    ##convert to data.frame
    temp.Previous<-as.data.frame(temp.Previous)

    ##set column names
    temp.ColNames<-sapply(1:length(temp.Repeated[,1]),function(x){
      paste(temp.Repeated[x,"Name"],"/",
            temp.Previous[temp.Previous[,"Dose"]==temp.Repeated[x,"Dose"],"Name"],
            sep="")
    })

    ##Calculate Recycling Ratio
    RecyclingRatio<-as.numeric(temp.Repeated[,"LxTx"])/as.numeric(temp.Previous[,"LxTx"])

    ##Just transform the matrix and add column names
    RecyclingRatio<-t(RecyclingRatio)
    colnames(RecyclingRatio)<-temp.ColNames

  }else{RecyclingRatio<-NA}


  # Calculate Recuperation Rate ---------------------------------------------


  ##Recuperation Rate
  if("R0" %in% LnLxTnTx[,"Name"]==TRUE){
    Recuperation<-round(LnLxTnTx[LnLxTnTx[,"Name"]=="R0","LxTx"]/
                          LnLxTnTx[LnLxTnTx[,"Name"]=="Natural","LxTx"],digits=4)
  }else{Recuperation<-NA}


  # Combine and Evaluate Rejection Criteria ---------------------------------

  RejectionCriteria <- data.frame(
    citeria = c(colnames(RecyclingRatio), "recuperation rate"),
    value = c(RecyclingRatio,Recuperation),
    threshold = c(
      rep(paste("+/-", rejection.criteria$recycling.ratio/100)
          ,length(RecyclingRatio)),
      paste("", rejection.criteria$recuperation.rate/100)
    ),
    status = c(

      if(is.na(RecyclingRatio)==FALSE){

        sapply(1:length(RecyclingRatio), function(x){
          if(abs(1-RecyclingRatio[x])>(rejection.criteria$recycling.ratio/100)){
            "FAILED"
          }else{"OK"}})}else{NA},

      if(is.na(Recuperation)==FALSE &
           Recuperation>rejection.criteria$recuperation.rate){"FAILED"}else{"OK"}

    ))

  ##============================================================================##
  ##PLOTTING
  ##============================================================================##

  # Plotting - Config -------------------------------------------------------

  ##grep plot parameter
  par.default <- par(no.readonly = TRUE)

  ##colours and double for plotting
  col <- get("col", pos = .LuminescenceEnv)

  col.doubled <- rep(col, each=2)

  layout(matrix(c(1,1,2,2,
                  1,1,2,2,
                  3,3,4,4,
                  3,3,4,4,
                  5,5,5,5),5,4,byrow=TRUE))

  par(oma=c(0,0,0,0), mar=c(4,4,3,3))

  ## 1 -> TL Lx
  ## 2 -> TL Tx
  ## 3 -> TL Lx Plateau
  ## 4 -> TL Tx Plateau
  ## 5 -> Legend

  ##recalculate signal.integral from channels to temperature
  signal.integral.temperature <- c(object@records[[TL.signal.ID[1]]]@data[signal.integral.min,1] :
                                     object@records[[TL.signal.ID[1]]]@data[signal.integral.max,1])


  ##warning if number of curves exceed colour values
  if(length(col)<length(TL.signal.ID/2)){
    cat("\n[analyse_SAR.TL.R] Warning: To many curves! Only the first",
        length(col),"curves are plotted!")
  }


  # # Plotting TL Lx Curves ----------------------------------------------------

  #open plot area LnLx
  plot(NA,NA,
       xlab="Temp. [\u00B0C]",
       ylab=paste("TL [a.u.]",sep=""),
       xlim=c(0.1,
              max(temp.sequence.structure[temp.sequence.structure[,"protocol.step"]=="SIGNAL","x.max"])),
       ylim=c(
         min(temp.sequence.structure[temp.sequence.structure[,"protocol.step"]=="SIGNAL","y.min"]),
         max(temp.sequence.structure[temp.sequence.structure[,"protocol.step"]=="SIGNAL","y.max"])),

       main=expression(paste(L[n],",",L[x]," curves",sep="")),
       log=log)


  ##plot curves
  sapply(seq(1,length(TL.signal.ID),by=2), function(x){


    lines(object@records[[TL.signal.ID[x]]]@data,col=col.doubled[x])

  })

  ##mark integration limits
  abline(v=min(signal.integral.temperature), lty=2, col="gray")
  abline(v=max(signal.integral.temperature), lty=2, col="gray")


  # Plotting TnTx Curves ----------------------------------------------------

  #open plot area TnTx
  plot(NA,NA,
       xlab="Temp. [\u00B0C]",
       ylab=paste("TL [a.u.]",sep=""),
       xlim=c(0.1,
              max(temp.sequence.structure[temp.sequence.structure[,"protocol.step"]=="SIGNAL","x.max"])),
       ylim=c(
         min(temp.sequence.structure[temp.sequence.structure[,"protocol.step"]=="SIGNAL","y.min"]),
         max(temp.sequence.structure[temp.sequence.structure[,"protocol.step"]=="SIGNAL","y.max"])),

       main=expression(paste(T[n],",",T[x]," curves",sep="")),
       log=log)


  ##plot curves
  sapply(seq(2,length(TL.signal.ID),by=2), function(x){


    lines(object@records[[TL.signal.ID[x]]]@data,col=col.doubled[x])

  })

  ##mark integration limits
  abline(v=min(signal.integral.temperature), lty=2, col="gray")
  abline(v=max(signal.integral.temperature), lty=2, col="gray")


  # Plotting Plateau Test LnLx -------------------------------------------------

  NTL.net.LnLx <- data.frame(object@records[[TL.signal.ID[1]]]@data[,1],
                             object@records[[TL.signal.ID[1]]]@data[,2]-
                               object@records[[TL.background.ID[1]]]@data[,2])

  Reg1.net.LnLx <- data.frame(object@records[[TL.signal.ID[3]]]@data[,1],
                              object@records[[TL.signal.ID[3]]]@data[,2]-
                                object@records[[TL.background.ID[3]]]@data[,2])


  TL.Plateau.LnLx <- data.frame(NTL.net.LnLx[,1], Reg1.net.LnLx[,2]/NTL.net.LnLx[,2])

  ##Plot Plateau Test
  plot(NA, NA,
       xlab = "Temp. [\u00B0C]",
       ylab = "TL [a.u.]",
       xlim = c(min(signal.integral.temperature)*0.9, max(signal.integral.temperature)*1.1),
       ylim = c(0, max(NTL.net.LnLx[,2])),
       main = expression(paste("Plateau test ",L[n],",",L[x]," curves",sep=""))
  )


  ##plot single curves
  lines(NTL.net.LnLx, col=col[1])
  lines(Reg1.net.LnLx, col=col[2])


  ##plot
  par(new=TRUE)
  plot(TL.Plateau.LnLx,
       axes=FALSE,
       xlab="",
       ylab="",
       ylim=c(0,
              quantile(TL.Plateau.LnLx[c(signal.integral.min:signal.integral.max),2],
                       probs = c(0.90), na.rm = TRUE)+3),
       col="darkgreen")
  axis(4)


  # Plotting Plateau Test TnTx -------------------------------------------------

  ##get NTL signal
  NTL.net.TnTx <- data.frame(object@records[[TL.signal.ID[2]]]@data[,1],
                             object@records[[TL.signal.ID[2]]]@data[,2]-
                               object@records[[TL.background.ID[2]]]@data[,2])

  ##get signal from the first regeneration point
  Reg1.net.TnTx <- data.frame(object@records[[TL.signal.ID[4]]]@data[,1],
                              object@records[[TL.signal.ID[4]]]@data[,2]-
                                object@records[[TL.background.ID[4]]]@data[,2])


  ##combine values
  TL.Plateau.TnTx <- data.frame(NTL.net.TnTx[,1], Reg1.net.TnTx[,2]/NTL.net.TnTx[,2])

  ##Plot Plateau Test
  plot(NA, NA,
       xlab = "Temp. [\u00B0C]",
       ylab = "TL [a.u.]",
       xlim = c(min(signal.integral.temperature)*0.9, max(signal.integral.temperature)*1.1),
       ylim = c(0, max(NTL.net.TnTx[,2])),
       main = expression(paste("plateau Test ",T[n],",",T[x]," curves",sep=""))
  )


  ##plot single curves
  lines(NTL.net.TnTx, col=col[1])
  lines(Reg1.net.TnTx, col=col[2])


  ##plot
  par(new=TRUE)
  plot(TL.Plateau.TnTx,
       axes=FALSE,
       xlab="",
       ylab="",
       ylim=c(0,
              quantile(TL.Plateau.TnTx[c(signal.integral.min:signal.integral.max),2],
                       probs = c(0.90), na.rm = TRUE)+3),
       col="darkgreen")
  axis(4)




  # Plotting Legend ----------------------------------------


  plot(c(1:(length(TL.signal.ID)/2)),
       rep(8,length(TL.signal.ID)/2),
       type = "p",
       axes=FALSE,
       xlab="",
       ylab="",
       pch=15,
       col=col[1:length(TL.signal.ID)],
       cex=2,
       ylim=c(0,10)
  )

  ##add text
  text(c(1:(length(TL.signal.ID)/2)),
       rep(4,length(TL.signal.ID)/2),
       paste(LnLxTnTx$Name,"\n(",LnLxTnTx$Dose,")", sep="")

  )

  ##add line
  abline(h=10,lwd=0.5)

  ##set failed text and mark De as failed
  if(length(grep("FAILED",RejectionCriteria$status))>0){

    mtext("[FAILED]", col="red")


  }

  ##reset par
  par(par.default)
  rm(par.default)

  # Plotting  GC  ----------------------------------------
  temp.sample <- data.frame(Dose=LnLxTnTx$Dose,
                            LxTx=LnLxTnTx$LxTx,
                            LxTx.Error=LnLxTnTx$LxTx*0.1,
                            TnTx=LnLxTnTx$TnTx
  )

  temp.GC <- get_RLum(plot_GrowthCurve(temp.sample,
                                               ...))[,c("De","De.Error")]

  ##add recjection status
  if(length(grep("FAILED",RejectionCriteria$status))>0){

    temp.GC <- data.frame(temp.GC, RC.Status="FAILED")

  }else{

    temp.GC <- data.frame(temp.GC, RC.Status="OK")

  }

  # Return Values -----------------------------------------------------------

  newRLumResults.analyse_SAR.TL <- set_RLum(
    class = "RLum.Results",
    data = list(
      De.values = temp.GC,
      LnLxTnTx.table = LnLxTnTx,
      rejection.criteria = RejectionCriteria))


  return(newRLumResults.analyse_SAR.TL)

}
