#' Analyse SAR CW-OSL measurements.
#'
#' The function analyses SAR CW-OSL curve data and provides a summary of the
#' measured data for every position. The output of the function is optimised
#' for SAR OSL measurements on quartz.
#'
#' The function works only for standard SAR protocol measurements introduced by
#' Murray and Wintle (2000) with CW-OSL curves. For the calculation of the
#' Lx/Tx value the function \link{calc_OSLLxTxRatio} is used. \cr\cr
#'
#' \bold{Provided rejection criteria}\cr\cr \sQuote{recyling ratio}: calculated
#' for every repeated regeneration dose point.\cr \sQuote{recuperation}:
#' recuperation rate calculated by comparing the Lx/Tx values of the zero
#' regeneration point with the Ln/Tn value (the Lx/Tx ratio of the natural
#' signal). For methodological background see Aitken and Smith (1988)\cr
#'
#' \sQuote{IRSL/BOSL}: the integrated counts (\code{signal.integral}) of an
#' IRSL curve are compared to the integrated counts of the first regenerated
#' dose point. It is assumed that IRSL curves got the same dose as the first
#' regenerated dose point. \strong{Note:} This is not the IR depletation ratio
#' described by Duller (2003).
#'
#' @param input.data \link{Risoe.BINfileData-class} (\bold{required}): input
#' data from a Risoe BIN file, produced by the function \link{read_BIN2R}.
#'
#' @param signal.integral \link{vector} (\bold{required}): channels used for
#' the signal integral, e.g. \code{signal.integral=c(1:2)}
#'
#' @param background.integral \link{vector} (\bold{required}): channels used
#' for the background integral, e.g. \code{background.integral=c(85:100)}
#'
#' @param position \link{vector} (optional): reader positions that want to be
#' analysed (e.g. \code{position=c(1:48)}. Empty positions are automatically
#' omitted. If no value is given all positions are analysed by default.
#'
#' @param run \link{vector} (optional): range of runs used for the analysis. If
#' no value is given the range of the runs in the sequence is deduced from the
#' Risoe.BINfileData object.
#'
#' @param set \link{vector} (optional): range of sets used for the analysis. If
#' no value is given the range of the sets in the sequence is deduced from the
#' \code{Risoe.BINfileData} object.
#'
#' @param dtype \code{\link{character}} (optional): allows to further limit the
#' curves by their data type (\code{DTYPE}), e.g., \code{dtype = c("Natural",
#' "Dose")} limits the curves to this two data types. By default all values are
#' allowed. See \link{Risoe.BINfileData-class} for allowed data types.
#'
#' @param keep.SEL \code{\link{logical}} (default): option allowing to use the
#' \code{SEL} element of the \link{Risoe.BINfileData-class} manually. NOTE: In
#' this case any limitation provided by \code{run}, \code{set} and \code{dtype}
#' are ignored!
#'
#' @param info.measurement \link{character} (with default): option to provide
#' information about the measurement on the plot output (e.g. name of the BIN
#' or BINX file).
#'
#' @param output.plot \link{logical} (with default): plot output
#' (\code{TRUE/FALSE})
#'
#' @param output.plot.single \link{logical} (with default): single plot output
#' (\code{TRUE/FALSE}) to allow for plotting the results in single plot
#' windows. Requires \code{output.plot = TRUE}.
#'
#' @param cex.global \link{numeric} (with default): global scaling factor.
#'
#' @param \dots further arguments that will be passed to the function
#' \code{\link{calc_OSLLxTxRatio}} (supported: \code{background.count.distribution} and \code{sigmab})
#' and can be used to adjust the plot. Supported" \code{mtext}, \code{log}
#'
#' @return A plot (optional) and \link{list} is returned containing the
#' following elements: \item{LnLxTnTx}{\link{data.frame} of all calculated
#' Lx/Tx values including signal, background counts and the dose points.}
#' \item{RejectionCriteria}{\link{data.frame} with values that might by used as
#' rejection criteria. NA is produced if no R0 dose point exists.}
#' \item{SARParameters}{\link{data.frame} of additional measurement parameters
#' obtained from the BIN file, e.g. preheat or read temperature (not valid for
#' all types of measurements).}
#'
#'
#' @note Rejection criteria are calculated but not considered during the
#' analysis to discard values.\cr\cr
#'
#' \bold{The analysis of IRSL data is not directly supported}. You may want to
#' consider using the functions \code{\link{analyse_SAR.CWOSL}} or
#' \code{\link{analyse_pIRIRSequence}} instead.\cr
#'
#' \bold{The development of this function will not be continued. We recommend
#' to use the function \link{analyse_SAR.CWOSL} or instead.}
#'
#'
#' @section Function version: 0.2.17
#'
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France), Margret C. Fuchs, HZDR, Freiberg (Germany)
#' @seealso \link{calc_OSLLxTxRatio}, \link{Risoe.BINfileData-class},
#' \link{read_BIN2R}
#'
#' and for further analysis \link{plot_GrowthCurve}
#'
#' @references Aitken, M.J. and Smith, B.W., 1988. Optical dating: recuperation
#' after bleaching. Quaternary Science Reviews 7, 387-393.
#'
#' Duller, G., 2003. Distinguishing quartz and feldspar in single grain
#' luminescence measurements. Radiation Measurements, 37 (2), 161-165.
#'
#' Murray, A.S. and Wintle, A.G., 2000. Luminescence dating of quartz using an
#' improved single-aliquot regenerative-dose protocol. Radiation Measurements
#' 32, 57-73.
#' @keywords datagen dplot
#'
#' @examples
#'
#'
#' ##load data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##analyse data
#' output <- Analyse_SAR.OSLdata(input.data = CWOSL.SAR.Data,
#'                               signal.integral = c(1:5),
#'                               background.integral = c(900:1000),
#'                               position = c(1:1),
#'                               output.plot = TRUE)
#'
#' ##combine results relevant for further analysis
#' output.SAR <- data.frame(Dose = output$LnLxTnTx[[1]]$Dose,
#'                          LxTx = output$LnLxTnTx[[1]]$LxTx,
#'                          LxTx.Error = output$LnLxTnTx[[1]]$LxTx.Error)
#' output.SAR
#'
#' @export
Analyse_SAR.OSLdata <- function(
  input.data,
  signal.integral,
  background.integral,
  position,
  run,
  set,
  dtype,
  keep.SEL = FALSE,
  info.measurement = "unkown measurement",
  output.plot = FALSE,
  output.plot.single = FALSE,
  cex.global = 1,
  ...
){

  ##============================================================================##
  ##CONFIG
  ##============================================================================##

  ##set colors gallery to provide more colors

  col <- get("col", pos = .LuminescenceEnv)

  ##============================================================================##
  ##ERROR HANDLING
  ##============================================================================##

  if(missing(input.data)==TRUE){stop("[Analyse_SAR.OSLdata] No input data given!")
  }else{sample.data<-input.data}

  if(missing(signal.integral)==TRUE){stop("[Analyse_SAR.OSLdata] No signal integral is given!")}
  if(missing(background.integral)==TRUE){stop("[Analyse_SAR.OSLdata] No background integral is given!")}

  ##set values for run and set if they are not defined by the user
  if(missing(position)==TRUE){position<-min(sample.data@METADATA[,"POSITION"]):max(sample.data@METADATA[,"POSITION"])}

  if(missing(run)==TRUE){run<-min(sample.data@METADATA[,"RUN"]):max(sample.data@METADATA[,"RUN"])}

  if(missing(set)==TRUE){set<-min(sample.data@METADATA[,"SET"]):max(sample.data@METADATA[,"SET"])}

  if(missing(dtype)){dtype <- c("Natural",
                                "N+dose",
                                "Bleach",
                                "Bleach+dose",
                                "Natural (Bleach)",
                                "N+dose (Bleach)",
                                "Dose",
                                "Background")}


  # Deal with extra arguments ----------------------------------------------------

  ##deal with addition arguments
  extraArgs <- list(...)

  background.count.distribution <-
    if ("background.count.distribution" %in% names(extraArgs)) {
      extraArgs$background.count.distribution
    } else
    {
      "non-poisson"
    }

  sigmab <- if("sigmab" %in% names(extraArgs)) {extraArgs$sigmab} else
  {NULL}

  ##============================================================================##
  ##CALCULATIONS
  ##============================================================================##


  ##loop over all positions
  for (i in position){

    ##checking if position is valid
    if(length(which(sample.data@METADATA["POSITION"]==i))>0){

      ##check if OSL curves are part of the data set
      if(nrow(sample.data@METADATA[sample.data@METADATA[,"LTYPE"]=="OSL",]) == 0){

        stop("[Analyse_SAR.OSLdata()] No 'OSL' curves found!")

      }

      if(!keep.SEL){
        ##select all OSL data depending on the run and set
        sample.data@METADATA[,"SEL"]<-FALSE
        sample.data@METADATA[sample.data@METADATA[,"LTYPE"]=="OSL" &
                               sample.data@METADATA[,"RUN"]%in%run==TRUE &
                               sample.data@METADATA[,"SET"]%in%set==TRUE &
                               sample.data@METADATA[,"DTYPE"]%in%dtype==TRUE, "SEL"] <- TRUE
      }

      ##grep all OSL curve IDs
      OSL.curveID<-sample.data@METADATA[sample.data@METADATA["SEL"]==TRUE &
                                          sample.data@METADATA["POSITION"]==i,"ID"]

      ##estimate LnLx.curveID and TnTx.curveID from records
      LnLx.curveID<-OSL.curveID[seq(1,length(OSL.curveID),by=2)]
      TnTx.curveID<-OSL.curveID[seq(2,length(OSL.curveID),by=2)]


      ##Provide Values For Growth Curve Fitting

      ##(1) get dose information
      Dose<-sapply(1:length(LnLx.curveID),function(x){
        Dose<-sample.data@METADATA[sample.data@METADATA["ID"]==LnLx.curveID[x],"IRR_TIME"]
      })

      ##(2) set LxTx curves
      LnLxTnTx.curves<-(sapply(1:length(LnLx.curveID),function(x){

        ##produce data.frames for Lx/Tx calculations
        Lx.HIGH<-sample.data@METADATA[sample.data@METADATA[,"ID"]==LnLx.curveID[x],"HIGH"]
        Lx.NPOINTS<-sample.data@METADATA[sample.data@METADATA[,"ID"]==LnLx.curveID[x],"NPOINTS"]
        Tx.HIGH<-sample.data@METADATA[sample.data@METADATA[,"ID"]==TnTx.curveID[x],"HIGH"]
        Tx.NPOINTS<-sample.data@METADATA[sample.data@METADATA[,"ID"]==TnTx.curveID[x],"NPOINTS"]

        Lx.curve<-data.frame(x=seq(Lx.HIGH/Lx.NPOINTS,Lx.HIGH,by=Lx.HIGH/Lx.NPOINTS),
                             y=unlist(sample.data@DATA[LnLx.curveID[x]]))
        Tx.curve<-data.frame(x=seq(Tx.HIGH/Tx.NPOINTS,Tx.HIGH,by=Tx.HIGH/Tx.NPOINTS),
                             y=unlist(sample.data@DATA[TnTx.curveID[x]]))

        return(list(Lx.curve,Tx.curve))
      }))

      ##(3) calculate Lx/Tx ratio
      LnLxTnTx <- get_RLum(
        merge_RLum(lapply(1:length(LnLxTnTx.curves[1, ]), function(k) {
          calc_OSLLxTxRatio(
            Lx.data = as.data.frame(LnLxTnTx.curves[1, k]),
            Tx.data = as.data.frame(LnLxTnTx.curves[2, k]),
            signal.integral = signal.integral,
            background.integral = background.integral,
            background.count.distribution = background.count.distribution,
            sigmab = sigmab
          )
        })))


      ##finally combine to data.frame including the record ID for further analysis
      LnLxTnTx <- cbind(LnLxTnTx,LnLx.curveID,TnTx.curveID)

      ##(4.1) set info concerning the kind of regeneration points

      ##generate unique dose id - this are also the # for the generated points
      temp.DoseID<-c(0:(length(Dose)-1))
      temp.DoseName<-paste("R",temp.DoseID,sep="")
      temp.DoseName<-cbind(Name=temp.DoseName,Dose)

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


      ##(5) Combine all values in a data.frame
      temp.LnLxTnTx<-data.frame(Name=temp.DoseName[,"Name"],
                                Dose=Dose,
                                Repeated=as.logical(temp.DoseName[,"Repeated"]))
      LnLxTnTx<-cbind(temp.LnLxTnTx,LnLxTnTx)
      LnLxTnTx[,"Name"]<-as.character(LnLxTnTx[,"Name"])

      ##(6) Calculate Recyling Ratio and Recuperation Rate

      ##(6.1)
      ##Calculate Recycling Ratio

      if(length(LnLxTnTx[LnLxTnTx[,"Repeated"]==TRUE,"Repeated"])>0){

        ##identify repeated doses
        temp.Repeated<-LnLxTnTx[LnLxTnTx[,"Repeated"]==TRUE,c("Name","Dose","LxTx")]

        ##find corresponding previous dose for the repeated dose
        temp.Previous<-t(sapply(1:length(temp.Repeated[,1]),function(x){
          LnLxTnTx[LnLxTnTx[,"Dose"]==temp.Repeated[x,"Dose"] &
                     LnLxTnTx[,"Repeated"]==FALSE,c("Name","Dose","LxTx")]
        }))

        ##convert to data.frame
        temp.Previous<-as.data.frame(temp.Previous)

        ##set column names
        temp.ColNames<-sapply(1:length(temp.Repeated[,1]),function(x){
          paste(temp.Repeated[x,"Name"],"/",
                temp.Previous[temp.Previous[,"Dose"]==temp.Repeated[x,"Dose"],"Name"]
                ,sep="")
        })

        ##Calculate Recycling Ratio
        RecyclingRatio<-as.numeric(temp.Repeated[,"LxTx"])/as.numeric(temp.Previous[,"LxTx"])

        ##Just transform the matrix and add column names
        RecyclingRatio<-t(RecyclingRatio)
        colnames(RecyclingRatio) <- unique(temp.ColNames)

      }else{RecyclingRatio<-NA}

      ##(6.2)
      ##Recuperation Rate

      if("R0" %in% LnLxTnTx[,"Name"]==TRUE){
        Recuperation<-round(LnLxTnTx[LnLxTnTx[,"Name"]=="R0","LxTx"]/LnLxTnTx[LnLxTnTx[,"Name"]=="Natural","LxTx"],digits=4)
      }else{Recuperation<-NA}


      ##(6.3) IRSL
      ##Print IRSL Curves if IRSL curve is set
      sample.data@METADATA[,"SEL"]<-FALSE
      sample.data@METADATA[sample.data@METADATA["LTYPE"]=="IRSL" &
                             sample.data@METADATA[,"RUN"]%in%run==TRUE &
                             sample.data@METADATA[,"SET"]%in%set==TRUE,"SEL"]<-TRUE


      ##get IRSL curve ID & ID for Reg1 again
      IRSL.curveID<-sample.data@METADATA[sample.data@METADATA["SEL"]==TRUE & sample.data@METADATA["POSITION"]==i,"ID"]

      ##if no IRSL curve the length of the object is 0
      if(length(IRSL.curveID)>0){

        ##chose an IRSL curve with a dose of the first regeneration point
        Reg1again.curveID<-LnLxTnTx[LnLxTnTx[,"Repeated"]==TRUE & LnLxTnTx[,"Dose"]==LnLxTnTx[2,"Dose"],"LnLx.curveID"]

        if(length(Reg1again.curveID)>0){

          ##BOSL/IRSL
          IRSL_BOSL<-round(sum(unlist(sample.data@DATA[IRSL.curveID])[signal.integral])
                           /sum(unlist(sample.data@DATA[Reg1again.curveID])[signal.integral]),digits=4)
        }else{IRSL_BOSL<-NA}
      }else{IRSL_BOSL<-NA}

      ##Combine the two values
      if(exists("RejectionCriteria")==FALSE){
        RejectionCriteria<-cbind(RecyclingRatio,Recuperation,IRSL_BOSL)
      }else{
        RejectionCriteria.temp<-cbind(RecyclingRatio,Recuperation,IRSL_BOSL)
        RejectionCriteria<-rbind(RejectionCriteria,RejectionCriteria.temp)
      }

      ##============================================================================##
      ##PLOTTING
      ##============================================================================##

      if(output.plot){

        ##set plot settings
        plot.settings <- list(
          mtext = sample.data@METADATA[sample.data@METADATA[,"ID"]==LnLx.curveID[1],"SAMPLE"],
          log = ""

        )

        ##modify arguments
        plot.settings <- modifyList(plot.settings, list(...))



        if(output.plot.single==FALSE){
          layout(matrix(c(1,2,1,2,3,4,3,5),4,2,byrow=TRUE))
        }
        ##warning if number of curves exceed colour values
        if(length(col)<length(LnLx.curveID)){
          cat("\n[Analyse_SAR.OSLdata()] Warning: To many curves! Only the first",
              length(col),"curves are plotted!")
        }

        ##==========================================================================
        ##plot Ln,Lx Curves


        ##get maximum value of LnLx curves
        LnLx.curveMax<-max(unlist(sample.data@DATA[LnLx.curveID]))

        ##get channel resolution (it should be the same for all values)
        HIGH<-sample.data@METADATA[sample.data@METADATA[,"ID"]==LnLx.curveID[1],"HIGH"]
        NPOINTS<-sample.data@METADATA[sample.data@METADATA[,"ID"]==LnLx.curveID[1],"NPOINTS"]

        xaxt.values<-seq(HIGH/NPOINTS,HIGH,by=HIGH/NPOINTS)

        ##open plot area LnLx
        plot(NA,NA,
             xlab="Time [s]",
             ylab=paste("OSL [cts/",HIGH/NPOINTS," s]",sep=""),
             xlim=c(HIGH/NPOINTS,HIGH),
             ylim=c(1,max(unlist(sample.data@DATA[LnLx.curveID]))),
             main=expression(paste(L[n],",",L[x]," curves",sep="")),
             log=plot.settings$log
        )
        ##plot curves and get legend values
        sapply(1:length(LnLx.curveID),function(x){
          yaxt.values<-unlist(sample.data@DATA[LnLx.curveID[x]])
          lines(xaxt.values,yaxt.values,col=col[x])
        })

        ##mark integration limits
        abline(v=xaxt.values[min(signal.integral)], lty=2, col="gray")
        abline(v=xaxt.values[max(signal.integral)], lty=2, col="gray")
        abline(v=xaxt.values[min(background.integral)], lty=2, col="gray")
        abline(v=xaxt.values[max(background.integral)], lty=2, col="gray")


        ##plot legend
        legend("topright",as.character(LnLxTnTx$Name),lty=c(rep(1,length(LnLx.curveID))),
               cex=0.8*cex.global,col=col, bg="gray")

        ##sample name
        mtext(side=3,plot.settings$mtext,cex=0.7*cex.global)

        ##========================================================================
        ##open plot area TnTx
        plot(NA,NA,
             xlab="Time [s]",
             ylab=paste("OSL [cts/",HIGH/NPOINTS," s]",sep=""),
             xlim=c(HIGH/NPOINTS,HIGH),
             ylim=c(1,max(unlist(sample.data@DATA[TnTx.curveID]))),
             main=expression(paste(T[n],",",T[x]," curves",sep="")),
             log=plot.settings$log
        )
        ##plot curves and get legend values
        sapply(1:length(TnTx.curveID),function(x){
          yaxt.values<-unlist(sample.data@DATA[TnTx.curveID[x]])
          lines(xaxt.values,yaxt.values,col=col[x])
        })

        ##mark integration limits
        abline(v=xaxt.values[min(signal.integral)], lty=2, col="gray")
        abline(v=xaxt.values[max(signal.integral)], lty=2, col="gray")
        abline(v=xaxt.values[min(background.integral)], lty=2, col="gray")
        abline(v=xaxt.values[max(background.integral)], lty=2, col="gray")


        ##plot legend
        legend("topright",as.character(LnLxTnTx$Name),lty=c(rep(1,length(TnTx.curveID))),
               cex=0.8*cex.global,col=col, bg="gray")

        ##sample name
        mtext(side=3,plot.settings$mtext,cex=0.7*cex.global)

        ##========================================================================
        ##Print TL curves for TnTx -
        sample.data@METADATA[,"SEL"]<-FALSE
        sample.data@METADATA[sample.data@METADATA["LTYPE"]=="TL","SEL"]<-TRUE

        ##check if TL any curves is measured within the sequence
        if(length(sample.data@METADATA[sample.data@METADATA[,"SEL"]==TRUE,1])>0){


          ##to ensure that the right TL curves are used the run and set number of the LnLx and TnTx curves are used
          LnLx.SET<-sapply(LnLx.curveID,function(x){sample.data@METADATA[sample.data@METADATA["ID"]==x,"SET"]})
          LnLx.RUN<-sapply(LnLx.curveID,function(x){sample.data@METADATA[sample.data@METADATA["ID"]==x,"RUN"]})
          TnTx.SET<-sapply(TnTx.curveID,function(x){sample.data@METADATA[sample.data@METADATA["ID"]==x,"SET"]})
          TnTx.RUN<-sapply(TnTx.curveID,function(x){sample.data@METADATA[sample.data@METADATA["ID"]==x,"RUN"]})

          ##get TL curve IDs in general considering the constraints
          TL.curveID<-sapply(1:length(TnTx.curveID),function(x){results<-
                                                                  sample.data@METADATA[sample.data@METADATA["SEL"]==TRUE & sample.data@METADATA["POSITION"]==i &
                                                                                         sample.data@METADATA["SET"]>=LnLx.SET[x] & sample.data@METADATA["RUN"]>=LnLx.RUN[x] &
                                                                                         sample.data@METADATA["SET"]<=TnTx.SET[x] & sample.data@METADATA["RUN"]<=TnTx.RUN[x],"ID"]})

          ##get maximum value of  TL curves
          TL.curveMax<-max(unlist(sample.data@DATA[TL.curveID]))

          ##get channel resolution (it should be the same for all values)
          HIGH<-unique(sample.data@METADATA[sample.data@METADATA["ID"]==TL.curveID[1],"HIGH"])
          NPOINTS<-unique(sample.data@METADATA[sample.data@METADATA["ID"]==TL.curveID[1],"NPOINTS"])
          xaxt.values<-seq(HIGH/NPOINTS,HIGH,by=HIGH/NPOINTS)

          ##get heating rate
          RATE<-unique(sample.data@METADATA[sample.data@METADATA["ID"]==TL.curveID[1],"RATE"])

          ##open plot area for TL curves
          plot(NA,NA,
               xlab="T [\u00B0C]",
               ylab=paste("TL [cts/",HIGH/NPOINTS," \u00B0C]",sep=""),
               xlim=c(HIGH/NPOINTS,HIGH),
               ylim=c(1,TL.curveMax),
               main="Cutheat - TL curves",
               sub=paste("(",RATE," K/s)",sep=""),
               log=if(plot.settings$log=="y" | plot.settings$log=="xy"){"y"}else{""}
          )

          ##plot curves and get legend values
          sapply(1:length(TL.curveID),function(x){
            yaxt.values<-unlist(sample.data@DATA[TL.curveID[x]])
            lines(xaxt.values,yaxt.values,col=col[x])
          })

          ##plot legend
          legend("topleft",as.character(LnLxTnTx$Name),lty=c(rep(1,length(TL.curveID))),
                 cex=0.8*cex.global,col=col, bg="white", bty="n")

          ##sample name
          mtext(side=3,plot.settings$mtext,cex=0.7*cex.global)


        }else{
          plot(NA,NA,xlim=c(0,100),ylim=c(0,100), main="Cutheat - TL curves")
          text(50,50,"no cutheat as TL curve detected")
        }

        ##======================================================================##
        ##Print IRSL Curves if IRSL curve is set

        if(is.na(IRSL_BOSL) == FALSE){
          ##get channel resolution (it should be the same for all values)
          HIGH<-unique(sample.data@METADATA[sample.data@METADATA["ID"]==IRSL.curveID ,"HIGH"])
          NPOINTS<-unique(sample.data@METADATA[sample.data@METADATA["ID"]==IRSL.curveID ,"NPOINTS"])
          xaxt.values<-seq(HIGH/NPOINTS,HIGH,by=HIGH/NPOINTS)

          ##open plot IRSL curve
          plot(NA,NA,
               xlab="Time [s]",
               ylab=paste("OSL and IRSL [cts/",HIGH/NPOINTS," s]",sep=""),
               xlim=c(0,HIGH),
               ylim=c(0,max(unlist(sample.data@DATA[Reg1again.curveID]))),
               main="IRSLT"
          )

          ##show integral limits
          abline(v=xaxt.values[min(signal.integral)], lty=2, col="gray")
          abline(v=xaxt.values[max(signal.integral)], lty=2, col="gray")

          ##print(length(sample.data@DATA[IRSL.curveID]))
          lines(xaxt.values,unlist(sample.data@DATA[IRSL.curveID]),col="red")
          lines(xaxt.values,unlist(sample.data@DATA[Reg1again.curveID[1]]),col="blue")

          ##legend
          legend("topright",c("R1 again","IRSL"),lty=c(1,1),col=c("blue","red"), bty="n")


          mtext(side=3,paste("IRSL/BOSL = ",IRSL_BOSL*100,"%",sep=""),
                cex=.8*cex.global
          )
        }

        if(((is.na(IRSL_BOSL)==TRUE) & length(IRSL.curveID)>0) |
             ((is.na(IRSL_BOSL)==FALSE) & length(IRSL.curveID)>0)){

          ##plot only IRSL curve
          plot(xaxt.values,unlist(sample.data@DATA[IRSL.curveID]),
               xlab="Time [s]",
               ylab=paste("IRSL [cts/",HIGH/NPOINTS," s]",sep=""),
               xlim=c(0,10),
               ylim=c(0,max(unlist(sample.data@DATA[IRSL.curveID]))),
               main="IRSL curve (10 s)",
               col="red",
               type="l"
          )
        }else{
          plot(NA,NA,xlim=c(0,10), ylim=c(0,10), main="IRSL curve")
          text(5,5,"no IRSL curve detected")
        }
        ##=========================================================================
        ##Plot header
        if(output.plot.single==TRUE){
          mtext(side=3,paste("ALQ Pos. ",i,sep=""))
        }else{
          mtext(side=3,paste("ALQ Pos. ",i,sep=""),outer=TRUE,line=-2.5)
        }

        ##Plot footer
        mtext(side=4,info.measurement,outer=TRUE,line=-1.5,cex=0.6*cex.global, col="blue")

        ##output on terminal for plot
        writeLines(paste("\n[Analyse_SAR.OSLdata()] >> Figure for position ",i," produced.",sep=""))

        ##reset mfrow
        par(mfrow=c(1,1))



      }#endif for output.plot
      ##preprate output of values
      ##==============================================================================

      ##Add LnLxTnTx values to the list
      if(exists("LnLxTnTx_List")==FALSE){LnLxTnTx_List<-list()}
      LnLxTnTx_List[[i]]<-LnLxTnTx
      rm(LnLxTnTx)


    }else{writeLines(paste("[Analyse_SAR.OSLdata()] >> Position ",i," is not valid and has been omitted!",sep=""))} #end if position checking

  }#end for loop

  ##============================================================================##
  ##OUTPUT OF FUNCTION
  ##============================================================================##

  ##get further information from the position used

  ##this is what you get from the Risoe file
  readTemp<-unique(sample.data@METADATA[sample.data@METADATA[,"POSITION"]==min(position) & sample.data@METADATA[,"LTYPE"]!="TL","TEMPERATURE"])

  cutheat<-unique(sample.data@METADATA[sample.data@METADATA[,"POSITION"]==min(position)  &
                                         sample.data@METADATA[,"LTYPE"]=="TL","HIGH"])
  if(length(cutheat)==0){cutheat=NA}

  systemID<-unique(sample.data@METADATA[sample.data@METADATA[,"POSITION"]==min(position),"SYSTEMID"])

  SARParameters<-data.frame(readTemp=readTemp,cutheat=cutheat,systemID=systemID)

  return(list(LnLxTnTx=LnLxTnTx_List,
              RejectionCriteria=RejectionCriteria,
              SARParameters=SARParameters))


}
