#' Plot single luminescence curves from a BIN file object
#'
#' Plots single luminescence curves from an object returned by the
#' \link{read_BIN2R} function.
#'
#' \bold{Nomenclature}\cr
#'
#' See \code{\link{Risoe.BINfileData-class}}
#'
#' \bold{curve.transformation}\cr
#'
#' This argument allows transforming continuous wave (CW) curves to pseudo
#' (linear) modulated curves. For the transformation, the functions of the
#' package are used.  Currently, it is not possible to pass further arguments
#' to the transformation functions. The argument works only for \code{ltype}
#' \code{OSL} and \code{IRSL}.\cr
#'
#' \bold{Irradiation time}\cr
#'
#' Plotting the irradiation time (s) or the given dose (Gy) requires that the
#' variable \code{IRR_TIME} has been set within the BIN-file.  This is normally
#' done by using the 'Run Info' option within the Sequence Editor or by editing
#' in R.
#'
#' @param BINfileData \link{Risoe.BINfileData-class} (\bold{required}):
#' requires an S4 object returned by the \link{read_BIN2R} function.
#' @param position \link{vector} (optional): option to limit the plotted curves
#' by position (e.g. \code{position = 1}, \code{position = c(1,3,5)}).
#' @param run \link{vector} (optional): option to limit the plotted curves by
#' run (e.g., \code{run = 1}, \code{run = c(1,3,5)}).
#' @param set \link{vector} (optional): option to limit the plotted curves by
#' set (e.g., \code{set = 1}, \code{set = c(1,3,5)}).
#' @param sorter \link{character} (with default): the plot output can be
#' ordered by "POSITION","SET" or "RUN". POSITION, SET and RUN are options
#' defined in the Risoe Sequence Editor.
#' @param ltype \link{character} (with default): option to limit the plotted
#' curves by the type of luminescence stimulation.  Allowed values:
#' \code{"IRSL"}, \code{"OSL"},\code{"TL"}, \code{"RIR"}, \code{"RBR"}
#' (corresponds to LM-OSL), \code{"RL"}.  All type of curves are plotted by
#' default.
#' @param curve.transformation \link{character} (optional): allows transforming
#' CW-OSL and CW-IRSL curves to pseudo-LM curves via transformation functions.
#' Allowed values are: \code{CW2pLM}, \code{CW2pLMi}, \code{CW2pHMi} and
#' \code{CW2pPMi}. See details.
#' @param dose_rate \link{numeric} (optional): dose rate of the irradition
#' source at the measurement date. If set, the given irradiation dose will be
#' shown in Gy.  See details.
#' @param temp.lab \link{character} (optional): option to allow for different
#' temperature units. If no value is set deg. C is chosen.
#' @param cex.global \link{numeric} (with default): global scaling factor.
#' @param \dots further undocumented plot arguments.
#' @return Returns a plot.
#' @note The function has been successfully tested for the Sequence Editor file
#' output version 3 and 4.
#' @section Function version: 0.4.1
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France),\cr Michael Dietze, GFZ Potsdam (Germany)
#' @seealso \code{\link{Risoe.BINfileData-class}},\code{\link{read_BIN2R}},
#' \code{\link{CW2pLM}}, \code{\link{CW2pLMi}}, \code{\link{CW2pPMi}},
#' \code{\link{CW2pHMi}}
#' @references Duller, G., 2007. Analyst. pp. 1-45.
#' @keywords dplot
#' @examples
#'
#'
#' ##load data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##plot all curves from the first position to the desktop
#' #pdf(file = "~/Desktop/CurveOutput.pdf", paper = "a4", height = 11, onefile = TRUE)
#'
#' ##example - load from *.bin file
#' #BINfile<- file.choose()
#' #BINfileData<-read_BIN2R(BINfile)
#'
#' #par(mfrow = c(4,3), oma = c(0.5,1,0.5,1))
#' #plot_Risoe.BINfileData(CWOSL.SAR.Data,position = 1)
#' #mtext(side = 4, BINfile, outer = TRUE, col = "blue", cex = .7)
#' #dev.off()
#'
#' @export
plot_Risoe.BINfileData<- function(
  BINfileData,
  position,
  run,
  set,
  sorter = "POSITION",
  ltype = c("IRSL","OSL","TL","RIR","RBR","RL"),
  curve.transformation,
  dose_rate,
  temp.lab,
  cex.global = 1,
  ...
){

  ##check if the object is of type Risoe.BINfileData
  if(class(BINfileData)!="Risoe.BINfileData"){stop("Wrong object! Object of type Risoe.BINfileData needed.")}

  temp<-BINfileData

  # Missing check ----------------------------------------------------------------

  ##set plot position if missing
  if(missing(position)==TRUE){position<-c(min(temp@METADATA[,"POSITION"]):max(temp@METADATA[,"POSITION"]))}
  if(missing(run)==TRUE){run<-c(min(temp@METADATA[,"RUN"]):max(temp@METADATA[,"RUN"]))}
  if(missing(set)==TRUE){set<-c(min(temp@METADATA[,"SET"]):max(temp@METADATA[,"SET"]))}

  ##temp.lab
  if(missing(temp.lab) == TRUE){temp.lab <- "\u00B0C"}


  ##fun
  extraArgs <- list(...) # read out additional arguments list
  fun       <- if("fun" %in% names(extraArgs)) {extraArgs$fun} else {FALSE}

  # Ordering --------------------------------------------------------------------

  ##(1) order by RUN, SET OR BY POSITION
  if(sorter=="RUN"){
    temp@METADATA<-temp@METADATA[order(temp@METADATA[,"RUN"]),]
  }else if(sorter=="SET"){
    temp@METADATA<-temp@METADATA[order(temp@METADATA[,"SET"]),]
  }else {
    temp@METADATA<-temp@METADATA[order(temp@METADATA[,"POSITION"]),]
  }



  # Select values for plotting ------------------------------------------------------------------

  ##(2) set SEL for selected position

  ##set all to FALSE
  temp@METADATA[,"SEL"]<-FALSE

  ##set TRUE
  temp@METADATA[(temp@METADATA[,"POSITION"] %in% position)==TRUE &
                  (temp@METADATA[,"RUN"] %in% run)==TRUE &
                  (temp@METADATA[,"SET"] %in% set)==TRUE &
                  (temp@METADATA[,"LTYPE"] %in% ltype)==TRUE,"SEL"]<-TRUE

  ##------------------------------------------------------------------------##
  ##PLOTTING
  ##------------------------------------------------------------------------##
  ##(3) plot curves
  for(i in 1:length(temp@METADATA[,"ID"])){

    ##print only if SEL == TRUE
    if(temp@METADATA[i,"SEL"]==TRUE)
    {

      ##find measured unit
      measured_unit<-if(temp@METADATA[i,"LTYPE"]=="TL"){" \u00B0C"}else{"s"}

      ##set x and y values
      values.x <- seq(temp@METADATA[i,"HIGH"]/temp@METADATA[i,"NPOINTS"],
                      temp@METADATA[i,"HIGH"],by=temp@METADATA[i,"HIGH"]/temp@METADATA[i,"NPOINTS"])
      values.y <- unlist(temp@DATA[temp@METADATA[i,"ID"]])
      values.xy <- data.frame(values.x, values.y)

      ##set curve transformation if wanted
      if((temp@METADATA[i,"LTYPE"] == "OSL" | temp@METADATA[i,"LTYPE"] == "IRSL") &
           missing(curve.transformation) == FALSE){

        if(curve.transformation=="CW2pLM"){

          values.xy <- CW2pLM(values.xy)

        }else if(curve.transformation=="CW2pLMi"){

          values.xy <- CW2pLMi(values.xy)[,1:2]

        }else if(curve.transformation=="CW2pHMi"){

          values.xy <- CW2pHMi(values.xy)[,1:2]

        }else if(curve.transformation=="CW2pPMi"){

          values.xy <- CW2pPMi(values.xy)[,1:2]

        }else{

          warning("Function for curve.transformation is unknown. No transformation is performed.")

        }

      }

      ##plot graph
      plot(values.xy,
           main=paste("pos=", temp@METADATA[i,"POSITION"],", run=", temp@METADATA[i,"RUN"],
                      ", set=", temp@METADATA[i,"SET"],sep=""
           ),
           type="l",
           ylab=paste(temp@METADATA[i,"LTYPE"]," [cts/",round(temp@METADATA[i,"HIGH"]/temp@METADATA[i,"NPOINTS"],digits=3)," ",
                      measured_unit,"]",sep=""),
           xlab=if(measured_unit=="\u00B0C"){paste("temp. [",temp.lab,"]",sep="")}else{"time [s]"},
           col=if(temp@METADATA[i,"LTYPE"]=="IRSL" | temp@METADATA[i,"LTYPE"]=="RIR"){"red"}
           else if(temp@METADATA[i,"LTYPE"]=="OSL" | temp@METADATA[i,"LTYPE"]=="RBR"){"blue"}
           else{"black"},
           sub=if(temp@METADATA[i,"LTYPE"]=="TL"){paste("(",temp@METADATA[i,"RATE"]," K/s)",sep="")}else{},
           lwd=1.2*cex.global,
           cex=0.9*cex.global
      )

      ##add mtext for temperature

      ##grep temperature (different for different verions)

      temperature<-if(temp@METADATA[i,"VERSION"]=="03"){temp@METADATA[i,"AN_TEMP"]}
      else{temp@METADATA[i,"TEMPERATURE"]}

      ##mtext
      mtext(side=3,
            if(temp@METADATA[i,"LTYPE"]=="TL"){paste("TL to ",temp@METADATA[i,"HIGH"], " ",temp.lab,sep="")}
            else{paste(temp@METADATA[i,"LTYPE"],"@",temperature," ",temp.lab ,sep="")},
            cex=0.9*cex.global)

      ##add mtext for irradiation
      mtext(side=4,cex=0.8*cex.global, line=0.5,
            if(temp@METADATA[i, "IRR_TIME"]!=0){

              if(missing("dose_rate")==TRUE){
                paste("dose = ",temp@METADATA[i, "IRR_TIME"], " s", sep="")
              }else{
                paste("dose = ",temp@METADATA[i, "IRR_TIME"]*dose_rate, " Gy", sep="")
              }
            }
      )#end mtext

    }#endif::selection

  }#endforloop

  if(fun==TRUE){sTeve()}

}
