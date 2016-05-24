#' @include get_Risoe.BINfileData.R set_Risoe.BINfileData.R
NULL

#' Class \code{"Risoe.BINfileData"}
#'
#' S4 class object for luminescence data in R. The object is produced as output
#' of the function \code{\link{read_BIN2R}}.
#'
#'
#' @name Risoe.BINfileData-class
#'
#' @docType class
#'
#' @slot METADATA Object of class "data.frame" containing the meta information for each curve.
#'
#' @slot DATA Object of class "list" containing numeric vector with count data.
#'
#' @slot .RESERVED Object of class "list" containing list of undocumented raw values for internal use only.
#'
#' @note
#'
#' \bold{Internal METADATA - object structure}
#'
#' \tabular{rllll}{
#' \bold{#} \tab \bold{Name} \tab \bold{Data Type} \tab \bold{V} \tab \bold{Description} \cr
#' [,1]  \tab ID  \tab \code{numeric} \tab RLum \tab Unique record ID (same ID as in slot \code{DATA})\cr
#' [,2]  \tab SEL \tab \code{logic} \tab RLum \tab Record selection, not part official BIN-format, triggered by TAG\cr
#' [,3]  \tab VERSION \tab \code{raw} \tab 03-07 \tab BIN-file version number \cr
#' [,4]  \tab LENGTH \tab \code{integer} \tab 03-07 \tab Length of this record\cr
#' [,5]  \tab PREVIOUS \tab \code{integer} \tab 03-07 \tab Length of previous record\cr
#' [,6]  \tab NPOINTS \tab \code{integer} \tab 03-07 \tab Number of data points in the record\cr
#' [,7]  \tab RUN \tab \code{integer} \tab 03-07 \tab Run number\cr
#' [,8]  \tab SET \tab \code{integer} \tab 03-07 \tab Set number\cr
#' [,9]  \tab POSITION \tab  \code{integer} \tab 03-07 \tab Position number\cr
#' [,10] \tab GRAIN \tab \code{integer} \tab 03-04 \tab Grain number\cr
#' [,11] \tab GRAINNUMBER \tab \code{integer} \tab 06-07 \tab Grain number\cr
#' [,12] \tab CURVENO \tab \code{integer} \tab 06-07 \tab Curve number\cr
#' [,13] \tab XCOORD \tab \code{integer} \tab 03-07 \tab X position of a single grain\cr
#' [,14] \tab YCOORD \tab \code{integer} \tab 03-07 \tab Y position of a single grain\cr
#' [,15] \tab SAMPLE \tab \code{factor} \tab 03-07 \tab Sample name\cr
#' [,16] \tab COMMENT \tab \code{factor} \tab 03-07 \tab Comment name\cr
#' [,17] \tab SYSTEMID \tab \code{integer} \tab 03-07 \tab Risoe system id\cr
#' [,18] \tab FNAME \tab \code{factor} \tab 06-07 \tab File name (*.bin/*.binx)\cr
#' [,19] \tab USER \tab \code{facotr} \tab 03-07 \tab User name\cr
#' [,20] \tab TIME \tab \code{character} \tab 03-07 \tab Data collection time (hh-mm-ss)\cr
#' [,21] \tab DATE \tab \code{factor} \tab 03-07 \tab Data collection date (ddmmyy)\cr
#' [,22] \tab DTYPE \tab \code{character} \tab 03-07 \tab Data type\cr
#' [,23] \tab BL_TIME \tab \code{numeric} \tab 03-07 \tab Bleaching time\cr
#' [,24] \tab BL_UNIT \tab \code{integer} \tab 03-07 \tab Bleaching unit (mJ, J, secs, mins, hrs)\cr
#' [,25] \tab NORM1 \tab \code{numeric} \tab 03-07 \tab Normalisation factor (1)\cr
#' [,26] \tab NORM2 \tab \code{numeric} \tab 03-07 \tab Normalisation factor (2)\cr
#' [,27] \tab NORM3 \tab \code{numeric} \tab 03-07 \tab Normalisation factor (3)\cr
#' [,28] \tab BG \tab \code{numeric} \tab 03-07 \tab Background level\cr
#' [,29] \tab SHIFT \tab \code{integer} \tab 03-07 \tab Number of channels to shift data\cr
#' [,30] \tab TAG \tab \code{integer} \tab 03-07 \tab Tag, triggers SEL\cr
#' [,31] \tab LTYPE \tab \code{character} \tab 03-07 \tab Luminescence type\cr
#' [,32] \tab LIGHTSOURCE \tab \code{character} \tab 03-07 \tab Light source\cr
#' [,33] \tab LPOWER \tab \code{numeric} \tab 03-07 \tab Optical stimulation power\cr
#' [,34] \tab LIGHTPOWER \tab \code{numeric} \tab 06-07 \tab Optical stimulation power\cr
#' [,35] \tab LOW \tab \code{numeric} \tab 03-07 \tab Low (temperature, time, wavelength)\cr
#' [,36] \tab HIGH \tab \code{numeric} \tab 03-07 \tab High (temperature, time, wavelength)\cr
#' [,37] \tab RATE \tab \code{numeric} \tab 03-07 \tab Rate (heating rate, scan rate)\cr
#' [,38] \tab TEMPERATURE \tab \code{integer} \tab 03-07 \tab Sample temperature\cr
#' [,39] \tab MEASTEMP \tab \code{integer} \tab 06-07 \tab Measured temperature\cr
#' [,40] \tab AN_TEMP \tab \code{numeric} \tab 03-07 \tab Annealing temperature\cr
#' [,41] \tab AN_TIME \tab \code{numeric} \tab 03-07 \tab Annealing time\cr
#' [,42] \tab TOLDELAY \tab \code{integer} \tab 03-07 \tab TOL 'delay' channels\cr
#' [,43] \tab TOLON \tab \code{integer} \tab 03-07 \tab TOL 'on' channels\cr
#' [,44] \tab TOLOFF \tab \code{integer} \tab 03-07 \tab TOL 'off' channels\cr
#' [,45] \tab IRR_TIME \tab \code{numeric} \tab 03-07 \tab Irradiation time\cr
#' [,46] \tab IRR_TYPE \tab \code{integer} \tab 03-07 \tab Irradiation type (alpha, beta or gamma)\cr
#' [,47] \tab IRR_UNIT \tab \code{integer} \tab 03-04 \tab Irradiation unit (Gy, Rads, secs, mins, hrs)\cr
#' [,48] \tab IRR_DOSERATE \tab \code{numeric} \tab 06-07 \tab Irradiation dose rate (Gy/s)\cr
#' [,49] \tab IRR_DOSERATEERR \tab \code{numeric} \tab 06-07 \tab Irradiation dose rate error (Gy/s)\cr
#' [,50] \tab TIMESINCEIRR \tab \code{integer} \tab 06-07 \tab Time since irradiation (s)\cr
#' [,51] \tab TIMETICK \tab \code{numeric} \tab 06-07 \tab Time tick for pulsing (s)\cr
#' [,52] \tab ONTIME \tab \code{integer} \tab 06-07 \tab On-time for pulsing (in time ticks)\cr
#' [,53] \tab STIMPERIOD \tab \code{integer} \tab 06-07 \tab Stimulation period (on+off in time ticks)\cr
#' [,54] \tab GATE_ENABLED \tab \code{raw} \tab 06-07 \tab PMT signal gating enabled\cr
#' [,55] \tab ENABLE_FLAGS \tab \code{raw} \tab 06-07 \tab PMT signal gating  enabled\cr
#' [,56] \tab GATE_START \tab \code{integer} \tab 06-07 \tab Start gating (in time ticks)\cr
#' [,57] \tab GATE_STOP \tab \code{ingeter} \tab 06-07 \tab Stop gating (in time ticks), 'Gateend' for version 04, here only GATE_STOP is used\cr
#' [,58] \tab PTENABLED \tab \code{raw} \tab 06-07 \tab Photon time enabled\cr
#' [,59] \tab DTENABLED \tab \code{raw} \tab 06-07 \tab PMT dead time correction enabled\cr
#' [,60] \tab DEADTIME \tab \code{numeric} \tab 06-07 \tab PMT dead time (s)\cr
#' [,61] \tab MAXLPOWER \tab \code{numeric} \tab 06-07 \tab Stimulation power to 100 percent (mW/cm^2)\cr
#' [,62] \tab XRF_ACQTIME \tab \code{numeric} \tab 06-07 \tab XRF acquisition time (s)\cr
#' [,63] \tab XRF_HV \tab \code{numeric} \tab 06-07 \tab XRF X-ray high voltage (V)\cr
#' [,64] \tab XRF_CURR \tab \code{integer} \tab 06-07 \tab XRF X-ray current (uA)\cr
#' [,65] \tab XRF_DEADTIMEF \tab \code{numeric} \tab 06-07 \tab XRF dead time fraction\cr
#' [,66] \tab SEQUENCE \tab \code{character} \tab 03-04 \tab Sequence name\cr
#' [,67] \tab DETECTOR_ID \tab \code{raw} \tab 07 \tab Detector ID\cr
#' [,68] \tab LOWERFILTER_ID \tab \code{integer} \tab 07 \tab Lower filter ID in reader\cr
#' [,69] \tab UPPERFILTER_ID \tab \code{integer} \tab 07 \tab Uper filter ID in reader\cr
#' [,70] \tab ENOISEFACTOR \tab \code{numeric} \tab 07 \tab Excess noise filter, usage unknown
#'
#' } V = BIN-file version (RLum means that it does not depend on a specific BIN
#' version)\cr
#'
#' Note that the \code{Risoe.BINfileData} object combines all values from
#' different versions from the BIN-file, reserved bits are skipped, however,
#' the function \code{\link{write_R2BIN}} reset arbitrary reserved bits. Invalid
#' values for a specific version are set to \code{NA}. Furthermore, the
#' internal R data types do not necessarily match the required data types for
#' the BIN-file data import! Data types are converted during data import.\cr
#'
#' \bold{LTYPE} values
#'
#' \tabular{rll}{ [,0] \tab TL \tab: Thermoluminescence \cr [,1] \tab OSL \tab:
#' Optically stimulated luminescence \cr [,2] \tab IRSL \tab: Infrared
#' stimulated luminescence \cr [,3] \tab M-IR \tab: Infrared monochromator
#' scan\cr [,4] \tab M-VIS \tab: Visible monochromator scan\cr [,5] \tab TOL
#' \tab: Thermo-optical luminescence \cr [,6] \tab TRPOSL \tab: Time Resolved
#' Pulsed OSL\cr [,7] \tab RIR \tab: Ramped IRSL\cr [,8] \tab RBR \tab: Ramped
#' (Blue) LEDs\cr [,9] \tab USER \tab: User defined\cr [,10] \tab POSL \tab:
#' Pulsed OSL \cr [,11] \tab SGOSL \tab: Single Grain OSL\cr [,12] \tab RL
#' \tab: Radio Luminescence \cr [,13] \tab XRF \tab: X-ray Fluorescence }
#'
#' \bold{DTYPE} values \tabular{rll}{ [,0] \tab 0 \tab Natural \cr [,1] \tab 1
#' \tab N+dose \cr [,2] \tab 2 \tab Bleach \cr [,3] \tab 3 \tab Bleach+dose \cr
#' [,4] \tab 4 \tab Natural (Bleach) \cr [,5] \tab 5 \tab N+dose (Bleach) \cr
#' [,6] \tab 6 \tab Dose \cr [,7] \tab 7 \tab Background }
#'
#' \bold{LIGHTSOURCE} values \tabular{rll}{ [,0] \tab 0 \tab Non \cr [,1] \tab
#' 1 \tab Lamp \cr [,2] \tab 2 \tab IR diodes/IR Laser \cr [,3] \tab 3 \tab
#' Calibration LED \cr [,4] \tab 4 \tab Blue Diodes \cr [,5] \tab 5 \tab White
#' lite \cr [,6] \tab 6 \tab Green laser (single grain) \cr [,7] \tab 7 \tab IR
#' laser (single grain) }
#'
#' (information on the BIN/BINX file format are kindly provided by Risoe, DTU
#' Nutech)
#'
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("Risoe.BINfileData", ...)}.
#'
#' @section Function version: 0.2.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso
#' \code{\link{plot_Risoe.BINfileData}}, \code{\link{read_BIN2R}},
#' \code{\link{write_R2BIN}},\code{\link{merge_Risoe.BINfileData}},
#' \code{\link{Risoe.BINfileData2RLum.Analysis}},
#' \code{\link{Risoe.BINfileData2RLum.Data.Curve}}
#'
#' @references Risoe DTU, 2013. The Sequence Editor User Manual - Feb 2013 and Risoe DTU, 2015. The
#' Sequence Editor User Manual - March 2015
#'
#' \code{http://www.nutech.dtu.dk/}
#'
#' @keywords classes
#'
#' @examples
#'
#' showClass("Risoe.BINfileData")
#'
#' @export
setClass("Risoe.BINfileData",
         slots = list(
           METADATA="data.frame",
           DATA = "list",
           .RESERVED = "list"
           )
         )

##set generic S4 function for object
#' @describeIn Risoe.BINfileData
#' Show structure of RLum and Risoe.BINfile class objects
#' @export
setMethod("show", signature(object = "Risoe.BINfileData"),
          function(object){

            version<-paste(unique(object@METADATA[,"VERSION"]), collapse = ", ")
            systemID<-paste(unique(object@METADATA[,"SYSTEMID"]), collapse = ", ")
            filename <- as.character(object@METADATA[1,"FNAME"])
            records.overall<-length(object@DATA)
            records.type<-table(object@METADATA[,"LTYPE"])
            user<-paste(unique(as.character(object@METADATA[,"USER"])), collapse = ", ")
            date<-paste(unique(as.character(object@METADATA[,"DATE"])), collapse = ", ")
            run.range<-range(object@METADATA[,"RUN"])
            set.range<-range(object@METADATA[,"SET"])
            grain.range <- range(object@METADATA[,"GRAIN"])
            pos.range<-range(object@METADATA[,"POSITION"])

            records.type.count <- sapply(1:length(records.type),
              function(x){paste(
              names(records.type)[x],"\t(n = ",records.type[x],")",sep="")
              })

            records.type.count <- paste(records.type.count,
                                        collapse="\n\t                      ")

            ##print
            cat("\n[Risoe.BINfileData object]")
            cat("\n\n\tBIN/BINX version     ", version)
            if(version>=6){
              cat("\n\tFile name:           ", filename)
            }
            cat("\n\tObject date:         ", date)
            cat("\n\tUser:                ", user)
            cat("\n\tSystem ID:           ", ifelse(systemID == 0,"0 (unknown)", systemID))
            cat("\n\tOverall records:     ", records.overall)
            cat("\n\tRecords type:        ", records.type.count)
            cat("\n\tPosition range:      ",pos.range[1],":",pos.range[2])
            cat("\n\tGrain range:         ",grain.range[1],":",grain.range[2])
            cat("\n\tRun range:           ",run.range[1],":",run.range[2])
            cat("\n\tSet range:           ",set.range[1],":",set.range[2])
          }#end function
          )#end setMethod


# constructor (set) method for object class -----------------------------------

#' @describeIn Risoe.BINfileData
#' The Risoe.BINfileData is normally produced as output of the function read_BIN2R.
#' This construction method is intended for internal usage only.
#'
#' @param METADATA Object of class "data.frame" containing the meta information
#' for each curve.
#'
#' @param DATA Object of class "list" containing numeric vector with count data.
#'
#' @param .RESERVED Object of class "list" containing list of undocumented raw
#' values for internal use only.
#' @export
setMethod("set_Risoe.BINfileData",
          signature = c(
            METADATA = "data.frame", DATA = "list", .RESERVED = "ANY"
          ),

          function(METADATA, DATA, .RESERVED) {
            if (missing(.RESERVED)) {
              .RESERVED <- list()
            }

            new(
              "Risoe.BINfileData",
              METADATA = METADATA,
              DATA = DATA,
              .RESERVED = .RESERVED
            )

          })


# accessor (get) method for object class -----------------------------------

#' @describeIn Risoe.BINfileData
#' Formal get-method for Risoe.BINfileData object. It does not allow accessing
#' the object directly, it is just showing a terminal message.
#'
#' @param object an object of class \code{\linkS4class{Risoe.BINfileData}}
#'
#' @param ... other arguments that might be passed
#'
#' @export
setMethod("get_Risoe.BINfileData",
          signature= "Risoe.BINfileData",
          definition = function(object, ...) {

            cat("[get_Risoe.BINfileData()] No direct access is provided for this object type. Use the function 'Risoe.BINfileData2RLum.Analysis' for object coercing.")

          })##end setMethod

##-------------------------------------------------------------------------------------------------##
##=================================================================================================##
