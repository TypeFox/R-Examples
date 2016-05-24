#' Analyse SAR CW-OSL measurements
#'
#' The function performs a SAR CW-OSL analysis on an
#' \code{\linkS4class{RLum.Analysis}} object including growth curve fitting.
#'
#' The function performs an analysis for a standard SAR protocol measurements
#' introduced by Murray and Wintle (2000) with CW-OSL curves. For the
#' calculation of the Lx/Tx value the function \link{calc_OSLLxTxRatio} is
#' used. For \bold{changing the way the Lx/Tx error is calculated} use the argument
#' \code{background.count.distribution} and \code{sigmab}, which will be passed to the function
#' \link{calc_OSLLxTxRatio}.\cr\cr
#'
#' \bold{Argument \code{object} is of type \code{list}}\cr\cr
#'
#' If the argument \code{object} is of type \code{\link{list}} containing \bold{only}
#' \code{\linkS4class{RLum.Analysis}} objects, the function re-calls itself as often as elements
#' are in the list. This is usefull if an entire measurement wanted to be analysed without
#' writing separate for-loops. To gain in full control of the parameters (e.g., \code{dose.points}) for
#' every aliquot (corresponding to one \code{\linkS4class{RLum.Analysis}} object in the list), in
#' this case the arguments can be provided as \code{\link{list}}. This \code{list} should
#' be of similar length as the \code{list} provided with the argument \code{object}, otherwise the function
#' will create an own list of the requested lenght. Function output will be just one single \code{\linkS4class{RLum.Results}} object.
#'
#' Please be careful when using this option. It may allow a fast an efficient data analysis, but
#' the function may also break with an unclear error message, due to wrong input data.\cr\cr
#'
#' \bold{Working with IRSL data}\cr\cr
#'
#' The function was originally designed to work just for 'OSL' curves,
#' following the principles of the SAR protocol. An IRSL measurement protocol
#' may follow this procedure, e.g., post-IR IRSL protocol (Thomsen et al.,
#' 2008). Therefore this functions has been enhanced to work with IRSL data,
#' however, the function is only capable of analysing curves that follow the
#' SAR protocol structure, i.e., to analyse a post-IR IRSL protocol, curve data
#' have to be pre-selected by the user to fit the standards of the SAR
#' protocol, i.e., Lx,Tx,Lx,Tx and so on. \cr
#'
#' Example: Imagine the measurement contains pIRIR50 and pIRIR225 IRSL curves.
#' Only one curve type can be analysed at the same time: The pIRIR50 curves or
#' the pIRIR225 curves.\cr\cr
#'
#' \bold{Supported rejection criteria}\cr\cr \sQuote{recycling.ratio}:
#' calculated for every repeated regeneration dose point.\cr
#'
#' \sQuote{recuperation.rate}: recuperation rate calculated by comparing the
#' Lx/Tx values of the zero regeneration point with the Ln/Tn value (the Lx/Tx
#' ratio of the natural signal). For methodological background see Aitken and
#' Smith (1988).\cr
#'
#' \sQuote{palaeodose.error}: set the allowed error for the De value, which per
#' default should not exceed 10\%.
#'
#' @param object \code{\linkS4class{RLum.Analysis}} (\bold{required}): input
#' object containing data for analysis, alternatively a \code{\link{list}} of
#' \code{\linkS4class{RLum.Analysis}} objects can be provided.
#'
#' @param signal.integral.min \code{\link{integer}} (\bold{required}): lower
#' bound of the signal integral. Can be a \code{\link{list}} of \code{\link{integer}s}, if \code{object} is
#' of type \code{\link{list}}. If the input is vector (e.g., \code{c(1,2)}) the 2nd value will be interpreted
#' as the minimum signal integral for the Tx curve.
#'
#' @param signal.integral.max \code{\link{integer}} (\bold{required}): upper
#' bound of the signal integral. Can be a \code{\link{list}} of \code{\link{integer}s}, if \code{object} is
#' of type \code{\link{list}}. If the input is vector (e.g., \code{c(1,2)}) the 2nd value will be interpreted
#' as the maximum signal integral for the Tx curve.
#'
#' @param background.integral.min \code{\link{integer}} (\bold{required}):
#' lower bound of the background integral. Can be a \code{\link{list}} of \code{\link{integer}s}, if \code{object} is
#' of type \code{\link{list}}. If the input is vector (e.g., \code{c(1,2)}) the 2nd value will be interpreted
#' as the minimum background integral for the Tx curve.
#'
#' @param background.integral.max \code{\link{integer}} (\bold{required}):
#' upper bound of the background integral. Can be a \code{\link{list}} of \code{\link{integer}s}, if \code{object} is
#' of type \code{\link{list}}. If the input is vector (e.g., \code{c(1,2)}) the 2nd value will be interpreted
#' as the maximum background integral for the Tx curve.
#'
#' @param rejection.criteria \code{\link{list}} (with default): provide a named list
#' and set rejection criteria in percentage for further calculation. Can be a \code{\link{list}} in
#' a \code{\link{list}}, if \code{object} is of type \code{\link{list}}
#'
#' Allowed #' arguments are \code{recycling.ratio}, \code{recuperation.rate},
#' \code{palaeodose.error} and \code{exceed.max.regpoint = TRUE/FALSE}.
#' Example: \code{rejection.criteria = list(recycling.ratio = 10)}.
#' Per default all numericla values are set to 10.
#'
#' @param dose.points \code{\link{numeric}} (optional): a numeric vector
#' containg the dose points values Using this argument overwrites dose point
#' values in the signal curves. Can be a \code{\link{list}} of \code{\link{numeric}} vectors,
#' if \code{object} is of type \code{\link{list}}
#'
#' @param mtext.outer \code{\link{character}} (optional): option to provide an
#' outer margin mtext. Can be a \code{\link{list}} of \code{\link{character}s},
#' if \code{object} is of type \code{\link{list}}
#'
#' @param plot \code{\link{logical}} (with default): enables or disables plot
#' output.
#'
#' @param plot.single \code{\link{logical}} (with default) or
#' \code{\link{numeric}} (optional): single plot output (\code{TRUE/FALSE}) to
#' allow for plotting the results in single plot windows. If a numerice vector
#' is provided the plots can be selected individually, i.e. \code{plot.single =
#' c(1,2,3,4)} will plot the TL and Lx, Tx curves but not the legend (5) or the
#' growth curve (6), (7) and (8) belong to rejection criteria plots. Requires
#' \code{plot = TRUE}.
#'
#' @param \dots further arguments that will be passed to the function
#' \code{\link{plot_GrowthCurve}} or \code{\link{calc_OSLLxTxRatio}}
#' (supported: \code{background.count.distribution} and \code{sigmab}). \bold{Please note} that
#' if you consider to use the early light subtraction method you should provide your own \code{sigmab}
#' value!
#'
#'
#' @return A plot (optional) and an \code{\linkS4class{RLum.Results}} object is
#' returned containing the following elements:
#' \item{De.values}{\link{data.frame} containing De-values, De-error and
#' further parameters} \item{LnLxTnTx.values}{\link{data.frame} of all
#' calculated Lx/Tx values including signal, background counts and the dose
#' points} \item{rejection.criteria}{\link{data.frame} with values that might
#' by used as rejection criteria. NA is produced if no R0 dose point exists.}
#' \item{Formula}{\link{formula} formula that have been used for the growth
#' curve fitting }\cr The output should be accessed using the function
#' \code{\link{get_RLum}}.
#'
#'
#' @note This function must not be mixed up with the function
#' \code{\link{Analyse_SAR.OSLdata}}, which works with
#' \link{Risoe.BINfileData-class} objects.\cr
#'
#' \bold{The function currently does only support 'OSL' or 'IRSL' data!}
#'
#' @section Function version: 0.7.1
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#'
#' @seealso \code{\link{calc_OSLLxTxRatio}}, \code{\link{plot_GrowthCurve}},
#' \code{\linkS4class{RLum.Analysis}}, \code{\linkS4class{RLum.Results}}
#' \code{\link{get_RLum}}
#'
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
#'
#' Thomsen, K.J., Murray, A.S., Jain, M., Boetter-Jensen, L., 2008. Laboratory
#' fading rates of various luminescence signals from feldspar-rich sediment
#' extracts. Radiation Measurements 43, 1474-1486.
#' doi:10.1016/j.radmeas.2008.06.002
#'
#' @keywords datagen plot
#'
#' @examples
#'
#' ##load data
#' ##ExampleData.BINfileData contains two BINfileData objects
#' ##CWOSL.SAR.Data and TL.SAR.Data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##transform the values from the first position in a RLum.Analysis object
#' object <- Risoe.BINfileData2RLum.Analysis(CWOSL.SAR.Data, pos=1)
#'
#' ##perform SAR analysis
#' results <- analyse_SAR.CWOSL(object,
#'                   signal.integral.min = 1,
#'                   signal.integral.max = 2,
#'                   background.integral.min = 900,
#'                   background.integral.max = 1000,
#'                   log = "x",
#'                   fit.method = "EXP")
#'
#' ##show De results
#' get_RLum(results)
#'
#' ##show LnTnLxTx table
#' get_RLum(results, data.object = "LnLxTnTx.table")
#'
#' @export
analyse_SAR.CWOSL<- function(
  object,
  signal.integral.min,
  signal.integral.max,
  background.integral.min,
  background.integral.max,
  rejection.criteria,
  dose.points = NULL,
  mtext.outer,
  plot = TRUE,
  plot.single = FALSE,
  ...
) {

# SELF CALL -----------------------------------------------------------------------------------
if(is.list(object)){

  ##make live easy
  if(missing("signal.integral.min")){
    signal.integral.min <- 1
    warning("[analyse_SAR.CWOSL()] 'signal.integral.min' missing, set to 1", call. = FALSE)
  }

  if(missing("signal.integral.max")){
    signal.integral.max <- 2
    warning("[analyse_SAR.CWOSL()] 'signal.integral.max' missing, set to 2", call. = FALSE)
  }

  ##now we have to extend everything to allow list of arguments ... this is just consequent
  signal.integral.min <- rep(list(signal.integral.min), length = length(object))
  signal.integral.max <- rep(list(signal.integral.max), length = length(object))
  background.integral.min <- rep(list(background.integral.min), length = length(object))
  background.integral.max <- rep(list(background.integral.max), length = length(object))


  if(!missing(rejection.criteria)){
    rejection.criteria <- rep(list(rejection.criteria), length = length(object))

  }


  if(!is.null(dose.points)){

    if(is(dose.points, "list")){
      dose.points <- rep(dose.points, length = length(object))

    }else{
      dose.points <- rep(list(dose.points), length = length(object))

    }

  }else{
    dose.points <- rep(list(NULL), length(object))

  }

  if(!missing(mtext.outer)){
    mtext.outer <- rep(as.list(mtext.outer), length = length(object))

  }else{
    mtext.outer <- rep(list(""), length = length(object))

  }

   ##run analysis
   temp <- lapply(1:length(object), function(x){

    analyse_SAR.CWOSL(object[[x]],
                      signal.integral.min = signal.integral.min[[x]],
                      signal.integral.max = signal.integral.max[[x]],
                      background.integral.min = background.integral.min[[x]],
                      background.integral.max = background.integral.max[[x]] ,
                      dose.points = dose.points[[x]],
                      mtext.outer = mtext.outer[[x]],
                      plot = plot,
                      plot.single = plot.single,
                      main = ifelse("main"%in% names(list(...)), list(...)$main, paste0("ALQ #",x)),
                      ...)

  })

  ##combine everything to one RLum.Results object as this as what was written ... only
  ##one object

  ##merge results and check if the output became NULL
  results <- merge_RLum(temp)

  ##DO NOT use invisible here, this will stop the function from stopping
  if(length(results) == 0){
    return(NULL)

  }else{
    return(results)

  }

}

# CONFIG  -----------------------------------------------------------------

  ##set error list, this allows to set error messages without breaking the function
  error.list <- list()

# General Integrity Checks ---------------------------------------------------

  ##GENERAL

    ##MISSING INPUT
    if(missing("object")){
      stop("[analyse_SAR.CWOSL()] No value set for 'object'!")
    }

    ##INPUT OBJECTS
    if(!is(object, "RLum.Analysis")){
      stop("[analyse_SAR.CWOSL()] Input object is not of type 'RLum.Analyis'!")
    }


    if(missing("signal.integral.min") & !is.list(object)){
      signal.integral.min <- 1
      warning("[analyse_SAR.CWOSL()] 'signal.integral.min' missing, set to 1", call. = FALSE)
    }

    if(missing("signal.integral.max") & !is.list(object)){
      signal.integral.min <- 2
      warning("[analyse_SAR.CWOSL()] 'signal.integral.max' missing, set to 2", call. = FALSE)
    }

    if(missing("background.integral.min")){
     stop("[analyse_SAR.CWOSL()] No value set for 'background.integral.min'!")
    }

    if(missing("background.integral.max")){
      stop("[analyse_SAR.CWOSL()] No value set for 'background.integral.max'!")
    }


      ##build signal and background integrals
      signal.integral <- c(signal.integral.min[1]:signal.integral.max[1])
      background.integral <- c(background.integral.min[1]:background.integral.max[1])

        ##account for the case that Lx and Tx integral differ
        if (length(signal.integral.min) == 2 &
            length(signal.integral.max) == 2) {
          signal.integral.Tx <-
            c(signal.integral.min[2]:signal.integral.max[2])

        }else{
          signal.integral.Tx <- NULL

        }

        if (length(background.integral.min) == 2 &
            length(background.integral.max) == 2) {
          background.integral.Tx <-
            c(background.integral.min[2]:background.integral.max[2])

        }else{
          background.integral.Tx <- NULL

        }

        ##Account for the case that the use did not provide everything ...
        if(is.null(signal.integral.Tx) & !is.null(background.integral.Tx)){
          signal.integral.Tx <- signal.integral

          warning("[analyse_SAR.CWOSL()] background integral for Tx curves set, but not for the signal integral; signal integral for Tx automatically set.")
        }

      if(!is.null(signal.integral.Tx) & is.null(background.integral.Tx)){
        background.integral.Tx <- background.integral

        warning("[analyse_SAR.CWOSL()] signal integral for Tx curves set, but not for the background integral; background integral for Tx automatically set.")
      }


    ##INTEGRAL LIMITS
    if(!is(signal.integral, "integer") | !is(background.integral, "integer")){
      stop("[analyse_SAR.CWOSL()] 'signal.integral' or 'background.integral' is not
           of type integer!")
    }


    ##CHECK IF DATA SET CONTAINS ANY OSL curve
    if(!TRUE%in%grepl("OSL", structure_RLum(object)$recordType) &&
       !TRUE%in%grepl("IRSL", structure_RLum(object)$recordType)){

      stop("[analyse_SAR.CWOSL()] No record of type 'OSL' or 'IRSL' are detected in the sequence
object!")

    }

    ##Check if any OSL curve is measured, if not set curve type on IRSL
    ##to allow further proceedings
    CWcurve.type  <- ifelse(!TRUE%in%grepl("OSL", structure_RLum(object)$recordType),
                            "IRSL","OSL")


# Rejection criteria ------------------------------------------------------

  #Set rejection criteria
  if(missing(rejection.criteria)){

    rejection.criteria <- list(
      recycling.ratio = 10, recuperation.rate = 10, palaeodose.error = 10, exceed.max.regpoint = FALSE)

  }else{

    ##recycling ratio
    temp.recycling.ratio <- if("recycling.ratio" %in% names(rejection.criteria)) {

      rejection.criteria$recycling.ratio

    } else {10}

    ##recuperation rate
    temp.recuperation.rate <- if("recuperation.rate" %in% names(rejection.criteria)) {

      rejection.criteria$recuperation.rate

    } else {10}

    ##paleaodose error
    temp.palaeodose.error <- if("palaeodose.error" %in% names(rejection.criteria)) {

      rejection.criteria$palaeodose.error

    } else {10}

    ##exceed.max.regpoint
    temp.exceed.max.regpoint <- if("exceed.max.regpoint" %in% names(rejection.criteria)) {

      rejection.criteria$exceed.max.regpoint

    } else {FALSE}



    ##combine
    rejection.criteria <- list(
      recycling.ratio = temp.recycling.ratio,
      recuperation.rate = temp.recuperation.rate,
      palaeodose.error = temp.palaeodose.error,
      exceed.max.regpoint = temp.exceed.max.regpoint)

    ##remove objects
    rm(temp.recycling.ratio,temp.recuperation.rate, temp.palaeodose.error )

  }

# Deal with extra arguments ----------------------------------------------------

  ##deal with addition arguments
  extraArgs <- list(...)

  main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
  {""}

  log <- if("log" %in% names(extraArgs)) {extraArgs$log} else
  {""}

  cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else
  {1}

  background.count.distribution <-
    if ("background.count.distribution" %in% names(extraArgs)) {
      extraArgs$background.count.distribution
    } else
    {
      "non-poisson"
    }

  sigmab <- if("sigmab" %in% names(extraArgs)) {extraArgs$sigmab} else
  {NULL}


# Protocol Integrity Checks --------------------------------------------------

  ##check overall structur of the object
  ##every SAR protocol has to have equal number of curves


  ##grep curve types from analysis value and remove unwanted information
  temp.ltype <- sapply(1:length(object@records), function(x) {

                ##export as global variable
                object@records[[x]]@recordType <<- gsub(" .*", "",
                                                        object@records[[x]]@recordType)

                object@records[[x]]@recordType

  })


  ##problem: FI lexsyg devices provide irradiation information in a separate curve
  if("irradiation"%in%temp.ltype){

    ##grep irraditation times
    temp.irradiation <- structure_RLum(object)
    temp.irradiation <- temp.irradiation[temp.irradiation$recordType == "irradiation",
                                         "x.max"]

    ##remove every 2nd entry (test dose) and add "0" dose for natural signal
    temp.Dose <- c(0,temp.irradiation)

    ##remove irradiation entries from file
    object <- set_RLum(
               class = "RLum.Analysis",
               records = get_RLum(object, recordType = c(CWcurve.type, "TL")),
               protocol = "SAR")

  }

  ##check if the wanted curves are a multiple of two
  ##gsub removes unwanted information from the curves
  if(table(temp.ltype)[CWcurve.type]%%2!=0){
    error.list[[1]] <- "[analyse_SAR.CWOSL()] Input OSL/IRSL curves are not a multiple of two."
  }

  ##check if the curve lengths differ
  temp.matrix.length <- unlist(sapply(1:length(object@records), function(x) {
                          if(object@records[[x]]@recordType==CWcurve.type){
                              length(object@records[[x]]@data[,1])
                          }
  }))

  if(length(unique(temp.matrix.length))!=1){
    error.list[[2]] <- "[analyse_SAR.CWOSL()] Input curves lengths differ."

  }

  ##just proceed if error list is empty
  if (length(error.list) == 0) {

    ##check background integral
    if (max(signal.integral) == min(signal.integral)) {
      signal.integral <-
        c(min(signal.integral) : (max(signal.integral) + 1))

      warning("[analyse_SAR.CWOSL()] integral signal limits cannot be equal, reset automatically!")

    }


    ##background integral should not longer than curve channel length
    if (max(background.integral) == min(background.integral)) {
      background.integral <-
        c((min(background.integral) - 1) : max(background.integral))

    }

    if (max(background.integral) > temp.matrix.length[1]) {
      background.integral <-
          c((temp.matrix.length[1] - length(background.integral)):temp.matrix.length[1])

      ##prevent that the background integral becomes negative
      if(min(background.integral) < max(signal.integral)){
        background.integral <- c((max(signal.integral) + 1):max(background.integral))

      }

      warning(
        "[analyse_SAR.CWOSL()] Background integral out of bounds. Set to: c(",
        min(background.integral),":", max(background.integral),")"
      )

    }

    ##Do the same for the Tx-if set
    if (!is.null(background.integral.Tx)) {

      if (max(background.integral.Tx) == min(background.integral.Tx)) {
        background.integral.Tx <-
          c((min(background.integral.Tx) - 1) : max(background.integral.Tx))

      }

      if (max(background.integral.Tx) > temp.matrix.length[2]) {
        background.integral.Tx <-
          c((temp.matrix.length[2] - length(background.integral.Tx)):temp.matrix.length[2])


        ##prevent that the background integral becomes negative
        if (min(background.integral.Tx) < max(signal.integral.Tx)) {
          background.integral.Tx <-
            c((max(signal.integral.Tx) + 1):max(background.integral.Tx))


        }

        warning(
          "Background integral for Tx out of bounds. Set to: c(",
          min(background.integral.Tx),
          ":",
          max(background.integral.Tx),
          ")"
        )

      }
    }


    # Grep Curves -------------------------------------------------------------

    ##grep relevant curves from RLum.Analyis object
    OSL.Curves.ID <-
      get_RLum(object, recordType = CWcurve.type, get.index = TRUE)

    ##separate curves by Lx and Tx (it makes it much easier)
    OSL.Curves.ID.Lx <-
      OSL.Curves.ID[seq(1,length(OSL.Curves.ID),by = 2)]
    OSL.Curves.ID.Tx <-
      OSL.Curves.ID[seq(2,length(OSL.Curves.ID),by = 2)]

    ##get index of TL curves
    TL.Curves.ID <-
      suppressWarnings(get_RLum(object, recordType = "TL", get.index = TRUE))

    ##separate TL curves
    TL.Curves.ID.Lx <-
      sapply(1:length(OSL.Curves.ID.Lx), function(x) {
        TL.Curves.ID[which(TL.Curves.ID == (OSL.Curves.ID.Lx[x] - 1))]
      })

    TL.Curves.ID.Tx <-
      sapply(1:length(OSL.Curves.ID.Tx), function(x) {
        TL.Curves.ID[which(TL.Curves.ID == (OSL.Curves.ID.Tx[x] - 1))]
      })


    # COMPONENT FITTING -------------------------------------------------------


    # for(x in seq(1,length(OSL.Curves.ID),by=2)){
    #
    #
    #   temp.fit.output <- fit_CWCurve(object@records[[OSL.Curves.ID[x]]],
    #                 n.components.max=3,
    #                 output.terminal = FALSE,
    #                 output.terminalAdvanced = FALSE,
    #                 plot = FALSE
    #
    #               )
    #   if(exists("fit.output") == FALSE){
    #
    #     fit.output <- get_RLum(temp.fit.output)
    #
    #   }else{
    #
    #     fit.output <- rbind(fit.output, get_RLum(temp.fit.output))
    #
    #   }
    #
    # }

    ##TODO

    # Calculate LnLxTnTx values  --------------------------------------------------

    ##calculate LxTx values using external function
    LnLxTnTx <- lapply(seq(1,length(OSL.Curves.ID),by = 2), function(x){
      temp.LnLxTnTx <- get_RLum(
        calc_OSLLxTxRatio(
          Lx.data = object@records[[OSL.Curves.ID[x]]]@data,
          Tx.data = object@records[[OSL.Curves.ID[x + 1]]]@data,
          signal.integral = signal.integral,
          signal.integral.Tx = signal.integral.Tx,
          background.integral = background.integral,
          background.integral.Tx = background.integral.Tx,
          background.count.distribution = background.count.distribution,
          sigmab = sigmab
        )
      )

        ##grep dose
        if (exists("temp.irradiation") == FALSE) {
          temp.Dose <- object@records[[OSL.Curves.ID[x]]]@info$IRR_TIME

          ##for the case that no information on the dose can be found
          if (is.null(temp.Dose)) {
            temp.Dose <- NA
          }

          temp.LnLxTnTx <-
            cbind(Dose = temp.Dose, temp.LnLxTnTx)

        }else{
          temp.LnLxTnTx <- cbind(Dose = temp.Dose[x], temp.LnLxTnTx)

        }
      })

    ##combine
    LnLxTnTx <- data.table::rbindlist(LnLxTnTx)

    # Set regeneration points -------------------------------------------------

    ##overwrite dose point manually
    if (!is.null(dose.points)) {
      if (length(dose.points) != length(LnLxTnTx$Dose)) {
        stop("[analyse_SAR.CWOSL()] length 'dose.points' differs from number of curves.")

      }

      LnLxTnTx$Dose <- dose.points

    }

    ##check whether we have dose points at all
    if (is.null(dose.points) & anyNA(LnLxTnTx$Dose)) {
      stop("[analyse_SAR.CWOSL()] 'dose.points' contains NA values or have not been set!")

    }

    #generate unique dose id - this are also the # for the generated points
    temp.DoseID <- c(0:(length(LnLxTnTx$Dose) - 1))
    temp.DoseName <- paste("R",temp.DoseID,sep = "")
    temp.DoseName <-
      cbind(Name = temp.DoseName,Dose = LnLxTnTx$Dose)


    ##set natural
    temp.DoseName[temp.DoseName[,"Name"] == "R0","Name"] <-
      "Natural"

    ##set R0
    temp.DoseName[temp.DoseName[,"Name"] != "Natural" &
                    temp.DoseName[,"Dose"] == 0,"Name"] <- "R0"

    ##correct numeration numeration of other dose points

    ##how many dose points do we have with 0?
    non.temp.zero.dose.number <- nrow(temp.DoseName[temp.DoseName[, "Dose"] != 0,])

    temp.DoseName[temp.DoseName[,"Name"] != "Natural" &
                    temp.DoseName[,"Name"] != "R0","Name"] <- paste("R",c(1:non.temp.zero.dose.number),sep =
                                                                      "")

    ##find duplicated doses (including 0 dose - which means the Natural)
    temp.DoseDuplicated <- duplicated(temp.DoseName[,"Dose"])

    ##combine temp.DoseName
    temp.DoseName <-
      cbind(temp.DoseName,Repeated = temp.DoseDuplicated)

    ##correct value for R0 (it is not really repeated)
    temp.DoseName[temp.DoseName[,"Dose"] == 0,"Repeated"] <- FALSE

    ##combine in the data frame
    temp.LnLxTnTx <- data.frame(Name = temp.DoseName[,"Name"],
                                Repeated = as.logical(temp.DoseName[,"Repeated"]))

    LnLxTnTx <- cbind(temp.LnLxTnTx,LnLxTnTx)
    LnLxTnTx[,"Name"] <- as.character(LnLxTnTx[,"Name"])

    # Calculate Recycling Ratio -----------------------------------------------

    ##Calculate Recycling Ratio

    if (length(LnLxTnTx[LnLxTnTx[,"Repeated"] == TRUE,"Repeated"]) > 0) {
      ##identify repeated doses
      temp.Repeated <-
        LnLxTnTx[LnLxTnTx[,"Repeated"] == TRUE,c("Name","Dose","LxTx")]

      ##find concering previous dose for the repeated dose
      temp.Previous <-
        t(sapply(1:length(temp.Repeated[,1]),function(x) {
          LnLxTnTx[LnLxTnTx[,"Dose"] == temp.Repeated[x,"Dose"] &
                     LnLxTnTx[,"Repeated"] == FALSE,c("Name","Dose","LxTx")]
        }))

      ##convert to data.frame
      temp.Previous <- as.data.frame(temp.Previous)

      ##set column names
      temp.ColNames <-
        sapply(1:length(temp.Repeated[,1]),function(x) {
          paste("Recycling ratio (", temp.Repeated[x,"Name"],"/",
                temp.Previous[temp.Previous[,"Dose"] == temp.Repeated[x,"Dose"],"Name"],
                ")",
                sep = "")
        })

      ##Calculate Recycling Ratio
      RecyclingRatio <-
        round(as.numeric(temp.Repeated[,"LxTx"]) / as.numeric(temp.Previous[,"LxTx"]),
              digits = 4)

      ##Just transform the matrix and add column names
      RecyclingRatio <- t(RecyclingRatio)
      colnames(RecyclingRatio) <- temp.ColNames

    }else{
      RecyclingRatio <- NA
    }



    # Calculate Recuperation Rate ---------------------------------------------


    ##Recuperation Rate (capable to handle multiple type of recuperation values)

    if (length(LnLxTnTx[LnLxTnTx[,"Name"] == "R0","Name"]) > 0) {
      Recuperation <-
        sapply(1:length(LnLxTnTx[LnLxTnTx[,"Name"] == "R0","Name"]),
               function(x) {
                 round(LnLxTnTx[LnLxTnTx[,"Name"] == "R0","LxTx"][x] /
                         LnLxTnTx[LnLxTnTx[,"Name"] == "Natural","LxTx"],
                       digits = 4)
               })
      ##Just transform the matrix and add column names
      Recuperation  <-  t(Recuperation)
      colnames(Recuperation)  <-
        unlist(strsplit(paste(
          "Recuperation rate",
          1:length(LnLxTnTx[LnLxTnTx[,"Name"] == "R0","Name"]), collapse = ";"
        ), ";"))

    }else{
      Recuperation <- NA
    }


    # Evaluate and Combine Rejection Criteria ---------------------------------


    temp.criteria <- c(ifelse(is.null(colnames(RecyclingRatio)), NA,colnames(RecyclingRatio)),
                       ifelse(is.null(colnames(Recuperation)), NA,colnames(Recuperation)))

    temp.value <- c(RecyclingRatio,Recuperation)
    temp.threshold <-
      c(rep(
        rejection.criteria$recycling.ratio / 100, length(RecyclingRatio)
      ),
      rep(
        paste("", rejection.criteria$recuperation.rate / 100),
        length(Recuperation)
      ))

    ##RecyclingRatio
    if (is.na(RecyclingRatio)[1] == FALSE) {
      temp.status.RecyclingRatio <-
        sapply(1:length(RecyclingRatio), function(x) {
          if (abs(1 - RecyclingRatio[x]) > (rejection.criteria$recycling.ratio / 100)) {
            "FAILED"
          }else{
            "OK"
          }
        })
    }else{
      temp.status.RecyclingRatio <- "OK"

    }

    ##Recuperation
    if (is.na(Recuperation)[1] == FALSE) {
      temp.status.Recuperation  <-
        sapply(1:length(Recuperation), function(x) {
          ifelse(Recuperation[x] > rejection.criteria$recuperation.rate,
                 "FAILED", "OK")


        })

    }else{
      temp.status.Recuperation <- "OK"

    }

    RejectionCriteria <- data.frame(
      Criteria = temp.criteria,
      Value = temp.value,
      Threshold = temp.threshold,
      Status = c(temp.status.RecyclingRatio,temp.status.Recuperation),
      stringsAsFactors = FALSE
    )

    ##============================================================================##
    ##PLOTTING
    ##============================================================================##

    if (plot == TRUE) {
      # Plotting - Config -------------------------------------------------------

      ##colours and double for plotting
      col <- get("col", pos = .LuminescenceEnv)

      if (plot.single[1] == FALSE) {
        ## read par settings
        par.default <- par(no.readonly = TRUE)

        layout(matrix(
          c(1,1,3,3,
            1,1,3,3,
            2,2,4,4,
            2,2,4,4,
            5,5,5,5),5,4,byrow = TRUE
        ))

        par(
          oma = c(0,0,0,0), mar = c(4,4,3,3), cex = cex * 0.6
        )

        ## 1 -> TL previous LnLx
        ## 2 -> LnLx
        ## 3 -> TL previous TnTx
        ## 4 -> TnTx
        ## 5 -> Legend

        ## set selected curves to allow plotting of all curves
        plot.single.sel <- c(1,2,3,4,5,6,7,8)

      }else{
        ##check for values in the single output of the function and convert
        if (!is(plot.single, "logical")) {
          if (!is(plot.single, "numeric")) {
            stop("[analyse_SAR.CWOSL()] Invalid data type for 'plot.single'.")
          }

          plot.single.sel  <- plot.single

        }else{
          plot.single.sel <- c(1,2,3,4,5,6,7,8)

        }

      }


      ##warning if number of curves exceed colour values
      if (length(col) < length(OSL.Curves.ID) / 2) {
        temp.message  <-
          paste(
            "\n[analyse_SAR.CWOSL()] To many curves! Only the first",
            length(col),"curves are plotted!"
          )
        warning(temp.message)
      }

      ##legend text
      legend.text <-
        paste(LnLxTnTx$Name,"\n(",LnLxTnTx$Dose,")", sep = "")


      ##get channel resolution (should be equal for all curves)
      resolution.OSLCurves <- round(object@records[[OSL.Curves.ID[1]]]@data[2,1] -
                                      object@records[[OSL.Curves.ID[1]]]@data[1,1],
                                    digits = 2)


      # Plotting TL Curves previous LnLx ----------------------------------------

      ##overall plot option selection for plot.single.sel
      if (1 %in% plot.single.sel) {
        ##check if TL curves are available
        if (length(TL.Curves.ID.Lx[[1]] > 0)) {
          ##It is just an approximation taken from the data
          resolution.TLCurves <-  round(mean(diff(
            round(object@records[[TL.Curves.ID.Lx[1]]]@data[,1], digits = 1)
          )), digits = 1)

          ylim.range <-
            sapply(seq(1,length(TL.Curves.ID.Lx),by = 1) ,function(x) {
              range(object@records[[TL.Curves.ID.Lx[x]]]@data[,2])

            })

          plot(
            NA,NA,
            xlab = "T [\u00B0C]",
            ylab = paste("TL [cts/",resolution.TLCurves," \u00B0C]",sep =
                           ""),
            xlim = c(object@records[[TL.Curves.ID.Lx[1]]]@data[1,1],
                     max(object@records[[TL.Curves.ID.Lx[1]]]@data[,1])),
            ylim = c(1,max(ylim.range)),
            main = main,
            log = if (log == "y" | log == "xy") {
              "y"
            }else{
              ""
            }
          )

          #provide curve information as mtext, to keep the space for the header
          mtext(side = 3,
                expression(paste(
                  "TL previous ", L[n],",",L[x]," curves",sep = ""
                )),
                cex = cex * 0.7)

          ##plot TL curves
          sapply(1:length(TL.Curves.ID.Lx) ,function(x) {
            lines(object@records[[TL.Curves.ID.Lx[x]]]@data,col = col[x])

          })



        }else{
          plot(
            NA,NA,xlim = c(0,1), ylim = c(0,1), main = "",
            axes = FALSE,
            ylab = "",
            xlab = ""
          )
          text(0.5,0.5, "No TL curve detected")

        }
      }#plot.single.sel

      # Plotting LnLx Curves ----------------------------------------------------

      ##overall plot option selection for plot.single.sel
      if (2 %in% plot.single.sel) {
        ylim.range <- sapply(1:length(OSL.Curves.ID.Lx) ,function(x) {
          range(object@records[[OSL.Curves.ID.Lx[x]]]@data[,2])
        })

        if((log == "x" | log == "xy") & object@records[[OSL.Curves.ID.Lx[[1]]]]@data[1,1] == 0){
          xlim <- c(object@records[[OSL.Curves.ID.Lx[1]]]@data[2,1],
                    max(object@records[[OSL.Curves.ID.Lx[1]]]@data[,1]) +
                      object@records[[OSL.Curves.ID.Lx[1]]]@data[2,1])


        }else{


        xlim  <- c(object@records[[OSL.Curves.ID.Lx[1]]]@data[1,1],
                   max(object@records[[OSL.Curves.ID.Lx[1]]]@data[,1]))

        }
        #open plot area LnLx
        plot(
          NA,NA,
          xlab = "Time [s]",
          ylab = paste(CWcurve.type," [cts/",resolution.OSLCurves," s]",sep =
                         ""),
          xlim = xlim,
          ylim = range(ylim.range),
          main = main,
          log = log
        )

        #provide curve information as mtext, to keep the space for the header
        mtext(side = 3, expression(paste(L[n],",",L[x]," curves",sep = "")),
              cex = cex * 0.7)

        ##plot curves
        sapply(1:length(OSL.Curves.ID.Lx), function(x) {

          if((log == "x" | log == "xy") & object@records[[OSL.Curves.ID.Lx[[x]]]]@data[1,1] == 0){
            object@records[[OSL.Curves.ID.Lx[[x]]]]@data[1,] <-
              object@records[[OSL.Curves.ID.Lx[[x]]]]@data[1,] +
              diff(c(object@records[[OSL.Curves.ID.Lx[[x]]]]@data[1,1],
                     object@records[[OSL.Curves.ID.Lx[[x]]]]@data[2,1]))

            warnings("[analyse_SAR.CWOSL()] curves shifted by one chanel for log-plot.")
          }

          lines(object@records[[OSL.Curves.ID.Lx[[x]]]]@data,col = col[x])

        })


        ##mark integration limit Lx curves
        abline(
          v = (object@records[[OSL.Curves.ID.Lx[1]]]@data[min(signal.integral),1]), lty =
            2, col = "gray"
        )
        abline(
          v = (object@records[[OSL.Curves.ID.Lx[1]]]@data[max(signal.integral),1]), lty =
            2, col = "gray"
        )
        abline(
          v = (object@records[[OSL.Curves.ID.Lx[1]]]@data[min(background.integral),1]), lty =
            2, col = "gray"
        )
        abline(
          v = (object@records[[OSL.Curves.ID.Lx[1]]]@data[max(background.integral),1]), lty =
            2, col = "gray"
        )

        ##mtext, implemented here, as a plot window has to be called first
        if (missing(mtext.outer)) {
          mtext.outer  <- ""
        }
        mtext(
          mtext.outer, side = 4, outer = TRUE, line = -1.7, cex = cex, col = "blue"
        )

      }# plot.single.sel

      # Plotting TL Curves previous TnTx ----------------------------------------

      ##overall plot option selection for plot.single.sel
      if (3 %in% plot.single.sel) {
        ##check if TL curves are available
        if (length(TL.Curves.ID.Tx[[1]] > 0)) {
          ##It is just an approximation taken from the data
          resolution.TLCurves <-  round(mean(diff(
            round(object@records[[TL.Curves.ID.Tx[1]]]@data[,1], digits = 1)
          )), digits = 1)


          ylim.range <- sapply(1:length(TL.Curves.ID.Tx) ,function(x) {
            range(object@records[[TL.Curves.ID.Tx[x]]]@data[,2])

          })



          plot(
            NA,NA,
            xlab = "T [\u00B0C]",
            ylab = paste("TL [cts/",resolution.TLCurves," \u00B0C]",sep = ""),
            xlim = c(object@records[[TL.Curves.ID.Tx[1]]]@data[1,1],
                     max(object@records[[TL.Curves.ID.Tx[1]]]@data[,1])),
            ylim = c(1,max(ylim.range)),
            main = main,
            log = if (log == "y" | log == "xy") {
              "y"
            }else{
              ""
            }
          )

          #provide curve information as mtext, to keep the space for the header
          mtext(side = 3,
                expression(paste(
                  "TL previous ", T[n],",",T[x]," curves",sep = ""
                )),
                cex = cex * 0.7)

          ##plot TL curves
          sapply(1:length(TL.Curves.ID.Tx) ,function(x) {
            lines(object@records[[TL.Curves.ID.Tx[x]]]@data,col = col[x])

          })



        }else{
          plot(
            NA,NA,xlim = c(0,1), ylim = c(0,1), main = "",
            axes = FALSE,
            ylab = "",
            xlab = ""
          )
          text(0.5,0.5, "No TL curve detected")

        }

      }#plot.single.sel

      # Plotting TnTx Curves ----------------------------------------------------

      ##overall plot option selection for plot.single.sel
      if (4 %in% plot.single.sel) {
        ylim.range <- sapply(1:length(OSL.Curves.ID.Tx) ,function(x) {
          range(object@records[[OSL.Curves.ID.Tx[x]]]@data[,2])

        })

        if((log == "x" | log == "xy") & object@records[[OSL.Curves.ID.Tx[[1]]]]@data[1,1] == 0){
          xlim <- c(object@records[[OSL.Curves.ID.Tx[1]]]@data[2,1],
                    max(object@records[[OSL.Curves.ID.Tx[1]]]@data[,1]) +
                      object@records[[OSL.Curves.ID.Tx[1]]]@data[2,1])


        }else{
          xlim <- c(object@records[[OSL.Curves.ID.Tx[1]]]@data[1,1],
                    max(object@records[[OSL.Curves.ID.Tx[1]]]@data[,1]))
        }

        #open plot area LnLx
        plot(
          NA,NA,
          xlab = "Time [s]",
          ylab = paste(CWcurve.type ," [cts/",resolution.OSLCurves," s]",sep =
                         ""),
          xlim = xlim,
          ylim = range(ylim.range),
          main = main,
          log = log
        )

        #provide curve information as mtext, to keep the space for the header
        mtext(side = 3,
              expression(paste(T[n],",",T[x]," curves",sep = "")),
              cex = cex * 0.7)

        ##plot curves and get legend values
        sapply(1:length(OSL.Curves.ID.Tx) ,function(x) {

          ##account for log-scale and 0 values
          if((log == "x" | log == "xy") & object@records[[OSL.Curves.ID.Tx[[x]]]]@data[1,1] == 0){
            object@records[[OSL.Curves.ID.Tx[[x]]]]@data[1,] <-
              object@records[[OSL.Curves.ID.Tx[[x]]]]@data[1,] +
                 diff(c(object@records[[OSL.Curves.ID.Tx[[x]]]]@data[1,1],
                      object@records[[OSL.Curves.ID.Tx[[x]]]]@data[2,1]))

            warnings("[analyse_SAR.CWOSL()] curves shifted by one chanel for log-plot.")

          }

          lines(object@records[[OSL.Curves.ID.Tx[[x]]]]@data,col = col[x])

        })

        ##mark integration limit Tx curves
        abline(
          v = (object@records[[OSL.Curves.ID.Tx[1]]]@data[min(signal.integral),1]), lty =
            2, col = "gray"
        )
        abline(
          v = (object@records[[OSL.Curves.ID.Tx[1]]]@data[max(signal.integral),1]), lty =
            2, col = "gray"
        )
        abline(
          v = (object@records[[OSL.Curves.ID.Tx[1]]]@data[min(background.integral),1]), lty =
            2, col = "gray"
        )
        abline(
          v = (object@records[[OSL.Curves.ID.Tx[1]]]@data[max(background.integral),1]), lty =
            2, col = "gray"
        )

      }# plot.single.sel

      # Plotting Legend ----------------------------------------

      ##overall plot option selection for plot.single.sel
      if (5 %in% plot.single.sel) {
        par.margin  <- par()$mar
        par.mai  <- par()$mai
        par(mar = c(1,1,1,1), mai = c(0,0,0,0))

        plot(
          c(1:(length(
            OSL.Curves.ID
          ) / 2)),
          rep(7,length(OSL.Curves.ID) / 2),
          type = "p",
          axes = FALSE,
          xlab = "",
          ylab = "",
          pch = 20,
          col = unique(col[1:length(OSL.Curves.ID)]),
          cex = 4 * cex,
          ylim = c(0,10)
        )

        ##add text
        text(c(1:(length(
          OSL.Curves.ID
        ) / 2)),
        rep(7,length(OSL.Curves.ID) / 2),
        legend.text,
        offset = 1,
        pos = 1)


        ##add line
        abline(h = 10,lwd = 0.5)

        ##set failed text and mark De as failed
        if (length(grep("FAILED", RejectionCriteria$status)) > 0) {
          mtext("[FAILED]", col = "red")


        }

        #reset margin
        par(mar = par.margin, mai = par.mai)

      }#plot.single.sel

      if (exists("par.default")) {
        par(par.default)

      }


    }##end plot == TRUE


    # Plotting  GC  ----------------------------------------

    temp.sample <- data.frame(
      Dose = LnLxTnTx$Dose,
      LxTx = LnLxTnTx$LxTx,
      LxTx.Error = LnLxTnTx$LxTx.Error,
      TnTx = LnLxTnTx$Net_TnTx
    )

    ##overall plot option selection for plot.single.sel
    if (plot == TRUE && 6 %in% plot.single.sel) {
      plot  <-  TRUE

    }else {
      plot  <- FALSE

    }

    ##Fit and plot growth curve
    temp.GC <- plot_GrowthCurve(temp.sample,
                                output.plot = plot,
                                ...)

    ##grep informaton on the fit object
    temp.GC.fit.Formula  <- get_RLum(temp.GC, "Formula")

    ##grep results
    temp.GC <- get_RLum(temp.GC)

    # Provide Rejection Criteria for Palaedose error --------------------------

    palaeodose.error.calculated <- ifelse(is.na(temp.GC[,1]) == FALSE,
                                          round(temp.GC[,2] / temp.GC[,1], digits = 5),
                                          NA)

    palaeodose.error.threshold <-
      rejection.criteria$palaeodose.error / 100

    if (is.na(palaeodose.error.calculated)) {
      palaeodose.error.status <- "FAILED"

    }else{
      palaeodose.error.status <- ifelse(
        palaeodose.error.calculated <= rejection.criteria$palaeodose.error,
        "OK", "FAILED"
      )

    }


    palaeodose.error.data.frame <- data.frame(
      Criteria = "Palaeodose error",
      Value = palaeodose.error.calculated,
      Threshold = palaeodose.error.threshold,
      Status =  palaeodose.error.status,
      stringsAsFactors = FALSE
    )

    ##add exceed.max.regpoint
    if (!is.na(temp.GC[,1])) {
      status.exceed.max.regpoint <-
        ifelse(max(LnLxTnTx$Dose) < temp.GC[,1], "FAILED", "OK")

    }else{
      status.exceed.max.regpoint <- "OK"

    }

    exceed.max.regpoint.data.frame <- data.frame(
      Criteria = "De > max. dose point",
      Value = as.numeric(temp.GC[,1]),
      Threshold = as.numeric(max(LnLxTnTx$Dose)),
      Status =  status.exceed.max.regpoint
    )


    ##add to RejectionCriteria data.frame
    RejectionCriteria <- rbind(RejectionCriteria,
                               palaeodose.error.data.frame,
                               exceed.max.regpoint.data.frame)


    ##add recjection status
    if (length(grep("FAILED",RejectionCriteria$Status)) > 0) {
      temp.GC <- data.frame(temp.GC, RC.Status = "FAILED")


    }else{
      temp.GC <- data.frame(temp.GC, RC.Status = "OK")


    }


    ##add information on the integration limits
    temp.GC.extened <-
      data.frame(
        signal.range = paste(min(signal.integral),":",
                             max(signal.integral)),
        background.range = paste(min(background.integral),":",
                                 max(background.integral)),
        signal.range.Tx = paste(min(ifelse(is.null(signal.integral.Tx),NA,signal.integral.Tx)),":",
                                max(ifelse(is.null(signal.integral.Tx),NA,signal.integral.Tx))),
        background.range.Tx = paste(min(ifelse(is.null(background.integral.Tx), NA,background.integral.Tx)) ,":",
                                    max(ifelse(is.null(background.integral.Tx), NA,background.integral.Tx))),
        stringsAsFactors = FALSE
      )



    # Set return Values -----------------------------------------------------------

    ##generate unique identifier
    UID <- .create_UID()

    temp.results.final <- set_RLum(
      class = "RLum.Results",
      data = list(
        De.values = as.data.frame(c(temp.GC, temp.GC.extened, UID = UID), stringsAsFactors = FALSE),
        LnLxTnTx.table = cbind(LnLxTnTx, UID = UID, stringsAsFactors = FALSE),
        rejection.criteria = cbind(RejectionCriteria, UID, stringsAsFactors = FALSE),
        Formula = temp.GC.fit.Formula,
        call = sys.call()
      )
    )


    # Plot graphical interpretation of rejection criteria -----------------------------------------

    if (plot == TRUE && 7 %in% plot.single.sel) {
      ##set graphical parameter
      if (!plot.single) {
        par(mfrow = c(1,2))
      }else{
        par(mfrow = c(1,1))
      }


      ##Rejection criteria
      temp.rejection.criteria <- get_RLum(temp.results.final,
                                          data.object = "rejection.criteria")

      temp.rc.reycling.ratio <- temp.rejection.criteria[grep("Recycling ratio",temp.rejection.criteria[,"Criteria"]),]

      temp.rc.recuperation.rate <- temp.rejection.criteria[grep("Recuperation rate",temp.rejection.criteria[,"Criteria"]),]

      temp.rc.palaedose.error <- temp.rejection.criteria[grep("Palaeodose error",temp.rejection.criteria[,"Criteria"]),]

      plot(
        NA,NA,
        xlim = c(-0.5,0.5),
        ylim = c(0,30),
        yaxt = "n", ylab = "",
        xaxt = "n", xlab = "",
        bty = "n",
        main = "Rejection criteria"
      )

      axis(
        side = 1, at = c(-0.2,-0.1,0,0.1,0.2), labels = c("- 0.2", "- 0.1","0/1","+ 0.1", "+ 0.2")
      )

      ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
      ##polygon for recycling ratio
      text(
        x = 0, y = 30.5, "Recycling ratio", pos = 1, srt = 0
      )
      polygon(
        x = c(
          -as.numeric(as.character(temp.rc.reycling.ratio$Threshold))[1],-as.numeric(as.character(temp.rc.reycling.ratio$Threshold))[1],
          as.numeric(as.character(temp.rc.reycling.ratio$Threshold))[1],
          as.numeric(as.character(temp.rc.reycling.ratio$Threshold))[1]
        ),
        y = c(21,29,29,21), col = "gray", border = NA
      )
      polygon(x = c(-0.3,-0.3,0.3,0.3) , y = c(21,29,29,21))


      ##consider possibility of multiple pIRIR signals and multiple recycling ratios
      if (nrow(temp.rc.recuperation.rate) > 0) {
        col.id  <- 1
        for (i in seq(1,nrow(temp.rc.recuperation.rate),
                      length(unique(temp.rc.recuperation.rate[,"Criteria"])))) {
          for (j in 0:length(unique(temp.rc.recuperation.rate[,"Criteria"]))) {
            points(
              temp.rc.reycling.ratio[i + j, "Value"] - 1,
              y = 25,
              pch = col.id,
              col = col.id,
              cex = 1.3 * cex
            )

          }
          col.id <- col.id + 1
        }
        rm(col.id)

        ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
        ##polygon for recuperation rate
        text(
          x = 0, y = 20.5, "Recuperation rate", pos = 1, srt = 0
        )
        polygon(
          x = c(
            0,
            0,
            as.numeric(as.character(
              temp.rc.recuperation.rate$Threshold
            ))[1],
            as.numeric(as.character(
              temp.rc.recuperation.rate$Threshold
            ))[1]
          ),
          y = c(11,19,19,11), col = "gray", border = NA
        )

        polygon(x = c(-0.3,-0.3,0.3,0.3) , y = c(11,19,19,11))
        polygon(
          x = c(-0.3,-0.3,0,0) , y = c(11,19,19,11), border = NA, density = 10, angle = 45
        )

        for (i in 1:nrow(temp.rc.recuperation.rate)) {
          points(
            temp.rc.palaedose.error[i, "Value"],
            y = 15,
            pch = i,
            col = i,
            cex = 1.3 * cex
          )

        }
      }

      ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
      ##polygon for palaeodose error
      text(
        x = 0, y = 10.5, "Palaeodose error", pos = 1, srt = 0
      )
      polygon(
        x = c(
          0,
          0,
          as.numeric(as.character(temp.rc.palaedose.error$Threshold))[1],
          as.numeric(as.character(temp.rc.palaedose.error$Threshold))[1]
        ),
        y = c(1,9,9,1), col = "gray", border = NA
      )
      polygon(x = c(-0.3,-0.3,0.3,0.3) , y = c(1,9,9,1))
      polygon(
        x = c(-0.3,-0.3,0,0) , y = c(1,9,9,1), border = NA, density = 10, angle = 45
      )


      for (i in 1:nrow(temp.rc.palaedose.error)) {
        points(
          temp.rc.palaedose.error[i, "Value"],
          y = 5,
          pch = i,
          col = i,
          cex = 1.3 * cex
        )
      }
    }


    if (plot == TRUE && 8 %in% plot.single.sel) {
      ##graphical represenation of IR-curve
      temp.IRSL <- suppressWarnings(get_RLum(object, recordType = "IRSL"))

      if(length(temp.IRSL) != 0){
        plot_RLum.Data.Curve(temp.IRSL, par.local = FALSE)

      }else{
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(x = c(1,1), y = c(1, 1), labels = "No IRSL curve detected!")

      }

    }


    ##It is doubled in this function, but the par settings need some more careful considerations ...
    if (exists("par.default")) {
      par(par.default)
      rm(par.default)
    }



    # Return --------------------------------------------------------------------------------------
    invisible(temp.results.final)

  }else{
    warning(paste0(
      "\n",
      paste(unlist(error.list), collapse = "\n"),"\n... >> nothing was done here!"
    ), call. = FALSE)
    invisible(NULL)

  }

}
