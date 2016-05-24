#' Calculate the Lx/Tx ratio for a given set of TL curves [beta version]
#'
#' Calculate Lx/Tx ratio for a given set of TL curves.
#'
#' -
#'
#' @param Lx.data.signal \code{\linkS4class{RLum.Data.Curve}} or
#' \code{\link{data.frame}} (\bold{required}): TL data (x =
#' temperature, y = counts) (TL signal)
#'
#' @param Lx.data.background \code{\linkS4class{RLum.Data.Curve}} or
#' \code{\link{data.frame}} (optional): TL data (x =
#' temperature, y = counts). If no data are provided no background subtraction
#' is performed.
#'
#' @param Tx.data.signal \code{\linkS4class{RLum.Data.Curve}} or
#' \code{\link{data.frame}} (\bold{required}): TL data (x =
#' temperature, y = counts) (TL test signal)
#'
#' @param Tx.data.background \code{\linkS4class{RLum.Data.Curve}} or
#' \code{\link{data.frame}} (optional): TL data (x =
#' temperature, y = counts). If no data are provided no background subtraction
#' is performed.
#'
#' @param signal.integral.min \code{\link{integer}} (\bold{required}): channel number
#' for the lower signal integral bound (e.g. \code{signal.integral.min = 100})
#'
#' @param signal.integral.max \code{\link{integer}} (\bold{required}): channel number
#' for the upper signal integral bound (e.g. \code{signal.integral.max = 200})
#'
#' @return Returns an S4 object of type \code{\linkS4class{RLum.Results}}.
#' Slot \code{data} contains a \link{list} with the following structure:\cr\cr
#' $ LxTx.table \cr .. $ LnLx \cr .. $ LnLx.BG \cr .. $ TnTx \cr .. $ TnTx.BG
#' \cr .. $ Net_LnLx \cr .. $ Net_LnLx.Error\cr
#'
#' @note \bold{This function has still BETA status!}
#'
#' @section Function version: 0.3.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France), Christoph Schmidt, University of Bayreuth (Germany)
#'
#' @seealso \code{\linkS4class{RLum.Results}}, \code{\link{analyse_SAR.TL}}
#'
#' @references -
#'
#' @keywords datagen
#'
#' @examples
#'
#'
#' ##load package example data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##convert Risoe.BINfileData into a curve object
#' temp <- Risoe.BINfileData2RLum.Analysis(TL.SAR.Data, pos = 3)
#'
#'
#' Lx.data.signal <- get_RLum(temp, record.id=1)
#' Lx.data.background <- get_RLum(temp, record.id=2)
#' Tx.data.signal <- get_RLum(temp, record.id=3)
#' Tx.data.background <- get_RLum(temp, record.id=4)
#' signal.integral.min <- 210
#' signal.integral.max <- 230
#'
#' output <- calc_TLLxTxRatio(Lx.data.signal,
#'                            Lx.data.background,
#'                            Tx.data.signal, Tx.data.background,
#'                            signal.integral.min, signal.integral.max)
#' get_RLum(output)
#'
#' @export
calc_TLLxTxRatio <- function(
  Lx.data.signal,
  Lx.data.background,
  Tx.data.signal,
  Tx.data.background,
  signal.integral.min,
  signal.integral.max
){


  ##--------------------------------------------------------------------------##
  ##(1) - a few integrity check

     ##check for MISSING objects
     if(missing(Lx.data.signal) == TRUE | missing(Tx.data.signal) == TRUE |
        missing(signal.integral.min) == TRUE |  missing(signal.integral.max) == TRUE){

       temp.missing <- paste(
                       c(if(missing(Lx.data.signal)){"Lx.data.signal"},
                         if(missing(Tx.data.signal)){"Tx.data.signal"},
                         if(missing(signal.integral.min)){"signal.integral.min"},
                         if(missing(signal.integral.max)){"signal.integral.max"}),
                       collapse = ", ")

          stop(paste("[calc_TLLxTxRatio()] Arguments are missing: ",temp.missing, ".", sep=""))

     }


     ##check DATA TYPE differences
     if(is(Lx.data.signal)[1]!=is(Tx.data.signal)[1]){
       stop("[calc_TLLxTxRatio()] Data type of Lx and Tx data differs!")}

     ##check for allowed data.types
     if(!is(Lx.data.signal, "data.frame") &
        !is(Lx.data.signal, "RLum.Data.Curve")){

       stop("[calc_TLLxTxRatio()] Input data type for not allowed. Allowed are 'RLum.Data.Curve' and 'data.frame'")

     }

  ##--------------------------------------------------------------------------##
  ## Type conversion (assuming that all input variables are of the same type)

  if(is(Lx.data.signal, "RLum.Data.Curve")){

    Lx.data.signal <- as(Lx.data.signal, "matrix")
    Tx.data.signal <- as(Tx.data.signal, "matrix")

    if(missing(Lx.data.background) == FALSE && is.null(Lx.data.background) == FALSE){

      Lx.data.background <- as(Lx.data.background, "matrix")

    }

    if(missing(Tx.data.background) == FALSE && is.null(Tx.data.background) == FALSE){

      Tx.data.background <- as(Tx.data.background, "matrix")

    }

  }

  ##(d) - check if Lx and Tx curves have the same channel length
     if(length(Lx.data.signal[,2])!=length(Tx.data.signal[,2])){

       stop("[calc_TLLxTxRatio()] Channel number of Lx and Tx data differs!")}


   ##(e) - check if signal integral is valid
   if(signal.integral.min < 1 | signal.integral.max > length(Lx.data.signal[,2])){
     stop("[calc_TLLxTxRatio()] Signal.integral is not valid!")}




#  Background Consideration --------------------------------------------------


   ##Lx.data
   if(missing(Lx.data.background)==FALSE){

     LnLx.BG <- sum(Lx.data.background[signal.integral.min:signal.integral.max,2])

    }else{

     LnLx.BG <- NA

    }

   ##Tx.data
      if(missing(Tx.data.background)==FALSE){

        TnTx.BG <- sum(Tx.data.background[signal.integral.min:signal.integral.max,2])

      }else{

        TnTx.BG <- NA

      }

# Calculate Lx/Tx values --------------------------------------------------

    LnLx <- sum(Lx.data.signal[signal.integral.min:signal.integral.max,2])
    TnTx <- sum(Tx.data.signal[signal.integral.min:signal.integral.max,2])


     ##calculate variance of background
     if(is.na(LnLx.BG) == FALSE & is.na(TnTx.BG) == FALSE){

       BG.Error <- sd(c(LnLx.BG, TnTx.BG))
     }


    if(is.na(LnLx.BG) == FALSE){

      net_LnLx <-  LnLx - LnLx.BG
      net_LnLx.Error <- abs(net_LnLx * BG.Error/LnLx.BG)

    }else{

      net_LnLx <- NA
      net_LnLx.Error <- NA

    }

    if(is.na(TnTx.BG) == FALSE){

         net_TnTx <-  TnTx - TnTx.BG
         net_TnTx.Error <- abs(net_TnTx * BG.Error/TnTx.BG)

    }else{

      net_TnTx <- NA
      net_TnTx.Error  <- NA

    }


    if(is.na(net_TnTx) == TRUE){

      LxTx <- LnLx/TnTx
      LxTx.Error <- NA

    }else{

      LxTx <- net_LnLx/net_TnTx
      LxTx.Error <- LxTx*((net_LnLx.Error/net_LnLx) + (net_TnTx.Error/net_TnTx))


    }



    ##COMBINE to a data.frame
    temp.results <- data.frame(LnLx,
                               LnLx.BG,
                               TnTx,
                               TnTx.BG,
                               net_LnLx,
                               net_LnLx.Error,
                               net_TnTx,
                               net_TnTx.Error,
                               LxTx,
                               LxTx.Error)

# Return values -----------------------------------------------------------

   newRLumResults.calc_TLLxTxRatio <- set_RLum(
     class = "RLum.Results",
     data=list(LxTx.table = temp.results))

   return(newRLumResults.calc_TLLxTxRatio)

}
