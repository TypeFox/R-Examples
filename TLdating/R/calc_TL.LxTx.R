#' calculation of the Lx/Tx matrix
#'
#' Internal function called by \link{analyse_TL.MAAD} and \link{analyse_TL.SAR}. \cr
#' This function separates the Lx matrix from the Tx matrix.
#' Then, it estimates the Lx/Tx matrix.
#' It also provides a name for each of the curves. \cr
#'
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): \code{TLum.Analysis} object
#'
#' @return
#'  The function provides an \linkS4class{TLum.Results} containing: \cr
#'  \describe{
#'    \item{\code{Temperatures}}{
#'      \link{numeric}: Vector with the temperature values.}
#'    \item{\code{Names}}{
#'      \link{character}: Vector with the curve names.}
#'    \item{\code{Datatype}}{
#'      \link{character}: Vector with the curve type.}
#'    \item{\code{Doses}}{
#'      \link{numeric}: Vector with the curve doses.}
#'    \item{\code{Testdoses}}{
#'      \link{numeric}: Vector with the curve test-doses.}
#'    \item{\code{Lx}}{
#'      Lx matrix.}
#'    \item{\code{Lx.error}}{
#'      Absolute error for the Lx matrix.}
#'    \item{\code{Tx}}{
#'      Tx matrix.}
#'    \item{\code{Tx.error}}{
#'      Absolute error for the Tx matrix}
#'    \item{\code{LxTx}}{
#'      Lx/Tx matrix.}
#'    \item{\code{LxTx.error}}{
#'      Absolute error for the Lx/Tx matrix.}
#'  }
#'
#' @references
#'  Aitken, M.J. (1985) Thermoluminescence Dating, Academic Press, London \cr
#'
#'  Murray & Wintle (2000). Luminescence dating of quartz using an improved single-aliquot regenerative-dose protocol. Radiation Measurements, Vol.32, No.1, p.57-73. \cr
#'
#' @seealso
#'  \link{analyse_TL.MAAD},
#'  \link{analyse_TL.SAR}.
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export calc_TL.LxTx

calc_TL.LxTx <- function(

  object

){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(object)){
    stop("[analyse_TL.MAAD] Error: Input object is missing.")

  }else if (!is(object,"TLum.Analysis")){
    stop("[analyse_TL.MAAD] Error: Input object is not of type 'TLum.Analysis'.")
  }

  #------------------------------------------------------------------------------

  nRecords <- length(object@records)

  nPoints <- length(object@records[[1]]@data)
  temperatures <- object@records[[1]]@temperatures
  Tmax <- max(temperatures)

  # ------------------------------------------------------------------------------
  # Value check

  # Check Tmax & nPoints
  for(i in 1: nRecords){
    temp.record <- object@records[[i]]

    temp.nPoints<- temp.record@metadata$NPOINTS
    temp.Tmax <- max(temp.record@temperatures)

    if(temp.nPoints != nPoints){
      stop("[analyse_TL.MAAD] Error: All the TL curves do not have the same number of points'.")
    }

    if(temp.Tmax != Tmax){
      stop("[analyse_TL.MAAD] Error: All the TL curves do not have the same maximum temperature'.")
    }
  }
  #------------------------------------------------------------------------------

  dTypes <- vector()
  doses <- vector()
  testdoses <- vector()
  names <- vector()

  Lx <- vector()
  Lx.error <- vector()

  Tx <- vector()
  Tx.error <- vector()


  #Test dose identification + Separation of Lx and Tx.
  for (i in 1:nRecords){

    temp.record <- object@records[[i]]
    temp.dType <- temp.record@metadata$DTYPE

    temp.dose <- temp.record@metadata$IRR_TIME

    temp.curve <- temp.record@data
    temp.curve.error <- temp.record@error

    # Check -------------------------------------------------
    if(length(temp.curve) != length(temp.curve.error)){
      stop("[analyse_TL.MAAD] Error: The signal and the Error vector have a different length.")
    }
    #-----------------------------------------

    if(temp.dType != "Testdose"){
      dTypes <- c(dTypes,temp.dType)
      doses <- c(doses,temp.dose)

      Lx <- cbind(Lx,temp.curve)
      Lx.error <- cbind(Lx.error,temp.curve.error)

    }else{
      testdoses <- c(testdoses,temp.dose)

      Tx <- cbind(Tx,temp.curve)
      Tx.error <- cbind(Tx.error,temp.curve.error)
    }
  }

  #------------------------------------------------------------------------------
  # Check: Lx & Tx length
  if(length(Tx) > 0){
    if(length(Lx) != length(Tx)){
      stop("[calc_TL.LxTx] Error: Lx and Tx do not have the same size.")
    }
  }else{
    warning("[calc_TL.LxTx] Warning: No Testdose.")
  }

  #------------------------------------------------------------------------------
  # determination of the name of the curve (Natural, R0, R1, R2,...)
  temp.r <- 1
  temp.a <- 1
  temp.b <- 1
  temp.i <- 1

  for(i in 1 : length(doses)){

    if(dTypes[i] == "Natural"){
      temp.name <- "N"

    }else if(dTypes[i] == "Dose"){
      temp.name <- paste("R", temp.r, sep="")
      temp.r <- temp.r + 1

    }else if(dTypes[i] == "N+dose"){
      temp.name <- paste("A", temp.a, sep="")
      temp.a <- temp.a + 1

    }else if(dTypes[i] == "Bleach"){
      temp.name <- "B0"

    }else if(dTypes[i] == "Bleach+dose"){
      temp.name <- paste("B", temp.b, sep="")
      temp.b <- temp.b + 1

    }else{
      temp.name <- paste("O", temp.i, sep="")
      temp.i <- temp.i + 1
    }

    names <- c(names, temp.name)
  }

  # Naming Lx
  names(doses) <- names

  colnames(Lx) <- names
  colnames(Lx.error) <- names

  # Naming Tx
  if(length(Tx) > 0){
    names(testdoses) <- names

    colnames(Tx) <- names
    colnames(Tx.error) <- names
  }
  #------------------------------------------------------------------------------

  #Lx/Tx
  if(length(Tx) > 0){
    LxTx <- Lx/Tx

    #error estimation
    Lx.error.r <- abs(Lx.error/Lx)
    Tx.error.r <- abs(Tx.error/Tx)
    LxTx.error.r <- sqrt(Lx.error.r^2 + Tx.error.r^2)
    LxTx.error <- abs(LxTx*LxTx.error.r)

  }else{
    LxTx <- Lx
    LxTx.error <- Lx.error

    warning("[calc_TL.LxTx] Warning: LxTx = Lx.")
  }

  #Replace NA, NaN and Inf by NA

  Lx[!is.finite(Lx)] <- NA
  Lx.error[!is.finite(Lx.error)] <- NA

  Tx[!is.finite(Tx)] <- NA
  Tx.error[!is.finite(Tx.error)] <- NA

  LxTx[!is.finite(LxTx)] <- NA
  LxTx.error[!is.finite(LxTx.error)] <- NA

  result <- list(Temperatures=temperatures,
                 Names=names,
                 Datatypes=dTypes,
                 Doses=doses,
                 Testdoses=testdoses,
                 LxTx=as.data.frame(LxTx),
                 LxTx.error=as.data.frame(LxTx.error),
                 Lx=as.data.frame(Lx),
                 Lx.error=as.data.frame(Lx.error),
                 Tx=as.data.frame(Tx),
                 Tx.error=as.data.frame(Tx.error)
                 )

  new.TLum.Results.calc_TL.LxTx <- set_TLum.Results(data = result)

  return (new.TLum.Results.calc_TL.LxTx)
}
