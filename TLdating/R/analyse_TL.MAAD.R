#' MAAD protocol for TL dating
#'
#' Function to estimate the ED in TL dating using the MAAD protocol. \cr
#' It provides an estimation of the palaeodose (Q) and/or the sublinearity correction (I).
#' The equivalent dose (ED) is estimated by the addition of Q and I. \cr
#' See details for more information.
#'
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the TL curves used for the ED estimation.
#' @param eval.Tmin
#'  \link{integer} (\bold{required}): Temperature ( °C) of the lower boundary for the signal integration.
#' @param eval.Tmax
#'  \link{integer} (\bold{required}): Temperature ( °C) of the upper boundary for the signal integration.
#' @param rejection.criteria
#'  \link{list} (with default): list containing the rejection criteria (in \%). See details.
#' @param fitting.parameters
#'  \link{list} (with default): list containing the fitting parameters. See details.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'
#' @details
#' This function estimates the equivalent dose for the thermoluminescence dating with the MAAD protocol.
#' It can provide an estimation of the palaeodose (Q) and the sublinearity correction (I) simultaniously or separately.
#' These are estimated using the growth curve approach (QC) (Aitken, 1985) and the dose plateau approach (DP).
#' Both approaches should provide a similar result. The equivalent dose is estimated by the addition of Q and I\cr
#' The Lx/Tx matrix is estimated using \link{calc_TL.LxTx}. \cr
#' The average TL curves for each dose step are estimate using \link{calc_TL.MAAD.average}. \cr
#' The plateau test values are estimated using \link{calc_TL.plateau}. \cr
#'
#' \bold{Rejection criteria} \cr
#' The rejection criteria are: \cr
#' \describe{
#'  \item{\code{testdose.error}}{
#'    \link{numeric}: Maximum error accepted on Tx (in \%).}
#'  \item{\code{paleodose.error}}{
#'    \link{numeric}: Maximum error accepted on Lx (in \%).}
#'  }
#'
#' \bold{Fitting parameters} \cr
#' The fitting parameters are:  \cr
#' \describe{
#'  \item{\code{method}}{
#'    \link{character}: Fitting method (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}).}
#'  \item{\code{fit.weighted}}{
#'    \link{logical}: If the fitting is weighted or not.}
#'  \item{\code{fit.use.slope}}{
#'    \link{logical}: If the slope of the Q growth curve is reused for the sublinearity correction.}
#'  \item{\code{fit.aDoses.min}}{
#'    \link{numeric}: Lowest additive dose used for the fitting.}
#'  \item{\code{fit.aDoses.max}}{
#'    \link{numeric}: Highest additive dose used for the fitting.}
#'  \item{\code{fit.rDoses.min}}{
#'    \link{numeric}: Lowest regenerative dose used for the fitting.}
#'  \item{\code{fit.rDoses.max}}{
#'    \link{numeric}: Highest regenerative dose used for the fitting.}
#' }
#' See also \link{calc_TL.MAAD.fit.Q} and \link{calc_TL.MAAD.fit.I}. \cr
#'
#' \bold{Plotting parameters} \cr
#' The plotting parameters are:  \cr
#' \describe{
#'  \item{\code{plot.Tmin}}{
#'    \link{numeric}: Lowest temperature plotted.}
#'  \item{\code{plot.Tmax}}{
#'    \link{numeric}: Highest temperature plotted.}
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#' See also \link{plot_TL.MAAD}. \cr
#'
#' @return
#'  The results are plotted using \link{plot_TL.MAAD}. \cr
#'
#'  The function also provides a \linkS4class{TLum.Results} containing: \cr
#'  \describe{
#'    \item{\code{De.GC}}{
#'      \link{list}: Results obtained with the dose plateau approach and their uncertainties
#'      (\code{De}, \code{De.error}, \code{Q}, \code{Q.error}, \code{I}, \code{I.error})}
#'    \item{\code{De.DP}}{
#'      \link{list}: Results obtained with the growth curve approach and their uncertainties
#'      (\code{De}, \code{De.error}, \code{Q}, \code{Q.error}, \code{I}, \code{I.error})}
#'    \item{\code{LnLxTnTx.table}}{
#'      \link{matrix}: Lx/Tx values}
#'    \item{\code{RC.Status}}{
#'    \link{character}: The acceptance result. }
#'  }
#'
#'
#' @references
#'  Aitken, M.J. (1985) Thermoluminescence Dating, Academic Press, London \cr
#'
#' @seealso
#'  \link{calc_TL.LxTx},
#'  \link{calc_TL.plateau},
#'  \link{calc_TL.MAAD.average},
#'  \link{calc_TL.MAAD.separate},
#'  \link{calc_TL.MAAD.fit.I},
#'  \link{calc_TL.MAAD.fit.Q},
#'  \link{analyse_TL.SAR}.
#'
#'
#' @author David Strebler, University of Cologne (Germany)
#'
#' @examples
#'
#' ##load data
#'
#' ##perform analysis
#'
#' @export analyse_TL.MAAD


analyse_TL.MAAD <- function(
  object,

  eval.Tmin,

  eval.Tmax,

  rejection.criteria = list(testdose.error = 10,
                            paleodose.error = 10),

  fitting.parameters=list(fit.method="LIN",
                          fit.weighted=FALSE,
                          fit.use.slope=FALSE,
                          fit.aDoses.min=0,
                          fit.aDoses.max=NA,
                          fit.rDoses.min=0,
                          fit.rDoses.max=NA),

  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA,
                           no.plot=FALSE)

  ){

  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------

  if (missing(object)){
    stop("[analyse_TL.SAR] Error: Input 'object' is missing.")

  }else if (!is(object,"TLum.Analysis")){
    stop("[analyse_TL.MAAD] Error: Input 'object' is not of type 'TLum.Analysis'.")
  }

  if (missing(eval.Tmin)){
    stop("[analyse_TL.MAAD] Error: Input 'eval.Tmin' is missing.")

  }else if(!is.numeric(eval.Tmin)){
    stop("[calc_TL.MAAD.average] Error: Input 'eval.Tmin' is not of type 'numeric'.")
  }

  if (missing(eval.Tmax)){
    stop("[analyse_TL.MAAD] Error: Input 'eval.Tmax' is missing.")

  }else if(!is.numeric(eval.Tmax)){
    stop("[calc_TL.MAAD.average] Error: Input 'eval.Tmax' is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------

  sample.name <- as.character(object@records[[1]]@metadata$SAMPLE)

  data <- calc_TL.LxTx(object)

  temperatures  <- get_TLum.Results(data, "Temperatures")

  nPoints <- length(temperatures)
  Tmax <- max(temperatures)
  Tstep <- Tmax/nPoints

  # ------------------------------------------------------------------------------
  # Value check

  if(eval.Tmin < 0){
    stop("[analyse_TL.MAAD] Error: eval.Tmin < 0.")
  }

  if(eval.Tmax > Tmax){
    stop(paste("[analyse_TL.MAAD] Error: eval.Tmax >", Tmax, "(Tmax)."))
  }

  if(eval.Tmin > eval.Tmax){
    stop("[analyse_TL.MAAD] Error: eval.Tmin > eval.Tmax.")
  }

  if(eval.Tmin < 0){
    stop("[analyse_TL.MAAD] Error: eval.Tmin < 0.")
  }
  # ------------------------------------------------------------------------------

  eval.min <- ceiling(eval.Tmin/Tstep)
  eval.max <- floor(eval.Tmax/Tstep)

  dTypes <- get_TLum.Results(data, "Datatypes")

  doses <- get_TLum.Results(data, "Doses")

  Lx <- as.matrix(get_TLum.Results(data, "Lx"))
  Lx.error <- as.matrix(get_TLum.Results(data, "Lx.error"))

  Tx <- as.matrix(get_TLum.Results(data, "Tx"))
  Tx.error <- as.matrix(get_TLum.Results(data, "Tx.error"))

  LxTx <- as.matrix(get_TLum.Results(data, "LxTx"))
  LxTx.error <- as.matrix(get_TLum.Results(data, "LxTx.error"))


  # ------------------------------------------------------------------------------
  # Separation of aLx and rLx
  # ------------------------------------------------------------------------------

  temp.data <- calc_TL.MAAD.separate(doses=doses,
                                     Lx=Lx,
                                     Lx.error=Lx.error,
                                     dTypes=dTypes)

  aDoses.f <- get_TLum.Results(temp.data, "aDoses")
  aNames.f <- get_TLum.Results(temp.data, "aNames")

  rDoses.f <- get_TLum.Results(temp.data, "rDoses")
  rNames.f <- get_TLum.Results(temp.data, "rNames")

  aLx.f <- as.matrix(get_TLum.Results(temp.data, "aLx"))
  aLx.f.error <- as.matrix(get_TLum.Results(temp.data, "aLx.error"))

  rLx.f <- as.matrix(get_TLum.Results(temp.data, "rLx"))
  rLx.f.error <- as.matrix(get_TLum.Results(temp.data, "rLx.error"))

  temp.data <- calc_TL.MAAD.separate(doses=doses,
                                     Lx=Tx,
                                     Lx.error=Tx.error,
                                     dTypes=dTypes)

  aTx.f <- as.matrix(get_TLum.Results(temp.data, "aLx"))
  aTx.f.error <- as.matrix(get_TLum.Results(temp.data, "aLx.error"))

  rTx.f <- as.matrix(get_TLum.Results(temp.data, "rLx"))
  rTx.f.error <- as.matrix(get_TLum.Results(temp.data, "rLx.error"))

  temp.data <- calc_TL.MAAD.separate(doses=doses,
                                     Lx=LxTx,
                                     Lx.error=LxTx.error,
                                     dTypes=dTypes)

  aLxTx.f <- as.matrix(get_TLum.Results(temp.data, "aLx"))
  aLxTx.f.error <- as.matrix(get_TLum.Results(temp.data, "aLx.error"))

  rLxTx.f <- as.matrix(get_TLum.Results(temp.data, "rLx"))
  rLxTx.f.error <- as.matrix(get_TLum.Results(temp.data, "rLx.error"))


  # ------------------------------------------------------------------------------
  # Average for aLx and rLx
  # ------------------------------------------------------------------------------

  # Additive doses
  aDoses.a <- vector()
  aNames.a <- vector()

  aLx.a <- vector()
  aLx.a.error <- vector()

  aTx.a <- vector()
  aTx.a.error <- vector()

  aLxTx.a  <- vector()
  aLxTx.a.error <- vector()

  if(length(aLx.f) > 0){
    temp.data <- calc_TL.MAAD.average(names=aNames.f,
                                      doses=aDoses.f,
                                      Lx=aLx.f,
                                      Lx.error=aLx.f.error)

    aDoses.a <- get_TLum.Results(temp.data, "doses")
    aNames.a <- get_TLum.Results(temp.data, "names")

    aLx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
    aLx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))

    if(length(aTx.f) > 0){


      temp.data <- calc_TL.MAAD.average(names=aNames.f,
                                        doses=aDoses.f,
                                        Lx=aTx.f,
                                        Lx.error=aTx.f.error)

      aTx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
      aTx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))
    }

    temp.data <- calc_TL.MAAD.average(names=aNames.f,
                                      doses=aDoses.f,
                                      Lx=aLxTx.f,
                                      Lx.error=aLxTx.f)

    aLxTx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
    aLxTx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))
  }

  #Regenerative doses
  rDoses.a <- vector()
  rNames.a <- vector()

  rLx.a <- vector()
  rLx.a.error <- vector()

  rTx.a <- vector()
  rTx.a.error <- vector()

  rLxTx.a  <- vector()
  rLxTx.a.error <- vector()

  if(length(rLx.f) > 0){
    temp.data <- calc_TL.MAAD.average(names=rNames.f,
                                    doses=rDoses.f,
                                    Lx=rLx.f,
                                    Lx.error=rLx.f.error)

    rDoses.a <- get_TLum.Results(temp.data, "doses")
    rNames.a <- get_TLum.Results(temp.data, "names")

    rLx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
    rLx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))

    if(length(rTx.f) > 0){
      temp.data <- calc_TL.MAAD.average(names=rNames.f,
                                        doses=rDoses.f,
                                        Lx=rTx.f,
                                        Lx.error=rTx.f.error)

      rTx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
      rTx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))
    }

    temp.data <- calc_TL.MAAD.average(names=rNames.f,
                                      doses=rDoses.f,
                                      Lx=rLxTx.f,
                                      Lx.error=rLxTx.f.error)

    rLxTx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
    rLxTx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))
  }

  # ------------------------------------------------------------------------------
  # Plateau test
  # ------------------------------------------------------------------------------

  #Plateau test (Additive step)
    #aLx
  if(length(aLx.a) > 0){
    aLx.a.plateau <- get_TLum.Results(calc_TL.plateau(Ln=aLx.a[,aNames.a == "N"],
                                                      Ln.error=aLx.a.error[,aNames.a == "N"],
                                                      Lx=aLx.a[,aNames.a != "N" & aDoses.a != 0],
                                                      Lx.error=aLx.a.error[,aNames.a != "N"& aDoses.a != 0]),
                                      "LnLx")
  }
    #atx
  if(length(aTx.a) > 0){
    aTx.a.plateau <- get_TLum.Results(calc_TL.plateau(Ln=aTx.a[,aNames.a == "N"],
                                                      Ln.error=aTx.a.error[,aNames.a == "N"],
                                                      Lx=aTx.a[,aNames.a != "N" & aDoses.a != 0],
                                                      Lx.error=aTx.a.error[,aNames.a != "N"& aDoses.a != 0]),
                                      "LnLx")
  }
    #aLxTx
  if(length(aLxTx.a) > 0){
    aLxTx.a.plateau <- get_TLum.Results(calc_TL.plateau(Ln=aLxTx.a[,aNames.a == "N"],
                                                        Ln.error=aLxTx.a.error[,aNames.a == "N"],
                                                        Lx=aLxTx.a[,aNames.a != "N" & aDoses.a != 0],
                                                        Lx.error=aLxTx.a.error[,aNames.a != "N"& aDoses.a != 0]),
                                        "LnLx")
  }

  #Plateau test (Additive step)
    #rLx
  if(length(rLx.a) > 0){
    if(length(aLx.a)>0){
      rLx.a.plateau <- get_TLum.Results(calc_TL.plateau(Ln=aLx.a[,aNames.a == "N"],
                                                        Ln.error=aLx.a.error[,aNames.a == "N"],
                                                        Lx=rLx.a[,rNames.a != "N" & rDoses.a != 0],
                                                        Lx.error=rLx.a.error[,rNames.a != "N"& rDoses.a != 0]),
                                        "LnLx")
    }
    else{
      rLx.a.plateau <- rLx.a
    }
  }
    #rTx
  if(length(rTx.a) > 0){
    if(length(aTx.a)>0){
      rTx.a.plateau <- get_TLum.Results(calc_TL.plateau(Ln=aTx.a[,aNames.a == "N"],
                                                        Ln.error=aTx.a.error[,aNames.a == "N"],
                                                        Lx=rTx.a[,rNames.a != "N" & rDoses.a != 0],
                                                        Lx.error=rTx.a.error[,rNames.a != "N"& rDoses.a != 0]),
                                        "LnLx")
    }else{
      rTx.a.plateau<-rTx.a
    }
  }
    #rLxTx
  if(length(rLxTx.a) > 0){
    if(length(aLxTx.a)>0){
      rLxTx.a.plateau <- get_TLum.Results(calc_TL.plateau(Ln=aLxTx.a[,aNames.a == "N"],
                                                          Ln.error=aLxTx.a.error[,aNames.a == "N"],
                                                          Lx=rLxTx.a[,rNames.a != "N" & rDoses.a != 0],
                                                          Lx.error=rLxTx.a.error[,rNames.a != "N"& rDoses.a != 0]),
                                          "LnLx")
    }else{
      rLxTx.a.plateau<-rLxTx.a
    }
  }

  #-------------------------------------------------------------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Selection criteria
  # ------------------------------------------------------------------------------

  fit.aDoses.min <- fitting.parameters$fit.aDoses.min
  fit.aDoses.max <- fitting.parameters$fit.aDoses.max

  fit.rDoses.min <- fitting.parameters$fit.rDoses.min
  fit.rDoses.max <- fitting.parameters$fit.rDoses.max

  paleodose.error <- rejection.criteria$paleodose.error
  testdose.error <- rejection.criteria$testdose.error

  # ------------------------------------------------------------------------------
  # Value check
  # additive doses
  if(length(aDoses.a) > 0){
    if(is.null(fit.aDoses.min) || is.na(fit.aDoses.min) || fit.aDoses.min < min(aDoses.a)){
      fit.aDoses.min <- min(aDoses.a[aDoses.a != 0])
    }

    if(is.null(fit.aDoses.max) || is.na(fit.aDoses.max) || fit.aDoses.max > max(aDoses.a)){
      fit.aDoses.max <- max(aDoses.a)
    }

    if(fit.aDoses.max < fit.aDoses.min){
      stop("[analyse_TL.SAR] Error: fit.aDoses.max < fit.aDoses.min")
    }
  }else{
    fit.aDoses.min <- NULL
    fit.aDoses.max <- NULL
  }

  # regenerative doses
  if(length(rDoses.a) > 0){
    if(is.null(fit.rDoses.min) || is.na(fit.rDoses.min) || fit.rDoses.min < min(rDoses.a)){
      fit.rDoses.min <- min(rDoses.a[rDoses.a != 0])
    }

    if(is.null(fit.rDoses.max) || is.na(fit.rDoses.max) || fit.rDoses.max > max(rDoses.a)){
      fit.rDoses.max <- max(rDoses.a)
    }

    if(fit.rDoses.max < fit.rDoses.min){
      stop("[analyse_TL.SAR] Error: fit.rDoses.max < fit.rDoses.min")
    }
  }else{
    fit.rDoses.min <- NULL
    fit.rDoses.max <- NULL
  }

  # rejection.criteria
  if(is.null(paleodose.error) || !is.finite(paleodose.error)){
    stop("[analyse_TL.SAR] Error: paleodose.error is missing.")
  }else if(paleodose.error<0){
    paleodose.error <- abs(paleodose.error)
    warning("[analyse_TL.SAR] Warning: paleodose.error is negative.")
  }

  if(is.null(testdose.error) || !is.finite(testdose.error)){
    stop("[analyse_TL.SAR] Error: paleodose.error is missing.")
  }else if(testdose.error<0){
    testdose.error <- abs(testdose.error)
    warning("[analyse_TL.SAR] Warning: testdose.error is negative.")
  }
  #---------------------------------------------------------
  Lx.error.criteria <- round(abs(paleodose.error/100),digits=3)
  Tx.error.criteria <- round(abs(testdose.error/100),digits=3)


  #Max paleodose error (addivite curve)
  if(length(aLx.f) > 0){
    temp.aLx.lim <- aLx.f[eval.min:eval.max, aDoses.f >= fit.aDoses.min & aDoses.f <= fit.aDoses.max]
    temp.aLx.lim.error <- aLx.f.error[eval.min:eval.max, aDoses.f >= fit.aDoses.min & aDoses.f <= fit.aDoses.max]
    temp.aLx.lim.error.r <- abs(temp.aLx.lim.error/temp.aLx.lim)

    temp.aLx.lim.error.r[!is.finite(temp.aLx.lim.error.r)] <- NA

    aLx.error.r.GC <- vector()

    for(i in 1:ncol(temp.aLx.lim.error.r)){
      temp.aLx.error.r.GC <- mean(temp.aLx.lim.error.r[,i],na.rm = TRUE)

      aLx.error.r.GC <- c(aLx.error.r.GC, temp.aLx.error.r.GC)
    }

    aLx.error.max <- max(aLx.error.r.GC, na.rm = TRUE)
  }else{
    aLx.error.max <- NA
  }

  #Max testdose error (additive curve)
  if(length(aTx.f) > 0){
    temp.aTx.lim <- aTx.f[eval.min:eval.max, aDoses.f >= fit.aDoses.min & aDoses.f <= fit.aDoses.max]
    temp.aTx.lim.error <- aTx.f.error[eval.min:eval.max, aDoses.f >= fit.aDoses.min & aDoses.f <= fit.aDoses.max]
    temp.aTx.lim.error.r <- abs(temp.aTx.lim.error/temp.aTx.lim)

    temp.aTx.lim.error.r[!is.finite(temp.aTx.lim.error.r)] <- NA

    aTx.error.r.GC <- vector()

    for(i in 1:ncol(temp.aTx.lim.error.r)){
      temp.aTx.error.r.GC <- mean(temp.aTx.lim.error.r[,i],na.rm = TRUE)

      aTx.error.r.GC <- c(aTx.error.r.GC, temp.aTx.error.r.GC)
    }

    aTx.error.max <- max(aTx.error.r.GC, na.rm = TRUE)
  }else{
    aTx.error.max <- NA
  }

  #Max paleodose error (regenerative curve)
  if(length(rLx.f) > 0){
    temp.rLx.lim <- rLx.f[eval.min:eval.max, rDoses.f >= fit.rDoses.min & rDoses.f <= fit.rDoses.max & rDoses.f != 0]
    temp.rLx.lim.error <- rLx.f.error[eval.min:eval.max, rDoses.f >= fit.rDoses.min & rDoses.f <= fit.rDoses.max & rDoses.f != 0]
    temp.rLx.lim.error.r <- abs(temp.rLx.lim.error/temp.rLx.lim)

    temp.rLx.lim.error.r[!is.finite(temp.rLx.lim.error.r)] <- NA

    rLx.error.r.GC <- vector()

    for(i in 1:ncol(temp.rLx.lim.error.r)){
      temp.rLx.error.r.GC <- mean(temp.rLx.lim.error.r[,i],na.rm = TRUE)

      rLx.error.r.GC <- c(rLx.error.r.GC, temp.rLx.error.r.GC)
    }

    rLx.error.max <- max(rLx.error.r.GC, na.rm = TRUE)
  }else{
    rLx.error.max <- NA
  }

  #Max testdose error (regenerative curve)
  if(length(rTx.f) > 0){
    temp.rTx.lim <- rTx.f[eval.min:eval.max,rDoses.a >= fit.rDoses.min & rDoses.a <= fit.rDoses.max]
    temp.rTx.lim.error <- rTx.f.error[eval.min:eval.max,rDoses.a >= fit.rDoses.min & rDoses.a <= fit.rDoses.max ]
    temp.rTx.lim.error.r <- abs(temp.rTx.lim.error/temp.rTx.lim)

    temp.rTx.lim.error.r[!is.finite(temp.rTx.lim.error.r)] <- NA

    rTx.error.r.GC <- vector()

    for(i in 1:ncol(temp.rTx.lim.error.r)){
      temp.rTx.error.r.GC <- mean(temp.rTx.lim.error.r[,i],na.rm = TRUE)

      rTx.error.r.GC <- c(rTx.error.r.GC, temp.rTx.error.r.GC)
    }

    rTx.error.max <- max(rTx.error.r.GC, na.rm = TRUE)
  }else{
    rTx.error.max <- NA
  }

  # Max paleodose error test (additive dose)
  if(length(aLx.f) > 0){
    if(aLx.error.max > Lx.error.criteria){
      test.aLx.error <- FALSE
    }else{
      test.aLx.error <- TRUE
    }
  }else{
    test.aLx.error <- TRUE
  }

  # Max testdose error test (additive dose)
  if(length(aTx.f) > 0){
    if(aTx.error.max > Tx.error.criteria){
      test.aTx.error <- FALSE
    }else{
      test.aTx.error <- TRUE
    }
  }else{
    test.aTx.error <- TRUE
  }

  # Max paleodose test (regenerative curve)
  if(length(rLx.f) > 0){
    if(rLx.error.max > Lx.error.criteria){
      test.rLx.error <- FALSE
    }else{
      test.rLx.error <- TRUE
    }
  }else{
    test.rLx.error <- TRUE
  }

  # Max testdose test (regenerative curve)
  if(length(rTx.f) > 0){
    if(rTx.error.max > Tx.error.criteria){
      test.rTx.error <- FALSE
    }else{
      test.rTx.error <- TRUE
    }
  }else{
    test.rTx.error <- TRUE
  }

  # Acceptance result
  if(test.aLx.error && test.aTx.error && test.rLx.error && test.rTx.error) {
    acceptance.result <- "OK"
  }else{
    acceptance.result <- "FAILED"
  }

  rejection.values <-list(aLx.error.max=aLx.error.max,
                          aTx.error.max=aTx.error.max,
                          rLx.error.max=rLx.error.max,
                          rTx.error.max=rTx.error.max,
                          test.aLx.error=test.aLx.error,
                          test.aTx.error=test.aTx.error,
                          test.rLx.error=test.rLx.error,
                          test.rTx.error=test.rTx.error)

  #----------------------------------------------------------------------------------------------------------------
  # D_e estimation
  #----------------------------------------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Palaeodose - Q
  # ------------------------------------------------------------------------------

  if(length(aLx.a) > 0){

    # Growth curve (GC)
    aLxTx.a.lim <- aLxTx.a[eval.min:eval.max,]
    aLxTx.a.error.lim <- aLxTx.a.error[eval.min:eval.max,]

    aLxTx.GC <- vector()
    aLxTx.GC.error <- vector()

    for(i in 1:length(aDoses.a)){

      temp.LxTx <- aLxTx.a.lim[,i]
      temp.LxTx.error <- aLxTx.a.error.lim[,i]
      temp.LxTx.w <- 1/(temp.LxTx.error^2)

      aLxTx.GC[i] <- sum(temp.LxTx.w*temp.LxTx,na.rm=TRUE)/sum(temp.LxTx.w,na.rm=TRUE)
      aLxTx.GC.error[i] <- 1/sqrt(sum(temp.LxTx.w))
    }

    # Q (GC)

      #Data
    temp.bool <- aDoses.a >= fit.aDoses.min & aDoses.a <= fit.aDoses.max
    temp.doses <- aDoses.a[aDoses.a >= fit.aDoses.min & aDoses.a <= fit.aDoses.max]
    temp.LxTx <- aLxTx.GC[aDoses.a >= fit.aDoses.min & aDoses.a <= fit.aDoses.max]
    temp.LxTx.error <- aLxTx.GC.error[aDoses.a >= fit.aDoses.min & aDoses.a <= fit.aDoses.max]


    temp.fit <- calc_TL.MAAD.fit.Q(LxTx = temp.LxTx,
                                LxTx.error = temp.LxTx.error,
                                doses = temp.doses,
                                fitting.parameters=fitting.parameters)

    GC.Q <- get_TLum.Results(temp.fit, "GC")

    Q.GC <- get_TLum.Results(temp.fit, "Q")
    Q.GC.error <- get_TLum.Results(temp.fit, "Q.error")
    Q.GC.slope <- get_TLum.Results(temp.fit, "summary")

    # Q (DP)
    Q.DP <- vector()
    Q.DP.error <- vector()
    Q.DP.slope <- list()

    for(i in 1:length(temperatures)){
      temp.aLxTx <- aLxTx.a[i,]
      temp.aLxTx.error <- aLxTx.a.error[i,]

      # Data
      temp.bool <- aDoses.a >= fit.aDoses.min & aDoses.a <= fit.aDoses.max

      temp.doses <- aDoses.a[temp.bool]
      temp.LxTx <- temp.aLxTx[temp.bool]
      temp.LxTx.error <- temp.aLxTx.error[temp.bool]

      # Regression
      temp.fit <- calc_TL.MAAD.fit.Q(LxTx = temp.LxTx,
                                   LxTx.error = temp.LxTx.error,
                                   doses = temp.doses,
                                   fitting.parameters=fitting.parameters)

      temp.GC <- get_TLum.Results(temp.fit, "GC")

      temp.Q <- get_TLum.Results(temp.fit, "Q")
      temp.Q.error <- get_TLum.Results(temp.fit, "Q.error")
      temp.slope <- get_TLum.Results(temp.fit, "summary")

      #save
      Q.DP <- c(Q.DP, temp.Q)
      Q.DP.error <- c(Q.DP.error, temp.Q.error)
      Q.DP.slope[[i]] <- temp.slope
    }

    # Q.DP average
    Q.DP.w <- 1/(Q.DP.error^2)

    Q.DP.lim <- Q.DP[eval.min:eval.max]
    Q.DP.lim.error <- Q.DP.error[eval.min:eval.max]
    Q.DP.lim.w <- 1/(Q.DP.lim.error^2)

    Q.DP.a <- sum(Q.DP.lim.w*Q.DP.lim,na.rm=TRUE)/sum(Q.DP.lim.w,na.rm=TRUE)
    Q.DP.a.error <- 1/sqrt(sum(Q.DP.lim.w, na.rm = TRUE))
  }else{
    aLxTx.GC <- vector()
    aLxTx.GC.error <- vector()

    GC.Q <- vector()

    Q.DP <- vector()
    Q.DP.error <- vector()
    Q.DP.a <- vector()
    Q.DP.a.error <- vector()
    Q.GC <- vector()
    Q.GC.error <- vector()

    Q.DP.slope <- list()
    Q.GC.slope <- list()

    }

  if(length(Q.DP) == 0){
    Q.DP.a <- 0
    Q.DP.a.error <- 0
    Q.GC <- 0
    Q.GC.error <- 0
  }

  # ------------------------------------------------------------------------------
  # Sublinearity correction - I
  # ------------------------------------------------------------------------------

  if(length(rLxTx.a) > 0){

    # Growth curve
    rLxTx.a.w <- 1/(rLxTx.a.error^2)

    rLxTx.a.lim <- rLxTx.a[eval.min:eval.max,]
    rLxTx.a.w.lim <- rLxTx.a.w[eval.min:eval.max,]

    GC.rLxTx <- vector()
    GC.rLxTx.error <- vector()

    for(i in 1:length(rDoses.a)){

      temp.LxTx <- rLxTx.a.lim[,i]
      temp.LxTx.w <- rLxTx.a.w.lim[,i]

      GC.rLxTx[i] <- sum(temp.LxTx.w*temp.LxTx, na.rm=TRUE)/sum(temp.LxTx.w, na.rm=TRUE)
      GC.rLxTx.error[i] <- 1/sqrt(sum(temp.LxTx.w))
    }

    #I (GC)
      # Data
    temp.bool <- rDoses.a >= fit.rDoses.min & rDoses.a <= fit.rDoses.max
    temp.doses <- rDoses.a[temp.bool]
    temp.LxTx <- GC.rLxTx[temp.bool]
    temp.LxTx.error <- GC.rLxTx.error[temp.bool]


    if(length(GC.Q)>0){
      temp.slope <- Q.GC.slope

      temp.fit <- calc_TL.MAAD.fit.I(LxTx = temp.LxTx,
                                   LxTx.error = temp.LxTx.error,
                                   doses = temp.doses,
                                   slope = temp.slope,
                                   fitting.parameters=fitting.parameters)
    }else{

      fitting.parameters$fit.use.slope =FALSE
      temp.fit <- calc_TL.MAAD.fit.I(LxTx = temp.LxTx,
                                   LxTx.error = temp.LxTx.error,
                                   doses = temp.doses,
                                   fitting.parameters=fitting.parameters)
    }


    GC.I <- get_TLum.Results(temp.fit, "GC")

    I.GC <- get_TLum.Results(temp.fit, "I")
    I.GC.error <- get_TLum.Results(temp.fit, "I.error")
    I.GC.slope <- get_TLum.Results(temp.fit, "summary")


    # I (DP)
    I.DP <- vector()
    I.DP.error <- vector()
    I.DP.slope <- list()

    for(i in 1:length(temperatures)){
      temp.rLxTx <- rLxTx.a[i,]
      temp.rLxTx.error <- rLxTx.a.error[i,]

      #temp.rLxTx.w <-rLxTx.a.w[i,]

        # selection of the doses used
      temp.dose <- rDoses.a[rDoses.a >= fit.rDoses.min & rDoses.a <= fit.rDoses.max]
      temp.LxTx <- temp.rLxTx[rDoses.a >= fit.rDoses.min & rDoses.a <= fit.rDoses.max]
      temp.LxTx.error <- temp.rLxTx.error[rDoses.a >= fit.rDoses.min & rDoses.a <= fit.rDoses.max]

      #Regression
      if(length(GC.Q)>0){
        temp.slope <- Q.DP.slope[[i]]

        temp.fit <- calc_TL.MAAD.fit.I(LxTx = temp.LxTx,
                                     LxTx.error = temp.LxTx.error,
                                     doses = temp.dose,
                                     slope = temp.slope,
                                     fitting.parameters=fitting.parameters)
      }else{

        fitting.parameters$fit.use.slope =FALSE
        temp.fit <- calc_TL.MAAD.fit.I(LxTx = temp.LxTx,
                                     LxTx.error = temp.LxTx.error,
                                     doses = temp.dose,
                                     fitting.parameters=fitting.parameters)
      }

      temp.GC <- get_TLum.Results(temp.fit, "GC")

      temp.I <- get_TLum.Results(temp.fit, "I")
      temp.I.error <- get_TLum.Results(temp.fit, "I.error")
      temp.slope <- get_TLum.Results(temp.fit, "summary")

      I.DP <- c(I.DP, temp.I)
      I.DP.error <- c(I.DP.error, temp.I.error)
      I.DP.slope[[i]] <- temp.slope
    }

    # I average
    I.DP.w <- 1/(I.DP.error^2)

    I.DP.lim <- I.DP[eval.min:eval.max]
    I.DP.lim.error <- I.DP.error[eval.min:eval.max]
    I.DP.lim.w <- I.DP.w[eval.min:eval.max]

    I.DP.a <- sum(I.DP.lim.w*I.DP.lim, na.rm=TRUE)/sum(I.DP.lim.w, na.rm=TRUE)
    I.DP.a.error <- 1/sqrt(sum(I.DP.lim.w, na.rm = TRUE))

  }else{
    GC.rLxTx <- vector()
    GC.rLxTx.error <- vector()

    I.DP <- vector()
    I.DP.error <- vector()
    I.DP.a <- vector()
    I.DP.a.error <- vector()
    I.GC <- vector()
    I.GC.error <- vector()
    I.DP.slope <- list()
    I.GC.slope <- list()
  }

  if(length(I.DP) == 0){
    I.DP.a <- 0
    I.DP.a.error <- 0
    I.GC <- 0
    I.GC.error <- 0
  }

  # ------------------------------------------------------------------------------
  # Equivalent dose - De (Q+I)
  # ------------------------------------------------------------------------------
  #De.GC
  De.GC <- sum(Q.GC, I.GC, na.rm = TRUE)
  De.GC.error <- sqrt(sum(Q.GC.error^2,I.GC.error^2, na.rm = TRUE))

  if(!is.finite(De.GC)){
    De.GC <- 0
  }

  if(!is.finite(De.GC.error)){
    De.GC.error <- De.GC
  }

  # De.DP
  De.DP <- sum(Q.DP.a, I.DP.a ,na.rm = TRUE)
  De.DP.error <- sqrt(sum(Q.DP.a.error^2,I.DP.a.error^2, na.rm = TRUE))


  if(!is.finite(De.DP)){
    De.DP <- 0
  }

  if(!is.finite(De.DP.error)){
    De.DP.error <- De.DP
  }


  #----------------------------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------------------------

  no.plot <- plotting.parameters$no.plot

  # ------------------------------------------------------------------------------
  #Check values
  if(is.null(no.plot) || is.na(no.plot) || !is.logical(no.plot)){
    no.plot <- FALSE
  }
  # ------------------------------------------------------------------------------

  if(!no.plot){
    plot_TL.MAAD(sample.name=sample.name,
                 fitting.parameters=fitting.parameters,
                 temperatures=temperatures,
                 eval.Tmin=eval.Tmin,
                 eval.Tmax=eval.Tmax,
                 aNames=aNames.a,
                 aDoses=aDoses.a,
                 aLx=aLx.a,
                 aTx=aTx.a,
                 aLxTx=aLxTx.a,
                 aLx.plateau=aLx.a.plateau,
                 aTx.plateau=aTx.a.plateau,
                 aLxTx.plateau=aLxTx.a.plateau,
                 rNames=rNames.a,
                 rDoses=rDoses.a,
                 rLx=rLx.a,
                 rTx=rTx.a,
                 rLxTx=rLxTx.a,
                 rLx.plateau=rLx.a.plateau,
                 rTx.plateau=rTx.a.plateau,
                 rLxTx.plateau=rLxTx.a.plateau,
                 DP.Q.line=Q.DP,
                 DP.Q.line.error=Q.DP.error,
                 GC.Q.slope=Q.GC.slope,
                 GC.Q.line=GC.Q,
                 GC.Q.LxTx=GC.rLxTx,
                 GC.Q.LxTx.error=GC.rLxTx.error,
                 DP.I.line=I.DP,
                 DP.I.line.error=I.DP.error,
                 GC.I.slope=I.GC.slope,
                 GC.I.line=GC.I,
                 GC.I.LxTx=aLxTx.GC,
                 GC.I.LxTx.error=aLxTx.GC.error,
                 Q.DP=Q.DP.a,
                 Q.DP.error=Q.DP.a.error,
                 Q.GC=Q.GC,
                 Q.GC.error=Q.GC.error,
                 I.DP=I.DP.a,
                 I.DP.error=I.DP.a.error,
                 I.GC=I.GC,
                 I.GC.error=I.GC.error,
                 De.GC=De.GC,
                 De.GC.error=De.GC.error,
                 De.DP=De.DP,
                 De.DP.error=De.DP.error,
                 rejection.values=rejection.values,
                 plotting.parameters=plotting.parameters
    )
  }


  #----------------------------------------------------------------------------------------------------------------
  #Export Results
  #----------------------------------------------------------------------------------------------------------------

  new.De.DP <- list(De = De.DP,
                    De.error=De.DP.error,
                    Q = Q.DP.a,
                    Q.error = Q.DP.a.error,
                    I = I.DP.a,
                    I.error = I.DP.a.error)

  new.De.GC <- list(De = De.GC,
                    De.error=De.GC.error,
                    Q = Q.GC,
                    Q.error = Q.GC.error,
                    I = I.GC,
                    I.error = I.GC.error)


  new.TLum.Results.analyse_TL.MAAD <- set_TLum.Results(data = list(DP = new.De.DP,
                                                                   GC = new.De.GC,
                                                                   LnLxTnTx.table = LxTx,
                                                                   RC.Status = acceptance.result)
                                                       )

  return(new.TLum.Results.analyse_TL.MAAD)
}
