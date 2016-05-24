#' SAR protocol for thermoluminescence dating
#'
#' This function calculates the equivalent dose (ED) using the SAR protocol. \cr
#' See details for more information.
#'
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the TL curves used for the ED calculation.
#' @param eval.Tmin
#'  \link{integer} (\bold{required}): Temperature (°C) of the lowest boundary for the signal integration.
#' @param eval.Tmax
#'  \link{integer} (\bold{required}): Temperature (°C) of the upper boundary for the signal integration.
#' @param rejection.criteria
#'  \link{list} (with default): list containing the rejection criteria (in \%). See details.
#' @param fitting.parameters
#'  \link{list} (with default): list containing the fitting parameters. See details.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'
#' @details
#' This function estimates the equivent dose in thermoluminescence dating using the SAR protocol.
#' The equivalent dose is estimated for each disc using the growth curve approaches (QC) (Aitken, 1985) and the dose plateau approach (DP).
#' Both approach should provide a similar result. \cr
#'
#' The Lx/Tx matrix is estimated using \link{calc_TL.LxTx}. \cr
#' The plateau test values are estimated using \link{calc_TL.plateau}. \cr
#'
#' \bold{Rejection criteria} \cr
#' The rejection criteria are: \cr
#' \describe{
#'  \item{\code{recycling.ratio}}{
#'    \link{numeric}: Maximum recycling ratio accepted (in \%).}
#'  \item{\code{recuperation.rate}}{
#'    \link{numeric}: Maximum recuperation rate accepted (in \%).}
#'  \item{\code{paleodose.error}}{
#'    \link{numeric}: Maximum error accepted on Lx (in \%).}
#'  \item{\code{testdose.error}}{
#'    \link{numeric}: Maximum error accepted on Tx (in \%).}
#'  }
#'
#' \bold{Fitting parameters} \cr
#' The fitting parameters are:  \cr
#' \describe{
#'  \item{\code{method}}{
#'    \link{character}: Fitting method (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}).}
#'  \item{\code{fit.weighted}}{
#'    \link{logical}: If the fitting is weighted or not.}
#'  \item{\code{fit.rDoses.min}}{
#'    \link{numeric}: lowest regenerative dose used for the fitting.}
#'  \item{\code{fit.rDoses.max}}{
#'    \link{numeric}: highest regenerative dose used for the fitting.}
#' }
#' See also \link{calc_TL.SAR.fit}. \cr
#'
#' \bold{Plotting parameters} \cr
#' The plotting parameters are:  \cr
#' \describe{
#'  \item{\code{plot.Tmin}}{
#'    \link{numeric}: lowest temperature plotted.}
#'  \item{\code{plot.Tmax}}{
#'    \link{numeric}: highest temperature plotted.}
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#' See also \link{plot_TL.SAR}. \cr
#'
#' @return
#'  The results are plotted using \link{plot_TL.SAR}. \cr
#'
#'  The function also provides an \linkS4class{TLum.Results} containing: \cr
#'  \describe{
#'    \item{\code{De.GC}}{
#'      \link{list}: Results obtained with the dose plateau approach and their uncertainties.
#'      (\code{De}, \code{De.error})}
#'    \item{\code{De.DP}}{
#'      \link{list}: Results obtained with the growth curve approach and their uncertainties.
#'      (\code{De}, \code{De.error})}
#'    \item{\code{LnLxTnTx.table}}{
#'      \link{matrix}: Lx/Tx values}
#'    \item{\code{RC.Status}}{
#'    \link{character}: Results of the rejection tests.
#'    }
#'  }
#'
#' @references
#'  Aitken, M.J. (1985) Thermoluminescence Dating, Academic Press, London \cr
#'
#'  Murray & Wintle (2000). Luminescence dating of quartz using an improved single-aliquot regenerative-dose protocol. Radiation Measurements, Vol.32, No.1, p.57-73. \cr
#'
#' @seealso
#'  \link{calc_TL.LxTx},
#'  \link{calc_TL.plateau},
#'  \link{calc_TL.SAR.fit},
#'  \link{analyse_TL.MAAD}.
#'
#' @author David Strebler, University of Cologne (Germany), \cr David Strebler
#'
#' @export analyse_TL.SAR


analyse_TL.SAR <- function(

  object,

  eval.Tmin,

  eval.Tmax,

  rejection.criteria = list(recycling.ratio = 10,
                            recuperation.rate = 10,
                            testdose.error = 10,
                            paleodose.error = 10
                            ),

  fitting.parameters=list(fit.method="LIN",
                          fit.weighted=FALSE,
                          fit.rDoses.min=NA,
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
    stop("[analyse_TL.SAR] Error: Input 'object' is not of type 'TLum.Analysis'.")
  }

  if (missing(eval.Tmin)){
    stop("[analyse_TL.SAR] Error: Input 'eval.Tmin' is missing.")

  }else if(!is.numeric(eval.Tmin)){
    stop("[calc_TL.MAAD.average] Error: Input 'eval.Tmin' is not of type 'numeric'.")
  }

  if (missing(eval.Tmax)){
    stop("[analyse_TL.SAR] Error: Input 'eval.Tmax' is missing.")

  }else if(!is.numeric(eval.Tmax)){
    stop("[calc_TL.MAAD.average] Error: Input 'eval.Tmax' is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------

  data <- calc_TL.LxTx(object=object)

  temperatures  <- get_TLum.Results(data, "Temperatures")
  Tmax <- max(temperatures)

  # ------------------------------------------------------------------------------
  # Value check
  if(eval.Tmin <= 0){
    stop("[analyse_TL.SAR] Error: eval.Tmin <= 0.")
  }

  if(eval.Tmax > Tmax){
    eval.Tmax <- Tmax
    warning(paste("[analyse_TL.SAR] Error: eval.Tmax >", Tmax, "(Tmax)."))
  }

  if(eval.Tmin > eval.Tmax){
    stop("[analyse_TL.SAR] Error: eval.Tmin > eval.Tmax.")
  }
  # ------------------------------------------------------------------------------

  sample.name <- as.character(object@records[[1]]@metadata$SAMPLE)
  sample.position <- object@records[[1]]@metadata$POSITION
  nPoints <- object@records[[1]]@metadata$NPOINTS
  Trate <- Tmax/nPoints

  eval.min <- ceiling(eval.Tmin/Trate)
  eval.max <- floor(eval.Tmax/Trate)

  rNames <- get_TLum.Results(data, "Names")
  rDoses <- get_TLum.Results(data, "Doses")
  temperatures  <- get_TLum.Results(data, "Temperatures")

  Lx <- as.matrix(get_TLum.Results(data, "Lx"))
  Lx.error <- as.matrix(get_TLum.Results(data, "Lx.error"))

  Tx <- as.matrix(get_TLum.Results(data, "Tx"))
  Tx.error <- as.matrix(get_TLum.Results(data, "Tx.error"))

  LxTx <- as.matrix(get_TLum.Results(data, "LxTx"))
  LxTx.error <- as.matrix(get_TLum.Results(data, "LxTx.error"))

  # ------------------------------------------------------------------------------
  # Value check
  if(length(Lx) == 0){
    stop("[analyse_TL.SAR] Error: No Lx.")
  }

  if(length(Tx) == 0){
    stop("[analyse_TL.SAR] Error: No Tx.")
  }

  if(length(LxTx) == 0){
    stop("[analyse_TL.SAR] Error: No LxTx.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # rLx, rTx & rLxTx
  # ------------------------------------------------------------------------------
  #rLx

  rLx<- Lx
  rLx.error <- Lx.error

  rLx.w <- 1/(rLx.error^2)

  rLx[!is.finite(rLx)] <- 0.0001
  rLx.error[!is.finite(rLx.error)] <- 0.0001
  rLx.w[!is.finite(rLx.w)] <- 0


  #rTx
  rTx<- Tx
  rTx.error <- Tx.error

  rTx.w <- 1/(rTx.error^2)

  rTx[!is.finite(rTx)] <- 0.0001
  rTx.error[!is.finite(rTx.error)] <- 0.0001
  rTx.w[!is.finite(rTx.w)] <- 0

  #rLxTx
  rLxTx<- LxTx
  rLxTx.error <- LxTx.error
  rLxTx.w <- 1/(rLxTx.error^2)

  rLxTx[!is.finite(rLxTx)] <- 0.0001
  rLxTx.error[!is.finite(rLxTx.error)] <- 0.0001
  rLxTx.w[!is.finite(rLxTx.w)] <- 0

  #-------------------------------------------------------------------------------------------------------------------------------------
  # Plateau test
  #-------------------------------------------------------------------------------------------------------------------------------------

  if(length(rLx) > 0){
    rLx.plateau <- get_TLum.Results(calc_TL.plateau(Ln=rLx[,rNames=="N"],
                                                    Ln.error=rLx.error[,rNames=="N"],
                                                    Lx=rLx[,rDoses!=0],
                                                    Lx.error=rLx.error[,rDoses!=0]),
                                    "LnLx")
  }


  if(length(rLx) > 0){
    rTx.plateau <- get_TLum.Results(calc_TL.plateau(Ln=rTx[,rNames=="N"],
                                                    Ln.error=rTx.error[,rNames=="N"],
                                                    Lx=rTx[,rDoses!=0],
                                                    Lx.error=rTx.error[,rDoses!=0]),
                                    "LnLx")
  }

  if(length(rLx) > 0){
    rLxTx.plateau <- get_TLum.Results(calc_TL.plateau(Ln=rLxTx[,rNames=="N"],
                                                      Ln.error=rLxTx.error[,rNames=="N"],
                                                      Lx=rLxTx[,rDoses!=0],
                                                      Lx.error=rLxTx.error[,rDoses!=0]),
                                      "LnLx")
  }

  #-------------------------------------------------------------------------------------------------------------------------------------
  # GC.rTx & GC.rLxTx (average)
  #-------------------------------------------------------------------------------------------------------------------------------------

  GC.rTx <- vector()
  GC.rTx.error <- vector()

  GC.rLxTx <- vector()
  GC.rLxTx.error <- vector()

  for(i in 1 : length(rDoses)){
    temp.name <- rNames[i]

    # GC.rTx
    temp.rTx <- rTx[eval.min:eval.max, i]
    temp.rTx.error <- rTx.error[eval.min:eval.max, i]
    temp.rTx.w <- rTx.w[eval.min:eval.max, i]

    temp.GC.rTx <- sum(temp.rTx.w*temp.rTx,na.rm=TRUE)/sum(temp.rTx.w,na.rm=TRUE)
    temp.GC.rTx.error <- 1/sqrt(sum(temp.rTx.w, na.rm = TRUE))

    names(temp.GC.rTx) <- temp.name
    names(temp.GC.rTx.error) <- temp.name

    # GC.rLxTx
    temp.rLxTx <- rLxTx[eval.min:eval.max, i]
    temp.rLxTx.error <- rLxTx.error[eval.min:eval.max, i]
    temp.rLxTx.w <- rLxTx.w[eval.min:eval.max, i]

    temp.GC.rLxTx <- sum(temp.rLxTx.w*temp.rLxTx, na.rm=TRUE)/sum(temp.rLxTx.w, na.rm=TRUE)
    temp.GC.rLxTx.error <- 1/sqrt(sum(temp.rLxTx.w, na.rm = TRUE))

    names(temp.GC.rLxTx) <- temp.name
    names(temp.GC.rLxTx.error) <- temp.name

    # Save
    GC.rTx <- c(GC.rTx, temp.GC.rTx)
    GC.rTx.error <- c(GC.rTx.error, temp.GC.rTx.error)

    GC.rLxTx <- c(GC.rLxTx, temp.GC.rLxTx)
    GC.rLxTx.error <- c(GC.rLxTx.error, temp.GC.rLxTx.error)
  }

  #-------------------------------------------------------------------------------------------------------------------------------------
  #Selection criteria
  #-------------------------------------------------------------------------------------------------------------------------------------

  recycling.ratio <- rejection.criteria$recycling.ratio
  recuperation.rate <- rejection.criteria$recuperation.rate
  testdose.error <- rejection.criteria$testdose.error
  paleodose.error <- rejection.criteria$paleodose.error

  fit.rDoses.min <- fitting.parameters$fit.rDoses.min
  fit.rDoses.max <- fitting.parameters$fit.rDoses.max

  # ------------------------------------------------------------------------------
  # Value check
  # regenerative doses
  if(is.null(fit.rDoses.min) || is.na(fit.rDoses.min) || fit.rDoses.min < min(rDoses)){
    fit.rDoses.min <- min(rDoses[rDoses != 0])
  }

  if(is.null(fit.rDoses.max) || is.na(fit.rDoses.max) || fit.rDoses.max > max(rDoses)){
    fit.rDoses.max <- max(rDoses)
  }

  if(fit.rDoses.max < fit.rDoses.min){
    stop("[analyse_TL.SAR] Error: fit.rDoses.max < fit.rDoses.min")
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
  # ------------------------------------------------------------------------------
  recycling.limit <- round(abs(recycling.ratio/100),digits=3)
  recuperation.limit <- round(abs(recuperation.rate/100),digits=3)
  Lx.error.criteria <- round(abs(paleodose.error/100),digits=3)
  Tx.error.criteria <- round(abs(testdose.error/100),digits=3)

  # recycling ratio
  recycling.ratio <- vector()

  rNames.duplicated <- rNames[duplicated(rDoses) & rDoses!=0]

  for(i in 1: length(rNames.duplicated)){
    temp.name.r <- rNames.duplicated[i]
    temp.dose <- rDoses[temp.name.r]
    temp.name.i <- rNames[rDoses==temp.dose & !duplicated(rDoses)]

    temp.name <- paste(temp.name.r, "/", temp.name.i,sep ="")

    temp.GC.rLxTx.i <- GC.rLxTx[rNames==temp.name.i]
    temp.GC.rLxTx.r <- GC.rLxTx[rNames==temp.name.r]

    temp.ratio <- temp.GC.rLxTx.r/temp.GC.rLxTx.i
    #temp.ratio.error
    names(temp.ratio) <- temp.name


    recycling.ratio <- c(recycling.ratio, temp.ratio)
    #recycling.ratio.error
  }

  # recuperation Rate
  recuperation.rate <- vector()

  rNames.zero <- rNames[rNames!="N" & rDoses==0]

  for(i in 1:length(rNames.zero)){
    temp.name.z <- rNames.zero[i]

    temp.name <- paste(temp.name.z,"/N", sep ="")

    temp.GC.rLxTx.n <- GC.rLxTx[rNames=="N"]
    temp.GC.rLxTx.z <- GC.rLxTx[rNames==temp.name.z]

    temp.rate <- temp.GC.rLxTx.z/temp.GC.rLxTx.n
    #temp.rate.error
    names(temp.rate) <- temp.name

    recuperation.rate <- c(recuperation.rate, temp.rate)
    #recuperation.rate.error
  }

  #Max paleodose error
  temp.rLx <- rLx[eval.min:eval.max, rDoses!=0]
  temp.rLx.error <- rLx.error[eval.min:eval.max, rDoses!=0]
  temp.rLx.error.r <- temp.rLx.error/temp.rLx

  rLx.error.r.GC <- vector()

  for(i in 1:ncol(temp.rLx.error.r)){
    temp.rLx.error.r.GC <- mean(temp.rLx.error.r[, i],na.rm = TRUE)

    rLx.error.r.GC <- c(rLx.error.r.GC, temp.rLx.error.r.GC)
  }

  Lx.error.max <- max(abs(rLx.error.r.GC), na.rm = TRUE)

  #Max testdose error
  temp.rTx <- rTx[eval.min:eval.max, rDoses!=0]
  temp.rTx.error <- rTx.error[eval.min:eval.max, rDoses!=0]
  temp.rTx.error.r <- temp.rTx.error/temp.rTx

  rTx.error.r.GC <- vector()

  for(i in 1:ncol(temp.rTx.error.r)){
    temp.rTx.error.r.GC <- mean(temp.rTx.error.r[, i],na.rm = TRUE)

    rTx.error.r.GC <- c(rTx.error.r.GC, temp.rTx.error.r.GC)
  }

  Tx.error.max <- max(abs(rTx.error.r.GC), na.rm = TRUE)


  # Recycling test
  test.recycling <- vector()

  for(i in 1:length(recycling.ratio)){

    if(abs(1 - recycling.ratio[i]) > recycling.limit){
      temp.test <- FALSE
    }else{
      temp.test <- TRUE
    }
    test.recycling <- c(test.recycling, temp.test)
  }

  if(FALSE %in% test.recycling){
    recycling.result <- "FAILED"
  }else{
    recycling.result <- "OK"
  }

  # Recuperation test
  test.recuperation <- vector()


  for(i in 1:length(recuperation.rate)){

    if(abs(recuperation.rate[i]) > recycling.limit){
      temp.test <- FALSE
    }else{
      temp.test <- TRUE
    }
    test.recuperation <- c(test.recuperation, temp.test)
  }

  if(FALSE %in% test.recuperation){
    recuperation.result <- "FAILED"
  }else{
    recuperation.result <- "OK"
  }

  # Paleodose error test
  if(Lx.error.max > Lx.error.criteria){
    test.Lx.error <- FALSE
    Lx.error.result <- "FAILED"
  }else{
    test.Lx.error <- TRUE
    Lx.error.result <- "OK"
  }

  # Testdose error test
  if(Tx.error.max > Tx.error.criteria){
    test.Tx.error <- FALSE
    Tx.error.result <- "FAILED"
  }else{
    test.Tx.error <- TRUE
    Tx.error.result <- "OK"
  }

  # Acceptance result
  if(!(FALSE %in% test.recycling)
     && !(FALSE %in% test.recuperation)
     && test.Lx.error
     && test.Tx.error) {
    acceptance.result <- "OK"
  }else{
    acceptance.result <- "FAILED"
  }

  rejection.values <- list(recycling.ratio=recycling.ratio,
                           recuperation.rate=recuperation.rate,
                           Lx.error.max=Lx.error.max,
                           Tx.error.max=Tx.error.max,
                           test.recycling=test.recycling,
                           test.recuperation=test.recuperation,
                           test.Lx.error=test.Lx.error,
                           test.Tx.error=test.Lx.error)


  # ------------------------------------------------------------------------------
  # testdose response
  # ------------------------------------------------------------------------------
  TxTn <- vector()

  for(i in 1 : length(rDoses)){
    temp.TxTn <- GC.rTx[i]/GC.rTx["N"]

    TxTn <- c(TxTn, temp.TxTn)
  }
  names(TxTn) <- rNames

  #-------------------------------------------------------------------------------------------------------------------------------------
  # De
  #-------------------------------------------------------------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Dose plateau approach
  # ------------------------------------------------------------------------------
  Q.DP <- vector()
  Q.DP.error <- vector()

  for(i in 1 : length(temperatures)){

    temp.LnTn <- rLxTx[i, rNames=="N"]
    temp.LnTn.error <- rLxTx.error[i, rNames=="N"]

    temp.rDoses <- rDoses[rNames!="N" & rDoses!=0]

    temp.rLxTx <- rLxTx[i, rNames!="N" & rDoses!=0]
    temp.rLxTx.error <- rLxTx.error[i, rNames!="N" & rDoses!=0]

    # Data
    temp.bool <- temp.rDoses >= fit.rDoses.min & temp.rDoses <= fit.rDoses.max

    temp.doses <- temp.rDoses[temp.bool]
    temp.LxTx <- temp.rLxTx[temp.bool]
    temp.LxTx.error <- temp.rLxTx.error[temp.bool]

    temp.fit <- calc_TL.SAR.fit(LnTn=temp.LnTn,
                                LnTn.error=temp.LnTn.error,
                                LxTx=temp.LxTx,
                                LxTx.error=temp.LxTx.error,
                                doses=temp.doses,
                                fitting.parameters=fitting.parameters)

    temp.Q.DP <- get_TLum.Results(temp.fit, "Q")
    temp.Q.DP.error <- get_TLum.Results(temp.fit, "Q.error")

    #save
    Q.DP <- c(Q.DP, temp.Q.DP)
    Q.DP.error <- c(Q.DP.error, temp.Q.DP.error)
  }
  Q.DP.w <- 1/(Q.DP.error^2)

  # Q.DP average
  if(length(Q.DP) > 0){
    Q.DP.lim <- Q.DP[eval.min:eval.max]
    Q.DP.lim.error <- Q.DP.error[eval.min:eval.max]
    Q.DP.lim.w <- 1/(Q.DP.lim.error^2)

    Q.DP.a <- sum(Q.DP.lim.w*Q.DP.lim,na.rm=TRUE)/sum(Q.DP.lim.w,na.rm=TRUE)
    Q.DP.a.error <- 1/sqrt(sum(Q.DP.lim.w, na.rm = TRUE))
  }else{
    Q.DP.a <- 0
    Q.DP.a.error <- 0
  }

  # ------------------------------------------------------------------------------
  # Growth curve approach
  # ------------------------------------------------------------------------------
  temp.LnTn <- GC.rLxTx[rNames == "N"]
  temp.LnTn.error <- GC.rLxTx.error[rNames == "N"]

  temp.rDoses <- rDoses[rNames!= "N" & rDoses!=0]

  temp.rLxTx <- GC.rLxTx[rNames != "N" & rDoses!=0]
  temp.rLxTx.error <- GC.rLxTx.error[rNames != "N" & rDoses!=0]

  # Data
  temp.bool <- temp.rDoses >= fit.rDoses.min & temp.rDoses <= fit.rDoses.max

  temp.doses <- temp.rDoses[temp.bool]
  temp.LxTx <- temp.rLxTx[temp.bool]
  temp.LxTx.error <- temp.rLxTx.error[temp.bool]


  temp.fit <- calc_TL.SAR.fit(LnTn=temp.LnTn,
                             LnTn.error=temp.LnTn.error,
                             LxTx=temp.LxTx,
                             LxTx.error=temp.LxTx.error,
                             doses=temp.doses,
                             fitting.parameters=fitting.parameters)

  GC.Q.line <- get_TLum.Results(temp.fit, "GC")
  Q.GC <- get_TLum.Results(temp.fit, "Q")
  Q.GC.error <- get_TLum.Results(temp.fit, "Q.error")
  Q.GC.slope <- get_TLum.Results(temp.fit, "summary")


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
    plot_TL.SAR(sample.name=sample.name,
                sample.position=sample.position,
                plotting.parameters=plotting.parameters,
                fitting.parameters=fitting.parameters,
                eval.Tmin=eval.Tmin,
                eval.Tmax=eval.Tmax,
                temperatures=temperatures,
                names=rNames,
                names.duplicated=rNames.duplicated,
                doses=rDoses,
                Lx=rLx,
                Tx=rTx,
                LxTx=rLxTx,
                Lx.plateau=rLx.plateau,
                Tx.plateau=rTx.plateau,
                LxTx.plateau=rLxTx.plateau,
                DP.Q.line=Q.DP,
                DP.Q.line.error=Q.DP.error,
                Q.DP = Q.DP.a,
                Q.DP.error = Q.DP.a.error,
                GC.Q.line=GC.Q.line,
                GC.Q.slope=Q.GC.slope,
                GC.Q.LxTx = GC.rLxTx,
                GC.Q.LxTx.error=GC.rLxTx.error,
                Q.GC=Q.GC,
                Q.GC.error=Q.GC.error,
                TxTn=TxTn,
                rejection.values=rejection.values
                )
  }

  #----------------------------------------------------------------------------------------------------------------
  # Export results
  #----------------------------------------------------------------------------------------------------------------

  new.De.DP <- data.frame(De=Q.DP.a,
                          De.error=Q.DP.a.error)

  new.De.GC <- data.frame(De=Q.GC,
                         De.error=Q.GC.error)

  rejection.results <- data.frame(criteria = c("Recycling ratio", "Recuperation rate", "Lx error", "Tx error"),
                                           value = c(recycling.ratio, recuperation.rate, Lx.error.max, Tx.error.max),
                                           threshold = c(recycling.limit, recuperation.limit, Lx.error.criteria, Tx.error.criteria),
                                           status = c(recycling.result, recuperation.result, Lx.error.result, Tx.error.result)
                                           )

  new.TLum.Results.analyse_TL.SAR <- set_TLum.Results(data = list(DP = new.De.DP,
                                                                  GC = new.De.GC,
                                                                  LnLxTnTx.table = rLxTx,
                                                                  rejection.criteria = rejection.results,
                                                                  RC.Status=acceptance.result))

  return(new.TLum.Results.analyse_TL.SAR)
}

