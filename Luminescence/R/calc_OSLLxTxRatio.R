#' Calculate Lx/Tx ratio for CW-OSL curves
#'
#' Calculate Lx/Tx ratios from a given set of CW-OSL curves assuming late light background subtraction.
#'
#' The integrity of the chosen values for the signal and background integral is
#' checked by the function; the signal integral limits have to be lower than
#' the background integral limits. If a \link{vector} is given as input instead
#' of a \link{data.frame}, an artificial \code{data.frame} is produced. The
#' error calculation is done according to Galbraith (2002).\cr
#'
#' \bold{sigmab}\cr
#'
#' The default value of \code{sigmab} is calculated assuming the background is
#' constant and \bold{would not} applicable when the background varies as,
#' e.g., as observed for the early light substraction method.
#'
#' \bold{background.count.distribution}\cr
#'
#' This argument allows selecting the distribution assumption that is used for
#' the error calculation. According to Galbraith (2002, 2014) the background
#' counts may be overdispersed (i.e. do not follow a poisson distribution,
#' which is assumed for the photomultiplier counts). In that case (might be the
#' normal case) it has to be accounted for the overdispersion by estimating
#' \eqn{\sigma^2} (i.e. the overdispersion value). Therefore the relative
#' standard error is calculated as:\cr\cr (a) \code{poisson}\cr
#' \deqn{rse(\mu_{S}) \approx \sqrt(Y_{0} + Y_{1}/k^2)/Y_{0} - Y_{1}/k} (b)
#' \code{non-poisson}\cr \deqn{rse(\mu_{S}) \approx \sqrt(Y_{0} + Y_{1}/k^2 +
#' \sigma^2(1+1/k))/Y_{0} - Y_{1}/k}
#'
#' \bold{Please note} that when using the early background subtraction method in
#' combination with the 'non-poisson' distribution argument, the corresponding Lx/Tx error
#' may considerably increase due to a high sigmab value.
#' Please check whether this is valid for your data set and  if necessary
#' consider to provide an own sigmab value using the corresponding argument \code{sigmab}.
#'
#' @param Lx.data \code{\linkS4class{RLum.Data.Curve}} or \link{data.frame}
#' (\bold{required}): requires a CW-OSL shine down curve (x = time, y = counts)
#'
#' @param Tx.data \code{\linkS4class{RLum.Data.Curve}} or \link{data.frame}
#' (optional): requires a CW-OSL shine down curve (x = time, y = counts). If no
#' input is given the Tx.data will be treated as \code{NA} and no Lx/Tx ratio
#' is calculated.
#'
#' @param signal.integral \code{\link{vector}} (\bold{required}): vector with the
#' limits for the signal integral.
#'
#' @param signal.integral.Tx \code{\link{vector}} (optional): vector with the
#' limits for the signal integral for the Tx curve. If nothing is provided the
#' value from \code{signal.integral} is used.
#'
#' @param background.integral \code{\link{vector}} (\bold{required}): vector with the
#' bounds for the background integral.
#'
#' @param background.integral.Tx \code{\link{vector}} (optional): vector with the
#' limits for the background integral for the Tx curve. If nothing is provided the
#' value from \code{background.integral} is used.
#'
#' @param background.count.distribution \code{\link{character}} (with default): Sets
#' the count distribution assumed for the error calculation. Possible arguments
#' \code{poisson} or \code{non-poisson}. See details for further information
#'
#' @param sigmab \link{numeric} (optional): Option to set a manual value for
#' the overdispersion (for LnTx and TnTx), used for the Lx/Tx error
#' calculation. The value should be provided as absolute squared count values,
#' e.g. \code{sigmab = c(300,300)}. Note: If only one value is provided this
#' value is taken for both (LnTx and TnTx) signals.
#'
#' @return Returns an S4 object of type \code{\linkS4class{RLum.Results}}.
#'
#' Slot \code{data} contains a \code{\link{list}} with the following structure:\cr
#' $LxTx.table (data.frame) \cr
#' .. $ LnLx \cr
#' .. $ LnLx.BG \cr
#' .. $ TnTx \cr
#' .. $ TnTx.BG \cr
#' .. $ Net_LnLx \cr
#' .. $ Net_LnLx.Error\cr
#' .. $ Net_TnTx.Error\cr
#' .. $ LxTx\cr
#' .. $ LxTx.Error \cr
#' $ calc.parameters (list) \cr
#' .. $ sigmab.LnTx\cr
#' .. $ sigmab.TnTx\cr
#' .. $ k \cr
#' $ call (original function call)\cr
#'
#' @note The results of this function have been cross-checked with the Analyst
#' (vers. 3.24b). Access to the results object via  \code{\link{get_RLum}}.\cr
#'
#' \bold{Caution:} If you are using early light subtraction (EBG), please either provide your
#' own \code{sigmab} value or use \code{background.count.distribution = "poisson"}.
#'
#'
#' @section Function version: 0.6.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\linkS4class{RLum.Data.Curve}},
#' \code{\link{Analyse_SAR.OSLdata}}, \code{\link{plot_GrowthCurve}},
#' \code{\link{analyse_SAR.CWOSL}}
#'
#'  @references Duller, G., 2007. Analyst.
#' \url{http://www.nutech.dtu.dk/english/~/media/Andre_Universitetsenheder/Nutech/Produkter\%20og\%20services/Dosimetri/radiation_measurement_instruments/tl_osl_reader/Manuals/analyst_manual_v3_22b.ashx}\cr
#'
#' Galbraith, R.F., 2002. A note on the variance of a background-corrected OSL
#' count. Ancient TL, 20 (2), 49-51.
#'
#' Galbraith, R.F., 2014. A further note on the variance of a
#' background-corrected OSL count. Ancient TL, 31 (2), 1-3.
#'
#' @keywords datagen
#'
#' @examples
#'
#' ##load data
#' data(ExampleData.LxTxOSLData, envir = environment())
#'
#' ##calculate Lx/Tx ratio
#' results <- calc_OSLLxTxRatio(Lx.data, Tx.data, signal.integral = c(1:2),
#'                              background.integral = c(85:100))
#'
#' ##get results object
#' get_RLum(results)
#'
#' @export
calc_OSLLxTxRatio <- function(
  Lx.data,
  Tx.data,
  signal.integral,
  signal.integral.Tx = NULL,
  background.integral,
  background.integral.Tx = NULL,
  background.count.distribution = "non-poisson",
  sigmab
){

  ##--------------------------------------------------------------------------##
  ##(1) - integrity checks


  if(missing(Tx.data) == FALSE){

    ##(a) - check data type
    if(is(Lx.data)[1]!=is(Tx.data)[1]){
      stop("[calc_OSLLxTxRatio()] Data type of Lx and Tx data differs!")
    }

    ##(b) - test if data.type is valid in general
    if(is(Lx.data)[1] == "RLum.Data.Curve"){

      Lx.data <- as(Lx.data, "data.frame")
      Tx.data <- as(Tx.data, "data.frame")


    }else{

      ##go further
      if((is(Lx.data)[1] != "data.frame" &
          is(Lx.data)[1] != "numeric") &
         is(Lx.data)[1] != "matrix"){
        stop("[calc_OSLLxTxRatio()] Data type error! Required types are data.frame or numeric vector.")
      }
    }

    ##(c) - convert vector to data.frame if nescessary
    if(is(Lx.data)[1] != "data.frame" &
       is(Lx.data)[1] != "matrix"){
      Lx.data <- data.frame(x=1:length(Lx.data),y=Lx.data)
      Tx.data <- data.frame(x=1:length(Tx.data),y=Tx.data)
    }

    ##(d) - check if Lx and Tx curves have the same channel length
    if(length(Lx.data[,2]) != length(Tx.data[,2])){
      stop("[calc_OSLLxTxRatio()] Channel number of Lx and Tx data differs!")}

  }else{

    Tx.data <- data.frame(x = NA,y = NA)

    ##support RLum.objects
    if(is(Lx.data)[1] == "RLum.Data.Curve"){
      Lx.data <- as(Lx.data, "data.frame")

    }

    ##check for matrix
    if(is(Lx.data)[1] == "matrix"){
      Lx.data <- as.data.frame(Lx.data)

    }

    ##no it should be a data.frame, if not, try to produce one
    if(is(Lx.data)[1]!="data.frame") {
      Lx.data <- data.frame(x = 1:length(Lx.data),y = Lx.data)
    }

  }#endif::missing Tx.data

  ##(e) - check if signal integral is valid
  if(min(signal.integral) < 1 | max(signal.integral>length(Lx.data[,2]))){
    stop("[calc_OSLLxTxRatio()] signal.integral is not valid!")}

  ##(f) - check if background integral is valid
  if(min(background.integral)<1 | max(background.integral>length(Lx.data[,2]))){
    stop(paste("[calc_OSLLxTxRatio()] background.integral is not valid! Max: ",length(Lx.data[,2]),sep=""))}

  ##(g) - check if signal and background integral overlapping
  if(min(background.integral)<=max(signal.integral)){
    stop("[calc_OSLLxTxRatio()] Overlapping of 'signal.integral' and 'background.integral' is not permitted!")}

  ##(h) - similar procedure for the Tx limits
  if(all(c(!is.null(signal.integral.Tx),!is.null(background.integral.Tx)))){

    if(min(signal.integral.Tx) < 1 | max(signal.integral.Tx>length(Tx.data[,2]))){
      stop("[calc_OSLLxTxRatio()] signal.integral.Tx is not valid!")}

    if(min(background.integral.Tx)<1 | max(background.integral.Tx>length(Tx.data[,2]))){
      stop(paste("[calc_OSLLxTxRatio()] background.integral.Tx is not valid! Max: ",length(Tx.data[,2]),sep=""))}

    if(min(background.integral.Tx)<=max(signal.integral.Tx)){
      stop("[calc_OSLLxTxRatio()] Overlapping of 'signal.integral.Tx' and 'background.integral.Tx' is not permitted!")}

  }
  else if(!all(c(is.null(signal.integral.Tx),is.null(background.integral.Tx)))){
    stop("[calc_OSLLxTxRatio()] You have to provide both: signal.integral.Tx and background.integral.Tx!")

  }else{
    signal.integral.Tx <- signal.integral
    background.integral.Tx <- background.integral

  }



  ##check sigmab
  if (!missing(sigmab)) {
    if (!is.null(sigmab)) {
      if (!is(sigmab, "numeric")) {
        stop("[calc_OSLLxTxRatio()] 'sigmab' has to be of type numeric.")
      }

      if (length(sigmab) > 2) {
        stop("[calc_OSLLxTxRatio()] Maximum allowed vector length for 'sigmab' is 2.")
      }
    }
  }


  ##--------------------------------------------------------------------------##
  ##(2) - read data and produce background subtracted values

  ## calculate k value - express the background as mutiple value from the number
  ## of signal integral channels, however, it can be < 1 also
  n <- length(signal.integral)
  m <- length(background.integral)
  k <- m/n

  n.Tx <- length(signal.integral.Tx)
  m.Tx <- length(background.integral.Tx)
  k.Tx <- m.Tx/n.Tx


  ##LnLx (comments are corresponding variables to Galbraith, 2002)
  Lx.curve <- Lx.data[,2]
  Lx.signal <- sum(Lx.curve[signal.integral])                #Y.0
  Lx.background <- sum(Lx.curve[background.integral])        #Y.1
  Lx.background <- Lx.background*1/k                         #mu.B
  LnLx <- Lx.signal - Lx.background

  ##TnTx
  Tx.curve <- ifelse(is.na(Tx.data[,1])==FALSE, Tx.data[,2], NA)
  Tx.signal <- sum(Tx.curve[signal.integral.Tx])
  Tx.background <- sum(Tx.curve[background.integral.Tx])*1/k.Tx
  TnTx <- (Tx.signal-Tx.background)

  ##--------------------------------------------------------------------------##
  ##(3)
  ## calculate Lx/Tx Errors according Galbraith (2002) and the personal
  ## communication of Galbraith (2014) via e-mail
  ## Nomenclature as stated in the articles

  ##(a)
  ## set Y.0 (sum OSL signal including the background) and
  ## Y.1 (total counts over m later channels)
  Y.0 <- Lx.signal
  Y.0_TnTx <- Tx.signal
  Y.1 <- sum(Lx.curve[background.integral])
  Y.1_TnTx <- sum(Tx.curve[background.integral.Tx])


  ##(b) estimate overdispersion (here called sigmab), see equation (4) in
  ## Galbraith (2002), Galbraith (2014)
  ## If else condition for the case that k < 2

  if(round(k,digits = 1) >= 2 & ((min(background.integral) + length(signal.integral)*(2+1)) <= length(Lx.curve))){

    ##(b)(1)(1)
    ## note that m = n*k = multiple of background.integral from signal.integral
    Y.i <- sapply(0:round(k,digits=0), function(i){
      sum(Lx.curve[
        (min(background.integral)+length(signal.integral)*i):
          (min(background.integral)+length(signal.integral)+length(signal.integral)*i)])
    })

    Y.i <- na.exclude(Y.i)
    sigmab.LnLx <- abs(var(Y.i) - mean(Y.i))  ##sigmab is denoted as sigma^2 = s.Y^2-Y.mean
    ##therefore here absolute values are given


  }else{

    ## provide warning if m is < 25, as suggested by Rex Galbraith
    ## low number of degree of freedom
    if (m < 25) {
      warning("Number of background integral channels is < 25. The calculation
              might be not reliable!")

    }

    sigmab.LnLx <- abs((var(Lx.curve[background.integral]) -
                          mean(Lx.curve[background.integral])) * n)

  }

  if (round(k.Tx, digits = 1) >= 2 &
      ((
        min(background.integral.Tx) + length(signal.integral.Tx) * (2 + 1)
      ) <= length(Tx.curve))) {
    ##(b)(1)(1)
    ## note that m.Tx = n.Tx*k.Tx = multiple of background.integral.Tx from signal.integral.Tx
    ## also for the TnTx signal
    Y.i_TnTx <- sapply(0:round(k.Tx, digits = 0), function(i) {
      sum(Tx.curve[(min(background.integral.Tx) + length(signal.integral.Tx) *
                      i):(
                        min(background.integral.Tx) + length(signal.integral.Tx) + length(signal.integral.Tx) *
                          i
                      )])
    })

    Y.i_TnTx <- na.exclude(Y.i_TnTx)
    sigmab.TnTx <- abs(var(Y.i_TnTx) - mean(Y.i_TnTx))

  } else{
    ## provide warning if m is < 25, as suggested by Rex Galbraith
    ## low number of degree of freedom
    if (m.Tx < 25) {
      warning(
        "Number of background integral channels for Tx is < 25. The calculation
        might be not reliable!"
      )

    }

    sigmab.TnTx <- abs((var(Tx.curve[background.integral.Tx]) -
                          mean(Tx.curve[background.integral.Tx])) * n.Tx)
  }


  ##account for a manually set sigmab value
  if (!missing(sigmab)) {
    if (!is.null(sigmab)) {
      if (length(sigmab) == 2) {
        sigmab.LnLx <- sigmab[1]
        sigmab.TnTx <- sigmab[2]

      }else{
        sigmab.LnLx <- sigmab[1]
        sigmab.TnTx <- sigmab[1]

      }
    }
  }

  ##(c)
  ## Calculate relative error of the background subtracted signal
  ## according to Galbratith (2002), equation (6) with changes
  ## from Galbraith (2014), equation 6
  ## Discussion with Rex Galbraith via e-mail (2014-02-27):
  ## Equation 6 is approriate to be implemented as standard

  if(background.count.distribution == "poisson"){

    ##(c.1) estimate relative standard error for assuming a poisson distribution
    LnLx.relError <- sqrt((Y.0 + Y.1/k^2))/(Y.0-Y.1/k)        ##  rse(mu.s)
    TnTx.relError <- sqrt((Y.0_TnTx + Y.1_TnTx/k^2))/(Y.0_TnTx-Y.1_TnTx/k)

  }else{

    ##(c.2) estimate relative standard error for a non-poisson distribution
    if(background.count.distribution != "non-poisson"){
      warning("Unknown method for background.count.distribution. A non-poisson distribution is assumed!")}

    LnLx.relError <- sqrt(Y.0 + Y.1/k^2 + sigmab.LnLx*(1+1/k))/
      (Y.0 - Y.1/k)
    TnTx.relError <- sqrt(Y.0_TnTx + Y.1_TnTx/k^2 + sigmab.TnTx*(1+1/k))/
      (Y.0_TnTx - Y.1_TnTx/k)

  }

  ##(d)
  ##calculate absolute standard error
  LnLx.Error <- abs(LnLx*LnLx.relError)
  TnTx.Error <- abs(TnTx*TnTx.relError)

  ##combine results
  LnLxTnTx <- cbind(
    Lx.signal,
    Lx.background,
    Tx.signal,
    Tx.background,
    LnLx,
    LnLx.Error,
    TnTx,
    TnTx.Error
  )

  ##--------------------------------------------------------------------------##
  ##(4) Calculate LxTx error according Galbraith (2014)

  #transform results in a data.frame
  LnLxTnTx <- as.data.frame((LnLxTnTx))

  #add col names
  colnames(LnLxTnTx)<-c("LnLx", "LnLx.BG",
                        "TnTx", "TnTx.BG",
                        "Net_LnLx", "Net_LnLx.Error",
                        "Net_TnTx", "Net_TnTx.Error")

  ##calculate Ln/Tx
  LxTx <- LnLxTnTx$Net_LnLx/LnLxTnTx$Net_TnTx

  ##calculate Ln/Tx error
  LxTx.relError <- sqrt(LnLx.relError^2 + TnTx.relError^2)
  LxTx.Error <- abs(LxTx * LxTx.relError)

  ##return combined values
  temp <- cbind(LnLxTnTx,LxTx,LxTx.Error)
  calc.parameters <- list(sigmab.LnLx = sigmab.LnLx,
                          sigmab.TnTx = sigmab.TnTx,
                          k = k)

  ##set results object
  temp.return <-
    set_RLum(
      class = "RLum.Results",
      data = list(
        LxTx.table = temp,
        calc.parameters = calc.parameters,
        call = sys.call())
    )

  invisible(temp.return)

}
