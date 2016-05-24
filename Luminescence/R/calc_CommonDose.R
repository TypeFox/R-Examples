#' Apply the (un-)logged common age model after Galbraith et al. (1999) to a
#' given De distribution
#'
#' Function to calculate the common dose of a De distribution.
#'
#' \bold{(Un-)logged model} \cr\cr When \code{log = TRUE} this function
#' calculates the weighted mean of logarithmic De values. Each of the estimates
#' is weighted by the inverse square of its relative standard error. The
#' weighted mean is then transformed back to the dose scale (Galbraith &
#' Roberts 2012, p. 14).\cr\cr The log transformation is not applicable if the
#' De estimates are close to zero or negative. In this case the un-logged model
#' can be applied instead (\code{log = FALSE}). The weighted mean is then
#' calculated using the un-logged estimates of De and their absolute standard
#' error (Galbraith & Roberts 2012, p. 14).
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De
#' \code{(data[,1])} and De error \code{(values[,2])}
#' @param sigmab \code{\link{numeric}} (with default): spread in De values
#' given as a fraction (e.g. 0.2). This value represents the expected
#' overdispersion in the data should the sample be well-bleached (Cunningham &
#' Walling 2012, p. 100).
#' @param log \code{\link{logical}} (with default): fit the (un-)logged common
#' age model to De data
#' @param \dots currently not used.
#' @return Returns a terminal output. In addition an
#' \code{\linkS4class{RLum.Results}} object is returned containing the
#' following element:
#'
#' \item{summary}{\link{data.frame} summary of all relevant model results.}
#' \item{data}{\link{data.frame} original input data} \item{args}{\link{list}
#' used arguments} \item{call}{\link{call} the function call}
#'
#' The output should be accessed using the function
#' \code{\link{get_RLum}}
#' @section Function version: 0.1
#' @author Christoph Burow, University of Cologne (Germany)
#' @seealso \code{\link{calc_CentralDose}}, \code{\link{calc_FiniteMixture}},
#' \code{\link{calc_FuchsLang2001}}, \code{\link{calc_MinDose}}
#' @references Galbraith, R.F. & Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks Radiation Measurements 4, 459-470.
#' \cr\cr Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley,
#' J.M., 1999. Optical dating of single grains of quartz from Jinmium rock
#' shelter, northern Australia. Part I: experimental design and statistical
#' models.  Archaeometry 41, 339-364. \cr\cr Galbraith, R.F. & Roberts, R.G.,
#' 2012. Statistical aspects of equivalent dose and error calculation and
#' display in OSL dating: An overview and some recommendations. Quaternary
#' Geochronology 11, 1-27. \cr\cr \bold{Further reading} \cr\cr Arnold, L.J. &
#' Roberts, R.G., 2009. Stochastic modelling of multi-grain equivalent dose
#' (De) distributions: Implications for OSL dating of sediment mixtures.
#' Quaternary Geochronology 4, 204-230. \cr\cr Bailey, R.M. & Arnold, L.J.,
#' 2006. Statistical modelling of single grain quartz De distributions and an
#' assessment of procedures for estimating burial dose. Quaternary Science
#' Reviews 25, 2475-2502. \cr\cr Cunningham, A.C. & Wallinga, J., 2012.
#' Realizing the potential of fluvial archives using robust OSL chronologies.
#' Quaternary Geochronology 12, 98-106. \cr\cr Rodnight, H., Duller, G.A.T.,
#' Wintle, A.G. & Tooth, S., 2006. Assessing the reproducibility and accuracy
#' of optical dating of fluvial deposits. Quaternary Geochronology 1,
#' 109-120.\cr\cr Rodnight, H., 2008. How many equivalent dose values are
#' needed to obtain a reproducible distribution?. Ancient TL 26, 3-10.
#' @examples
#'
#' ## load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ## apply the common dose model
#' calc_CommonDose(ExampleData.DeValues$CA1)
#'
#' @export
calc_CommonDose <- function(
  data,
  sigmab,
  log=TRUE,
  ...
) {
  
  ##============================================================================##
  ## CONSISTENCY CHECK OF INPUT DATA
  ##============================================================================##
  
  if(missing(data)==FALSE){
    if(is(data, "data.frame") == FALSE & is(data,"RLum.Results") == FALSE){
      stop("[calc_CentralDose] Error: 'data' object has to be of type
           'data.frame' or 'RLum.Results'!")
    }else{
      if(is(data, "RLum.Results") == TRUE){
        data <- get_RLum(data, signature(object = "De.values"))
      }
    }
  }
  try(colnames(data)<- c("ED","ED_Error"), silent = TRUE)
  if(colnames(data[1])!="ED"||colnames(data[2])!="ED_Error") {
    cat(paste("Columns must be named 'ED' and 'ED_Error'"), fill = FALSE)
    stop(domain=NA)
  }
  if(!missing(sigmab)) {
    if(sigmab <0 | sigmab >1) {
      cat(paste("sigmab needs to be given as a fraction between",
                "0 and 1 (e.g. 0.2)"), fill = FALSE)
      stop(domain=NA)
    }
  }
  
  
  ##============================================================================##
  ## ADDITIONAL ARGUMENTS
  ##============================================================================##
  settings <- list(verbose = TRUE)
  settings <- modifyList(settings, list(...))
  
  ##============================================================================##
  ## CALCULATIONS
  ##============================================================================##
  
  # set default value of sigmab
  if (missing(sigmab)) sigmab<- 0
  
  # calculate  yu = log(ED) and su = se(logED)
  if (log) {
    yu<- log(data$ED)
    su<- sqrt( (data$ED_Error/data$ED)^2 + sigmab^2 )
  }
  else {
    yu<- data$ED
    su<- sqrt((data$ED_Error)^2 + sigmab^2)
  }
  
  # calculate weights
  wu<- 1/su^2
  delta<- sum(wu*yu)/sum(wu)
  n<- length(yu)
  
  #standard error
  sedelta<- 1/sqrt(sum(wu))
  if (!log) {
    sedelta<- sedelta/delta
  }
  
  if (log){
    delta<- exp(delta)
  }
  
  ##============================================================================##
  ## TERMINAL OUTPUT
  ##============================================================================##
  
  if (settings$verbose) {
    cat("\n [calc_CommonDose]")
    cat(paste("\n\n----------- meta data --------------"))
    cat(paste("\n n:                      ",n))
    cat(paste("\n log:                    ",if(log==TRUE){"TRUE"}else{"FALSE"}))
    cat(paste("\n----------- dose estimate ----------"))
    cat(paste("\n common dose:            ", round(delta,2)))
    cat(paste("\n SE:                     ", round(delta*sedelta, 2)))
    cat(paste("\n rel. SE [%]:            ", round(sedelta*100,2)))
    cat(paste("\n------------------------------------\n\n"))
  }
  
  ##============================================================================##
  ## RETURN VALUES
  ##============================================================================##
  
  summary<- data.frame(de=delta,
                       de_err=delta*sedelta)
  
  call<- sys.call()
  args<- list(log=log, sigmab=sigmab)
  
  newRLumResults.calc_CommonDose<- set_RLum(
    class = "RLum.Results",
    data = list(summary = summary,
                data = data,
                args = args,
                call = call))
  
  invisible(newRLumResults.calc_CommonDose)
  
}
