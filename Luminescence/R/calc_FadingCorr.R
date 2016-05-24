#' Apply a fading correction according to Huntley & Lamothe (2001) for a given
#' g-value.
#'
#' This function runs the iterations that are needed to calculate the corrected
#' age including the error for a given g-value according to Huntley & Lamothe
#' (2001).
#'
#' The error of the fading-corrected age is determined using a Monte Carlo
#' simulation approach. Solving of the equation is realised using
#' \code{\link{uniroot}}. Large values for \code{n.MCruns} will significantly
#' increase the computation time.\cr
#'
#' \bold{\code{n.MCruns = 'auto'}}
#'
#' The error estimation based on a stochastic process, i.e. for a small number of MC runs the calculated
#' error varies considerably every time the function is called, even with the same input values.
#' The argument option \code{n.MCruns = 'auto'} tries to find a stable value for the standard error, i.e.
#' the standard deviation of values calculated during the MC runs (\code{age.corr.MC}),
#' within a given precision (2 digits) by increasing the number of MC runs stepwise and
#' calculating the corresponding error.
#'
#' If the determined error does not differ from the 9 values calculated previously
#' within a precision of (here) 3 digits the calculation is stopped as it is assumed that the error
#' is stable. Please note that (a) the duration depends on the input values as well as on
#' the provided computation ressources and it may take a while, (b) the length (size) of the output
#' vector \code{age.corr.MC}, where all the single values produced during the MC runs are stored,
#' equals the number of MC runs (here termed observations).
#'
#' To avoid an endless loop the calculation is stopped if the number of observations exceeds 10^7.
#' This limitation can be overwritten by setting the number of MC runs manually,
#' e.g. \code{n.MCruns = 10000001}. Note: For this case the function is not checking whether the calculated
#' error is stable.\cr
#'
#'
#' \bold{\code{seed}}
#'
#' This option allows to recreate previously calculated results by setting the seed
#' for the R random number generator (see \code{\link{set.seed}} for details). This option
#' should not be mixed up with the option \bold{\code{n.MCruns = 'auto'}}. The results may
#' appear similar, but they are not comparable!
#'
#'
#' @param g_value \link{vector} (\bold{required}): g-value and error obtained
#' from separate fading measurements (see example)
#'
#' @param tc \link{numeric} (\bold{required}): time in seconds (time between
#' irradiation and the prompt measurement, cf. Huntely & Lamothe 2001)
#'
#' @param age.faded \link{numeric} \link{vector} (\bold{required}): uncorrected
#' age with error in ka (see example)
#'
#' @param n.MCruns \link{integer} (with default): number of Monte Carlo
#' simulation runs for error estimation. If \code{n.MCruns = 'auto'} is used the function
#' tries to find a 'stable' error for the age. Note: This may take a while!
#'
#' @param seed \link{integer} (optional): sets the seed for the random number generator
#' in R using \code{\link{set.seed}}
#'
#' @param txtProgressBar \link{logical} (with default): enables or disables
#' \code{\link{txtProgressBar}}
#'
#'
#' @return Returns an S4 object of type \code{\linkS4class{RLum.Results}}. Slot
#' \code{data} contains a \code{\link{list}} with the following structure:\cr
#' $ age.corr (data.frame) \cr
#' .. $ age \cr
#' .. $ age.error \cr
#' .. $ age.faded \cr
#' .. $ age.faded.error \cr
#' .. $ g_value \cr
#' .. $ g_value.error \cr
#' .. $ tc \cr
#' .. $ n.MCruns \cr
#' .. $ observations \cr
#' .. $ seed \cr
#' $ age.corr.MC (numeric)\cr
#'
#' \code{Age.corr.MC} contain all possible ages from the Monte Carlo (error)
#' simulation.
#'
#'
#' @note The upper age limit is set to 500 ka!
#'
#'
#' @section Function version: 0.3.0
#'
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#'
#'
#' @seealso \code{\linkS4class{RLum.Results}}, \code{\link{get_RLum}},
#' \code{\link{uniroot}}
#'
#'
#' @references Huntley, D.J., Lamothe, M., 2001. Ubiquity of anomalous fading
#' in K-feldspars and the measurement and correction for it in optical dating.
#' Canadian Journal of Earth Sciences, 38, 1093-1106.
#'
#'
#' @keywords datagen
#'
#'
#' @examples
#' results <- calc_FadingCorr(g_value = c(3.3,0.03), tc = 752,
#'                 age.faded = c(100,10),
#'                 n.MCruns=100)
#'
#' get_RLum(results)
#'
#' @export
calc_FadingCorr <- function(
  g_value,
  tc,
  age.faded,
  n.MCruns = 10000,
  seed,
  txtProgressBar = TRUE
){
  
  ##============================================================================##
  ##DEFINE FUNCTION
  ##============================================================================##
  
  f <- function(x, af,kappa,tc){1-kappa*(log(x/tc)-1) - (af/x)}
  
  ##============================================================================##
  ##CALCULATION
  ##============================================================================##
  
  ##calculate kappa
  kappa <- g_value[1] / log(10) / 100
  
  ##transform tc in ka years
  tc <- tc / 60 / 60 / 24 / 365 / 1000
  
  ##calculate mean value
  temp <- uniroot(f, c(0.1,500), tol = 0.001, tc = tc, af = age.faded[1], kappa = kappa,
                  check.conv = TRUE)
  
  ##--------------------------------------------------------------------------##
  ##Monte Carlo simulation for error estimation
  tempMC.sd.recent <- NA
  tempMC.sd.count <- 1:10
  counter <- 1
  
  ##show some progression bar of the process
  if (n.MCruns == 'auto') {
    n.MCruns.i <- 10000
    
    cat("\n[calc_FadingCorr()] ... trying to find stable error value ...")
    if (txtProgressBar) {
      cat("\n -------------------------------------------------------------\n")
      cat(paste0("   ",paste0("(",0:9,")", collapse = "   "), "\n"))
    }
  }else{
    n.MCruns.i <- n.MCruns
    
  }
  
  
  
  # Start loop  ---------------------------------------------------------------------------------
  
  ##set object and preallocate memory
  tempMC <- vector("numeric", length = 1e+07)
  tempMC[] <- NA
  i <- 1
  j <- n.MCruns.i
  
  while(length(unique(tempMC.sd.count))>1 | j > 1e+07){
    
    ##set previous
    if(!is.na(tempMC.sd.recent)){
      tempMC.sd.count[counter] <- tempMC.sd.recent
      
    }
    
    ##set seed
    if (!missing(seed)) {
      seed.set <- seed
      set.seed(seed)
      
    }else{
      seed.set <- "not set"
      
    }
    
    ##pre-allocate memory
    g_valueMC <- vector("numeric", length = n.MCruns.i)
    age.fadeMC <- vector("numeric", length = n.MCruns.i)
    kappaMC <- vector("numeric", length = n.MCruns.i)
    
    ##set-values
    g_valueMC <- rnorm(n.MCruns.i,mean = g_value[1],sd = g_value[2])
    age.fadedMC <- rnorm(n.MCruns.i,mean = age.faded[1],sd = age.faded[2])
    kappaMC <- g_valueMC / log(10) / 100
    
    
    ##calculate for all values
    tempMC.i <- lapply(1:length(age.fadedMC),function(x) {
      uniroot(
        f,
        c(0.1,500),
        tol = 0.001,
        tc = tc,
        af = age.fadedMC[[x]],
        kappa = kappaMC[[x]],
        check.conv = TRUE
      )$root
      
    })
    
    ##write values in vector
    tempMC[i:j] <- unlist(tempMC.i)
    i <- j + 1
    j <- j + n.MCruns.i
    
    ##stop here if a fixed value is set
    if(n.MCruns != 'auto'){
      break
    }
    
    ##set recent
    tempMC.sd.recent <- round(sd(tempMC, na.rm = TRUE), digits = 3)
    
    if (counter %% 10 == 0) {
      counter <- 1
      
    }else{
      counter <- counter + 1
      
    }
    
    ##show progress in terminal
    if (txtProgressBar) {
      text <- rep("CHECK",10)
      if (counter %% 2 == 0) {
        text[1:length(unique(tempMC.sd.count))] <- "-----"
      }else{
        text[1:length(unique(tempMC.sd.count))] <- " CAL "
      }
      
      
      
      cat(paste("\r ",paste(rev(text), collapse = " ")))
    }
    
  }
  
  ##--------------------------------------------------------------------------##
  
  ##remove all NA values from tempMC
  tempMC <- na.exclude(tempMC)
  
  ##obtain corrected age
  age.corr <- data.frame(age = round(temp$root, digits = 2),
                         age.error = round(sd(tempMC),digits = 2),
                         age.faded = age.faded[1],
                         age.faded.error = age.faded[2],
                         g_value = g_value[1],
                         g_value.error = g_value[2],
                         tc = tc,
                         n.MCruns = n.MCruns,
                         observations = length(tempMC),
                         seed = seed.set)
  
  ##============================================================================##
  ##OUTPUT VISUAL
  ##============================================================================##
  
  cat("\n\n[calc_FadingCorr()]\n")
  cat("\n\t Fading correction according to Huntley & Lamothe (2001):\n")
  cat(paste("\n\t Age (faded): ",age.faded[1]," ka \u00b1 ",
            age.faded[2]," ka",sep=""))
  cat(paste("\n\t g-value: ",g_value[1], "%/decade \u00b1 ",
            g_value[2]," %/decade",sep=""))
  cat(paste("\n\t tc: ",format(tc, digits = 4, scientific = TRUE), " ka",sep=""))
  cat(paste("\n\t kappa: ",mean(kappa),sep=""))
  cat(paste("\n\t seed: ",seed.set))
  cat(paste("\n\t n.MCruns: ",n.MCruns))
  cat(paste("\n\t observations: ",
            format(length(tempMC), digits = 2, scientific =TRUE),sep=""))
  cat("\n\n\t ----------------------------------")
  cat(paste("\n\t Age (corr.): ",age.corr[1]," ka \u00b1 ",age.corr[2]," ka",sep=""))
  cat("\n\t ----------------------------------\n")
  
  ##============================================================================##
  ##OUTPUT RLUM
  ##============================================================================##
  
  temp.RLum.Results <- set_RLum(
    class = "RLum.Results",
    data = list(age.corr = age.corr,
                age.corr.MC = tempMC))
  
  return(temp.RLum.Results)
  
}
