#' Bootstrap correction to obtain desired failure probability
#'
#' Simulation-based iterative procedure to correct for possible bias with respect to the
#' failure probability alpha
#'
#' For a desired failure probability alpha, the corresponding parameter which is to use
#' when determining a high-risk zone is found in an iterative procedure. The simulation procedure
#' is the same as in \code{\link[highriskzone]{eval_method}}. In every iteration,
#' the number of high-risk zones with at least one unobserved event located outside is
#' compared with the desired failure probability. If necessary, the value of \code{cutoff} is
#' increased or decreased. The final value \code{alphastar} can than be used in
#' \code{\link[highriskzone]{det_hrz}}.
#'
#' If there are restriction areas in the observation window, use \code{\link[highriskzone]{bootcor_restr}}
#' instead.
#'
#' @param ppdata Observed spatial point process of class ppp.
#' @param cutoff Desired failure probability alpha, which is the probability of having
#'                 unobserved events outside the high-risk zone.
#' @param numit Number of iterations to perform (per tested value for cutoff). Default value is 1000.
#' @param tol Tolerance: acceptable difference between the desired failure probability and the fraction of
#'             high-risk zones not covering all events. Default value is 0.02.
#' @param nxprob Probability of having unobserved events.
#'                Default value is 0.1.
#' @param intens (optional) estimated intensity of the observed process (object of class "im",
#'                see \code{\link[spatstat]{density.ppp}}). If not given,
#'                it will be estimated.
#' @param covmatrix  (optional) Covariance matrix of the kernel of a normal distribution, only meaningful
#'                    if no intensity is given. If not given, it will be estimated.
#' @param simulate The type of simulation, can be one of \code{"thinning", "intens"} or \code{"clintens"}
#' @param radiusClust (optional) radius of the circles around the parent points in which the cluster
#'                    points are located. Only used for \code{simulate = "clintens"}.
#' @param clustering a value >= 1 which describes the amount of clustering; the
#'          adjusted estimated intensity of the observed pattern is divided by
#'          this value; it also is the parameter of the Poisson distribution
#'          for the number of points per cluster. Only used for \code{simulate = "clintens"}.
#' @param verbose logical. Should information on tested values/progress be printed?
#' @importFrom stats quantile
#' @export
#' @return An object of class bootcorr, which consists of a list of the final value for alpha (\code{alphastar})
#'         and a data.frame \code{course} containing information on the simulation course, e.g. the tested values.
#' @references Monia Mahling, Michael \enc{H?hle}{Hoehle} & Helmut \enc{K?chenhoff}{Kuechenhoff} (2013),
#' \emph{Determining high-risk zones for unexploded World War II bombs by using point process methodology.}
#' Journal of the Royal Statistical Society, Series C 62(2), 181-199.
#' @references Monia Mahling (2013),
#' \emph{Determining high-risk zones by using spatial point process methodology.}
#' Ph.D. thesis, Cuvillier Verlag \enc{G?ttingen}{Goettingen},
#' available online: http://edoc.ub.uni-muenchen.de/15886/
#' Chapter 6
#' @seealso \code{\link[highriskzone]{det_hrz}}, \code{\link[highriskzone]{eval_method}}, \code{\link[highriskzone]{bootcor_restr}}
#' @examples
#' \dontrun{
#' data(craterB)
#' set.seed(4321)
#'
#' bc <- bootcor(ppdata=craterB, cutoff=0.2, numit=100, tol=0.02, nxprob=0.1)
#' bc
#' summary(bc)
#' plot(bc)
#'
#' 
#' hrzbc <- det_hrz(craterB, type = "intens", criterion = "indirect", 
#' cutoff = bc$alphastar, nxprob = 0.1)
#' }

bootcor <- function(ppdata, cutoff, numit = 1000, tol=0.02,
                        nxprob = 0.1, intens = NULL,
                        covmatrix = NULL, simulate="intens", radiusClust=NULL, clustering=5, verbose=TRUE) {

  #check if input arguments have correct values
  roundnumit <- round(numit)
  if ( roundnumit != numit ) {
    warning("numit must be a natural number. It is now rounded to: ", roundnumit)
    numit <- roundnumit
  }
  match.arg(simulate, choices=c("thinning", "intens", "clintens"))


  #here the intensity is being estimated
  if ( simulate == "intens" ) {

    origintens <- est_intens(ppdata, covmatrix=covmatrix)
    intensSim <- origintens$intensest
    intensSim$v <- (1/(1 - nxprob))*origintens$intensest$v

  } else if ( simulate == "clintens" ) {

    if( is.null(radiusClust) ) {
      radiusClust <- quantile(nndist(ppdata), p=0.7, type=8)
    }
    intensSim <- det_nsintens(ppdata=ppdata, radius=radiusClust)
  }

  result <- matrix(data=NA, nrow=0, ncol=6)
  
  numout <- 0
  i <- 1
  k <- 1
  alphastar <- cutoff


  while(i <= numit){

    if ( simulate == "thinning" ) {
      thinned <- thin(full=ppdata, nxprob=nxprob)
      observed <- thinned$observed
      unobserved <- thinned$unobserved
    }
    if ( simulate == "intens" ) {
      thinned <- sim_intens(ppdata, intensSim, nxprob)
      observed <- thinned$observed
      unobserved <- thinned$unobserved
    }
    if ( simulate == "clintens" ) {
      ppsim <- sim_nsprocess(ppdata=ppdata, intens=intensSim, radius=radiusClust,
                             clustering=clustering, thinning=nxprob)
      thinned <- thin(full=ppsim, nxprob=nxprob)
      observed <- thinned$observed
      unobserved <- thinned$unobserved
    }


    if (is.null(intens)){
      estim <- est_intens(observed, covmatrix=covmatrix)
      intens <- estim$intensest
      covmatrix <- estim$covmatrix
    }

    resultdetHRZ <- det_hrz(ppdata=observed, type="intens", criterion="indirect",
                              cutoff=alphastar, intens=intens,
                              nxprob=nxprob, covmatrix=covmatrix)
    resultevalHRZ <- eval_hrz(hrz=resultdetHRZ$zone, unobspp=unobserved, obspp=observed)

    if(resultevalHRZ$numbermiss > 0){
      numout <- numout + 1
    }


    poutmin <- numout/numit
    poutmax <- (numout + numit - i)/numit

    resstep <- c(k, i, alphastar, numout, poutmin, poutmax)
    result <- rbind(result, resstep)

    if(poutmin > cutoff + tol){
      alphastar <- alphastar * i/(numit + 1)
      if(verbose){
        cat("Decrease alphastar to ", alphastar, " after ", i, " iterations with numout=", numout , "\n", sep="")
      }
      i <- 0
      numout <- 0
      k <- k + 1
    }

    if(poutmax < cutoff - tol){
      alphastar <- alphastar * (1 + (numit - i + 1)/numit)
      if(verbose){
        cat("Increase alphastar to ", alphastar, " after ", i, " iterations with numout=", numout , "\n", sep="")
      }
      i <- 0
      numout <- 0
      k <- k + 1
    }

    i <- i + 1

  }
  
  resultdf <- as.data.frame(result, row.names = as.character(1:dim(result)[1]))
  colnames(resultdf) <- c("k", "i", "alphastar", "numout", "poutmin", "poutmax")

  res <- list(alphastar=alphastar, course=resultdf)
  
  class(res) <- "bootcorr"
  return(res)

}

#' Bootstrap correction to obtain desired failure probability
#'
#' Simulation-based iterative procedure to correct for possible bias with respect to the
#' failure probability alpha
#'
#' For a desired failure probability alpha, the corresponding parameter which is to use
#' when determining a high-risk zone is found in an iterative procedure. The simulation procedure
#' is the same as in \code{\link[highriskzone]{eval_method}}. In every iteration,
#' the number of high-risk zones with at least one unobserved event located outside is
#' compared with the desired failure probability. If necessary, the value of \code{cutoff} is
#' increased or decreased. The final value \code{alphastar} can than be used in
#' \code{\link[highriskzone]{det_hrz}}.
#'
#' If there are restriction areas in the observation window, use \code{\link[highriskzone]{bootcor_restr}}
#' instead.
#'
#' @param ppdata Observed spatial point process of class ppp.
#' @param cutoff Desired failure probability alpha, which is the probability of having
#'                 unobserved events outside the high-risk zone.
#' @param numit Number of iterations to perform (per tested value for cutoff). Default value is 1000.
#' @param tol Tolerance: acceptable difference between the desired failure probability and the fraction of
#'             high-risk zones not covering all events. Default value is 0.02.
#' @param nxprob Probability of having unobserved events.
#'                Default value is 0.1.
#' @param intens (optional) estimated intensity of the observed process (object of class "im",
#'                see \code{\link[spatstat]{density.ppp}}). If not given,
#'                it will be estimated.
#' @param covmatrix  (optional) Covariance matrix of the kernel of a normal distribution, only meaningful
#'                    if no intensity is given. If not given, it will be estimated.
#' @param simulate The type of simulation, can be one of \code{"thinning", "intens"} or \code{"clintens"}
#' @param radiusClust (optional) radius of the circles around the parent points in which the cluster
#'                    points are located. Only used for \code{simulate = "clintens"}.
#' @param clustering a value >= 1 which describes the amount of clustering; the
#'          adjusted estimated intensity of the observed pattern is divided by
#'          this value; it also is the parameter of the Poisson distribution
#'          for the number of points per cluster. Only used for \code{simulate = "clintens"}.
#' @param verbose logical. Should information on tested values/progress be printed?
#' @export
#' @return An object of class bootcorr, which consists of a list of the final value for alpha (\code{alphastar})
#'         and a data.frame \code{course} containing information on the simulation course, e.g. the tested values.
#' @references Monia Mahling, Michael \enc{H?hle}{Hoehle} & Helmut \enc{K?chenhoff}{Kuechenhoff} (2013),
#' \emph{Determining high-risk zones for unexploded World War II bombs by using point process methodology.}
#' Journal of the Royal Statistical Society, Series C 62(2), 181-199.
#' @references Monia Mahling (2013),
#' \emph{Determining high-risk zones by using spatial point process methodology.}
#' Ph.D. thesis, Cuvillier Verlag \enc{G?ttingen}{Goettingen},
#' available online: http://edoc.ub.uni-muenchen.de/15886/
#' Chapter 6
#' @seealso \code{\link[highriskzone]{det_hrz}}, \code{\link[highriskzone]{eval_method}}, \code{\link[highriskzone]{bootcor_restr}}
#' @examples
#' \dontrun{
#' data(craterB)
#' set.seed(4321)
#'
#' bc <- bootcor(ppdata=craterB, cutoff=0.2, numit=100, tol=0.02, nxprob=0.1)
#' bc
#' summary(bc)
#' plot(bc)
#'
#' 
#' hrzbc <- det_hrz(craterB, type = "intens", criterion = "indirect", 
#' cutoff = bc$alphastar, nxprob = 0.1)
#' }
bootcorr <- bootcor
