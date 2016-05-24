### Thinning

#' Thinning of the observations (for evaluating the method)
#'
#' The thinning is done by drawing independently from a Bernoulli distribution. This function is needed for functions eval_method, sim_clintens, sim_intens
#' 
#' @export
#' @param full  all observations of the point pattern
#' @param nxprob  probability of having unobserved events 
#' @return A list of observed and unobserved point patterns. Both of class ppp.
#' @importFrom stats rbinom
#' @seealso \code{\link[stats]{rbinom}}, \code{\link[spatstat]{ppp}}
#' @examples
#'  data(craterB)
#'  thdata <- thin(craterB, nxprob=0.1)
#'  thdata
#'  plot(thdata$observed); points(thdata$unobserved, col=4)

thin <- function(full, nxprob){
  win <- full$window
  
  indicator <- rbinom(n=full$n, size=1, prob=nxprob)
  
  observed <- full[indicator == 0]
  unobserved <- full[indicator == 1]
  
  return(list(observed=observed, unobserved=unobserved))
}


### for intens


#' Simulation on given intensity
#' 
#' Generation of a random point pattern using the inhomogeneous Poisson process (if lambda is not constant)
#' and thinning of this data, to obtain "observed" and "unobserved" events.
#'
#' @param ppdata  Observed spatial point process of class ppp
#' @param intensSim  Intensity to use for the simulation
#' @param nxprob  Probability of having unobserved events
#' @return A list of of observed and unobserved point patterns (see \code{\link{thin}})
#' @seealso \code{\link{thin}}, \code{\link[spatstat]{rpoispp}}



sim_intens <- function(ppdata, intensSim, nxprob) {
  ppsim <- rpoispp(lambda=intensSim)
  ppsim$window <- ppdata$window
  
  # Warnings  
  if ( max(ppdata$window$xrange[2] - ppdata$window$xrange[1], ppdata$window$yrange[2] - ppdata$window$yrange[1])
       /spatstat::spatstat.options()$npixel > 4000 ) {
    warning( "For given size of ppdata$window the number of pixels are too low. Set higher by using spatstat::spatstat.options(npixel).")
  }
  
  
  ok <- inside.owin(ppsim$x, ppsim$y, ppsim$window)
  nout <- sum(!ok)
  if ( nout > 0 ) {
    warning(paste(nout, ngettext(nout, "point was", "points were"), 
                  "rejected as lying outside the specified window. Set number of pixels higher."))
    rr <- ripras(ppsim$x, ppsim$y)
    bb <- bounding.box.xy(ppsim$x, ppsim$y)
    bb <- boundingbox(rr, bb, ppsim$window)
    rejectwindow <- if (!is.null(rr)) 
      rebound.owin(rr, bb)
    else bb
    rejects <- ppp(ppsim$x[!ok], ppsim$y[!ok], window = rejectwindow, 
                   check = FALSE)
    ppsim$x <- ppsim$x[ok]
    ppsim$y <- ppsim$y[ok]
    ppsim$n <- length(ppsim$x)
  }
  #
  
  
  thinned <- thin(full=ppsim, nxprob=nxprob)
  thinned
}








### for clintens


#' Determination of the intensity for the Neyman Scott simulation.
#'
#' Used in function sim_nsppp.
#'
#' @param ppdata  observed point pattern whose estimated intensity (adjusted for
#'          thinning and divided by "clustering") is used for simulating the
#'          parent process
#' @param radius  radius of the circles around the parent points in which the cluster
#'          points are located
#' @return A pixel image (object of class "im"). See \code{\link[spatstat]{density.ppp}}.  
#' @seealso \code{\link[spatstat]{density.ppp}}, \code{\link[spatstat]{boundingbox}}, 
#'          \code{\link[spatstat]{owin}}, \code{\link[ks]{Hscv}}




det_nsintens <- function(ppdata, radius){
  
  pppsim <- ppdata
  frameWin <- boundingbox(pppsim$window)
  dilatedWin <- owin(frameWin$xrange + 1.2*c(-radius, radius), frameWin$yrange + 1.2*c(-radius, radius))
  pppsim$window <- dilatedWin
  
  covmatrixsim <- Hscv(cbind(pppsim$x, pppsim$y))
  intensest <- density(pppsim, varcov=covmatrixsim, positive=TRUE)
  
  return(intensest)
}


#' Simulation of the Neyman-Scott process.  
#'
#' Simulation of the Neyman-Scott process. Only applicable if the intensity was estimated
#' for an appropriately enlarged window.
#' More details in \code{sim_nsppp}.
#'
#' @param ppdata  observed point pattern whose estimated intensity (adjusted for
#'          thinning and divided by "clustering") is used for simulating the
#'          parent process
#' @param intens  estimated intensity
#' @param radius  radius of the circles around the parent points in which the cluster
#'          points are located (Maximum radius of a random cluster)
#' @param clustering  a value larger or equal 1 which describes the amount of clustering; the
#'          adjusted estimated intensity of the observed pattern is divided by
#'          this value; it is also the parameter of the Poisson distribution
#'          for the number of points per cluster
#' @param thinning  constant thinning probability (in case the observed pattern is a
#'          thinned version of a full pattern); usually equal to the probability of having
#'          unobserved events
#' @importFrom stats rpois
#' @return The simulated point pattern (an object of class "ppp").
#'        Additionally, some intermediate results of the simulation are returned as 
#'        attributes of this point pattern: see \code{\link[spatstat]{rNeymanScott}}.



sim_nsprocess <- function(ppdata, intens, radius, clustering=5, thinning=0){
  
  # there are no warnings needed here, because rNeymanScott is doing all the work
  
  intenssim <- intens
  intenssim$v <- (1/((1 - thinning)*clustering))*intens$v
  
  nclust <- function(x0, y0, radius, lambdaclust) {
    n <- rpois(1, lambdaclust)
    return(runifdisc(n, radius, centre=c(x0, y0)))
  }
#  resultppp <- rNeymanScott(kappa=intenssim, rmax=radius, rcluster=nclust, radius=radius, 
#                            lambdaclust=clustering, win=ppdata$window)
  resultppp <- rNeymanScott(kappa=intenssim, expand=radius, rcluster=nclust, radius=radius, 
                            lambdaclust=clustering, win=ppdata$window)

    return(resultppp)
}


#' Generation of a realisation of a Neyman-Scott process
#'
#' This algorithm generates a realisation of a Neyman-Scott process whose
#' expected number of points equals the number of observations in a given
#' pattern. 
#'
#' First, the algorithm generates a Poisson point process (see \code{\link[spatstat]{rpoispp}} for
#' details) of parent points with intensity kappa, which is a pixel image
#' object of class "im" (see \code{\link[spatstat]{im.object}}).\cr
#' This pixel image is derived from the observed pattern using \code{\link[spatstat]{density.ppp}}.
#' The bandwidth is not chosen in advance.\cr
#' If only a thinned version of the original pattern has been observed,
#' this can be taken into account using the parameter \code{thinning}.
#' Usually, not the estimated intensity itself is used for simulating the
#' parent process, but its values are divided by a constant named "clustering".\cr
#' Second, each parent point is replaced by a random cluster of points, created
#' by calling the function \code{\link[spatstat]{runifdisc}}. Each cluster consists of a Poisson
#' distributed number of points (with \code{clustering} being the expected number of
#' points in each cluster) which are located in a disc of a given \code{radius}.
#' These clusters are combined to yield a single point pattern which is
#' then returned as the result.\cr
#' The estimation of the intensity (on an adequate window) and the
#' simulation of the Neyman-Scott process are performed seperately,
#' so the intensity does not need to be reestimated in every iteration.\cr
#' The resulting process is a \enc{Mat?rn}{Matern} process whose parent process is an
#' inhomogeneous Poisson point process.
#'
#' @param ppdata  observed point pattern, whose estimated intensity (adjusted for
#'          thinning and divided by "clustering") is used for simulating the
#'          parent process
#' @param radius  radius of the circles around the parent points in which the cluster
#'          points are located (Maximum radius of a random cluster)
#' @param clustering  a value larger or equal 1 which describes the amount of clustering; the
#'          adjusted estimated intensity of the observed pattern is divided by
#'          this value; it is also the parameter of the Poisson distribution
#'          for the number of points per cluster
#' @param thinning  constant thinning probability (in case the observed pattern is a
#'          thinned version of a full pattern); usually equal to the probability of having
#'          unobserved events
#' @export  
#' @return The simulated point pattern (an object of class "ppp").
#'        Additionally, some intermediate results of the simulation are returned as 
#'        attributes of this point pattern: see \code{\link[spatstat]{rNeymanScott}}.
#' @seealso \code{\link[spatstat]{rNeymanScott}}, \code{\link[spatstat]{rThomas}}, 
#'          \code{\link[spatstat]{rMatClust}}
#' @examples      
#'  data(craterA)
#'  data(craterB)
#'  set.seed(100)
#'  sim_pp1 <- sim_nsppp(craterA, radius=300, clustering=15, thinning=0.1)
#'  sim_pp2 <- sim_nsppp(craterB, radius=300, clustering=15, thinning=0.1)
#'  op <- par(mfrow = c(1, 2))
#'  plot(sim_pp1, main = "simulated cluster process 1")
#'  plot(sim_pp2, main = "simulated cluster process 2")
#'  par(op)

sim_nsppp <- function ( ppdata, radius, clustering=5, thinning=0 ) {
  intensest <- det_nsintens( ppdata = ppdata, radius = radius)
  resultppp <- sim_nsprocess( ppdata = ppdata, intens = intensest, radius = radius, clustering = clustering, thinning = thinning)
  return(resultppp)
}
