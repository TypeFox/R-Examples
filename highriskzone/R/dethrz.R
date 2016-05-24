#' Determination of the high-risk zone.
#'
#' \code{det_hrz} determines the high-risk zone through the method of fixed radius 
#' (type = "dist" and criterion = "direct"), the quantile-based method (type = "dist" and 
#' criterion = "area"/"indirect") and the intensity-based method (type =  "intens").
#'
#' There are different methods implemented to determine a high-risk zone.
#' \describe{
#' \item{ Method of fixed radius }{
#'        In this method, the high-risk zone is determined by drawing a circle around each
#'        observed event with a fixed radius. This method will be used when \code{type = "dist"}
#'        and \code{criterion = "direct"}. \code{cutoff} then is the radius.
#'        }
#' \item{ Quantile-based method }{
#'        This method is a development of the above. Here the radius is not fixed. It uses
#'        the distance of every observed event to the nearest other event, which is calculated by the
#'        nearest-neighbour distance. The radius is assessed by the p-quantile of the empirical
#'        distribution function of the nearest-neighbour distance. This method will be used when
#'        \code{type = "dist"} and \code{criterion = "indirect"} or \code{"area"}. If 
#'        \code{criterion = "indirect"}, then \code{cutoff} is the quantile that should be used.
#'        If \code{criterion = "area"} then \code{cutoff} is the area that the  high-risk zone 
#'        has to have at the end and from that the quantile/the radii are determined. When the
#'        calculation is done via the area, it can not really be classified to the quantile-based
#'        method. It is rather a third "distance-based" method.
#'        }
#' \item{ Intensity-based method }{
#'        The first step of this method is to estimate the intensity of the observed events.
#'        The high-risk zone is then the field in which the estimated intensity exceeds a 
#'        certain  value. This value is called threshold c.
#'        The method will be used when \code{type = "intens"}. There are three different ways to
#'        get to a high-risk zone: 
#'        \enumerate{
#'           \item Fixing the threshold c: \code{criterion = "direct"} 
#'           \item Fixing the area of the high-risk zone: \code{criterion = "area"}
#'           \item Fixing the failure probability alpha, which is the probability of having 
#'                 unobserved events outside the high-risk zone: \code{criterion = "indirect"}
#'                 Here, the point process is assumed to be an inhomogeneous Poisson process.
#'          }
#'        For further information see Mahling et al. (2013) (References).
#'        }
#' }
#'
#' If there are restriction areas in the observation window, use \code{\link[highriskzone]{det_hrz_restr}}
#' instead.
#'
#' @param ppdata  Observed spatial point process of class ppp.
#' @param type  Method to use, can be one of \code{"dist"} (method of fixed radius or quantile-based method), or
#'              \code{"intens"} (intensity-based method)
#' @param criterion  criterion to limit the high-risk zone, can be one of 
#'        \code{"area"} (giving size of hrz), \code{"indirect"} (giving quantile/alpha depending on type),
#'        or \code{"direct"} (giving radius/threshold c depending on type) 
#' @param cutoff  Value of criterion (area, radius, quantile, alpha or threshold). 
#'                 Depending on criterion and type: 
#'                 If criterion = "direct", cutoff is the threshold. If criterion = "indirect", cutoff is
#'                 the quantile for the quantile-based method and the failure probability alpha for the 
#'                 intensity-base method. If criterion = "area", cutoff is the area the high-risk zone should
#'                 have.
#' @param distancemap  (optional) distance map: distance of every pixel to the nearest observation 
#'                     of the point pattern; only needed for \code{type="dist"}. If not given, 
#'                     it will be computed by \code{\link[spatstat]{distmap}}.
#' @param intens  (optional) estimated intensity of the observed process (object of class "im"), 
#'                only needed for type="intens". If not given,
#'                it will be estimated using \code{\link[spatstat]{density.ppp}}.
#' @param nxprob  Probability of having unobserved events.
#'                Default value is 0.1.
#' @param covmatrix  (optional) Covariance matrix of the kernel of a normal distribution, only needed for 
#'                    \code{type="intens"} if no intensity is given. If not given, it will be estimated
#'                    using \code{\link[ks]{Hscv}}. 
#' @export  
#' @aliases highriskzone.object highriskzone
#' @return An object of class "\code{highriskzone}", which is a list of
#'    \item{ typehrz, criterion, cutoff }{ see arguments}
#'    \item{ zone }{ Determined high-risk zone: Object of class "owin" based on a binary mask. 
#'                   See \code{\link[spatstat]{owin}}. }
#'    \item{ threshold }{ determined threshold. If criterion="area", it is either the distance (if type="dist")
#' or the threshold c (for type="intens"). If criterion="indirect", it is either the quantile of the
#' nearest-neighbour distance which is used as radius (if type="dist") or the threshold c (for type="intens"). If criterion="direct",
#' it equals the cutoff for both types.}
#'    \item{ calccutoff }{ determined cutoff-value. For type="dist" and criterion="area", this is the
#' quantile of the nearest-neighbour distance. For type="intens" and criterion="area", it is the failure
#' probability alpha. For all other criterions it is NA.}
#'    \item{ covmatrix }{ If not given (and \code{type="intens"}), it is estimated. See \code{\link[ks]{Hscv}}.}
#' @references Monia Mahling, Michael \enc{H?hle}{Hoehle} & Helmut \enc{K?chenhoff}{Kuechenhoff} (2013),
#' \emph{Determining high-risk zones for unexploded World War II bombs by using point process methodology.}
#' Journal of the Royal Statistical Society, Series C 62(2), 181-199.
#' @references Monia Mahling (2013),
#' \emph{Determining high-risk zones by using spatial point process methodology.}
#' Ph.D. thesis, Cuvillier Verlag \enc{G?ttingen}{Goettingen},
#' available online: http://edoc.ub.uni-muenchen.de/15886/
#' @seealso \code{\link[spatstat]{distmap}}, \code{\link[spatstat]{eval.im}}, \code{\link[spatstat]{owin}},
#'                \code{\link{eval_method}}, \code{\link[highriskzone]{det_hrz_restr}}
#' @examples
#'  data(craterA)
#'  spatstat::spatstat.options(npixel=400)
#' ## type: dist
#' hrzd1 <- det_hrz(craterA, type = "dist", criterion = "area", cutoff = 1000000, nxprob = 0.1)
#' hrzd2 <- det_hrz(craterA, type = "dist", criterion = "indirect", cutoff = 0.9, nxprob = 0.1)
#' hrzd3 <- det_hrz(craterA, type = "dist", criterion = "direct", cutoff = 100, nxprob = 0.1)
#' 
#' op <- par(mfrow = c(2, 2))
#' plot(craterA)
#' plot(hrzd1, zonecol = 2, win = craterA$window, plotwindow = TRUE)
#' plot(hrzd2, zonecol = 3,  win = craterA$window, plotwindow = TRUE)
#' plot(hrzd3, zonecol = 4,  win = craterA$window, plotwindow = TRUE)
#' par(op)
#' 
#' \dontrun{
#' # or first calculate the distancemap and use it:
#' distm <- distmap(craterA)
#' hrzd <- det_hrz(craterA, type = "dist", criterion = "direct", cutoff = 100,
#'                 distancemap = distm, nxprob = 0.1)
#'                 
#' ## type: intens
#' hrzi1 <- det_hrz(craterA, type = "intens", criterion = "area", cutoff = 1000000, nxprob = 0.1)
#' hrzi2 <- det_hrz(craterA, type = "intens", criterion = "indirect", cutoff = 0.1, nxprob = 0.1)
#' hrzi3 <- det_hrz(craterA, type = "intens", criterion = "direct", cutoff = 0.0001, nxprob = 0.1)
#' }
#'                  
#' ## More detailed examples on http://highriskzone.r-forge.r-project.org/



det_hrz <- function(ppdata, type, criterion, cutoff, distancemap=NULL, intens=NULL, nxprob=0.1, covmatrix=NULL){
  
  win <- ppdata$window
  calccutoff <- NA
  
  
  
  # set the right values
  if(type=="dist"){
    
    if(is.null(distancemap)){distancemap <- distmap(ppdata)}
    
    if(criterion=="area"){
      
      res_det_radius <- det_radius(ppdata=ppdata, distancemap=distancemap, areahrz=cutoff, win=win)
      threshold <- res_det_radius$thresh
      calccutoff <- res_det_radius$cutoffdist
      
    }else{
      
      threshold <- ifelse(criterion=="indirect", quantile(nndist(ppdata), p=cutoff, type=8), cutoff)
      
    }
    
    HRZimage <- eval.im(distancemap < threshold)
    
  }
  
  if(type=="intens"){
    
    if(is.null(intens)){
      
      estim <- est_intens(ppdata, covmatrix=covmatrix)
      intens <- estim$intensest
      covmatrix <- estim$covmatrix
      
    }
    
    if(criterion=="area"){
      
      res_det_thresholdfromarea <- det_thresholdfromarea(win=win, intens=intens, nxprob=nxprob, areahrz=cutoff)
      threshold <- res_det_thresholdfromarea$thresh
      calccutoff <- res_det_thresholdfromarea$calccutoff
      
    }else{
      
      threshold <- ifelse(criterion=="indirect", det_threshold(intens=intens, alpha=cutoff, nxprob=nxprob), cutoff)
      
    }
    
    HRZimage <- eval.im(intens > threshold)
    
  }
  
  Rwindow <- owin(xrange=win$xrange, yrange=win$yrange, mask=as.matrix(HRZimage))
  result <- list(typehrz=type, criterion=criterion, cutoff=cutoff, zone=Rwindow, threshold=threshold, calccutoff=calccutoff, covmatrix=covmatrix)
  class(result) <- "highriskzone"
  return(result)
  
}





#' Checks the arguments of det_hrz
#' 
#' For each argument it is checked if it is of a correct value or class.
#' 
#' @inheritParams det_hrz
#' @seealso \code{\link{det_hrz}}


check_det_hrz_input <- function(ppdata, type, criterion, cutoff, distancemap, intens, nxprob, covmatrix){
  
  #errors
  #check if arguments have correct values
  if ( !is.ppp(ppdata) ) {
    stop("data is not of class ppp")
  }
  
  match.arg( type, choices=c("intens", "dist") )  
  
  match.arg( criterion, choices=c("area", "indirect", "direct") )
  
  
  #this is for every criterion
  stopifnot( cutoff > 0 | length(cutoff) != 1 ) 
  #special: alpha and the quantile can only be in [0,1]
  if ( criterion == "indirect" ) stopifnot( cutoff < 1 ) 
  
  if ( !is.im(distancemap) & !is.null(distancemap) ) stop("distancemap must be of class im (see distmap) or NULL.")
  
  if ( !is.im(intens) & !is.null(intens) ) stop("wrong input for intens. The intensity must be of class im (see density.ppp) or NULL.")
  
  if ( nxprob > 0 & nxprob < 1 ) stop("nxprob is a probability and therefore it has to be in the interval [0, 1]")
  
  if ( !is.null(covmatrix) && (!isSymmetric(covmatrix) | !all(eigen(covmatrix)$values > 0)) ){
    stop("covmatrix has to be symmetric and positive semidefinit")
  }
  
}





       