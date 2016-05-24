################################################################################
# Teilfunktionen Bestimmung HRZ                                                #
################################################################################

#' Determination of the area of a high-risk zone using the nearest-neighbour distance.
#'
#' Used in function det_radius.
#'
#' @param cutoffval  distance used as radius of the discs
#' @param distancemap  distance map (object of class "im", see \code{\link[spatstat]{distmap}}): 
#'                     distance of every location in the observation window to the nearest event
#' @param win  observation window of class owin 
#' @return A numerical value giving the area of the window.
#' @seealso \code{\link[spatstat]{eval.im}}, \code{\link[spatstat]{owin}}, \code{\link[spatstat]{area.owin}}


#- fr?her area_dist
det_nnarea <- function(cutoffval, distancemap, win){
  HRZimage <- eval.im(distancemap < cutoffval)
  safetyregion <- owin(xrange=win$xrange, yrange=win$yrange, mask=as.matrix(HRZimage))
  safetyregion$m[is.na(safetyregion$m)] <- FALSE
  arearegion <- area.owin(safetyregion)
  return(arearegion)
}


#' Determination of the nearest-neighbour distance which results in a high-risk zone with desired area
#'
#' Used in function det_hrz.
#'
#' @param ppdata  observed spatial point pattern of class ppp.
#' @param distancemap  distance map (object of class "im", see \code{\link[spatstat]{distmap}}): 
#'                     distance of every location in the observation window to the nearest event
#' @param areahrz  given area of the high-risk zone
#' @param win  observation window of class owin 
#' @importFrom stats uniroot
#' @return A list of 
#'     \item{ cutoffdist }{ quantile of the nearest-neighbour distance }
#'     \item{ thresh }{ distance }
#' @seealso \code{\link{det_nnarea}}, \code{\link[stats]{quantile}}, \code{\link[stats]{uniroot}}

#- fr?her dist_area 
det_radius <- function(ppdata, distancemap, areahrz, win){
  
  f <- function(cutoffval){
    det_nnarea(cutoffval, distancemap=distancemap, win=win) - areahrz
  }
  
  ### neu
  mindist <- min(distancemap$v, na.rm=TRUE)
  maxdist <- max(distancemap$v, na.rm=TRUE)
  ### neu
  if(sign(f(mindist)) != sign(f(maxdist))){
    thres <- uniroot(f, lower=mindist, upper=maxdist)
    thresh <- thres$root
  }else{
    thresh <- ifelse(f(mindist) > 0, mindist, maxdist)
  }
  
  quantval <- function(quant){
    quantile(nndist(ppdata), p=quant, type=8) - thresh
  }
  
  ### neu
  if(sign(quantval(0)) != sign(quantval(1))){
    cutoffval <- uniroot(quantval, lower=0, upper=1)
    cutoffdist <- cutoffval$root
  }else{     
    cutoffdist <- NA
  }
  
  result <- list(cutoffdist=cutoffdist, thresh=thresh) 
  return(result)
  
}



#################################################################################
# Direktere Methode zur Bestimmung der Wahrscheinlichkeit fuer mindestens einen
# uebersehenen Blindgaenger, bei der die Poissonverteilung der Anzahl der
# Blindgaenger in einem bestimmten Gebiet verwendet wird


#' calculation of alpha (failure probability), when having the threshold c
#'
#' This function is used for the intensity-based method. It determines the
#' probability to have at least one unobserved event outside the high-risk zone. 
#' A Poisson distribution is used for the number of unobserved events in a certain area or field.
#' Used in functions det_threshold, det_thresholdfromarea.
#'
#' @param intens  estimated intensity of the observed process (object of class "im", see \code{\link[spatstat]{density.ppp}})
#' @param threshold  threshold c: The high-risk zone is the field in which the estimated intensity 
#'                   exceeds this value.
#' @param nxprob  probability of having unobserved events   
#' @return value of alpha


#- fr?her prob_pois
#* (pMiss in alpha umbenennen)
det_alpha <- function(intens, threshold, nxprob=0.1) {
  
  #Area of one pixel in the image
  areaPixel <- intens$xstep * intens$ystep
  
  intensmatrix <- as.matrix(intens)
  intensmatrix[is.na(intensmatrix)] <- 0
  
  # Determine safety region based on threshold - consider those
  # points which are located outside the safety region
  WwR <- (intensmatrix <= threshold)
  
  # Integral of the intensity function over area inside window W, but outside
  # the safety region
  intens.intWwR <- sum(intensmatrix[WwR==1] * areaPixel)
  
  # We know N_Z(W\R) \sim \Po(q / (1 - q) * \int_{W\R} \lambda(t) dt)
  alpha <- 1 - exp(-(nxprob / (1 - nxprob) * intens.intWwR))
  
  return(alpha)
}



#' Calculation of the threshold c, when having failure probability alpha.
#'
#' The high-risk zone is the field in which the estimated intensity 
#' exceeds the threshold c, which is determined here, having the failure probability
#' alpha.
#' This function is for the intensity-based method. Used in function det_hrz.
#'
#' @param intens  estimated intensity of the observed process (object of class "im", see \code{\link[spatstat]{density.ppp}})
#' @param alpha   failure probability: probability to have at least one unobserved event
#'                outside the high-risk zone
#' @param nxprob  probability of having unobserved events
#' @importFrom stats uniroot
#' @return value of the threshold c
#' @seealso \code{\link{det_alpha}}, \code{\link[stats]{uniroot}}


#- fr?her: intens_indirect
det_threshold <- function(intens, alpha=1e-5, nxprob=0.1){
  
  f <- function(logthreshold) {
    det_alpha(intens=intens, threshold=exp(logthreshold), nxprob=nxprob) - alpha
  }
  
  thres <- uniroot(f, lower=-100, upper=log(max(range(intens))))
  #thres <- uniroot(f, lower=-100, upper=max(0,log(max(range(intens.hat)))))   obere Zeile funktioniert evtl. bei Markierung nicht, wenn kein Punkt mehr in Teilprozess
  threshold <- exp(thres$root)
  
  return(threshold)
}




#' Calculation of the area of the high-risk zone. 
#'
#' This function is used for the intensity-based method. Calculation of the
#' area of the high-risk zone given the observation window,
#' the intensity matrix and the threshold c. Used in function 
#' det_thresholdfromarea.
#'
#' @param win  observation window
#' @param intensmatrix  matrix of the estimated intensity of the observed process (\code{as.matrix(intens)})
#' @param threshold  threshold c: The high-risk zone is the field in which the estimated intensity 
#'                   exceeds this value
#' @return A numerical value giving the area of the high-risk zone.
#' @seealso \code{\link[spatstat]{owin}}, \code{\link[spatstat]{area.owin}}


#- fr?her intens_areahrz 
det_area <- function(win, intensmatrix, threshold){
  
  R <- intensmatrix > threshold
  
  safetyregion <- owin(xrange=win$xrange, yrange=win$yrange, mask=R)
  
  safetyregion$m[is.na(safetyregion$m)] <- FALSE
  arearegion <- area.owin(safetyregion)
  return(arearegion)
}



#' Determination of alpha and the threshold c which results
#' in a high-risk zone with desired area. 
#'
#' This function is used for the intensity-based method. Used in function det_hrz.
#'
#' @param intens  estimated intensity of the observed process (object of class "im", see \code{\link[spatstat]{density.ppp}}) 
#' @param areahrz  area of the high-risk zone 
#' @param win  observation window
#' @param nxprob  probability of having unbserved events  
#' @importFrom stats uniroot
#' @return A list of
#'   \item{ threshold }{ Value of the threshold c. The high-risk zone is the field in which the estimated intensity 
#'                   exceeds this value }
#'   \item{ calccutoff }{ failure probability alpha for given area; probability to have at least unobserved event outside the high-risk zone }
#' @seealso \code{\link{det_area}}, \code{\link{det_alpha}}

#- fr?her: intens_area  
det_thresholdfromarea <- function(intens, areahrz, win, nxprob=0.1){
  
  intensmatrix <- as.matrix(intens)
  
  f <- function(logthreshold){
    det_area(win=win, intensmatrix=intensmatrix, exp(logthreshold)) - areahrz
  }
  
  ### neu
  minint <- -100000
  maxint <- log(max(range(intens)))
  ### neu
  if(sign(f(minint)) != sign(f(maxint))){
    thres <- uniroot(f, lower=minint, upper=maxint)
    thresh <- exp(thres$root)
    
    fal <- function(alpha) {
      det_alpha(intens, threshold=thresh, nxprob=nxprob) - alpha
    }
    
    al <- uniroot(fal, lower=-1, upper=2)
    alpha <- al$root
    
  }else{
    thresh <- ifelse(f(minint) < 0, exp(minint), exp(maxint)) 
    alpha <- NA
  }
  
  result <- list(threshold=thresh, calccutoff=alpha)
  return(result)
}
