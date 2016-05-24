#' Evaluation of the high-risk zone.  
#'
#' Evaluation of the high-risk zone, which is only possible with simulated or thinned data or if the locations
#' of the unobserved events have been revealed..
#'
#' @param hrz  High-risk zone of class owin based on a binary mask (see \code{\link[spatstat]{area.owin}})
#' @param unobspp  Unobserved spatial point process 
#' @param obspp  Observed spatial point process
#' @importFrom methods is
#' @export  
#' @return An object of class "\code{hrzeval}", which is a list of
#'    \item{ numbermiss }{ number of unobserved events outside the high-risk zone }
#'    \item{ numberunobserved }{ number of events in the unobserved point pattern }
#'    \item{ missingfrac }{ fraction of unobserved events outside the high-risk zone (numbermiss/numberunobserved) }
#'    \item{ arearegion }{ area of the high-risk zone }
#'    \item{ numberobs }{ number of events in the observed point pattern }
#'    \item{ out }{ subset of the unobserved events, which are outside the high-risk zone }
#'    \item{ insd }{ subset of the unobserved events, which are inside the high-risk zone }
#' @seealso \code{\link[spatstat]{inside.owin}}, \code{\link[spatstat]{area.owin}} 
#' @examples
#'  data(craterB)
#'  # thin data
#'  set.seed(100)
#'  thdata <- thin(craterB, nxprob=0.1)
#'  
#'  # determine hrz for the "observed events"
#'  hrz <- det_hrz(thdata$observed, type = "dist", criterion = "area", cutoff = 1500000, nxprob = 0.1)
#'  
#'  # evaluate the hrz
#'  evaluation <- eval_hrz(hrz = hrz$zone, unobspp = thdata$unobserved, obspp = thdata$observed)
#'  evaluation$missingfrac
#'  
#'  op <- par(mar=c(1, 4, 1, 6) , xpd=TRUE)
#'  plot(evaluation, hrz = hrz, obspp = thdata$observed, plothrz = TRUE, plotobs = TRUE, 
#'  insidecol = "magenta", outsidecol = "magenta", obscol = "blue", insidepch = 1, 
#'  outsidepch = 19, main = "Evaluation visualized")
#'  legend(2400, 2456.4061, c("observed", "unobs inside", "unobs outside"), 
#'  col = c("blue", "magenta", "magenta"), yjust=1, pch=c(1, 1, 19), cex=0.8)
#'  par(op)


eval_hrz <- function(hrz, unobspp, obspp=NULL){
  
  #check if input arguments have correct values
  if ( !is.owin(hrz) & !is(hrz, "highriskzone") ) stop("hrz must be of class owin or of class highriskzone")
  if ( !is.ppp(unobspp) ) stop("unobspp must be of class ppp")
  if ( !is.ppp(obspp) & !is.null(obspp) ) stop("obspp must be of class ppp or NULL")
  
  #get the zone from the highrsikzone object
  if ( is(hrz, "highriskzone") ) hrz <- hrz$zone
  
  
  hrz$m[is.na(hrz$m)] <- FALSE
  
  #check if the unobserved events are inside the high-risk zone
  isin <- inside.owin(unobspp$x, unobspp$y, hrz)
  
  insd <- subset(unobspp, isin)
  out <- subset(unobspp, !isin)
  
  
  numbermiss <- sum(!isin)
  numberunobserved <- unobspp$n      
  missingfrac <- numbermiss/numberunobserved
  arearegion <- area.owin(hrz)
  
  numberobs <- ifelse(is.null(obspp) || all(is.na(obspp)), NA, obspp$n)
  
  result <- list(numbermiss=numbermiss, numberunobserved=numberunobserved, missingfrac=missingfrac,
                 arearegion=arearegion, numberobs=numberobs, out=out, insd=insd)
  class(result) <- "hrzeval"
  return(result)
}



