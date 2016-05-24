# Evaluierung der Prozeduren

#' Evaluation of the procedures determining the high-risk zone. 
#'
#' Evaluates the performance of the three methods:
#' \itemize{
#'    \item Method of fixed radius
#'    \item Quantile-based method
#'    \item Intensity-based method
#' }
#' For further details on the methods, see \code{\link{det_hrz}} or the paper of Mahling et al. (2013)(References). \cr
#' There are three ways to simulate data for the evaluation.
#' 
#' The three simulation types are:
#' \describe{
#' \item{ Data-based simulation }{
#'        Here a given data set is used. The data set is thinned as explained below.
#'        Note that this method is very different from the others, since it is using the real data.
#'        }
#' \item{ Simulation of an inhomogeneous Poisson process }{
#'        Here, an inhomogeneous Poisson process is simulated and then that data is thinned.
#'        }
#' \item{ Simulation of a Neyman-Scott process }{
#'        Here a Neyman-Scott process is simulated (see \code{\link{sim_nsppp}}, \code{\link[spatstat]{rNeymanScott}})
#'        and this data is then also thinned.
#'        }
#' }
#' Thinning:\cr
#' Let \eqn{X} be the spatial point process, which is the location of all events and let \eqn{Y}
#' be a subset of \eqn{X} describing the observed process. The process of unobserved events
#' then is \ifelse{latex}{\eqn{Z = X \setminus Y}}{Z = X \ Y} , meaning that \eqn{Z} and \eqn{Y} are disjoint and together 
#' forming \eqn{X}.\cr
#' Since \eqn{Z} is not known, in this function an observed or simulated spatial point pattern 
#' \code{ppdata} is taken as the full pattern (which we denote by \ifelse{latex}{\eqn{\tilde X}}{X'}) comprising the 
#' observed events \ifelse{latex}{\eqn{\tilde Y}}{Y'} as well as the unobserved \ifelse{latex}{\eqn{\tilde Z}}{Z'}.
#' Each event in \ifelse{latex}{\eqn{\tilde X}}{X'} is assigned to one of the two processes \ifelse{latex}{\eqn{\tilde Y}}{Y'} or 
#' \ifelse{latex}{\eqn{\tilde Z}}{Z'} by drawing independent Bernoulli random numbers. \cr
#' The resulting process of observed events \ifelse{latex}{\eqn{\tilde Y}}{Y'} is used to determine the high-risk zone. 
#' Knowing now the unobserved process, it can be seen how many events are outside and inside the
#' high-risk zone.\cr
#' 
#' \code{type} and \code{criterion} may be vectors in this function.
#'
#' @inheritParams det_hrz
#' @param numit Number of iterations
#' @param simulate The type of simulation, can be one of \code{"thinning", "intens"} or \code{"clintens"}
#' @param radiusClust (Optional) radius of the circles around the parent points in which the cluster
#'                    points are located. Only used for \code{simulate = "clintens"}. 
#' @param clustering  a value >= 1 which describes the amount of clustering; the
#'          adjusted estimated intensity of the observed pattern is divided by
#'          this value; it is also the parameter of the Poisson distribution
#'          for the number of points per cluster. Only used for \code{simulate = "clintens"}.
#' @param pbar  logical. Should progress bar be printed?
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export  
#' @return A \code{data.frame} with variables 
#'    \item{ Iteration }{ Iterationstep of the result }
#'    \item{ Type, Criterion, Cutoff, nxprob }{ see arguments }
#'    \item{ threshold }{ determined threshold. If criterion="area", it is either the distance (if type="dist")
#' or the threshold c (for type="intens"). If criterion="indirect", it is either the quantile of the
#' nearest-neighbour distance which is used as radius (if type="dist") or the threshold c (for type="intens"). If criterion="direct",
#' it equals the cutoff for both types.}
#'    \item{ calccutoff }{ determined cutoff-value. For type="dist" and criterion="area", this is the
#' quantile of the nearest-neighbour distance. For type="intens" and criterion="area", it is the failure
#' probability alpha. For all other criterions it is NA.}
#'    \item{ covmatrix11, covmatrix12, covmatrix21, covmatrix22 }{ values in the covariance matrix. 
#' covmatrix11 and covmatrix22 are the diagonal elements (variances). }
#'    \item{ numbermiss }{ number of unobserved points outside the high-risk zone }
#'    \item{ numberunobserved }{ number of observations in the unobserved point pattern \ifelse{latex}{\eqn{\tilde Z}}{Z'} }
#'    \item{ missingfrac }{ fraction of unobserved events outside the high-risk zone (numbermiss/numberunobserved) }
#'    \item{ arearegion }{ area of the high-risk zone }
#'    \item{ numberobs }{ number of observations in the observed point pattern \ifelse{latex}{\eqn{\tilde Y}}{Y'} }
#'    
#' @seealso \code{\link{det_hrz}}, \code{\link[spatstat]{rNeymanScott}}, \code{\link{thin}}, \code{\link{sim_nsppp}}, \code{\link{sim_intens}}
#' @examples    
#'  data(craterB)
#'  
#'  # the input values are mainly the same as in det_hrz, so for more example ideas, 
#'  # see the documentation of det_hrz.
#'  evalm <- eval_method(craterB, type = c("dist", "intens"), criterion = c("area", "area"), 
#'                       cutoff = c(1500000, 1500000), nxprob = 0.1, numit = 10, 
#'                       simulate = "clintens", radiusClust = 300, 
#'                       clustering = 15, pbar = FALSE)
#'  evalm_d <- subset(evalm, evalm$Type == "dist")
#'  evalm_i <- subset(evalm, evalm$Type == "intens")
#'  
#'  # pout:  fraction of high-risk zones that leave at least one unobserved event uncovered
#'  # pmiss:  Mean fraction of unobserved events outside the high-risk zone
#'  data.frame(pmiss_d = mean(evalm_d$missingfrac),
#'             pmiss_i = mean(evalm_i$missingfrac),
#'             pout_d = ( sum(evalm_d$numbermiss > 0) / nrow(evalm_d) ), 
#'             pout_i = ( sum(evalm_i$numbermiss > 0) / nrow(evalm_i) ))


eval_method <- function(ppdata, type, criterion, cutoff, numit = 100, 
                        nxprob = 0.1, distancemap = NULL, intens = NULL, 
                        covmatrix = NULL, simulate, radiusClust=NULL, clustering=5,
                        pbar = TRUE) {
  
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
  
  # for the output later
  if ( simulate == "clintens" ) {
    result <- matrix(data=NA, nrow=0, ncol=18)
    #colnames(result) <- c("Iteration", "Type", "Criterion", "Cutoff", "nxprob", "threshold", 
    #                      "calccutoff", "covmatrix11", "covmatrix12", "covmatrix21", "covmatrix22", 
    #                      "numbermiss", "numberunobserved", "missingfrac", "arearegion", "numberobserved", 
    #                      "clusterrad", "clustering")
  } else {
    result <- matrix( data=NA, nrow=0, ncol=16 )
    #colnames(result) <- c("Iteration", "Type", "Criterion", "Cutoff", "nxprob", "threshold", 
    #                      "calccutoff", "covmatrix11", "covmatrix12", "covmatrix21", "covmatrix22", 
    #                      "numbermiss", "numberunobserved", "missingfrac", "arearegion", "numberobserved")    
  }
  
  
  nproc <- length(type)
  
  # Progress Bar to see how far the function is
  if ( pbar == TRUE ) pb <- txtProgressBar(min = 0, max = numit, style = 3)
  for ( i in 1:numit ) {
    
    ## simulations:
    if ( simulate == "thinning" ) {
      thinned <- thin(full=ppdata, nxprob=nxprob)   #thinned statt thinnedbombs
      observed <- thinned$observed                 #observed statt obsbombs
      unobserved <- thinned$unobserved             #unobserved statt uxo
    }
    if ( simulate == "intens" ) {
      thinned <- sim_intens(ppdata, intensSim, nxprob)
      observed <- thinned$observed
      unobserved <- thinned$unobserved
    }
    if ( simulate == "clintens" ) {
      ppsim <- sim_nsprocess(ppdata=ppdata, intens=intensSim, radius=radiusClust, 
                             clustering=clustering, thinning=nxprob)
      thinned <- thin(full=ppsim, nxprob=nxprob)   #thinned statt thinnedbombs
      observed <- thinned$observed                 #observed statt obsbombs
      unobserved <- thinned$unobserved             #unobserved statt uxo
    }
    
    
    
    if (is.null(distancemap) && "dist" %in% type){
      distancemap <- distmap(observed)
    }
    
    if (is.null(intens) && "intens" %in% type){
      estim <- est_intens(observed, covmatrix=covmatrix)
      intensest <- estim$intensest
      covmatrix <- estim$covmatrix
    }
    
    for (j in 1:nproc){
      
      typej <- type[j]
      criterionj <- criterion[j]
      cutoffj <- cutoff[j]
      
      resultdetHRZ <- det_hrz(ppdata=observed, type=typej, criterion=criterionj, 
                              cutoff=cutoffj, distancemap=distancemap, intens=intensest, 
                              nxprob=nxprob, covmatrix=covmatrix)
      resultevalHRZ <- eval_hrz(hrz=resultdetHRZ$zone, unobspp=unobserved, obspp=observed)
      
      covmatrix11 <- ifelse(all(is.null(resultdetHRZ$covmatrix)), NA, resultdetHRZ$covmatrix[1, 1])
      covmatrix12 <- ifelse(all(is.null(resultdetHRZ$covmatrix)), NA, resultdetHRZ$covmatrix[1, 2])
      covmatrix21 <- ifelse(all(is.null(resultdetHRZ$covmatrix)), NA, resultdetHRZ$covmatrix[2, 1])
      covmatrix22 <- ifelse(all(is.null(resultdetHRZ$covmatrix)), NA, resultdetHRZ$covmatrix[2, 2])
      
      if ( simulate == "clintens" ) {
        resstep <- c(i, typej, criterionj, cutoffj, nxprob, resultdetHRZ$threshold,
                     resultdetHRZ$calccutoff, covmatrix11, covmatrix12,
                     covmatrix21, covmatrix22, resultevalHRZ$numbermiss,
                     resultevalHRZ$numberunobserved, resultevalHRZ$missingfrac, 
                     resultevalHRZ$arearegion, resultevalHRZ$numberobs,
                     radiusClust, clustering)
      } else {
        resstep <- c(i, typej, criterionj, cutoffj, nxprob, resultdetHRZ$threshold,
                     resultdetHRZ$calccutoff, covmatrix11, covmatrix12,
                     covmatrix21, covmatrix22, resultevalHRZ$numbermiss,
                     resultevalHRZ$numberunobserved, resultevalHRZ$missingfrac, 
                     resultevalHRZ$arearegion, resultevalHRZ$numberobs)
      }
      
      result <- rbind(result, resstep)
      
    }
    if ( pbar == TRUE ) setTxtProgressBar(pb, i)
  }
  #resultdf <- as.data.frame(result, row.names = as.character(1:(nproc*numit)))
  resultdf <- data.frame(as.ordered(result[ , 1]), 
                         apply(result[ , 2:3], 2, as.character), 
                         apply(result[ , 4:NCOL(result)], 2, as.numeric), 
                         row.names = as.character(1:(nproc*numit)))
  
  
  if ( simulate == "clintens" ) {
    colnames(resultdf) <- c("Iteration", "Type", "Criterion", "Cutoff", "nxprob", "threshold", 
                         "calccutoff", "covmatrix11", "covmatrix12", "covmatrix21", "covmatrix22", 
                         "numbermiss", "numberunobserved", "missingfrac", "arearegion", "numberobserved", 
                         "clusterrad", "clustering")
  } else {
    colnames(resultdf) <- c("Iteration", "Type", "Criterion", "Cutoff", "nxprob", "threshold", 
                         "calccutoff", "covmatrix11", "covmatrix12", "covmatrix21", "covmatrix22", 
                         "numbermiss", "numberunobserved", "missingfrac", "arearegion", "numberobserved")    
  }
  
  resultdf
  
}



