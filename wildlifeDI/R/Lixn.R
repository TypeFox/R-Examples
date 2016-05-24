# ---- roxygen documentation ----
#
#' @title Minta's Spatial-temporal interaction statistics
#'
#' @description
#' The function \code{Lixn} measures dynamic interaction between two animals following
#' the methods outlined by Minta (1992).
#'
#' @details
#' The function \code{Lixn} can be used to calculate the Minta (1992) measures of dynamic
#' interaction between two animals. The Minta statistic tests how the two animals simultaneiously utilize
#' an area shared between the two individuals. Three coefficients are produced \eqn{L_{AA}}, \eqn{L_{BB}}, 
#' and \eqn{L_{ixn}}. Each of these statistics are based on a contingency table that compares the observed 
#' frequency of those fixes that are simultaneous and within/outside the shared area to expectations based on 
#' area overlap proportions (if \code{method="spatial"}) or expectations derived from all fixes (if 
#' \code{method="frequency"}) -- see Minta (1992) for more details. A Chi-squared statistic can then
#' be used to examine the significance between the observed and expected use of the shared area.
#' \cr\cr
#' Minta (1992) suggests the following interpretations of the coefficients. When \eqn{L_{AA}}
#' is near 0, the first animal's use of the shared area is random (or as expected). When
#' \eqn{L_{AA} > 0} it signifies spatial attraction to the shared area, or greater than
#' expected use. When \eqn{L_{AA} < 0} it signifies spatial avoidance of the shared area, or
#' less than expected use. Interpretation of \eqn{L_{BB}} is the same as for \eqn{L_{AA}} with
#' respect to the second animal. \eqn{L_{ixn}} tells us far more about the nature of the
#' interaction between the two individuals. As \eqn{L_{ixn}} nears 0, both animals use the
#' shared area randomly, with regards to the other animal. If \eqn{L_{ixn} > 0} the animals
#' use the shared area more \emph{simultaneously}, whereas if \eqn{L_{ixn} < 0} it is an
#' indication of \emph{solitary} use, or avoidance. This is why \eqn{L_{ixn}} is termed the temporal
#' interaction coefficient. A Chi-squared test can be used to identify the significance
#' of the \eqn{L_{AA}}, \eqn{L_{BB}}, and \eqn{L_{ixn}} values.
#' \cr\cr
#' NOTEs: 
#' \cr
#' 1. With modern telemetry datasets, where home ranges are readily estimated, choosing \code{method = 'spatial'}
#' is most appropriate. 
#' \cr
#' 2. When the home ranges do not overlap the Lixn statistic is not defined and the function returns a 
#' string of NA's.
#' \cr
#' 3. When one home range completely encloses another the Lixn statistic is not defined and the function returns
#' a string of NA's and \code{'ContainsB'} (or \code{'ContainsB'}) under the p.IXN result.
#' \cr
#' 4. Further to points 2 and 3, the Lixn statistic is not appropriate in situations where the overlap area is 
#' either very large or very small relative to either home range (i.e., a situation with almost complete enclosure 
#' or virtually no overlap). Thus, it is advised that \code{Lixn} be used only in situations where there are 
#' suitable marginal areas for areaA, areaB, and areaAB -- see Minta (1992).
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param method method for computing the marginal distribution from which expected
#'    values are computed. If \code{method = "spatial"}, the marginal values are calculated based on areas
#'    of the shared and unshared portions of the home ranges. If \code{method = "frequency"}, the marginal 
#'    values are calculated based on the number of all fixes within the shared and unshared portions of 
#'    the home ranges -- see Details.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param hr1 (-- required if method = 'spatial') home range polygon associated with \code{traj1}. Must
#'    be an object that coerces to class \code{SpatialPolygons*}.
#' @param hr2 (-- required if method = 'spatial') same as \code{hr1}, but for \code{traj2}.
#' @param OZ (-- required if method = 'frequency') shared area polygon associated with spatial use overlap
#'    between \code{traj1} and \code{traj2}. Must be an object that coerces to class \code{SpatialPolygons*}.
#'
#' @return
#' This function returns a list of objects representing the calculated values from the
#' Minta statistic and associated \emph{p}-values from the Chi-squared test.
#' \itemize{
#' \item pTable -- contingency table showing marginal probabilities of expected use 
#'    based on the selecton of the \code{method} parameter.
#' \item nTable -- contingency table showing observed frequency of use of the shared 
#' area based on simultaneous fixes.
#' \item oTable -- the odds for each cell in the contingency table.
#' \item Laa -- the calculated value of the \eqn{L_{AA}} statistic
#' \item p.AA -- the associated \emph{p}-value
#' \item Lbb -- the calculated value of the \eqn{L_{BB}} statistic
#' \item p.BB -- the associated \emph{p}-value
#' \item Lixn -- the calculated value of the \eqn{L_{ixn}} statistic
#' \item p.IXN -- the associated \emph{p}-value
#' }
#'
#' @references
#' Minta, S.C. (1992) Tests of spatial and temporal interaction among animals.
#' \emph{Ecological Applications}, \bold{2}: 178-188
#'
#' @keywords indices
#' @seealso GetSimultaneous
#' @examples
#' \dontrun{
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' library(adehabitatHR)
#' library(sp)
#' #use minimum convex polygon for demonstration...
#' hr37 <- mcp(SpatialPoints(ld(deer37)[,1:2]))
#' hr38 <- mcp(SpatialPoints(ld(deer38)[,1:2]))
#' #tc = 7.5 minutes, dc = 50 meters
#' Lixn(deer37, deer38,  method='spatial', tc=7.5*60, hr1=hr37, hr2=hr38)
#' }
#' @export
#
# ---- End of roxygen documentation ----
Lixn <- function(traj1,traj2,method="spatial",tc=0,hr1,hr2,OZ=NULL){
  output <- NULL
  #Get simultaneous fixes
  trajs <- GetSimultaneous(traj1,traj2,tc)
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  n <- nrow(tr1)
  pts1 <- SpatialPoints(tr1[,1:2])
  pts2 <- SpatialPoints(tr2[,1:2])
  
  #---- method specific computations ----
  if (method == 'spatial'){

    #could check that hr1 and hr2 exist and are spatial polygons here.
    if (gContains(hr1,hr2) == TRUE){output <- list(pTable=NA,nTable=NA,oTable=NA,Laa=NA,p.AA=NA,Lbb=NA,p.BB=NA,Lixn=NA,p.IXN="ContainsA")}
    else if (gContains(hr2,hr1) == TRUE){output <- list(pTable=NA,nTable=NA,oTable=NA,Laa=NA,p.AA=NA,Lbb=NA,p.BB=NA,Lixn=NA,p.IXN="ContainsB")}
    else {
      areaA <- gDifference(hr1,hr2)
      areaB <- gDifference(hr2,hr1)
      areaAB <- gIntersection(hr1,hr2)
      
      if (is.null(areaAB) == FALSE){
        
        #get intersected polygon for simultaneous fixes
        A1.int <- gCovers(areaA,pts1,byid=T)
        AB1.int <- gCovers(areaAB,pts1,byid=T)
        B2.int <- gCovers(areaB,pts2,byid=T)
        AB2.int <- gCovers(areaAB,pts2,byid=T)
        
        #check that the intersection vectors are same length before computing marginal values
        l.vec <- c(length(A1.int),length(B2.int),length(AB1.int),length(AB2.int))
        if (diff(range(l.vec)) == 0){
          #compute marginal values for simultaneous fixes
          n11 <- sum(AB1.int*AB2.int)
          n22 <- sum(A1.int*B2.int)
          n12 <- sum(A1.int*AB2.int)
          n21 <- sum(AB1.int*B2.int)
        } else {n11 <- 0; n12 <- 0; n21 <- 0; n22 <- 0}
        
        #compute expected values -- spatial
        a <- gArea(hr1)
        b <- gArea(hr2)
        ab <- gArea(areaAB)
        p11 <- (ab^2)/(a*b)
        p12 <- (1 - (ab/a))*(ab/b)
        p21 <- (ab/a)*(1 - (ab/b))
        p22 <- (1-(ab/a))*(1-(ab/b))
 
      } else {output <- list(pTable=NA,nTable=NA,oTable=NA,Laa=NA,p.AA=NA,Lbb=NA,p.BB=NA,Lixn=NA,p.IXN=NA)}
    
    }

  } else if (method == 'frequency'){

    areaAB <- OZ
    areaA <- gDifference(gBoundary(pts1),areaAB)
    areaB <- gDifference(gBoundary(pts2),areaAB)

    if (is.null(areaAB) == FALSE){
      
      #get intersected polygon for simultaneous fixes
      A1.int <- gCovers(areaA,pts1,byid=T)
      AB1.int <- gCovers(areaAB,pts1,byid=T)
      B2.int <- gCovers(areaB,pts2,byid=T)
      AB2.int <- gCovers(areaAB,pts2,byid=T)
    
      #check that the intersection vectors are same length before computing marginal values
      l.vec <- c(length(A1.int),length(B2.int),length(AB1.int),length(AB2.int))
      if (diff(range(l.vec)) == 0){
        #compute marginal values for simultaneous fixes
        n11 <- sum(AB1.int*AB2.int)
        n22 <- sum(A1.int*B2.int)
        n12 <- sum(A1.int*AB2.int)
        n21 <- sum(AB1.int*B2.int)
      } else {n11 <- 0; n12 <- 0; n21 <- 0; n22 <- 0}
    

      #Compute expected values -- frequency
      t1 <- ld(traj1)
      t2 <- ld(traj2)
      r <- nrow(t1)
      s <- nrow(t2)
      p1 <- SpatialPoints(t1[,1:2])
      p2 <- SpatialPoints(t2[,1:2])
      #get pts that intersect areaAB for ALL fixes
      AB1.int <- gCovers(areaAB,p1,byid=T)
      AB2.int <- gCovers(areaAB,p2,byid=T)
      rAB <- length(which(AB1.int == T))
      sAB <- length(which(AB2.int == T))
      p11 <- (rAB*sAB)/(r*s)
      p12 <- (1 - (rAB/r))*(sAB/s)
      p21 <- (rAB/r)*(1-(sAB/s))
      p22 <- (1 - (rAB/r))*(1 - (sAB/s))
    
    } else {output <- list(pTable=NA,nTable=NA,oTable=NA,Laa=NA,p.AA=NA,Lbb=NA,p.BB=NA,Lixn=NA,p.IXN=NA)}  

  } else {stop(paste("The method - ",method,", is not recognized. Please try again.",sep=""))}

  #---- end of method specific computations ---
  
  #Compute chi-square results if appropriate
  if (is.null(output)){
    #compute summary statistics 
    w <- (n11*n22)/(n12*n21)
    L <- log(w)
    se.L <- sqrt((1/n11) + (1/n12) + (1/n21) + (1/n22))
    #Chi-square of marginal frequency values
    n.1 <- n11 + n21
    n.2 <- n12 + n22
    n1. <- n11 + n12
    n2. <- n21 + n22
    
    # Intrinsic Hypothesis
    #------------ NOT RETURNED -----------------------
    #chiINT <- (n*((n11*n22 - n12*n21)^2))/(as.numeric(n1.)*as.numeric(n2.)*as.numeric(n.1)*as.numeric(n.2))
    #phiINT <- sqrt(chiINT/n)
    #--------------------------------------------------------------------
    
    #compute summary and chi-values for Extrinsic Hypothesis.
    p.1 <- p11+p21
    p.2 <- p12+p22
    p1. <- p11+p12
    p2. <- p21+p22
    #see email from Eric Howe
    #Not returned, not sure value of Chi.tot
    chi.tot <- (((n11-p11*n)^2)/p11*n) + (((n12-p12*n)^2)/p12*n) + (((n21-p21*n)^2)/(p21*n)) + (((n22-p22*n)^2)/(p22*n))
    Laa <- log((p.2*n.1)/(p.1*n.2))
    chiAA <- ((n.1-p.1*n)^2)/(p.2*p.1*n)
    Lbb <- log((p2.*n1.)/(p1.*n2.))
    chiBB <- ((n1.-p1.*n)^2)/(p2.*p1.*n)
    Lixn <- log(((n11/p11)+(n22/p22))/((n12/p12)+(n21/p21)))
    chiIXN <- (((n11/p11) + (n22/p22) - (n12/p12) - (n21/p21))^2)/(n*((1/p11) + (1/p12) + (1/p21) + (1/p22)))
    #odds for each cell
    o11 <- n11/(p11*n)
    o12 <- n12/(p12*n)
    o21 <- n21/(p21*n)
    o22 <- n22/(p22*n)
    
    #compute chi-square p-values for easier interpretation.
    #Intrinsic Hypothesis - Not Returned
    #p.INT <- 1 - pchisq(chiINT,df=1)
    #p.tot <- 1 - pchisq(chi.tot,df=3)  #something off, see email from Eric Howe, how to interpret?
    
    #Extrinsic Hypothesis
    p.AA <- 1 - pchisq(chiAA,df=1)
    p.BB <- 1 - pchisq(chiBB,df=1)
    p.IXN <- 1 - pchisq(chiIXN,df=1)      #Fixed Typo Here...
    #create an output data-frame
    pTable <- matrix(c(p11,p12,p21,p22),ncol=2,byrow=T,dimnames=list(c("B","b"),c("A","b")))
    nTable <- matrix(c(n11,n12,n21,n22),ncol=2,byrow=T,dimnames=list(c("B","b"),c("A","b")))
    oTable <- matrix(c(o11,o12,o21,o22),ncol=2,byrow=T,dimnames=list(c("B","b"),c("A","b")))
    
    #Not Returning IntHyp as not easily interpreted and ExtHyp is better
    #IntHyp <- list(n=n,L=L,se.L=se.L,p.INT=p.INT,phi.INT=phiINT)
    #ExtHyp <- list(Laa=Laa,p.AA=p.AA,Lbb=Lbb,p.BB=p.BB,Lixn=Lixn,p.IXN=p.IXN)
    output <- list(pTable=pTable,nTable=nTable,oTable=oTable,Laa=Laa,p.AA=p.AA,Lbb=Lbb,p.BB=p.BB,Lixn=Lixn,p.IXN=p.IXN) 
  }

  return(output)
}
#====================End of Minta Function =====================================
