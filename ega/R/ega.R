#' Clarke and Parkes (Consensus) error grid analysis
#'
#' @name ega-package
#' @import ggplot2
#' @docType package
NULL


#' 5072 paired reference and test glucose values.
#'
#' A dataset containing 5072 paired reference method and test method
#' glucose values (in mg/dL).
#'
#' @format A data frame with 5072 rows and 2 variables:
#' \describe{
#'   \item{ref}{Reference method glucose value, in mg/dL}
#'   \item{test}{Test method glucose value, in mg/dL}
#' }
#' @source The data is from a modified clinical dataset.
"glucose_data"


#' @export
#' @title Assign Clarke error grid zones to paired glucose values
#' @description \code{referenceVals} and \code{testVals} are assumed to contain
#' paired glucose values from a reference method and a test method,
#' respectively. The discrepancy between the two values is used to place the
#' pair into a Clarke error grid zone according to the criteria described
#' in the original paper by Clarke et. al. (see reference below).
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @return A character vector is returned, with each element being one of
#' \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, or \code{"E"}.
#' @examples
#' zones <- getClarkeZones(glucose_data$ref, glucose_data$test)
#'
#' # counts
#' table(zones)
#'
#' # percentages
#' table(zones)/length(zones)*100
#'
#' @references
#' Clarke, W. L., D. Cox, L. A. Gonder-Frederick, W. Carter, and S. L. Pohl.
#' "Evaluating Clinical Accuracy of Systems for Self-Monitoring of Blood
#' Glucose." Diabetes Care 10, no. 5 (September 1, 1987): 622-28.
getClarkeZones <- function(referenceVals, testVals) {

  zones <- vector(mode="character", length = length(referenceVals))

  bias <- testVals-referenceVals

  # absolute relative error = abs(bias)/reference*100
  are <- abs(bias)/referenceVals*100

  eq1 <- (7/5)*(referenceVals-130)
  eq2 <- referenceVals+110

  # zone D: ref < 70 and (test > 70 and test < 180) or
  #   ref > 240 and (test > 70 and test < 180)
  test_ok <- testVals > 70 & testVals < 180
  zoneD <- (referenceVals < 70 & test_ok) | (referenceVals > 240 & test_ok)

  zones[zoneD] <- "D"

  # assign A after D, since part of A will overwrite D

  # zone A: are <= 20  or (ref < 58.3 and test < 70)
  zoneA <- (are <= 20) | (referenceVals < 58.3 & testVals < 70)

  zones[zoneA] <- "A"

  # zone E: (ref <= 70 and test >= 180) or (ref >=180 and test <=70)
  zoneE <- (referenceVals <= 70 & testVals >= 180) |
    (referenceVals >= 180 & testVals <= 70)

  zones[zoneE] <- "E"



  # zone C: (ref >= 130 and ref <= 180 and test < eq1) or
  #   (ref > 70 and ref > 180 and ref > eq2)
  zoneC <- (referenceVals >= 130 & referenceVals <= 180 & testVals < eq1) |
    (referenceVals > 70 & testVals > 180 & testVals > eq2)

  zones[zoneC] <- "C"

  # the rest are zone B
  zones <- replace(zones, zones=="", "B")

  return(zones)

}


#' @export
#' @title Assign Parkes (Consensus) error grid zones to paired glucose values
#' @description \code{referenceVals} and \code{testVals} are assumed to contain
#' paired glucose values from a reference method and a test method,
#' respectively. The discrepancy between the two values, as well as the
#' type of error grid desired (Type 1 or Type 2 diabetes), is used to place the
#' pair into a Parkes (Consensus) error grid zone, according to the
#' criteria described in the second reference below.
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @param type An integer (1 or 2) specifying whether to obtain zones for Type 1
#' or Type 2 diabetes. Defaults to 1.
#' @return A character vector is returned, with each element being one of
#' \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, or \code{"E"}.
#' @examples
#' zones <- getParkesZones(glucose_data$ref, glucose_data$test)
#'
#' # counts
#' table(zones)
#'
#' # percentages
#' table(zones)/length(zones)*100
#' @references
#' Parkes, J. L., S. L. Slatin, S. Pardo, and B.H. Ginsberg. "A New Consensus
#' Error Grid to Evaluate the Clinical Significance of Inaccuracies in the
#' Measurement of Blood Glucose." Diabetes Care 23, no. 8 (August 2000):
#' 1143-48
#'
#' Pfutzner, Andreas, David C. Klonoff, Scott Pardo, and Joan L. Parkes.
#' "Technical Aspects of the Parkes Error Grid." Journal of Diabetes Science
#' and Technology 7, no. 5 (September 2013): 1275-81
getParkesZones <- function(referenceVals, testVals, type=1) {

  if (type != 1 & type != 2)
    stop("'type' must be 1 or 2.")

  zones <- vector(mode="character", length = length(referenceVals))

  # B lower lines, Type 1:
  # 50/0->50/30->170/145->385/300->5000/4495.45
  bl_xy <- list(c(50,0), c(50,30), c(170,145), c(385,300), c(5000,4495.45))

  # B lower lines, Type 2:
  # 50/0->50/30->90/80->330/230->550/450
  if (type == 2)
    bl_xy <- list(c(50,0), c(50,30), c(90,80), c(330,230), c(5000,4900))

  # B lower line equations
  bl_lines <- getLineEqs(bl_xy, referenceVals)

  bl1 <- bl_lines[[1]]
  bl2 <- bl_lines[[2]]
  bl3 <- bl_lines[[3]]
  bl4 <- bl_lines[[4]]

  # B upper lines, Type 1:
  # 0/50->30/50->140/170->280/380->5000/5729.3
  bu_xy <- list(c(0,50), c(30,50), c(140,170), c(280,380), c(5000,5729.3))

  # B upper lines, Type 2:
  # 0/50->30/50->230/330->440/550
  if (type == 2)
    bu_xy <- list(c(0,50), c(30,50), c(230,330), c(5000,5327.14))

  # B upper line equations
  bu_lines <- getLineEqs(bu_xy, referenceVals)

  bu1 <- bu_lines[[1]]
  bu2 <- bu_lines[[2]]
  bu3 <- bu_lines[[3]]

  if (type == 1)
    bu4 <- bu_lines[[4]]

  # C lower lines. Type 1:
  # 120/0->120/30->260/130->5000/2091.38
  cl_xy <- list(c(120,0), c(120,30), c(260,130), c(5000,2091.38))
  first_c <- 120

  # C lower lines, Type 2:
  # 90/0->260/130->550/250
  if (type == 2) {
    cl_xy <- list(c(90,0), c(260,130), c(5000,2091.38))
    first_c <- 90
  }

  # C lower line equations:
  cl_lines <- getLineEqs(cl_xy, referenceVals)

  cl1 <- cl_lines[[1]]
  cl2 <- cl_lines[[2]]

  if (type == 1)
    cl3 <- cl_lines[[3]]

  # C upper lines, Type 1:
  # 0/60->30/60->50/80->70/110->5000/11526.84
  cu_xy <- list(c(0,60), c(30,60), c(50,80), c(70,110), c(5000,11526.84))

  # C upper lines, Type 2:
  # 0/60->30/60->280/550
  if (type == 2)
    cu_xy <- list(c(0,60), c(30,60), c(5000,9801.2))

  # C upper line equations:
  cu_lines <- getLineEqs(cu_xy, referenceVals)

  cu1 <- cu_lines[[1]]
  cu2 <- cu_lines[[2]]

  if (type == 1) {
    cu3 <- cu_lines[[3]]
    cu4 <- cu_lines[[4]]
  }

  # D lower lines, Type 1:
  # 250/0->250/40->5000/1781.67
  dl_xy <- list(c(250,0), c(250,40), c(5000,1781.67))

  # D lower lines, Type 2:
  # 250/0->250/40->410/110->550/160
  if (type == 2)
    dl_xy <- list(c(250,0), c(250,40), c(410,110), c(5000,1749.28))

  # D lower line equations:
  dl_lines <- getLineEqs(dl_xy, referenceVals)

  dl1 <- dl_lines[[1]]
  dl2 <- dl_lines[[2]]

  if (type == 2)
    dl3 <- dl_lines[[3]]

  # D upper lines, Type 1:
  # 0/100->25/100->50/125->80/215->5000/36841.67
  du_xy <- list(c(0,100), c(25,100), c(50,125), c(80,215), c(5000,36841.67))
  first_d <- 100

  # D upper lines, Type 2:
  # 0/80->25/80->35/90->125/550
  if (type == 2) {
    du_xy <- list(c(0,80), c(25,80), c(35,90), c(5000,25466.67))
    first_d <- 80
  }

  # D upper line equations
  du_lines <- getLineEqs(du_xy, referenceVals)

  du1 <- du_lines[[1]]
  du2 <- du_lines[[2]]
  du3 <- du_lines[[3]]

  if (type == 1)
    du4 <- du_lines[[4]]

  # E lines, Type 1:
  # 0/150->35/155->5000/130900
  eu_xy <- list(c(0,150), c(35,155), c(5000,130900))

  # E lines, Type 2:
  # 0/200->35/200->50/550
  if (type == 2)
    eu_xy <- list(c(0,200), c(35,200), c(5000,116050))

  # E line equations
  eu_lines <- getLineEqs(eu_xy, referenceVals)

  eu1 <- eu_lines[[1]]
  eu2 <- eu_lines[[2]]

  # zone B lower
  zoneB_lower <- (between_y(bl_xy[1:2], testVals) & referenceVals >= 50) |
    (between_y(bl_xy[2:3], testVals) & testVals < bl2) |
    (between_y(bl_xy[3:4], testVals) & testVals < bl3) |
    (between_y(bl_xy[4:5], testVals) & testVals < bl4)

  # zone B upper
  zoneB_upper <- (between_x(bu_xy[1:2], referenceVals) & testVals > 50) |
    (between_y(bu_xy[2:3], testVals) & testVals > bu2) |
    (between_y(bu_xy[3:4], testVals) & testVals > bu3)

  if (type == 1) {
    zoneB_upper <- zoneB_upper |
      (between_y(bu_xy[4:5], testVals) & testVals > bu4)
  }

  # zone C lower
  if (type == 1) {

    zoneC_lower <- (between_y(cl_xy[1:2], testVals) & referenceVals >= first_c) |
      (between_y(cl_xy[2:3], testVals) & testVals < cl2) |
      (between_y(cl_xy[3:4], testVals) & testVals < cl3)

  }
  else {

    zoneC_lower <- (between_y(cl_xy[1:2], testVals) & testVals < cl1) |
      (between_y(cl_xy[2:3], testVals) & testVals < cl2)

  }

  # zone C upper
  zoneC_upper <- (between_x(cu_xy[1:2], referenceVals) & testVals > 60) |
    (between_y(cu_xy[2:3], testVals) & testVals > cu2)

  if (type == 1) {
    zoneC_upper <- zoneC_upper |
      (between_y(cu_xy[3:4], testVals) & testVals > cu3) |
      (between_y(cu_xy[4:5], testVals) & testVals > cu4)
  }

  # zone D lower
  zoneD_lower <- (between_y(dl_xy[1:2], testVals) & referenceVals >= 250) |
    (between_y(dl_xy[2:3], testVals) & testVals < dl2)

  if (type == 2) {
    zoneD_lower <- zoneD_lower |
      (between_y(dl_xy[3:4], testVals) & testVals < dl3)
  }

  # zone D upper
  zoneD_upper <- (between_x(du_xy[1:2], referenceVals) & testVals > first_d) |
    (between_y(du_xy[2:3], testVals) & testVals > du2) |
    (between_y(du_xy[3:4], testVals) & testVals > du3)

  if (type == 1) {
    zoneD_upper <- zoneD_upper |
      (between_y(du_xy[4:5], testVals) & testVals > du4)
  }

  # zone E
  zoneE_upper <- (between_x(eu_xy[1:2], referenceVals) & testVals > eu1) |
    (between_y(eu_xy[2:3], testVals) & testVals > eu2)

  zones[zoneB_lower] <- "B"
  zones[zoneC_lower] <- "C"
  zones[zoneD_lower] <- "D"

  zones[zoneB_upper] <- "B"
  zones[zoneC_upper] <- "C"
  zones[zoneD_upper] <- "D"

  zones[zoneE_upper] <- "E"

  # the rest are zone A
  zones <- replace(zones, zones=="", "A")

  zones

}


getLineEqs <- function(xylist, xvals) {

  num_lines <- length(xylist)-1
  line_eqs <- vector(mode="list", length=num_lines)

  for (i in 1:num_lines) {

    xy1 <- xylist[[i]]
    xy2 <- xylist[[i+1]]

    x1 <- xy1[1]
    y1 <- xy1[2]
    x2 <- xy2[1]
    y2 <- xy2[2]

    m <- (y2-y1)/(x2-x1)
    b <- (-x1*m)+y1

    line_eqs[[i]] <- m*xvals+b

  }

  line_eqs

}


between_y <- function(xylist, yvals) {

  xy1 <- xylist[[1]]
  xy2 <- xylist[[2]]

  yvals > xy1[2] & yvals <= xy2[2]

}

between_x <- function(xylist, yvals) {

  xy1 <- xylist[[1]]
  xy2 <- xylist[[2]]

  yvals > xy1[1] & yvals <= xy2[1]

}
