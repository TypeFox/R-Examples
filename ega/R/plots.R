#' @export
#' @title Plot a Clarke Error Grid
#' @description The function uses \code{\link[ggplot2]{ggplot}} to draw the
#' Clarke error grid lines according to the criteria described in the
#' original publication by Clarke et. al. (see reference below). If zones
#' have not already been assigned via the \code{zones} parameter, the
#' function \code{\link{getClarkeZones}} is called first. The values in
#' \code{referenceVals} and \code{testVals} are then superimposed as a scatter
#' plot. Some basic plot parameters can be specified as arguments, but the
#' return value can also be stored and modified further before plotting
#' (see examples and vignette).
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @param title The main plot title. Defaults to "Clarke Error Grid".
#' @param xlab The x-axis label. Defaults to "Reference Glucose
#' Concentration (mg/dL)".
#' @param ylab The y-axis label. Defaults to "Test Glucose Concentration
#' (mg/dL)".
#' @param linesize The size to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is 0.5.
#' @param linetype The type of line to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "solid".
#' @param linecolor The color of the zone lines. The acceptable values are the
#' same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "black".
#' @param linealpha The alpha (transparency) level to be used when drawing
#' the zone lines. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_segment}}. The default is 0.6.
#' @param pointsize The size to be used when plotting the glucose data points.
#' The acceptable values are the same as for \code{\link[ggplot2]{geom_point}}.
#' The default is 2.
#' @param pointalpha The alpha (transparency) level to be used when plotting
#' the glucose data points. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_point}}. The default is 1.
#' @param zones An optional character vector specifying the Clarke zones
#' for each paired value. If this is not supplied, \code{\link{getClarkeZones}}
#' will be called to generate zone labels.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the return
#' value is not assigned, a plot is drawn.
#' @examples
#' library(ggplot2)
#'
#' # default
#' plotClarkeGrid(glucose_data$ref, glucose_data$test)
#'
#' # with options
#' plotClarkeGrid(glucose_data$ref, glucose_data$test,
#'               pointsize=1.5,
#'               pointalpha=0.6,
#'               linetype="dashed")
#'
#' # store return value and modify
#' ceg <- plotClarkeGrid(glucose_data$ref, glucose_data$test)
#'
#' ceg + theme_gray() +
#'    theme(plot.title = element_text(size = rel(2), colour = "blue"))
#' @seealso
#' \code{\link{getClarkeZones}} \code{\link[ggplot2]{ggplot}}
#' @references
#' Clarke, W. L., D. Cox, L. A. Gonder-Frederick, W. Carter, and S. L. Pohl.
#' "Evaluating Clinical Accuracy of Systems for Self-Monitoring of Blood
#' Glucose." Diabetes Care 10, no. 5 (September 1, 1987): 622-28.
plotClarkeGrid <- function(referenceVals, testVals,
                           title="Clarke Error Grid",
                           xlab="Reference Glucose Concentration (mg/dL)",
                           ylab="Test Glucose Concentration (mg/dL)",
                           linesize=0.5,
                           linetype="solid",
                           linecolor="black",
                           linealpha=0.6,
                           pointsize=2,
                           pointalpha=1,
                           zones=NA) {

  # use default zone assignment if none is provided
  if (is.na(zones))
    zones <- getClarkeZones(referenceVals, testVals)

  tolerance <- 0.2

  # create a df for ggplot
  data <- data.frame("ref"=referenceVals, "test"=testVals, "zones"=zones)

  zoneA <- subset(data, zones=="A")
  maxA <- max(zoneA[["ref"]])

  # calculate line for zone A upper
  # A upper starts at (58 1/3, 70)
  x1 <- 58.3
  y1 <- 70
  x2 <- maxA
  y2 <- maxA+tolerance*maxA
  slope_Au <- (y2-y1)/(x2-x1) # rise over run!
  intercept_Au <- (-58.3*slope_Au)+70

  # choose large x,y end values so line goes off the chart
  xend_Au <- 5000
  yend_Au <- slope_Au*5000+intercept_Au # y=mx+b

  # calculate line for zone A lower
  x1 <- 70
  y1 <- 56
  x2 <- maxA
  y2 <- maxA-tolerance*maxA
  slope_Al <- (y2-y1)/(x2-x1)
  intercept_Al <- (-70*slope_Al)+56

  # choose large x,y end values so line goes off the chart
  xend_Al <- 5000
  yend_Al <- slope_Al*5000+intercept_Al

  # zone C upper
  x1 <- 70
  y1 <- 180
  slope_Cu <- 1
  intercept_Cu <- (-70*slope_Cu)+180

  # intersection of C upper and A upper
  # C upper: y = ax+c
  # A upper: y = bx+d

  # intersection: (d-c)/(a-b), (ad-bc)/(a-b)
  xend_Cu <- (intercept_Au-intercept_Cu)/(slope_Cu-slope_Au)
  yend_Cu <- ((slope_Cu*intercept_Au)-(slope_Au*intercept_Cu))/(slope_Cu-slope_Au)

  eu_coords <- list(c(0,180), c(70,180), c(70,Inf))
  el_coords <- list(c(180,0), c(180,70), c(Inf,70))
  dr_coords <- list(c(240,70), c(240,180), c(Inf,180))
  dl_coords <- list(c(0,70), c(58.3,70), c(70,84), c(70,180), c(0,180))
  au_coords <- list(c(58.3,70), c(xend_Au,yend_Au))
  al_coords <- list(c(70,0), c(70,56), c(xend_Al,yend_Al))
  cl_coords <- list(c(130,0), c(180,70))
  cu_coords <- list(c(70,180), c(xend_Cu,yend_Cu))

  labels_list <- list(c(240,230,"A"), c(120,185,"B"),
                      c(350,230,"B"), c(120,300,"C"),
                      c(163,20,"C"), c(35,130,"D"),
                      c(350,130,"D"), c(35,300,"E"),
                      c(350,35,"E"))

  # to appease CRAN
  ref <- test <- NULL

  ceg <- ggplot(data, aes(x=ref, y=test)) +

    # label the clinially relevant levels
    scale_x_continuous(breaks=c(70, 100, 150, 180, 240, 300, 350, 400, 450,
                                         500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                       expand = c(0,0)) +

    scale_y_continuous(breaks=c(70, 100, 150, 180, 250, 300, 350, 400, 450,
                                         500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                       expand=c(0, 0)) +

    # color by location type
    geom_point(aes(color=zones), size=pointsize, alpha=pointalpha) +

    # draw zone lines
    annotate_lines(eu_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(el_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(dr_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(dl_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(au_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(al_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(cl_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(cu_coords, linesize, linetype, linecolor, linealpha) +

    # now add the zone text labels
    annotate_labels(labels_list) +

    # dummy values to expand right and top
    annotate("text", x = 0, y = 350, size=6, label = "") +
    annotate("text", x = 400, y = 0, size=6, label = "") +

    theme_bw() +
    theme(legend.position="none") +
    ggtitle(title) +
    xlab(xlab) +
    ylab(ylab)

  ceg

}


#' @export
#' @title Plot a Parkes (Consensus) Error Grid
#' @description The function uses \code{\link[ggplot2]{ggplot}} to draw the
#' Parkes (consensus) error grid lines according to the criteria described in
#' the publications listed in the References section (see below). If zones
#' have not already been assigned via the \code{zones} parameter, the
#' function \code{\link{getParkesZones}} is called first. The values in
#' \code{referenceVals} and \code{testVals} are then superimposed as a scatter
#' plot. Some basic plot parameters can be specified as arguments, but the
#' return value can also be stored and modified further before plotting
#' (see examples and vignette).
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @param type An integer (1 or 2) specifying whether to plot the grid for Type 1
#' or Type 2 diabetes. Defaults to 1.
#' @param title The main plot title. Defaults to "Parkes (Consensus) Error Grid
#' for Type [type] Diabetes".
#' @param xlab The x-axis label. Defaults to "Reference Glucose
#' Concentration (mg/dL)".
#' @param ylab The y-axis label. Defaults to "Test Glucose Concentration
#' (mg/dL)".
#' @param linesize The size to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is 0.5.
#' @param linetype The type of line to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "solid".
#' @param linecolor The color of the zone lines. The acceptable values are the
#' same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "black".
#' @param linealpha The alpha (transparency) level to be used when drawing
#' the zone lines. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_segment}}. The default is 0.6.
#' @param pointsize The size to be used when plotting the glucose data points.
#' The acceptable values are the same as for \code{\link[ggplot2]{geom_point}}.
#' The default is 2.
#' @param pointalpha The alpha (transparency) level to be used when plotting
#' the glucose data points. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_point}}. The default is 1.
#' @param zones An optional character vector specifying the Clarke zones
#' for each paired value. If this is not supplied, \code{\link{getClarkeZones}}
#' will be called to generate zone labels.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the return
#' value is not assigned, a plot is drawn.
#' @examples
#' library(ggplot2)
#'
#' # default
#' plotParkesGrid(glucose_data$ref, glucose_data$test)
#'
#' # with options
#' plotParkesGrid(glucose_data$ref, glucose_data$test,
#'               pointsize=2,
#'               pointalpha=0.5,
#'               linesize=2,
#'               linealpha=0.3,
#'               linetype="dotdash")
#'
#' # store return value and modify
#' peg <- plotParkesGrid(glucose_data$ref, glucose_data$test, type=2)
#'
#' peg + theme_gray() +
#'    theme(plot.title = element_text(size = rel(2), colour = "red"))
#' @seealso
#' \code{\link{getParkesZones}} \code{\link[ggplot2]{ggplot}}
#' @references
#' Parkes, J. L., S. L. Slatin, S. Pardo, and B.H. Ginsberg. "A New Consensus
#' Error Grid to Evaluate the Clinical Significance of Inaccuracies in the
#' Measurement of Blood Glucose." Diabetes Care 23, no. 8 (August 2000):
#' 1143-48
#'
#' Pfutzner, Andreas, David C. Klonoff, Scott Pardo, and Joan L. Parkes.
#' "Technical Aspects of the Parkes Error Grid." Journal of Diabetes Science
#' and Technology 7, no. 5 (September 2013): 1275-81
plotParkesGrid <- function(referenceVals, testVals, type=1,
                           title=NA,
                           xlab="Reference Glucose Concentration (mg/dL)",
                           ylab="Test Glucose Concetration (mg/dL)",
                           linesize=0.5,
                           linetype="solid",
                           linecolor="black",
                           linealpha=0.6,
                           pointsize=2,
                           pointalpha=1,
                           zones=NA) {

  if (type != 1 & type != 2)
    stop("'type' must be 1 or 2.")

  if (is.na(title))
    title <- paste("Parkes (Consensus) Error Grid for Type", type,
                   "Diabetes")

  # use default zone assignment if none is provided
  if (is.na(zones))
    zones <- getParkesZones(referenceVals, testVals, type)

  # create a df for ggplot
  data <- data.frame("ref"=referenceVals, "test"=testVals, "zones"=zones)

  # zone B upper, Type 1
  # 0/50->30/50->140/170->280/380->5000/5729.3
  bu_coords <- list(c(0,50), c(30, 50), c(140, 170),
                    c(280, 380), c(5000, 5729.3))

  # zone B upper, Type 2
  # 0/50->30/50->230/330->440/550
  if (type == 2)
    bu_coords <- list(c(0,50), c(30,50), c(230,330), c(5000,5327.14))

  # zone B lower, Type 1
  # 50/0->50/30->170/145->385/300->5000/4495.45
  bl_coords <- list(c(50,0), c(50,30), c(170,145),
                    c(385,300), c(5000, 4495.45))

  # zone B lower, Type 2
  # 50/0->50/30->90/80->330/230->550/450
  if (type == 2)
    bl_coords <- list(c(50,0), c(50,30), c(90,80), c(330,230), c(5000,4900))

  # zone C upper, Type 1
  # 0/60->30/60->50/80->70/110->5000/11526.84
  cu_coords <- list(c(0,60), c(30,60), c(50,80),
                    c(70,110), c(5000, 11526.84))

  # zone C upper, Type 2
  # 0/60->30/60->280/550
  if (type == 2)
    cu_coords <- list(c(0,60), c(30,60), c(5000,9801.2))

  # zone C lower, Type 1
  # 120/0->120/30->260/130->5000/2091.38
  cl_coords <- list(c(120,0), c(120,30), c(260,130),
                    c(5000, 2091.38))

  # zone C lower, Type 2
  # 90/0->260/130->550/250
  if (type == 2)
    cl_coords <- list(c(90,0), c(260,130), c(5000,2091.38))

  # zone D upper, Type 1
  # 0/100->25/100->50/125->80/215->5000/36841.67
  du_coords <- list(c(0,100), c(25,100), c(50,125),
                    c(80,215), c(5000,36841.67))

  # zone D upper, Type 2
  # 0/80->25/80->35/90->125/550
  if (type == 2)
    du_coords <- list(c(0,80), c(25,80), c(35,90), c(5000,25466.67))

  # zone D lower, Type 1
  # 250/0->250/40->5000/1781.67
  dl_coords <- list(c(250,0), c(250,40), c(5000,1781.67))

  # zone D lower, Type 2
  # 250/0->250/40->410/110->550/160
  if (type == 2)
    dl_coords <- list(c(250,0), c(250,40), c(410,110), c(5000,1749.28))

  # zone E upper, Type 1
  # 0/150->35/155->5000/130900
  eu_coords <- list(c(0,150), c(35,155), c(5000,130900))

  # zone E upper, Type 2
  # 0/200->35/200->50/550
  if (type == 2)
    eu_coords <- list(c(0,200), c(35,200), c(5000,116050))

  label_list <- list(c(320,320,"A"), c(220,360,"B"),
                     c(385,235,"B"), c(140,375,"C"),
                     c(405,145,"C"), c(415,50,"D"),
                     c(75,383,"D"), c(21,383,"E"))

  # to appease CRAN
  ref <- test <- NULL

  peg <- ggplot(data, aes(x=ref, y=test)) +

    # label the clinially relevant levels
    scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400,
                                         450, 500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                                expand = c(0,0)) +

    scale_y_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400,
                                         450, 500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                                expand=c(0, 0)) +

    # color by zone
    geom_point(aes(color=zones), size=pointsize, alpha=pointalpha) +

    # draw zone lines
    annotate_lines(bu_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(bl_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(cu_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(cl_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(du_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(dl_coords, linesize, linetype, linecolor, linealpha) +
    annotate_lines(eu_coords, linesize, linetype, linecolor, linealpha) +

    # zone text labels
    annotate_labels(label_list, size=6) +

    # dummy values to force minimum size
    annotate("text", x = 0, y = 550, size=6, label = "") +
    annotate("text", x = 550, y = 0, size=6, label = "") +

    theme_bw() +
    theme(legend.position="none") +
    ggtitle(title) +
    xlab(xlab) +
    ylab(ylab)

  peg

}


annotate_lines <- function(xy_list, size=0.5, type="solid",
                           color="black", alpha=0.6) {

  lines <- length(xy_list)-1

  segs <- vector(mode="list", length=lines)

  for (i in 1:lines) {

    xy1 <- xy_list[[i]]
    xy2 <- xy_list[[i+1]]

    segs[[i]] <- annotate("segment", x=xy1[1], y=xy1[2], xend=xy2[1],
                        yend=xy2[2], size=size, linetype=type,
                        color=color, alpha=alpha)

  }

  segs

}


annotate_labels <- function(label_list, size=6) {

  lbls <- vector(mode="list", length=length(label_list))

  for (i in 1:length(lbls)) {

    lbl <- label_list[[i]]
    x <- as.numeric(lbl[1])
    y <- as.numeric(lbl[2])
    txt <- lbl[3]

    lbls[[i]] <- annotate("text", x=x, y=y, size=size, label=txt)

  }

  lbls

}
