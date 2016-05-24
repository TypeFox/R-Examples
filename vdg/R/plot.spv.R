#' Plot VDGs or FDS plots
#' 
#' Produce Variance Dispersion Graphs and/or Fraction of Design Space plots for
#' experimental designs. There are methods for the S3 classes \code{spv}, 
#' \code{spvlist}, \code{spvforlist} and \code{spvlistforlist} -- see 
#' \code{\link{spv}}.
#' 
#' @aliases plot.spv plot.spvlist plot.spvforlist plot.spvlistforlist
#' @param x an object of type \code{spv} for a single experimental design or an
#' object of type \code{spvlist} for multiple designs.
#' @param which either a numeric vector of integers or a character vector
#' indicating which plots to produce. The possible plots are: 
#' \describe{
#' \item{\code{1} or \code{"fds"}}{A (variance ratio) FDS plot}
#' \item{\code{2} or \code{"vdgsim"}}{A VDG with only the simulated prediction variance points plotted}
#' \item{\code{3} or \code{"vdgquantile"}}{A VDG with only the quantile regression lines corresponding to \code{tau} shown}
#' \item{\code{4} or \code{"vdgboth"}}{A combination of \code{2} and \code{3}}
#' \item{\code{5} or \code{"boxplots"}}{Parallel boxplots of the prediction variance}
#' }
#' @param np scalar; the number of points to use for calculating the fraction of design space criterion.
#' @param alpha the alpha transparency coefficient for the plots
#' @param points.colour colour for points in scatterplot of SPV against the radius
#' @param points.size size for points in scatterplot of SPV against the radius
#' @param tau the tau parameter for \code{\link[quantreg]{rq}} (quantile
#' regression)
#' @param radii either a numeric vector containing the radii to use for
#' calculating the mean spherical SPV over the spherical design space, or an integer
#' (length one vector) giving the number of radii to use for calculationg
#' the mean spherical SPV. If missing, the mean spherical SPV is not used.
#' @param hexbin logical indicating whether hexagonal binning should be used to display 
#' density instead of colour transparency
#' @param bins argument passed to \code{\link{stat_binhex}} to determine the 
#' number of hexagons used for binning.
#' @param VRFDS logical indicating whether to construct a variance ratio FDS plot or not (only for class \code{spvlist}). The
#' first design is used as reference design in case of \code{VRFDS} is \code{TRUE}
#' @param df degrees-of-freedom parameter passed to \code{\link{bs}}
#' @param lines.size line size passed to \code{\link[ggplot2]{geom_line}}
#' @param origin numeric vector specifying the origin of the design space
#' @param method optional; passed to \code{\link[proxy]{dist}} to overwrite
#' defaults of "Euclidean" for spherical regions or "supremum" for cubiodal
#' regions
#' @param arrange Logical indicating whether to return a single graphical object arranging the 
#' resulting plots in a single plot window via \code{\link{grid.arrange}}, or whether to return the
#' list of graphical objects containing the plots.
#' @param \dots additional arguments passed to \code{\link[proxy]{dist}}
#' @return Returns a list of \code{\link{ggplot}} graphical objects (or grobs) with names corresponding
#' to the character version of \code{which}. These plot objects can be manipulated by changing plot 
#' aesthetics and theme elements.
#' @keywords hplot
#' @author Pieter C. Schoonees
#' @method plot spv
#' @export
#' @import ggplot2
#' @import quantreg
#' @importFrom proxy dist
#' @importFrom splines bs
#' @importFrom gridExtra grid.arrange
#' @examples
#' 
#' # Single design (class 'spv')
#' # Larger n should be used in actual cases
#' library(rsm)
#' bbd3 <- as.data.frame(bbd(3)[,3:5])
#' colnames(bbd3) <- paste0("x", 1:3)
#' quad.3f <- formula(~ x1*x2*x3 - x1:x2:x3 + I(x1^2) + I(x2^2) + I(x3^2))
#' set.seed(1234)
#' out <- spv(n = 1000, design = bbd3, type = "spherical", formula = quad.3f)
#' out
#' plot(out)
#' 
#' # List of designs (class 'spvlist')
#' \dontrun{
#' library(Vdgraph)
#' data(SCDH5); data(SCDDL5)
#' des.list <- list(SCDH5 = SCDH5, SCDDL5 = SCDDL5)
#' quad.5f <- formula(~ x1 + x2 + x3 + x4 + x5 + x1:x2 + x1:x3 + x1:x4 + x1:x5
#'                     + x2:x3 + x2:x4 + x2:x5 + x3:x4 + x3:x5 + x4:x5 
#'                    + I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2))
#' out2 <- spv(n = 500, design = des.list, type = "spherical", formula = quad.5f)
#' out2
#' plot(out2)
#' }
#' 
#' # List of formulae (class 'spvforlist')
#' \dontrun{
#' fact3 <- expand.grid(x1 = c(-1,1), x2 = c(-1, 1), x3 = c(-1,1))
#' lin.3f <- formula(~ x1 + x2 + x3)
#' int.3f <- formula(~ (x1+x2+x3)^2)
#' set.seed(4312)
#' out3 <- spv(n = 500, design = fact3, type = "cuboidal", 
#'              formula = list(linear = lin.3f, interaction = int.3f))
#' out3
#' plot(out3)
#' }
#' 
#' # List of formulae and designs (class 'spvlistforlist')
#' \dontrun{
#' fact3.n <- rbind(fact3, 0, 0, 0)
#' set.seed(4312)
#' out4 <- spv(n = 200, design = list(factorial = fact3, factorial.with.cntr = fact3.n), 
#'              type = "cuboidal", formula = list(linear = lin.3f, interaction = int.3f))
#' out4
#' plot(out4)
#' }
plot.spv <- function (x, which = c("fds", "vdgsim", "vdgquantile", "vdgboth", "boxplots"), 
                      np = 50, alpha = 7/sqrt(length(x$spv)), points.colour = "#39BEB1",
                      points.size = 2, tau = c(0.05, 0.95), radii = 21, hexbin = FALSE, bins = 80, 
                      df = 10, lines.size = 1, origin = rep(0, ncol(x$sample)), method, 
                      arrange = FALSE, ...) {
  # Avoid global variable notes for R CMD check and ggplot2
  Radius <- SPV <- Fraction <- Location <- NULL
  
  # Handle which depending on whether it is numeric or character (gets transformed to numeric)
  pnms <- c("fds", "vdgsim", "vdgquantile", "vdgboth", "boxplots")
  show <- rep(FALSE, 5)
  if (is.character(which)) {
    which <- match.arg(which, several.ok = TRUE)
    which <- sort(match(which, pnms))
  } 
  if (!is.numeric(which)) stop("Argument 'which' is of incorrect type.")
  show[which] <- TRUE
  type <- x$type
  
  if (x$at && show[1L]){
    show[1L] <- FALSE
    which <- which[!(which %in% 1L)]
    message("Plot 1 = 'fds' cannot be produced: 'at' is TRUE (inaccurate FDS plot)")
  }
  
  add.meanspv <- x$type == "spherical" && !is.null(radii)
  
  if (is.null(tau) & !add.meanspv){
    if(any(3L:4L %in% which)) 
      message("Plots 3 = 'vdgquantile' and/or 4 = 'vdgboth' cannot be produced: 
              'tau' is NULL and mean SPV not requested/possible")
    show[3L:4L] <- FALSE
    which <- which[!(which %in% 3L:4L)]
  }
  
  if (!x$at && show[5L]){
    show[5L] <- FALSE
    which <- which[!(which %in% 5L)]
    message("Plot 5 = 'boxplots' cannot be produced: 'at' is FALSE")
  }
  
  pnms <- pnms[show]
  
  if (missing(method)) method <- switch(type, spherical = "Euclidean", cuboidal = "supremum", 
                                       lhs = "supremum", mlhs = "supremum", slhs = "supremum",  
                                       rslhs = "supremum")
  xvec <- proxy::dist(x$sample, matrix(origin, nrow = 1, ncol = ncol(x$sample)), method = method, ...)
  method <- attr(xvec, "method")
  xvec <- as.numeric(xvec)
  
  if (add.meanspv){
      if(length(radii) == 1) radii <- seq(from = 0, to = max(xvec), length.out = radii)
      mspv <- meanspv(formula = x$formula, radii = radii, FtF.inv = x$FtF.inv, 
                      n = ifelse(x$unscaled, 1, x$ndes))
      tmp3 <- as.data.frame(mspv)
      tmp3$Location <- "Mean"
  }
  
  if (any(show[-1L])) tmp1 <- data.frame(Radius = xvec, SPV = x$spv)
  
  if (show[1L]){
    maxmin <- range(x$spv)  
    pts <- 0:np/np
    tmp2 <- data.frame(Fraction = pts, SPV = quantile(x$spv, probs = pts, type = 1))
  }
  
  if (show[3L] || show[4L]){
    if (!is.null(tau)){
      pts <- seq(from = min(tmp1$Radius), to = max(tmp1$Radius), length = np)
      fits <- lapply(tau, function(x, data) quantreg::rq(SPV ~ bs(Radius, df = df), tau = x, 
                                                         data = data), data = tmp1)
      newdf <- data.frame(Radius = rep(pts, length(tau)), SPV = as.numeric(
        sapply(fits, predict, newdata = data.frame(Radius = pts))), 
        Location = rep(paste("tau =", tau), each = np))
      if (exists("tmp3", inherits = FALSE)) tmp3 <- rbind(tmp3, newdf)
      else tmp3 <- newdf
    }
  }

  if (exists("tmp3", inherits = FALSE)) tmp3$Location <- as.factor(tmp3$Location)
  
  if (show[1L]){
    plot1 <- ggplot(tmp2, aes(x = Fraction, y = SPV)) + ggtitle("Fraction of Design Space Plot") + 
      xlab("Fraction of Design Space") + 
      geom_line(size = lines.size, colour = points.colour) +
      theme(plot.title = element_text(vjust = 1))
  }
  
  if (show[2L]) {
    plot2 <- ggplot(tmp1, aes(x = Radius, y = SPV)) + ggtitle("Variance Dispersion Graph") +
      geom_point(alpha = alpha, colour = points.colour, size = points.size) + 
      theme(plot.title = element_text(vjust = 1)) + 
      xlab(paste0("Distance to Origin (", method,")"))
    if (hexbin) {
      plot2 <- plot2 + geom_hex(bins = bins) + 
        scale_fill_gradientn(colours = rev(topo.colors(5)[-(4:5)]), name = "Frequency", na.value = NA)
    }
  }
  
  if (show[3L]) {
    plot3 <- ggplot(tmp1, aes(x = Radius, y = SPV)) + ggtitle("Variance Dispersion Graph") +
      theme(plot.title = element_text(vjust = 1), legend.text.align = 0.5) 
    if(hexbin){
      plot3 <- plot3 + geom_hex(bins = bins) + 
        scale_fill_gradientn(colours = rev(topo.colors(5)[-(4:5)]), name = "Frequency", na.value = NA)
    } 
    plot3 <- plot3 + geom_line(mapping = aes(x = Radius, y = SPV, linetype = Location), data = tmp3, 
                               size = lines.size) + xlab(paste0("Distance to Origin (", method,")"))
  }
  
  if (show[4L]) {
    plot4 <- ggplot(tmp1, aes(x = Radius, y = SPV)) + ggtitle("Variance Dispersion Graph") +
      geom_point(alpha = alpha, colour = points.colour, size = points.size) + 
      theme(plot.title = element_text(vjust = 1), legend.text.align = 0.5) 
    if (hexbin) plot4 <- plot4 + geom_hex(bins = bins) + 
      scale_fill_gradientn(colours = rev(topo.colors(5)[-(4:5)]), name = "Frequency", na.value = NA)
    plot4 <- plot4 + geom_line(mapping = aes(x = Radius, y = SPV, linetype = Location), data = tmp3, 
                size = lines.size) + xlab(paste0("Distance to Origin (", method,")"))
  }
  
  if (show[5L]) {
    plot5 <- ggplot(tmp1, aes(x = Radius, y = SPV)) + ggtitle("Boxplots") +
      geom_boxplot(aes(x = as.factor(round(Radius, getOption("digits"))))) + 
      xlab(paste0("Distance to Origin (", method,")")) +
      theme(plot.title = element_text(vjust = 1), legend.text.align = 0.5)
  }
  
if (length(which)) {
    out <- mget(paste0("plot", which))
    names(out) <- pnms
    if(x$unscaled) out <- lapply(out, '+', ylab("Unscaled Prediction Variance (UPV)"))
    else out <- lapply(out, '+', ylab("Scaled Prediction Variance (SPV)"))
    if (arrange) do.call(gridExtra::grid.arrange, out)
    else return(out)
  }
}
