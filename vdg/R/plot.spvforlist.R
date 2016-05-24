#' @rdname plot.spv
#' @method plot spvforlist
#' @export
plot.spvforlist <- function (x, which = c("fds", "vdgsim", "vdgquantile", "vdgboth", "boxplots"), 
                             np = 50, alpha = 7/sqrt(length(x[[1]]$spv)), points.colour = "#39BEB1",
                             points.size = 2, tau = c(0.05, 0.95), radii = 21, hexbin = FALSE, bins = 80,
                             df = 10, lines.size = 1, origin = rep(0, ncol(x[[1]]$sample)), method, 
                             arrange = FALSE, ...) {
  
  # Avoid global variable notes for R CMD check and ggplot2
  Radius <- SPV <- Formula <- Fraction <- Location <- NULL
  
  # Handle which depending on whether it is numeric or character (gets transformed to numeric)
  pnms <- c("fds", "vdgsim", "vdgquantile", "vdgboth", "boxplots")
  show <- rep(FALSE, 5)
  if (is.character(which)) {
    which <- match.arg(which, several.ok = TRUE)
    which <- sort(match(which, pnms))
  } 
  if (!is.numeric(which)) stop("Argument 'which' is of incorrect type.")
  show[which] <- TRUE
  type <- x[[1]]$type
  
  if (x[[1]]$at && show[1L]){
    show[1L] <- FALSE
    which <- which[!(which %in% 1L)]
    message("Plot 1 = 'fds' cannot be produced: 'at' is TRUE (inaccurate FDS plot)")
  }
  
  add.meanspv <- x[[1]]$type == "spherical" && !is.null(radii)
  
  if (is.null(tau) & !add.meanspv){
    if (any(3L:4L %in% which)) 
      message("Plots 3 = 'vdgboth' and/or 4 cannot be produced: 'tau' is NULL and mean SPV not requested/possible")
    show[3L:4L] <- FALSE
    which <- which[!(which %in% 3L:4L)]
  }
  
  if (!x[[1]]$at && show[5L]){
    show[5L] <- FALSE
    which <- which[!(which %in% 5L)]
    message("Plot 5 = 'boxplots' cannot be produced: 'at' is FALSE")
  }
  
  pnms <- pnms[show]
  
  if (missing(method)) method <- switch(type, spherical = "Euclidean", cuboidal = "supremum", 
                                       lhs = "supremum", mlhs = "supremum", slhs = "supremum",  rslhs = "supremum")
  xvec <- proxy::dist(x[[1]]$sample, matrix(origin, nrow = 1, ncol = ncol(x[[1]]$sample)), method = method, ...)
  method <- attr(xvec, "method")
  xvec <- as.numeric(xvec)
  ndes <- length(x)
  nspv <- length(xvec)
  ntau <- length(tau)
  spvmat <- sapply(x, '[[', "spv")
  nms <- names(x)
  if (is.null(nms)) nms <- seq_along(x)
  
  if (add.meanspv){
    if (length(radii) == 1) radii <- seq(from = 0, to = max(xvec), length.out = radii)
    mspv <- lapply(x, function(y) as.data.frame(meanspv(formula = y$formula, radii = radii, FtF.inv = y$FtF.inv, n = ifelse(y$unscaled, 1, y$ndes))))
    tmp3 <- do.call(rbind, mspv)
    tmp3$Formula <- rep(names(x), each = length(radii))
    tmp3$Location <- "Mean"
  }
  
  if (any(show[-1L])){
    tmp1 <- data.frame(Radius = rep(xvec, ndes), SPV = as.numeric(spvmat), 
                       Formula = factor(rep(nms, each = nspv)))
  } 
  
  if (show[1L]){
    maxmin <- range(sapply(x, '[[', "spv")) 
    pts <- 0:np/np
    quantmat <- sapply(x, function(xx) quantile(xx$spv, probs = pts, type = 1))
    tmp2 <- data.frame(Fraction = rep(pts, ndes), SPV = as.numeric(quantmat),
                       Formula = factor(rep(colnames(quantmat), each = length(pts))))
  }
  
  if (show[3L] || show[4L]){
    if (!is.null(tau)){
      pts <- seq(from = min(tmp1$Radius), to = max(tmp1$Radius), length = np)
      aggfun <- function(x){
        fits <- lapply(tau, function(y) quantreg::rq(SPV ~ bs(Radius, df = df), tau = y, data = x))
        sapply(fits, predict, newdata = data.frame(Radius = pts))
      }
      preds <- sapply(split(tmp1, tmp1$Formula), aggfun)
      newdf <- data.frame(Radius = rep(pts, ntau*ndes), SPV = as.numeric(preds), 
                           Location = rep(rep(paste("tau =", tau), each = np), ndes),
                           Formula = factor(rep(colnames(preds), each = ntau*np)))
      if (exists("tmp3", inherits = FALSE)) tmp3 <- rbind(tmp3, newdf)
      else tmp3 <- newdf
    }
  }
  
  if (exists("tmp3", inherits = FALSE)) tmp3$Location <- as.factor(tmp3$Location)
  
  if (show[1L]) {
    plot1 <- ggplot(tmp2, aes(x = Fraction, y = SPV, colour = Formula)) + ggtitle("Fraction of Design Space Plot") + 
      xlab("Fraction of Design Space") + geom_line(size = lines.size, ...) +
      theme(plot.title = element_text(vjust = 1))
  }
  
  if (show[2L]) {
    plot2 <- ggplot(tmp1, aes(x = Radius, y = SPV, order = Formula)) + 
        ggtitle("Variance Dispersion Graph") +
        geom_point(alpha = alpha, colour = points.colour, size = points.size) + 
        theme(plot.title = element_text(vjust = 1)) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1))) + xlab(paste0("Distance to Origin (", method,")")) +
        facet_wrap(~ Formula)
    if (hexbin){
      plot2 <- plot2 + geom_hex(bins = bins) + 
        scale_fill_gradientn(colours = rev(topo.colors(5)[-(4:5)]), name = "Frequency", na.value = NA)
    }
  }
  
  if (show[3L]) {
      plot3 <- ggplot(tmp3, aes(x = Radius, y = SPV)) + ggtitle("Variance Dispersion Graph") +  
        theme(plot.title = element_text(vjust = 1)) +  
        geom_line(aes(group = interaction(Location, Formula), linetype = Formula, colour = Location), 
                  size = lines.size)  + theme(legend.text.align = 0.5) + xlab(paste0("Distance to Origin (", method,")"))
  }
  
  if (show[4L]) {
    plot4 <- ggplot(tmp1, aes(x = Radius, y = SPV, order = Formula)) + 
        ggtitle("Variance Dispersion Graph") +
        geom_point(alpha = alpha, colour = points.colour, size = points.size) + 
        theme(plot.title = element_text(vjust = 1)) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1))) + xlab(paste0("Distance to Origin (", method,")"))
    if (hexbin){
      plot4 <- plot4 + geom_hex(bins = bins) + 
        scale_fill_gradientn(colours = rev(topo.colors(5)[-(4:5)]), name = "Frequency", na.value = NA)
    }
    plot4 <- plot4 + geom_line(data = tmp3, aes(x = Radius, y = SPV, linetype = Location, 
                                                order = Formula), size = lines.size) + facet_wrap(~ Formula)
  }
  
  if (show[5L]) {
    if (!exists("tmp1", inherits = FALSE)){
      tmp1 <- data.frame(Radius = rep(xvec, ndes), SPV = as.numeric(spvmat), 
                         Formula = factor(rep(nms, each = nspv)))
    }
    plot5 <- ggplot(tmp1, aes(x = Radius, y = SPV)) + ggtitle("Boxplots") +
      theme(plot.title = element_text(vjust = 1), 
            legend.text.align = 0.5) +
      geom_boxplot(aes(x = as.factor(round(Radius, getOption("digits"))))) + xlab(paste0("Distance to Origin (", method,")")) +
      facet_wrap(~ Formula)
  }
  
  if (length(which)) {
    out <- mget(paste0("plot", which))
    names(out) <- pnms
    if (x[[1]]$unscaled) out <- lapply(out, '+', ylab("Unscaled Prediction Variance (UPV)"))
    else out <- lapply(out, '+', ylab("Scaled Prediction Variance (SPV)"))
    if (arrange) do.call(gridExtra::grid.arrange, out)
    else return(out)
  }
}