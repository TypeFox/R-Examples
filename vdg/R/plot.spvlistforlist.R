#' @rdname plot.spv
#' @method plot spvlistforlist
#' @export
plot.spvlistforlist <- function (x, which =c("fds", "vdgsim", "vdgquantile", "vdgboth", "boxplots"), 
                                 np = 50, alpha = 7/sqrt(length(x[[1]][[1]]$spv)), points.colour = "#39BEB1",
                                 points.size = 2, tau = c(0.05, 0.95),  radii = 21, hexbin = FALSE, bins = 80, 
                                 df = 10, lines.size = 1, origin = rep(0, ncol(x[[1]][[1]]$sample)), method, 
                                 arrange = FALSE, ...) {
  
  # Avoid global variable notes for R CMD check and ggplot2
  Radius <- SPV <- Design <- Formula <- Fraction <- Location <- NULL
  
  # Handle which depending on whether it is numeric or character (gets transformed to numeric)
  pnms <- c("fds", "vdgsim", "vdgquantile", "vdgboth", "boxplots")
  show <- rep(FALSE, 5)
  if (is.character(which)) {
    which <- match.arg(which, several.ok = TRUE)
    which <- sort(match(which, pnms))
  } 
  if (!is.numeric(which)) stop("Argument 'which' is of incorrect type.")
  show[which] <- TRUE
  type <- x[[1]][[1]]$type
  
  if (x[[1]][[1]]$at && show[1L]){
    show[1L] <- FALSE
    which <- which[!(which %in% 1L)]
    message("Plot 1 = 'fds' cannot be produced: 'at' is TRUE (inaccurate FDS plot)")
  }
  
  add.meanspv <- x[[1]][[1]]$type == "spherical" && !is.null(radii)
  
  if (is.null(tau) & !add.meanspv){
    if (any(3L:4L %in% which)) 
      message("Plots 3 = 'vdgquantile' and/or 4 = 'vdgboth' cannot be produced: 'tau' is NULL and mean SPV not requested/possible")
    show[3L:4L] <- FALSE
    which <- which[!(which %in% 3:4)]
  }
  
  if (!x[[1]][[1]]$at && show[5L]){
    show[5L] <- FALSE
    which <- which[!(which %in% 5L)]
    message("Plot 5 = 'boxplots' cannot be produced: 'at' is FALSE")
  }
  
  pnms <- pnms[show]
  
  if (missing(method)) method <- switch(type, spherical = "Euclidean", cuboidal = "supremum", 
                                       lhs = "supremum", mlhs = "supremum", slhs = "supremum",  rslhs = "supremum")
  xvec <- proxy::dist(x[[1]][[1]]$sample, matrix(origin, nrow = 1, ncol = ncol(x[[1]][[1]]$sample)), method = method, ...)
  method <- attr(xvec, "method")
  xvec <- as.numeric(xvec)
  nfor <- length(x)
  ndes <- length(x[[1]])
  nspv <- length(xvec)
  ntau <- length(tau)
  spvmat <- do.call(cbind, lapply(x, function(y) do.call(cbind, lapply(y, "[[", "spv"))))
  fornms <- names(x)
  desnms <- names(x[[1]])
  names(spvmat) <- paste(rep(desnms, nfor), rep(fornms, each = ndes), sep = ".")
  
  if (add.meanspv){
    if (length(radii) == 1) radii <- seq(from = 0, to = max(xvec), length.out = radii)
    mspv <- lapply(x, function(y) lapply(y, function(z) as.data.frame(meanspv(formula = z$formula, radii = radii, 
                                                                              FtF.inv = z$FtF.inv, n = ifelse(z$unscaled, 1, z$ndes)))))
    tmp3 <- do.call(rbind, lapply(mspv, function(x) do.call(rbind, x)))
    tmp3$Formula <- rep(fornms, each = ndes * length(radii))
    tmp3$Design <- rep(rep(desnms, each = length(radii)), nfor)
    tmp3$Location <- "Mean"
  }
  
  if (any(show[-1L])){
    tmp1 <- data.frame(Radius = rep(xvec, ndes*nfor), SPV = as.numeric(spvmat), 
                       Design = factor(rep(rep(desnms, each = nspv), nfor)), Formula = factor(rep(fornms, each = ndes*nspv)))
  } 
  
  if (show[1L]){
    maxmin <- range(spvmat) 
    pts <- 0:np/np
    quantmat <- apply(spvmat, 2, function(xx) quantile(xx, probs = pts, type = 1))
    tmp2 <- data.frame(Fraction = rep(pts, ndes*nfor), SPV = as.numeric(quantmat),
                       Formula = factor(rep(fornms, each = ndes*(np + 1))), 
                       Design = rep(rep(desnms, each = np + 1), nfor))
  }
  
  if (show[3L] || show[4L]){
    if (!is.null(tau)){
      if (!exists("tmp1", inherits = FALSE)){
        tmp1 <- data.frame(Radius = rep(xvec, ndes*nfor), SPV = as.numeric(spvmat), 
                           Design = factor(rep(rep(desnms, each = nspv), nfor)), Formula = factor(rep(fornms, each = nfor*nspv)))
      }
      pts <- seq(from = min(tmp1$Radius), to = max(tmp1$Radius), length = np)
      aggfun <- function(x){
        fits <- lapply(tau, function(y) quantreg::rq(SPV ~ bs(Radius, df = df), tau = y, data = x))
        sapply(fits, predict, newdata = data.frame(Radius = pts))
      }
      preds <- sapply(split(tmp1, list(tmp1$Formula, tmp1$Design)), aggfun)
      newdf <- data.frame(Radius = rep(pts, ntau*ndes*nfor), SPV = as.numeric(preds), 
                         Location = rep(rep(paste("tau =", tau), each = np), nfor*ndes),
                         Formula = factor(rep(rep(sort(fornms), each = ntau*np), ndes)), 
                         Design = factor(rep(sort(desnms), each = ntau*np*nfor)))
      if (exists("tmp3", inherits = FALSE)) tmp3 <- rbind(tmp3, newdf)
      else tmp3 <- newdf
    }
  }
    
  if (exists("tmp3", inherits = FALSE)) tmp3$Location <- as.factor(tmp3$Location)
  
  if (show[1L]){
    plot1 <- ggplot(tmp2, aes(x = Fraction, y = SPV, colour = Design)) + ggtitle("Fraction of Design Space Plot") + 
      xlab("Fraction of Design Space") + geom_line(size = lines.size) +
      theme(plot.title = element_text(vjust = 1)) + facet_wrap(~ Formula)
  }
  
  if (show[2L]) {
    plot2 <- ggplot(tmp1, aes(x = Radius, y = SPV, order = Formula)) + 
      ggtitle("Variance Dispersion Graph") +
      geom_point(alpha = alpha, colour = points.colour, size = points.size) + 
      theme(plot.title = element_text(vjust = 1)) + 
      guides(colour = guide_legend(override.aes = list(alpha = 1))) + facet_wrap(~ Design + Formula) + 
      xlab(paste0("Distance to Origin (", method,")"))
    if (hexbin){
      plot2 <- plot2 + geom_hex(bins = bins) + 
        scale_fill_gradientn(colours = rev(topo.colors(5)[-(4:5)]), name = "Frequency", na.value = NA)
    }
  }
  
  if (show[3L]) {  
      plot3 <- ggplot(tmp3, aes(x = Radius, y = SPV, colour = Design)) + ggtitle("Variance Dispersion Graph") +  
        theme(plot.title = element_text(vjust = 1)) +  
        geom_line(aes(group = interaction(Location, Formula, Design), linetype = Location), 
                  size = lines.size)  + theme(legend.text.align = 0.5) + xlab(paste0("Distance to Origin (", method,")")) +
        facet_wrap(~ Formula)
  }
  
  if (show[4L]) {
    plot4 <- ggplot(tmp1, aes(x = Radius, y = SPV, order = Formula)) + 
        ggtitle("Variance Dispersion Graph") +
        geom_point(alpha = alpha, colour = points.colour, size = points.size) + 
        theme(plot.title = element_text(vjust = 1)) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1))) + facet_wrap(~ Design + Formula) + 
        xlab(paste0("Distance to Origin (", method,")"))
    if (hexbin){
      plot4 <- plot4 + geom_hex(bins = bins) + 
        scale_fill_gradientn(colours = rev(topo.colors(5)[-(4:5)]), name = "Frequency", na.value = NA)
    }
    plot4 <- plot4 + geom_line(data = tmp3, aes(group = interaction(Location, Design, Formula), linetype = Location, 
                          order = Design), size = lines.size, colour = 1)  + 
                    theme(legend.text.align = 0.5) + xlab(paste0("Distance to Origin (", method,")")) + facet_wrap(~ Design + Formula)
  }
  
  if (show[5L]) {
    plot5 <- ggplot(tmp1, aes(x = Radius, y = SPV)) + ggtitle("Boxplots") +
      theme(plot.title = element_text(vjust = 1), 
            legend.text.align = 0.5) +
      geom_boxplot(mapping = aes(group = Radius)) + xlab(paste0("Distance to Origin (", method,")")) +
      facet_wrap(~ Design + Formula)
  }
  
  if (length(which)) {
    out <- mget(paste0("plot", which))
    names(out) <- pnms
    if (x[[1]][[1]]$unscaled) out <- lapply(out, '+', ylab("Unscaled Prediction Variance (UPV)"))
    else out <- lapply(out, '+', ylab("Scaled Prediction Variance (SPV)"))
    if (arrange) do.call(gridExtra::grid.arrange, out)
    else return(out)
  }
}