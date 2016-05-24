#' Visualization of \code{frfast} objects
#' @description Useful for drawing the estimated regression function, 
#' first and second derivative (for each factor's level).
#' @param x \code{frfast} object.
#' @param y NULL.
#' @param fac Vector which determines the level to take into account
#' in the plot. By default is \code{NULL}.
#' @param der Number or vector which determines any inference process. 
#' By default \code{der} is \code{NULL}. If this term is \code{0}, the plot 
#' shows the initial estimate. If it is \code{1} or \code{2},
#' it is designed for the first or second derivative, respectively.
#' @param points Draw the original data into the plot. By default it is
#' \code{TRUE}.
#' @param xlab A title for the \code{x} axis. 
#' @param ylab A title for the \code{y} axis. 
#' @param ylim The \code{y} limits of the plot.
#' @param main An overall title for the plot.
#' @param col A specification for the default plotting color.
#' @param CIcol A specification for the default confidence intervals
#' plotting color (for the fill).
#' @param CIlinecol A specification for the default confidence intervals
#' plotting color (for the edge).
#' @param pcol  A specification for the points color.
#' @param abline Draw an horizontal line into the plot of the second derivative 
#' of the model.
#' @param ablinecol The color to be used for \code{abline}.
#' @param lty The line type. Line types can either be specified as an integer
#' (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash, 
#' 5 = longdash, 6 = twodash).  See details in \code{\link{par}}.
#' @param CIlty The line type for confidence intervals. Line types can either 
#' be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed,
#' 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#' @param lwd The line width, a positive number, defaulting to 1.  
#' See details in \code{\link{par}}.
#' @param CIlwd The line width for confidence intervals, a positive number, 
#' defaulting to 1.
#' @param cex A numerical value giving the amount by which plotting symbols
#' should be magnified relative to the default. See details in \code{\link{par}}.
#' @param alpha Alpha transparency for overlapping elements expressed 
#' as a fraction between 0 (complete transparency) and 1 (complete opacity).
#' @param \ldots Other options.
#' 
#'@return Simply produce a plot.
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#' library(npregfast)
#' data(barnacle)
#' 
#' # Nonparametric regression without interactions
#' fit <- frfast(DW ~ RC, data = barnacle, nboot = 100) 
#' plot(fit)
#' plot(fit, der = 0)
#' plot(fit, der = 0, points = FALSE)
#' #plot(fit, der = 1, col = "red", CIcol = "blue")
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100) 
#' plot(fit2)
#' plot(fit2, der = 0, fac = "lens")
#' #plot(fit2, der = c(0,1), fac = c("barca","lens"))
#' 
#' @importFrom graphics lines par plot
#' @import ggplot2
#' @export




plot.frfast <- function(x = model, y, fac = NULL, der = NULL, points = TRUE, 
                        xlab = model$name[2], ylab = model$name[1], ylim = NULL,
                        main = NULL, col = "black", CIcol = "black", 
                        CIlinecol = "transparent", pcol = "grey80",  abline = TRUE, 
                        ablinecol = "red", lty = 1, CIlty = 2, lwd = 1, 
                        CIlwd = 1, cex = 1.4, alpha = 0.2,...) {
  
  
 
  
  model <- x
  nf <- model$nf
  fi <- length(fac)
  co <- length(der)
  facini <- fac
  
  # Control
  if ((nf == 1) & (fi >= 1)) 
    stop("Argument \"fac\" not suported. 
         There is not factor in the model.")
  
  if(!is.null(fac) &   !identical(rep(TRUE,fi), fac %in% model$label)) {
    stop("\"",paste(fac[!fac %in% model$label]),"\" is not a factor's level.")
  }
  
  if (sum(der > 2) >= 1) 
    stop("",paste(der[which(der > 2)])," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  
  if (missing(model)) 
    stop("Argument \"model\" is missing, with no default. 
         Must be a frfast object.")
  
  
  
  
  
  
  
  #  if (fi == 0 & nf > 1 | fi > 1 | co == 0 &
  #nf > 1 | co > 1 | co == 0 & nf == 1) {
  
  if (fi == 0) fi <- nf
  if (co == 0) co <- 3
  jnf <- c()
  # par(mfrow = c(fi, co))
  
  if (is.null(fac)) {
    jnf <- c(1:nf)
    fac <- model$label
  } else {
    for (i in 1:length(fac)) {
      jnf[i] <- which(model$label == fac[i])
    }
  }
  
  if (is.null(der)) {
    jder <- c(1:3)
  } else {
    jder <- der + 1
  }
  
  p <- list()
  cont <- 0
  for (j in jnf) {
    for (i in jder) {
      cont <- cont +1
      if (i == 1) {
        ylab2 <- ylab
        rgo <-  max(model$ydata[model$fmod == model$numlabel[j]], na.rm = TRUE) -
          min(model$ydata[model$fmod == model$numlabel[j]], na.rm = TRUE)
        ylim2 <- c(min(model$ydata[model$fmod == model$numlabel[j]], na.rm = TRUE) - (rgo*0.05), 
                   max(model$ydata[model$fmod == model$numlabel[j]], na.rm = TRUE) + (rgo*0.05))
      } else {
        rgo <-  max(model$p[, der = i, fac = j], na.rm = TRUE) -
          min(model$p[, der = i, fac = j], na.rm = TRUE)
        ylim2 <- c(min(model$p[, der = i, fac = j], na.rm = TRUE)-(rgo*0.05), 
                   max(model$p[, der = i, fac = j], na.rm = TRUE)+(rgo*0.05))
      }
      if (i == 2)  ylab2 <- "First derivative"
      if (i == 3)  ylab2 <- "Second derivative"
      
      if (is.null(main)){
        if(i == jder[1] & model$nf > 1){
          title <- paste("Level", model$label[j])
        }else{
          title <- NULL
        }
      }else{
        # if(length(main) == 1) main <- rep()
        if(i == jder[1]){
          title <- main[j]
        }else{
          title <- NULL
        }
      }
      
      if (is.null(ylim)) ylim <- ylim2  #### ver esto!!!!!
      
      #plot
      data_bin <- data.frame(x = model$x,
                          pl = model$pl[, der = i, fac = j],
                          pu = model$pu[, der = i, fac = j],
                          p = model$p[, der = i, fac = j])
      
      data_ori <- data.frame(xdata = model$xdata[model$fmod == model$numlabel[j]],
                          ydata = model$ydata[model$fmod == model$numlabel[j]])
      
       if ((points == TRUE) & (i == 1)) {
        points_layer <- geom_point(data = data_ori, aes_string(x = "xdata",
                                       y = "ydata"),
                                  colour = pcol, size = cex)
      }else{
        points_layer <- NULL
      }
      
      if ((i == 3) & (abline == TRUE)) {
        abline_layer <- geom_hline(yintercept = 0, colour = ablinecol)
      }else{
        abline_layer <- NULL
      }
      
      p[[cont]] <- ggplot() +
        points_layer +
        geom_ribbon(data = data_bin, aes_string(x = "x", 
                                         ymin = "pl", 
                                         ymax = "pu"), 
                    alpha = alpha, fill = CIcol, linetype = lty,
                    size = CIlwd, col = CIlinecol) +
        geom_line(data = data_bin, aes_string(x = "x", 
                                       y = "p"), 
                  size = lwd, colour = col, linetype = lty, na.rm=TRUE) +
        abline_layer +
        coord_cartesian(ylim = ylim) +
        # ylim(ylim) +
        ylab(ylab2) +
        xlab(xlab) +
        ggtitle(title)
        
  
      
      ylim <- NULL # hay que poner el ylim nulo para el siguiente plot
    }
  }
  
  if(fi*co == 1) {
    suppressWarnings(print(p[[1]]))
  }else{
  
    # NOTE: This ugly hack is here because of a bug in gridExtra which calls
    # a ggplot2 function directly instead of namespacing it.  The bug is fixed
    # in the gridExtra GitHub version, but not on CRAN. Hopefully gridExtra
    # will submit the fix to CRAN and I can remove this ugliness.
    # https://github.com/baptiste/gridextra/issues/5
    # Thanks to Dean Attali!!
    if (!"package:ggplot2" %in% search()) {
      suppressPackageStartupMessages(attachNamespace("ggplot2"))
      on.exit(detach("package:ggplot2"))
    }  
    
    
    
   args.list <- c(p, fi, co)
   names(args.list) <- c(c(rep("", fi*co), "nrow", "ncol"))
   suppressWarnings(do.call(gridExtra::grid.arrange, args.list))
  }


} 