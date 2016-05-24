#' Visualization of the differences between the estimated curves 
#' for two factor's levels
#' @description Useful for drawing the differences between the estimation of 
#' curves (initial estimate, first or second derivative) for  two factor's levels.
#' Missing values of factor's levels is not allowed. 
#'@param model Parametric or nonparametric regression out 
#' obtained by \code{\link{frfast}} function.
#'@param level2 Second factor's level at which to perform the 
#'differences between curves.
#'@param level1 First factor's level at which to perform the 
#'differences between curves.
#' @param der Number or vector which determines any inference process. 
#' By default \code{der} is \code{NULL}. If this term is \code{0}, the plot 
#' shows the differences between estimated regression functions. If it is 
#' \code{1} or \code{2}, it is designed for the first or second derivative, 
#' respectively.
#'@param est.include Draws the estimates of the model. 
#'By default it is \code{FALSE}.
#' @param xlab A title for the x axis. 
#' @param ylab A title for the y axis. 
#' @param ylim The \code{y} limits of the plot.
#' @param main An overall title for the plot.
#' @param col A specification for the default plotting color.
#' @param CIcol A specification for the default confidence intervals
#' plotting color.
#' @param CIlinecol A specification for the default confidence intervals
#' plotting color (for the edge).
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
#' @param alpha Alpha transparency for overlapping elements expressed 
#' as a fraction between 0 (complete transparency) and 1 (complete opacity).
#' @param \ldots Other options.
#' @return Simply produce a plot.
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' @examples
#' library(npregfast)
#' data(barnacle)
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100) 
#' plotdiff(fit2, level2 = "lens", level1 = "barca")
#' plotdiff(fit2, level2 = "lens", level1 = "barca", der = 1, col = "blue", CIcol = "grey")
#' plotdiff(fit2, "lens", "barca", der = c(0, 1), ylim = c(-0.05, 0.05))
#' 
#' 
#' 
#' @importFrom graphics lines par plot
#' @import ggplot2
#' @export


plotdiff <- function(model, level2, level1, der = NULL, est.include = FALSE,
                     xlab = model$name[2], ylab = model$name[1], ylim = NULL,
                     main = NULL, col = "black", CIcol = "black",
                     CIlinecol = "transparent", abline = TRUE, ablinecol = "red",
                     lty = 1, CIlty = 2, lwd = 1, CIlwd = 1.5, 
                     alpha = 0.2, ...) {
  nf <- model$nf
  # co=length(der)
  jnf <- c()
  jnf[1] <- which(model$label == level1)  #'B' plot.diff(ajus,'A','B');plot.diff(ajus,1,2) 
  jnf[2] <- which(model$label == level2)  #'A'
  # if(length(der)==0) {jder=c(1:3)}else{jder=der+1}
  
  ## Argumentos control
  if (missing(model)) 
    stop("Argument \"model\" is missing, with no default. 
         Must be a frfast object.")
  
  if (missing(level1) & missing(level2)) 
    stop("Argument 'level1' and/or 'level2' are missing, with no default")
  
  if (level1 == level2) 
    stop("Argument 'level1' and 'level2' are not different")
  
  if(!isTRUE(level1 %in% model$label)) {
    stop("\"",paste(level1),"\" is not a factor's level.")
  }
  
  if(!isTRUE(level2 %in% model$label)) {
    stop("\"",paste(level2),"\" is not a factor's level.")
  }
  
  if (sum(der > 2) >= 1) 
    stop("",paste(der[which(der > 2)])," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  
  if ((nf == 1)) 
    stop("Function \"plot.diff\" not suported. 
         There is not factor in the model.")
  
  
  p <-list()
  
  if (est.include == FALSE) {
    if (is.null(der)) der <- c(0, 1, 2)
    
    der <- der + 1
    
   #par(mfrow = c(1, length(der)))
    cont = 0
    for (i in der) {
      cont = cont + 1
      if (i == 1) ylab2 <- ylab
      if (i == 2) ylab2 <- "First derivative"
      if (i == 3) ylab2 <- "Second derivative"
      if (sum(model$diff[, der = i, jnf[2], jnf[1]], na.rm = T) == 0) { # para ver si -1* o no
        
        if (is.null(ylim)) {
          rgo <- max(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T) -
                       min(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T)
                       
          ylim <- c(min(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T) - 
                      (rgo * 0.05),
                    max(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T) + 
                      (rgo * 0.05))
        }
        if (is.null(main)){
          if(i == der[1] ){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == der[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        
        data_bin <- data.frame(x = model$x,
                               p = -1 * model$diff[, der = i, jnf[1], jnf[2]],
                               pl = -1 * model$diffl[, der = i, jnf[1], jnf[2]],
                               pu = -1 * model$diffu[, der = i, jnf[1], jnf[2]])
        
        if (abline == TRUE){
          abline_layer <- geom_hline(yintercept = 0, colour = ablinecol)
        }else{
          abline_layer <- NULL
        }
        

        p[[cont]] <- ggplot() +
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
        
        
      } else {
        
        if (is.null(ylim)) {
          rgo <- max(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T) -
            min(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T)
          ylim <- c(min(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T) - 
                      (rgo * 0.05),
                    max(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T) +
                      (rgo * 0.05))
        }
        if (is.null(main)){
          if(i == der[1] ){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == der[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        
        data_bin <- data.frame(x = model$x,
                               p = model$diff[, der = i, jnf[2], jnf[1]],
                               pl = model$diffl[, der = i, jnf[2], jnf[1]],
                               pu = model$diffu[, der = i, jnf[2], jnf[1]])
        
        
        
        if (abline == TRUE){
          abline_layer <- geom_hline(yintercept = 0, colour = ablinecol)
        }else{
          abline_layer <- NULL
        }
        
        
        
        p[[cont]] <- ggplot() +
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

      }
      ylim <- NULL # hay que poner el ylim nulo para el siguiente plot
    }
    
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
    
    
    
    args.list <- c(p, 1, length(der))
    names(args.list) <- c(c(rep("", length(der)), "nrow", "ncol"))
    suppressWarnings(do.call(gridExtra::grid.arrange, args.list))
    
    
  } else {  # est.include = TRUE
    
    if (is.null(der))  der <- c(0, 1, 2)
    jder <- der + 1  #if(length(der)==0) {jder=c(1:3)}else{jder=der+1}
   # par(mfrow = c(nf + 1, length(der)))
    cont = 0
    for (i in jder) {
      cont = cont + 1
      if (i == 1) 
        ylab2 <- ylab
      if (i == 2) 
        ylab2 <- "First derivative"
      if (i == 3) 
        ylab2 <- "Second derivative"
      if (sum(model$diff[, der = i, jnf[2], jnf[1]], na.rm = T) == 0) {
        
        if (is.null(ylim)) {
          rgo <- max(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T) -
            min(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T)
          ylim <- c(min(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T) - 
                      (rgo * 0.05),
                    max(-1 * (model$diff[, der = i, jnf[1], jnf[2]]), na.rm = T) +
                      (rgo * 0.05))
        }

        
        
        if (is.null(main)){
          if(i == jder[1] ){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == jder[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        

        
        data_bin <- data.frame(x = model$x,
                               p = -1 * model$diff[, der = i, jnf[1], jnf[2]],
                               pl = -1 * model$diffl[, der = i, jnf[1], jnf[2]],
                               pu = -1 * model$diffu[, der = i, jnf[1], jnf[2]])
        
        
        
        if (abline == TRUE){
          abline_layer <- geom_hline(yintercept = 0, colour = ablinecol)
        }else{
          abline_layer <- NULL
        }
        
        
        
        p[[cont]] <- ggplot() +
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


        
        for (j in length(jnf):1) {
          cont = cont + 1
          if (jnf[j] == jnf[2]) {
            title <- paste("Level", model$label[jnf[2]])
          } 
          if  (jnf[j] == jnf[1]) {
            title <- paste("Level", model$label[jnf[1]])
          }
          
          rgo <- max(model$p[, der = i, jnf[j]], na.rm = T) -
            min(model$p[, der = i, jnf[j]], na.rm = T)
          ylim <- c(min(model$p[, der = i, jnf[j]], na.rm = T) - (rgo * 0.05), 
                    max(model$p[, der = i, jnf[j]], na.rm = T) + (rgo * 0.05))
          
          data_bin <- data.frame(x = model$x,
                                 pl = model$pl[, der = i, fac = jnf[j]],
                                 pu = model$pu[, der = i, fac = jnf[j]],
                                 p = model$p[, der = i, fac = jnf[j]])
          
          p[[cont]] <- ggplot() +
            geom_ribbon(data = data_bin, aes_string(x = "x", 
                                             ymin = "pl", 
                                             ymax = "pu"), 
                        alpha = alpha, fill = CIcol, linetype = lty,
                        size = CIlwd, col = CIlinecol) +
            geom_line(data = data_bin, aes_string(x = "x", 
                                           y = "p"), 
                      size = lwd, colour = col, linetype = lty, na.rm=TRUE) +
            coord_cartesian(ylim = 
                              ylim) +
            # ylim(ylim) +
            ylab(ylab2) +
            xlab(xlab) +
            ggtitle(title)
        }
      
        ylim <- NULL # hay que poner el ylim nulo para el siguiente plot
        
       
        
      } else {
        
        
        if (is.null(main)){
          if(i == jder[1]){
            title <- "Differences"
          }else{
            title <- NULL
          }
        }else{
          if(i == jder[1]){
            title <- main
          }else{
            title <- NULL
          }
        }
        
        
        
        if (is.null(ylim)) {
          rgo <- max(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T) -
            min(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T)
          ylim <- c(min(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T) - 
                      (rgo * 0.05),
                    max(1 * (model$diff[, der = i, jnf[2], jnf[1]]), na.rm = T) +
                      (rgo * 0.05))
        }
        
        
        data_bin <- data.frame(x = model$x,
                               p = model$diff[, der = i, jnf[2], jnf[1]],
                               pl = model$diffl[, der = i, jnf[2], jnf[1]],
                               pu = model$diffu[, der = i, jnf[2], jnf[1]])
        
        
        
        if (abline == TRUE){
          abline_layer <- geom_hline(yintercept = 0, colour = ablinecol)
        }else{
          abline_layer <- NULL
        }
        
        
        
        p[[cont]] <- ggplot() +
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

        
        for (j in length(jnf):1) {
          cont = cont + 1
          if (jnf[j] == jnf[2]) {
            title <- paste("Level", model$label[jnf[2]])
          } 
          if (jnf[j] == jnf[1]) {
            title <- paste("Level", model$label[jnf[1]])
          } 
          
          rgo <- max(model$p[, der = i, jnf[j]], na.rm = T) -
            min(model$p[, der = i, jnf[j]], na.rm = T)
          ylim <- c(min(model$p[, der = i, jnf[j]], na.rm = T) - (rgo * 0.05), 
                    max(model$p[, der = i, jnf[j]], na.rm = T) + (rgo * 0.05))
          
          data_bin <- data.frame(x = model$x,
                                 pl = model$pl[, der = i, fac = jnf[j]],
                                 pu = model$pu[, der = i, fac = jnf[j]],
                                 p = model$p[, der = i, fac = jnf[j]])
          
          p[[cont]] <- ggplot() +
            geom_ribbon(data = data_bin, aes_string(x = "x", 
                                             ymin = "pl", 
                                             ymax = "pu"), 
                        alpha = alpha, fill = CIcol, linetype = lty,
                        size = CIlwd, col = CIlinecol) +
            geom_line(data = data_bin, aes_string(x = "x", 
                                           y = "p"), 
                      size = lwd, colour = col, linetype = lty, na.rm=TRUE) +
            coord_cartesian(ylim = ylim) +
            # ylim(ylim) +
            ylab(ylab2) +
            xlab(xlab) +
            ggtitle(title)
          
          
         
        }
        ylim <- NULL # hay que poner el ylim nulo para el siguiente plot
        }
    }
    
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
    
    
      args.list <- c(p, nf + 1, length(der))
      names(args.list) <- c(c(rep("", (nf + 1) * length(der)), "nrow", "ncol"))
      suppressWarnings(do.call(gridExtra::grid.arrange, args.list))
      
    }
  }
 