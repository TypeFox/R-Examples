##############################################################################
#                             PLOT METHODS
##############################################################################

#-------------------------------------------------------------------#
#' Plot method for correlograms and variograms

#' @param x Result of correlogram or variogram analysis.
#' @param var Variable to plot for multiple analyses with \code{\link{eco.correlog}}
#' (see examples).
#' @seealso  \code{\link{eco.correlog}} \code{\link{eco.cormantel}}  \code{\link{eco.variogram}}
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' 
#' # single Moran's I correlogram analysis
#' moran.example <- eco.correlog(Z=eco$P[,1], eco$XY, smax=10, size=1000)
#' plot(moran.example)
#' 
#' # multiple Moran's I correlogram analysis
#' moran.example2 <- eco.correlog(Z=eco$P, eco$XY, smax=10, size=1000)
#' plot(moran.example2, var ="P2")
#' plot(moran.example2, var ="P3")
#' 
#' corm <- eco.cormantel(M = dist(eco$P), size=1000,smax=7, XY = eco$XY,
#' nsim = 99)
#' plot(corm)
#' 
#' variog <- eco.variogram(Z = eco$P[, 2],XY =  eco$XY)
#' plot(variog)
#' }
#' 
#' @rdname eco.correlog-methods
#' 
#' @aliases plot,eco.correlog-method
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @exportMethod plot 


setMethod("plot", "eco.correlog", function(x, var = NULL, 
																				 xlabel = NULL, 
																				 ylabel = NULL, 
																				 title = NULL, 
																				 legend = TRUE,
																				 background = c("grey", "white"),
																				 errorbar = FALSE,
																				 intervals = TRUE,
																				 significant = FALSE,
																				 xlim = NULL,
																				 ylim = NULL,
																				 axis.size = 14,
																				 title.size = 16,
                                         meanplot = TRUE,
                                         nsim = 999) {
	
  # tricky solution for global binding problems during check. Set the 
  # variables as NULL
  d.mean <- obs <- uppr <- lwr <- max.sd <- min.sd <- null.uppr <- null.lwr <- NULL
	
	if(length(x@OUT) == 1) {
	  var2 <- 1
	  plot.method <- "uniplot"
	} else if(is.null(var)) {
	  plot.method <- "multiplot"
	} else {
	  var2 <- which(names(x@OUT) %in% var)
	  plot.method <- "uniplot"
	}
	
  randtest <- (x@TEST)[1]
  if(length(randtest) == 0) { 
    randtest <- "none"
  }
	
	method <- (x@METHOD)[1]
	
	method2 <- pmatch(method, c("Mantel test",
															"Partial Mantel test", 
															"Moran's I", 
															"Geary's, C",
															"bivariate Moran's Ixy",
															"empirical variogram",
															"Kinship"
															))
	if(length(method2) == 0) {
		stop("invalid input to plot")
	}

	if(plot.method == "uniplot") {
	
		datos <- as.data.frame(x@OUT[[var2]])
	
	
	########## title,  x and y labels ############
		
		if(is.null(title)) {
		  title <- names(x@OUT[var2])
		}
		
		if(is.null(xlabel)) {
		  xlabel <- "Great circle distance"
		}
		
		if(is.null(ylabel)) {
		  ylabel <- method
		}
	
	########## x and y axes ############
		
	if(!is.null(xlim)) {
		xlim <- ggplot2::scale_x_continuous(limits = xlim)
	}
	
	if(!is.null(ylim)) {
		ylim <- ggplot2::scale_y_continuous(limits = ylim)
	}
	
	localenv <- environment()
	
	
	########## basic plot ############
	
	z <- ggplot2::ggplot(datos, environment = localenv) + 
		ggplot2::geom_line(ggplot2::aes(x = d.mean, y = obs)) + 
		ggplot2::xlab(xlabel) + 
		ggplot2::ylab(ylabel) + 
		ggplot2::labs(title = title) 
	
	
	######## background ########
	
	theme <- match.arg(background)
	if(theme == "grey") {
	  
	  themecol <-  ggplot2::theme_grey()
	} else {
	  themecol <- ggplot2::theme_bw()
	}
	
	z <- z + themecol

	######## permutation and bootstrap cases ########
	
	if(randtest == "permutation") {
	  #labeling S and NS points
	  pval2 <- as.numeric(datos$p.val < 0.05)
	  if(sum(pval2) == 0) {
	    labelp <- "NS"
	  } else if (sum(pval2) == length(pval2)) {
	    labelp <- "P<0.05"
	  } else {
	    labelp <- c("P < 0.05", "NS")
	  }
	  pval2[pval2 == 0] <- "#F8766D"
	  pval2[pval2 == 1] <- "#00B0F6"
	  datos$pval2 <- pval2
	  puntos <- ggplot2::geom_point(ggplot2::aes(x = d.mean, y = obs, colour = pval2), 
	                      size = 5)
	  z <- z + puntos
	}
	
	if(randtest == "bootstrap")  {
	  intervalo <- ggplot2::geom_ribbon(ggplot2::aes(x = d.mean, ymax = uppr, 
	                                                 ymin = lwr), 
	                                    fill = "blue",
	                                    alpha = 0.2)
	  
	  puntos <- ggplot2::geom_point(ggplot2::aes(x = d.mean, y = obs, colour = "blue"), 
	                                size = 5)
	  z <- z + puntos + intervalo
	  legend <- FALSE
	}
	
	if(randtest == "none" & method == "empirical variogram") {
	  z <- z + ggplot2::ylab("Semivariance") + 
	    ggplot2::geom_point(ggplot2::aes(x = d.mean,y = obs), 
	                        colour = "#F8766D", size =5)
	  
	  legend <- FALSE
	}
	
 ########## error bar ############
	
if(errorbar) {
  datos$sd <- 1.96 * datos$sd.jack
  datos$min.sd <- datos$mean.jack - datos$sd
  datos$max.sd <- datos$mean.jack + datos$sd
  z <- z + ggplot2::geom_errorbar(data =datos,
                                  ggplot2::aes(x = d.mean, 
                                               ymax = max.sd , 
                                               ymin= min.sd))
}

#######Kinship additionals #######################################

if(method == "Kinship" & intervals) {
  
  z <- z +  ggplot2::geom_line(ggplot2::aes(x = d.mean, y = null.uppr), 
                               directions = "hv", 
                               linetype = 2, 
                               colour = "red") +  
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = null.lwr),
                       direction = "hv", 
                       linetype = 2, 
                       colour = "red") + 
    ggplot2::geom_ribbon(ggplot2::aes(x = d.mean, ymin = null.lwr, 
                                      ymax = null.uppr), 
                         fill = 90,
                         alpha = 0.05) 
}

####### legend (permutation case) ##########################
	
	if(legend) {
		z <- z +  ggplot2::scale_colour_discrete(name  ="P value",
																							labels= labelp) + 
			ggplot2::theme(axis.text = ggplot2::element_text(size = axis.size), 
										 axis.title = ggplot2::element_text(size = title.size, face = "bold"), 
										 legend.position = "right",  
		plot.title = ggplot2::element_text(size = title.size, face = "bold"))
	}
	
	if(!legend) {
		z <- z + ggplot2::theme(legend.position="none", axis.text = ggplot2::element_text(size = axis.size), 
									 axis.title = ggplot2::element_text(size = title.size, face = "bold"))
	}

 ########## xlim, ylim ############
	
	if(!is.null(xlim)) {
		z <- z + xlim
	}
	if(!is.null(ylim)) {
		z <- z + ylim
	}

	print(z)
	######################
	

	} else if(plot.method == "multiplot") {
	  if(is.null(title)) {
	    title <- ""
	  }
	  
	  man <- int.multiplot(x, significant = significant, 
                         plotit = FALSE, nsim = nsim)
	 z <- plot(man, meanplot = meanplot, background = background)
   
cat("\n", "show significant =", significant, "\n\n")
cat("significant loci:","\n", man$significant.loci.names, "\n\n")
z$data <- man

cat("meanplot =", meanplot, "\n\n")
if(meanplot) {
  cat("nsim =", nsim, "(number of simulation)", "\n\n")
  cat("mean distance, observed value and jacknife CI (lwr-uppr)", "\n")
  print(round(man$mean.correlogram, 4))
}

	}

 invisible(z)
})


#-------------------------------------------------------------------#
#' plot eco.lsa
#' @rdname eco.lsa-methods
#' @aliases plot,eco.lsa-method
#' @keywords internal

setMethod("plot", "eco.lsa", function(x, ...) {
  if(x@TEST == "permutation") {
  eco.rankplot(x, ...)
  } else if(x@TEST == "bootstrap") {
    eco.forestplot(x, ...)
  }
})


#-------------------------------------------------------------------#
#' plot eco.IBD
#' @keywords internal
#' @rdname eco.IBD 
#' @aliases plot,eco.IBD-method


setMethod("plot", "eco.IBD", function(x, ...) {
    callNextMethod(...)
})


#-------------------------------------------------------------------#
#' int.multiplot method. Graphical processing of multiple correlograms 
#' @rdname int.multiplot
#' @keywords internal

int.multiplot<- function(correlog, 
                         significant = TRUE,
                         plotit = TRUE,
                         ...) {
  
  
  # tricky solution for global binding problems during check. Set the 
  # variables as NULL
  d.mean <- obs <- uppr <- lwr <- max.sd <- min.sd <- null.uppr <- null.lwr <- value <- variable <-  NULL
  
  data <-correlog@IN$Z
  N <- length(correlog@OUT)
  loci.names <- colnames(correlog@IN$Z)
  
  if(significant) {
    sign <- numeric()
    j <- 1
    for(i in 1:N) {
      if(any(correlog@OUT[[i]][, 3] < 0.05)) {
        sign[j] <- i
        j <- j + 1
      }
    }
    
    if(length(sign) == 0) {
      stop("no significant correlograms")
    }
  } else {
    sign <- 1:N
  }
  
  #mean correlogram
  mean.correlog <- matrix(0, ncol = length(sign), nrow = length(correlog@CARDINAL))
  l.sign <- 1:length(sign)
  for(j in l.sign) {
    mat <- correlog@OUT[[sign[j]]][,2]
    mean.correlog[, j] <- mat
  }
  colnames(mean.correlog) <- colnames(data)[sign]
  rownames(mean.correlog) <- rownames(correlog@OUT[[1]])
  mean.correlog <- as.data.frame(mean.correlog)
  
  intervalos <- matrix(0, nrow = length(correlog@CARDINAL), ncol = 2)
  
  #ci for mean correlogram
  intervalos <- int.jackknife(t(mean.correlog), mean)
  intervalos <- data.frame(cbind(correlog@OUT[[1]][,1],
                                 intervalos$obs, t(intervalos$CI)))
  colnames(intervalos) <- c("d.mean", "obs","lwr", "uppr")
  rownames(intervalos) <- rownames(correlog@OUT[[1]])
  intervalos <- as.data.frame(intervalos)
  
  
  salida <- list(correlogram.alleles = mean.correlog,
                 mean.correlogram = intervalos,
                 method = correlog@METHOD,
                 significant.loci.names = colnames(data)[sign],
                 significant.loci.number = sign)
  
  class(salida) <- "int.multiplot"
  
  invisible(salida)
  
}


#-------------------------------------------------------------------#
#' plot int.multiplot
#' @rdname int.multiplot
#' @aliases plot,int.multiplot
#' @keywords internal

setMethod("plot", "int.multiplot", function(x, 
                                        xlabel = NULL, 
                                        ylabel = NULL, 
                                        title = NULL, 
                                        background = c("grey", "white"),
                                        meanplot = TRUE,
                                        multiplot = TRUE) { 
  
  
  # tricky solution for global binding problems during check. Set the 
  # variables as NULL
  d.mean <- obs <- uppr <- lwr <- max.sd <- min.sd <- null.uppr <- null.lwr <- value <- variable <-  NULL
  
  datos <- data.frame(x$mean.correlogram[,1], x$correlogram.alleles)
  intervalos <- x$mean.correlogram
  colnames(datos)[1] <- "d.mean"
  
  theme <- match.arg(background)
  
  if(theme == "grey") {
    themecol <-  ggplot2::theme_grey()
  } else {
    themecol <- ggplot2::theme_bw()
  }
  
  if(is.null(xlabel)) {
    xlabel <- "Great circle distance"
  } 
  
  if(is.null(ylabel)) {
    ylabel <- x$method
  }
  
  if(is.null(title)) {
    title <- ""
  }
  
  
  localenv <- environment()
  
  #PLOT FOR VARIABLES
  
  test.data <- reshape2::melt(datos, id ="d.mean")
  
  multi.correlog <- ggplot2::ggplot(data=test.data, environment = localenv)+
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = value, colour = variable), size = 1)+
    themecol+
    ggplot2::ylab(ylabel)+
    ggplot2::xlab(xlabel)+
    ggplot2::labs(title = title) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                   axis.title = ggplot2::element_text(size = 14, face = "bold"), 
                   legend.position = "right")
  
  
  if(meanplot) {
    multi.correlog <- multi.correlog + 
      ggplot2::geom_point(data = intervalos, ggplot2::aes(x = d.mean, y = obs), size = 5)+
      ggplot2::geom_errorbar(data = intervalos, ggplot2::aes(x = d.mean, ymax = uppr, ymin= lwr), size=1)+
      ggplot2::geom_line(data = intervalos, ggplot2::aes(x = d.mean, y = obs), size= 1.8)
  }
  
  ##############
  
  #PLOT FOR MEAN
  
  intervalos <- x$mean.correlogram
  
  mean.correlog <- ggplot2::ggplot(intervalos, environment = localenv) + 
    ggplot2::geom_point(ggplot2::aes(x = d.mean, y = obs), 
                        size = 5)+  
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = obs)) + 
    ggplot2::geom_errorbar(ggplot2::aes(x = d.mean, ymax = uppr, ymin= lwr))+
    ggplot2::xlab(xlabel) + 
    ggplot2::ylab(ylabel) + 
    ggplot2::labs(title = title) +
    themecol+
    #ggplot2::labs(title = title) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                   axis.title = ggplot2::element_text(size = 14, face = "bold"), 
                   plot.title = ggplot2::element_text(size = 16, face = "bold")) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                   axis.title = ggplot2::element_text(size = 14, face = "bold"), 
                   legend.position = "right")
  
  if(meanplot) {
    print(mean.correlog)
  } 
  if(multiplot) {
    print(multi.correlog)
  }
  
  salida <- list(mean.correlog = mean.correlog, multi.correlog = multi.correlog)
  
  invisible(salida)
})

