#' Forestplot graphs
#' 
#' @description This program generates a forest plot 
#' for the confidence interval of each individual of the
#' input data (as row number)
#' and the corresponding observed value of the statistic.
#' 
#' @param input Matrix/data frame, with three columns in the 
#' folowing order: observed value, lower and upper values of the 
#' confidence interval. 
#' @param xlabel Optional label for x axis.
#' @param ylabel Optional label for y axis.
#' @param titlelabel Optional title label.
#' @param legendlabel Optional legend label.
#' 
#' @examples
#' \dontrun{
#' 
#' require(ggplot2)
#' # simulated confidence intervals for the null hypotesis of a variable "a"
#' set.seed(8)
#' a<-runif(10, -2, 2)
#' infer <- runif(10, -1, 1)
#' super <- runif(10, -1, 1)
#' infer2 <- pmin(infer, super)
#' super2 <- pmax(infer, super)
#' data <- data.frame(a, infer2, super2)
#' forest <- eco.forestplot(data)
#' forest
#' 
#' # the forestplot method support the use of ggplot2 syntax
#' forest <- forest + theme_bw() + theme(legend.position="none")
#' forest
#' }
#' 
#' 
#' @rdname forestplot-methods
#' 
#' @aliases eco.forestplot,generic-method
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @exportMethod eco.forestplot

setGeneric("eco.forestplot", 
					 function(input, 
					 				 xlabel = NULL,
					 				 ylabel = NULL,
					 				 titlelabel = NULL,
					 				 legendlabel = NULL) {
	standardGeneric("eco.forestplot")
})



#' @rdname forestplot-methods
#' @aliases eco.forestplot,eco.lsa-method
#' @exportMethod eco.forestplot

setMethod("eco.forestplot", 
					"eco.lsa",
					function(input,
									 xlabel,
									 ylabel,
									 titlelabel,
									 legendlabel) {
	
					  
	if(input@TEST != "bootstrap") {
	  stop("this method is available for eco.lsa with bootstrap test")
	}
	
	datos <- input@OUT
	
	ind <- rownames(datos)
	obs <- datos$obs
	lwr <- datos$lwr
	uppr <- datos$uppr
	method <- input@METHOD
	
	data.select <- data.frame(ind, obs, lwr, uppr)


	
	if(is.null(xlabel)) {
		xlabel <- method
	} 
	
	if(is.null(ylabel)) {
		ylabel <- "Individual"
	} 
	
	if(is.null(titlelabel)) {
		titlelabel <- " "
	}
	
	if(is.null(legendlabel)) {
		legendlabel <- paste("  ", method)
	}
	
	p <- ggplot2::ggplot(data.select, ggplot2::aes(x = c(1:length(ind)), 
																								 y = obs, 
																								 ymin = lwr,
																								 ymax = uppr)) + 
		ggplot2::geom_pointrange(ggplot2::aes(colour = obs), size=0.9) +
	  ggplot2::geom_point(size=2.5,  shape=1) +
	  ggplot2::geom_point(size=2.7,  shape=1) +

		#ggplot2::geom_line(ggplot2::aes(x= XVALUE, y = YVALUE), lwd = 1, alpha = 0.4) + 
		ggplot2::coord_flip() + 
		ggplot2::scale_colour_gradient(legendlabel, 
																	 low = "green", 
																	 high = "red") +
		ggplot2::ylab(xlabel) +
		ggplot2::xlab(ylabel) +
		ggplot2::labs(title = titlelabel)+
		ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
									 axis.title = ggplot2::element_text(size = 14, 
									 																	 face = "bold"), 
									 plot.title = ggplot2::element_text(size = 16, 
									 																	 face = "bold"))
	attr(p, "data") <- data.select
	p
	
})


#' @rdname forestplot-methods
#' @aliases eco.forestplot,dataframeORmatrix-method
#' @exportMethod eco.forestplot

setMethod("eco.forestplot", 
					"dataframeORmatrix",
					function(input,
									 xlabel,
									 ylabel,
									 titlelabel,
									 legendlabel) {
						
						
						datos <- input
						
						ind <- rownames(datos)
						obs <- datos[, 1]
						lwr <- datos[, 2]
						uppr <- datos[, 3]
						
						data.select <- data.frame(ind, obs, lwr, uppr)
						
						if(is.null(xlabel)) {
							xlabel <- "value"
						} 
						
						if(is.null(ylabel)) {
							ylabel <- "Individual"
						} 
						
						if(is.null(titlelabel)) {
							titlelabel <- ""
						}
						
						if(is.null(legendlabel)) {
							legendlabel <- ""
						}
						
						p <- ggplot2::ggplot(data.select, ggplot2::aes(x = c(1:length(ind)), 
						                                               y = obs, 
						                                               ymin = lwr,
						                                               ymax = uppr)) + 
						  ggplot2::geom_pointrange(ggplot2::aes(colour = obs), size=0.9) +
						  ggplot2::scale_colour_gradient(low = "blue", 
						                                 high = "red") +
						  ggplot2::geom_point(size=2.5,  shape=1) +
						  ggplot2::geom_point(size=2.7,  shape=1) +
							ggplot2::coord_flip() +  
							ggplot2::ylab(xlabel) +
							ggplot2::xlab(ylabel) +
							ggplot2::labs(title = titlelabel)+
							ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
														 axis.title = ggplot2::element_text(size = 14, 
														 																	 face = "bold"), 
														 plot.title = ggplot2::element_text(size = 16, 
														 																	 face = "bold"))
						
						attr(p, "data") <- data.select
						p

					})

