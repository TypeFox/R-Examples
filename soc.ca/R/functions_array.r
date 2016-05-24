#' Array of maps
#' 
#' This function takes a list of map objects and arranges them into an array.
#' @param x a list of objects created by one of the mapping functions in the
#'   soc.ca package or any other ggplot2 plot
#' @param ncol the number of columns the plots are arranged into
#' @param title the main title of the array
#' @param fixed.coord if TRUE the limits of all plots are set to the same as the
#'   largest plot
#' @param padding the distance between the most extreme position and the axis
#'   limit
#' @examples
#' \dontrun{
#' example(soc.ca)
#' map.array(list(map.ind(result), map.mod(result)), ncol = 2)
#' }
#' @export
map.array   <- function(x, ncol = 1, title = "", fixed.coord = TRUE, padding = 0.15){

if (identical(fixed.coord, TRUE))  x  <- fix.coords(x, padding = padding)

do.call(grid.arrange, c(x, ncol = ncol, top = title))
}


#' Ellipse array
#' 
#' Create seperate maps with ellipses for each level in a factor arranged in an
#' array.
#' @param object a soc.ca class object
#' @param variable a factor of the same length as the data.frame used to create
#'   object
#' @param dim the dimensions in the order they are to be plotted. The first 
#'   number defines the horizontal axis and the second number defines the 
#'   vertical axis.
#' @param draw.ellipses if TRUE ellipses are drawn
#' @param ncol the number of columns the plots are arranged into
#' @param titles a vector of the same length as the number of levels in
#'   variable. These are the titles given to each subplot
#' @param main.title the main title for all the plots
#' @param ... sends any further arguments to \link{map.select} and \link{map.ellipse}.
#' @examples
#' \dontrun{
#' example(soc.ca)
#' map.ellipse.array(result, active[, 1])
#' }
#' @export
map.ellipse.array <- function(object, variable, dim = c(1,2),
                              draw.ellipses = TRUE, ncol = 2,
                              titles = levels(variable), main.title = "",  ...){
  
  var.levels      <- levels(variable)
  list.of.maps    <- list()
    
  for (i in 1:length(var.levels)){
  ind.ind         <- which(variable == var.levels[i]) 
  ellipse.ind     <- variable == var.levels[i]
  ellipse.ind[ellipse.ind == FALSE]  <- NA 
  
  p               <- map.select(object, dim = dim, list.ind = ind.ind, map.title = titles[i], label = FALSE, ...)
  
  if(identical(draw.ellipses, TRUE)) p  <- map.ellipse(object, p, ellipse.ind, ellipse.label = FALSE, label.size = 4, ...)
  list.of.maps[[i]]  <- p
  }
  
  # Standardize the coordinates
  list.of.maps <- fix.coords(list.of.maps)
  # Plot the maps
  map.array(list.of.maps, ncol = ncol, title = main.title)
  return(invisible(list.of.maps))
}

#################################################################################################
#' Array of several CSA maps
#' 
#' Creates an array of Class Specific Mulitple Correspondence analysises
#' 
#' @param object a \link{soc.ca} result object
#' @param variable a factor with the same order and length as those used for the active modalities in object
#' @param dim indicates what dimensions to map and in which order to plot them
#' @param ncol the number of columns the maps are arranged into
#' @param titles a vector of the same length as the number of levels in \code{variable}. These are the titles given to each subplot
#' @param main.title the main title for all the maps
#' @param FUN the mapping function used for the plots; \link{map.active}, \link{map.ctr}, \link{map.ind}, \link{map.select} or \link{map.sup}
#' @param fixed.coord if TRUE the limits of all plots are set to the same as the largest plot
#' @param ... sends any further arguments to the mapping functions
#' @export  
#' @examples
#' \dontrun{
#' example(soc.csa)
#' map.csa.all(result, active[, 1])
#' map.csa.all(result, active[, 1], FUN = map.ctr, ctr.dim = 1)
#' }

map.csa.all   <- function(object, variable, dim = c(1,2), ncol = 2, FUN = map.ind,
                          fixed.coord = TRUE, main.title = "", titles = levels(variable),...){
  
  csa.list     <- csa.all(object, variable)
  csa.results  <- csa.list$results
  plot.list    <- lapply(csa.results, FUN, dim = dim, ...) 
  
  # Map titles
  if (length(titles) > 1){
  for (i in 1:length(titles)){
  plot.list[[i]]  <- plot.list[[i]] + ggtitle(titles[i])
  }
  }
  
  # Coord standardization
  if (identical(fixed.coord, TRUE)) plot.list  <- fix.coords(plot.list)  
  map.array(plot.list, ncol = 2, title = main.title)
}

#######################################################################3
### Standardize coordinates

fix.coords <- function(plot.list, padding = 0.15){

minimums.x    <- vector(, length = length(plot.list)) 
minimums.y    <- vector(, length = length(plot.list))
maximums.x    <- vector(, length = length(plot.list))
maximums.y    <- vector(, length = length(plot.list))

for ( i in 1:length(plot.list)){
p                <- plot.list[[i]]
minimums.x[i]    <- p$ca.scales$lim.min.x
minimums.y[i]    <- p$ca.scales$lim.min.y
maximums.x[i]    <- p$ca.scales$lim.max.x
maximums.y[i]    <- p$ca.scales$lim.max.y
}

# Standardize the coordinates
xlim         <- c(max(maximums.x) + padding, min(minimums.x) - padding)
ylim         <- c(max(maximums.y) + padding, min(minimums.y) - padding) 

for ( i in 1:length(plot.list)){
  plot.list[[i]]  <- plot.list[[i]] + coord_fixed(xlim = xlim, ylim = ylim)
}
return(plot.list)
}

#' Map the coordinates of the individuals in a CSA and its MCA
#'
#' @param csa.object a result object created by the \link{soc.csa} function
#' @param mca.dim the dimension from the original MCA
#' @param csa.dim the dimension from the CSA
#' @param smooth if TRUE a line is added to the plot
#' @param method the method used by ggplot to set the line see \link{geom_smooth}
#' @seealso \link{soc.csa}, \link{map.csa.all}, link{map.csa.mca.array}
#' @export
#' @examples
#' example(soc.csa)
#' csa.res  <- soc.csa(result, class.age)
#' map.csa.mca(csa.res, mca.dim = 2, csa.dim = 1)
map.csa.mca <- function(csa.object, mca.dim = 1, csa.dim = 1, smooth = TRUE, method = "auto"){
  mca.res           <- csa.object$original.result
  class.indicator   <- csa.object$original.class.indicator
  mca.coord         <- mca.res$coord.ind[class.indicator, mca.dim]
  mca               <- mca.coord
  csa               <- csa.object$coord.ind[, csa.dim]
  ggdata            <- data.frame(mca = mca, csa = csa)
  titles            <- paste(c("CSA Dim:", "MCA Dim:"), c(csa.dim, mca.dim))
  
  p                 <- ggplot(ggdata, aes(x = mca, y = csa))
  if(smooth == TRUE) p   <- p +  geom_smooth(fill = "grey90", color = "red", method = method)
  p                 <- p + geom_point(shape = 21, size = 3, alpha = 0.8) + theme_min() + xlab(titles[2]) + ylab(titles[1])
  
  p$ca.scales       <- breaksandscales(ggdata)
  p
}

#' CSA-MCA array
#' 
#' Create an array of \link{map.csa.mca} maps
#' 
#' @param csa.object a result object created by the \link{soc.csa} function
#' @param ndim the number of dimensions to include in the array, starting from 1
#' @param fixed.coord if TRUE the limits of all plots are set to the same as the largest plot
#' @param ... for further arguments see \link{map.csa.mca}
#' @export
#' @examples
#' example(soc.csa)
#' csa.res <- soc.csa(result, class.age)
#' map.csa.mca.array(csa.res, ndim = 3)

map.csa.mca.array <- function(csa.object, ndim = 3, fixed.coord = TRUE, ...){
  
  plot.list  <- list()
  count      <- 1
  for(j in 1:ndim){
    for (i in 1:ndim){
      plot.list[[count]]  <- map.csa.mca(csa.object, mca.dim = i, csa.dim = j, ...)
      count  <- count + 1
    }
  }
  
  suppressMessages(print(map.array(plot.list, ncol = ndim, fixed.coord = fixed.coord)))
}
