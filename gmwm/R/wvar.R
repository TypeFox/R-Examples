# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' @title Wavelet Variance
#' @description Calculates the (MODWT) wavelet variance
#' @param x         A \code{vector} with dimensions N x 1, or a \code{lts} object, or a \code{gts} object, or a \code{imu} object. 
#' @param decomp    A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition.
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @param robust    A \code{boolean} that triggers the use of the robust estimate.
#' @param eff       A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param alpha     A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
#' @param freq      A \code{numeric} that provides the rate of samples.
#' @param from.unit A \code{string} indicating the unit which the data is converted from.
#' @param to.unit   A \code{string} indicating the unit which the data is converted to.
#' @param ... Further arguments passed to or from other methods.
#' @return A \code{list} with the structure:
#' \describe{
#'   \item{"variance"}{Wavelet Variance}
#'   \item{"ci_low"}{Lower CI}
#'   \item{"ci_high"}{Upper CI}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"alpha"}{p value used for CI}
#'   \item{"unit"}{String representation of the unit}
#' }
#' @details 
#' If `nlevels` is not specified, it is set to floor(log2(length(x)))
#' @author JJB
#' @rdname wvar
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' # Default
#' wvar(x)
#' # Robust
#' wvar(x, robust = TRUE, eff=0.3)
#' # 90% confidence interval
#' wvar(x, alpha = 0.10)
#' 
#' # IMU Object
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyros = 1:3, accels = 4:6, freq = 100)
#' df = wvar.imu(test)
#' }
#' @export
wvar = function(x, ...) {
  UseMethod("wvar")
}

#' @rdname wvar
#' @export
wvar.lts = function(x, decomp = "modwt", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  warning('`lts` object is detected. This function can only operate on the combined process.')
  freq = attr(x, 'freq')
  unit = attr(x, 'unit')
  x = x[,ncol(x)]

  wvar.default(x, decomp, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.gts = function(x, decomp="modwt", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  freq = attr(x, 'freq')
  unit = attr(x, 'unit')
  x = x[,1]
  
  wvar.default(x, decomp, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.default = function(x, decomp = "modwt", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, freq = 1, from.unit = NULL, to.unit = NULL, ...){
  if(is.null(x)){
    stop("`x` must contain a value")
  }else if((is.data.frame(x) || is.matrix(x))){
    if(ncol(x) > 1) stop("There must be only one column of data supplied.")
  }
  
  if(is.null(nlevels)){
    nlevels = floor(log2(length(x)))
  }

  # check freq
  if(!is(freq,"numeric") || length(freq) != 1){ stop("'freq' must be one numeric number.") }
  if(freq <= 0) { stop("'freq' must be larger than 0.") }
  
  # check unit
  all.units = c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')
  if( (!is.null(from.unit) && !from.unit %in% all.units) || (!is.null(to.unit) && !to.unit %in% all.units) ){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
  }
  
  obj =  .Call('gmwm_modwt_wvar_cpp', PACKAGE = 'gmwm',
               signal=x, nlevels=nlevels, robust=robust, eff=eff, alpha=alpha, 
               ci_type="eta3", strWavelet="haar", decomp = decomp)

  scales = .Call('gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)/freq
  
  # NO unit conversion
  if( is.null(from.unit) && is.null(to.unit)==F ){
    warning("'from.unit' is NULL. Unit conversion was not done.")
  }
  
  # unit conversion
  if (!is.null(from.unit)){
    if (!is.null(to.unit)){
      convert.obj = unitConversion(scales, from.unit = from.unit, to.unit = to.unit)
      
      if (convert.obj$converted) {
        # YES unit conversion
        scales = convert.obj$x
        message(paste0('Unit of object is converted from ', from.unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  if(!is.null(from.unit) && !is.null(to.unit)){ 
    unit = to.unit
  }else{
    unit = from.unit}
  
  create_wvar(obj, decomp, robust, eff, alpha, scales, unit)
}

#' @rdname wvar
#' @export
wvar.imu = function(x, decomp = "modwt", nlevels = NULL, alpha = 0.05, robust = F, eff = 0.6, to.unit = NULL, ...){

  if(!is.imu(x)){
    stop("`wvar.imu()` requires an IMU Object")
  }
  
  mlevels = floor(log2(nrow(x)))
  
  if(is.null(nlevels)){
    nlevels = mlevels
  }
  if(nlevels > mlevels){
    stop("`nlevels` must be less than ", mlevels,", which is the max number of levels.")
  }
  
  # freq conversion
  x.freq = attr(x, 'freq')
  scales = .Call('gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)/x.freq
  
  # NO unit conversion
  from.unit = attr(x, 'unit')
  if( is.null(from.unit) && is.null(to.unit)==F ){
    warning("The unit of the object is NULL. Unit conversion was not done.")
  }
  
  # unit conversion
  if (!is.null(from.unit)){
    if (!is.null(to.unit)){
      # start_end = c(start, end)
      obj = unitConversion(scales, from.unit = from.unit, to.unit = to.unit)
      
      if (obj$converted) {
        # YES unit conversion
        scales = obj$x
        
        message(paste0('Unit of object is converted from ', from.unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  obj.list = .Call('gmwm_batch_modwt_wvar_cpp', PACKAGE = 'gmwm', 
                   x, nlevels, robust, eff, alpha, ci_type="eta3", strWavelet="haar", decomp)

  total.len = nlevels*ncol(x)
  
  # Initialize empty data frame with right number of rows
  obj = data.frame(WV = numeric(total.len),
                   scales = numeric(total.len),
                   low = numeric(total.len),
                   high = numeric(total.len),
                   axis = character(total.len),
                   sensor = character(total.len), stringsAsFactors=FALSE)
  
  # Axis length
  x.axis = attr(x, 'axis')
  naxis = length(x.axis)
  
  # The correct unit
  if(!is.null(from.unit) && !is.null(to.unit)){ 
    unit = to.unit
  }else{
    unit = from.unit}
  
  # Put data into data frame
  t = 1
  for (i in 1:ncol(x)){
    # Cast for Analytical IMU Results
    obj.list[[i]] = create_wvar(obj.list[[i]], decomp, robust, eff, alpha, scales, unit)

    # Cast for Graphing IMU Results
    x.num.sensor = attr(x, 'num.sensor')
    obj[t:(t+nlevels-1),] = data.frame(WV = obj.list[[i]]$variance,
                                       scales = scales,
                                       low = obj.list[[i]]$ci_low,
                                       high = obj.list[[i]]$ci_high,
                                       axis = x.axis[(i-1)%%naxis+1], 
                                       sensor = if(i <= x.num.sensor[1]){"Gyroscope"}else{"Accelerometer"},
                                       stringsAsFactors=FALSE)
    t = t + nlevels
  }
  
  out = structure(list(dataobj=obj.list, plotobj=obj, unit = unit), class="wvar.imu")
  
  out
}


#' @title Create a Wvar object
#' @description Structures elements into a WVar object
#' @param obj    A \code{matrix} with dimensions N x 3, that contains the wavelet variance, low ci, hi ci.
#' @param decomp A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition
#' @param robust A \code{boolean} that triggers the use of the robust estimate.
#' @param eff    A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param alpha  A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
#' @param scales A \code{vec} that contains the amount of decomposition done at each level.
#' @param unit   A \code{string} that contains the unit expression of the frequency.
#' @return A \code{list} with the structure:
#' \describe{
#'   \item{"variance"}{Wavelet Variance}
#'   \item{"ci_low"}{Lower CI}
#'   \item{"ci_high"}{Upper CI}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"alpha"}{p value used for CI}
#'   \item{"unit"}{String representation of the unit}
#' }
#' @keywords internal
create_wvar = function(obj, decomp, robust, eff, alpha, scales, unit){
  structure(list(variance = obj[,1],
                       ci_low = obj[,2], 
                       ci_high = obj[,3], 
                       robust = robust, 
                       eff = eff,
                       alpha = alpha,
                       scales = scales,
                       decomp = decomp,
                       unit = unit), class = "wvar")
}

#' @title Print Wavelet Variances
#' @description Displays the summary table of wavelet variance.
#' @method print wvar
#' @export
#' @keywords internal
#' @param x A \code{wvar} object.
#' @param ... further arguments passed to or from other methods.
#' @author JJB
#' @return Summary table
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = wvar(x)
#' print( out )
print.wvar = function(x, ...){
  mat = matrix(unlist(x[1:3]),ncol=3,byrow=F)
  colnames(mat) = c("Variance", "Low CI", "High CI")
  rownames(mat) = x$scales
  print(mat)
}

#' @title Summary of Wavelet Variances
#' @description Displays the summary table of wavelet variance in addition to CI values and supplied efficiency.
#' @method summary wvar
#' @export
#' @keywords internal
#' @param object A \code{wvar} object.
#' @return Summary table and other properties of the object.
#' @param ... additional arguments affecting the summary produced.
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(x)
#' summary( out )
summary.wvar = function(object, ...){
  name = if(object$robust){
    "robust" 
  }else{
    "classical"
  }
  cat("Results of the wavelet variance calculation using the ",name, " method.\n",sep="")
  if(object$robust){
    cat("Robust was created using efficiency=",object$eff,"\n",sep="")
  }
  
  cat("The confidence interval was generated using (1-",object$p*2,")*100 \n",sep="")
  
  print(object)
}


#' @title Wrapper to ggplot Wavelet Variances Graph
#' @description Creates the wavelet variance graph
#' @method plot wvar
#' @export
#' @param x A \code{wvar} object.
#' @template CommonParams
#' @return A ggplot2 graph containing the wavelet variances.
#' @note Parameter line.type, line.color, point.size, point.shape, legend.label must contain 2 elements.
#' @author JJB, Wenchao
#' @seealso \code{\link{autoplot.wvar}}
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(x)
#' plot( out )
plot.wvar = function(x, transparence = 0.1, background = 'white', bw = F, 
                     CI.color = "#003C7D", line.type = c('solid','dotted'), line.color = c('#003C7D', '#999999'),
                     point.size = c(5,0), point.shape = c(20,46),
                     title = NA, title.size= 15, 
                     axis.label.size = 13, axis.tick.size = 11, 
                     axis.x.label = expression(paste("Scale ", tau)),
                     axis.y.label = expression(paste("Wavelet Variance ", nu)),
                     legend.title = '',  legend.label = NULL,
                     legend.key.size = 1, legend.title.size = 13, 
                     legend.text.size = 13, ...){
  autoplot.wvar(x, transparence = transparence, background = background, bw = bw, 
                CI.color = CI.color, line.type = line.type, line.color = line.color,
                point.size = point.size, point.shape = point.shape,
                title = title, title.size= title.size, 
                axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                axis.x.label = axis.x.label,
                axis.y.label = axis.y.label,
                legend.title = legend.title,  legend.label = legend.label,
                legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                legend.text.size = legend.text.size )
}

#' @title Graph Wavelet Variances
#' @description Creates the wavelet variance graph
#' @method autoplot wvar
#' @export
#' @keywords internal
#' @param object A \code{wvar} object.
#' @template CommonParams
#' @return A ggplot2 graph containing the wavelet variances.
#' @note Parameter line.type, line.color, point.size, point.shape, legend.label must contain 2 elements.
#' @author JJB, Wenchao
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(x)
#' autoplot( out )
autoplot.wvar = function(object, transparence = 0.1, background = 'white', bw = F, 
                         CI.color = "#003C7D", line.type = c('solid','dotted'), line.color = c('#003C7D', '#999999'),
                         point.size = c(5,0), point.shape = c(20,46),
                         title = NA, title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)),
                         legend.title = '',  legend.label =  NULL,
                         legend.key.size = 1, legend.title.size = 13, 
                         legend.text.size = 13, ...){
  
  .x=low=high=trans_breaks=trans_format=math_format=value=variable=NULL
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  #check parameter
  params = c('line.type', 'line.color', 'point.size', 'point.shape', 'legend.label')
  requireLength = c(2, 2, 2, 2, 2)
  legend.label.default = c(bquote("Empirical WV"~hat(nu)), bquote("CI("*hat(nu)*", "*.(1 - object$alpha)*")" )) 
  #legend.label.default = c(expression(paste("Empirical WV ", hat(nu))), expression(paste("CI(", hat(nu)," ,", 1 - object$alpha, ")" )) )
  default = list(c('solid','dotted'), c('#003C7D', '#999999'),  c(5, 0), c(20,46), legend.label.default)
  nullIsFine = c(rep(F,4), T)
  for (i in 1:length(params)){
    one_param = params[i]
    if( length(get(one_param))!=requireLength[i]){
      isNull = is.null(get(one_param))
      if(isNull && nullIsFine[i]){}else{
        warning(paste('Parameter', one_param, 'requires', requireLength[i],'elements,','but', length(get(one_param)),
                      'is supplied.','Default setting is used.'))
      }
      assign(one_param, default[[i]])
    }
  }
  
  if(bw){
    line.color = c("#000000", "#404040")
    CI.color = "grey50"
  }
  
  #process parameter (insert some values)
  params = params[-5];from = 2; to = 3; times = 1;
  for(i in 1:length(params)){
    real_param = get(params[i])
    target = real_param[from]
    stuff = rep(target, times)
    one_param = params[i]
    
    assign(one_param, c(real_param, stuff))
  }

  #other parameter
  breaks = c('var', 'low')
  legend.color = c(NA, alpha(CI.color, transparence) )
  legend.linetype = c(line.type[1], 'blank')
  legend.pointshape = c(point.shape[1], NA)
  
  WV = data.frame(var = object$variance, low = object$ci_low, high = object$ci_high, scale = object$scales)
  melt.wv = melt(WV, id.vars = 'scale')
  p = ggplot() + geom_line(data = melt.wv, mapping = aes(x = scale, y = value, color = variable, linetype = variable)) +
    geom_point(data = melt.wv, mapping =aes(x = scale, y = value, color = variable, size = variable, shape = variable)) +
    
    scale_linetype_manual(name = legend.title, values = c(line.type), breaks = breaks, labels = legend.label ) +
    scale_shape_manual(name = legend.title, values = c(point.shape), breaks = breaks, labels = legend.label)+
    scale_size_manual(name = legend.title, values = c(point.size), breaks = breaks, labels = legend.label) +
    scale_color_manual(name = legend.title,values = c(line.color), breaks = breaks, labels = legend.label) +
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  #if(CI){
  p = p + geom_ribbon(data = WV, mapping = aes(ymin = low , ymax = high, x = scale, y = NULL), alpha = transparence, fill = CI.color, show.legend = T) +
    guides(colour = guide_legend(override.aes = list(fill = legend.color, linetype = legend.linetype, shape = legend.pointshape)))
  #}
  if( background == 'white'||bw){
    p = p + theme_bw() 
  }
  
  #decide where to place the legend
  legendPlace = placeLegend(WV$var[1], WV$low[ length(WV$low) ], WV$high[ length(WV$high)])  
  p = p +
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      legend.justification=legendPlace[1:2], legend.position=legendPlace[3:4],
      legend.text.align = 0)  
  
  if (is.na(title)){
    name = if(object$robust){ "Robust"} else{ "Classic" }
    p = p +
      ggtitle(paste0("Haar Wavelet Variance Representation for ", name, " Calculation"))
  }
  p
}

#' @title Detail Implementation to Compare Wavelet Variances
#' @description Compare the estimates given by the classical and robust methods of calculating the wavelet variance.
#' @method autoplot wvarComp
#' @export
#' @keywords internal
#' @param object A \code{data frame} that contains data in order to plot
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE)
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of confidence interval. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines for wavelet variance and the edge of confidence interval, respectively. Length of vector must equal to 2. 
#' @param point.size A \code{vector} of \code{integer} that indicates the size of point
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of point
#' @param title A \code{string} that indicates the title of the graph 
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param axis.x.label A \code{string} that indicates the label on x axis
#' @param axis.y.label A \code{string} that indicates the label on y axis
#' @param facet.label.size An \code{integer} that indicates the size of facet label
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend
#' @param legend.title.size An \code{integer} that indicates the size of title on legend.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param nrow An \code{integer} that indicates number of rows
#' @param ... Additional parameters
#' @author JJB, Wenchao
#' @seealso \code{\link{compare.wvar}}
autoplot.wvarComp = function(object, split = TRUE, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                             CI.color = NULL, line.type = NULL, point.size = NULL, point.shape = NULL,
                             title = "Haar Wavelet Variance Representation", title.size= 15, 
                             axis.label.size = 13, axis.tick.size = 11, 
                             axis.x.label = expression(paste("Scale ", tau)),
                             axis.y.label = expression(paste("Wavelet Variance ", nu)),
                             facet.label.size = 13,
                             legend.label = NULL,
                             legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                             legend.text.size = 13, nrow = 1, ...){
  
  scales=low=high=WV=emp=theo=trans_breaks=trans_format=math_format=.x=dataset=NULL
  
  #p = ggplot(data = obj, mapping = aes(x = scales, y = WV)) + geom_line(mapping = aes(color = dataset), linetype = line.type[1]) + geom_point(mapping = aes(color = dataset), size = point.size, shape = point.shape) +
  p = ggplot(data = object, mapping = aes(x = scales, y = WV)) + geom_line(mapping = aes(color = dataset), linetype = line.type[1]) + geom_point(mapping = aes(color = dataset, size = dataset, shape = dataset) ) +
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) 
  if (!is.null(line.color)){
    #legend.label should work. Not work here. But it is changed when creating 'obj' (in wrapper function)
    p = p + scale_color_manual(name = legend.title, values = line.color , labels = legend.label)
  } 
      
  p = p + scale_size_manual(name = legend.title, values = point.size , labels = legend.label) +
      scale_shape_manual(name = legend.title, values = point.shape , labels = legend.label)

  if(CI){
    p = p + 
      geom_line(mapping = aes(y = low, color = dataset), linetype = line.type[2]) + geom_line(mapping = aes(y = high, color = dataset), linetype = line.type[2]) + 
      geom_ribbon(mapping = aes(ymin = low, ymax = high, fill = dataset), alpha = transparence) 
    if(!is.null(CI.color)){
      p = p + scale_fill_manual(name = legend.title, values = alpha(CI.color, transparence) , labels = legend.label)
    }
  }
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  if (split){
    p = p + facet_wrap(~dataset,nrow = nrow) + theme(legend.position="none")
  }else{
    if(is.null(line.color)){
      p = p + scale_colour_hue(name = legend.title)
    }
    if(is.null(CI.color)){
      p = p + scale_fill_discrete(name = legend.title)
    }
  }
  
  p = p +  xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      strip.text = element_text(size = facet.label.size)) 
  
  p
}


#' @title Compare Wavelet Variances
#' @description Compare the estimates given by the classical and robust methods of calculating the wavelet variance.
#' @param ... Any number of \code{wvar} objects can be passed in.
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE)
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param auto.label.wvar A \code{boolean} that indicates whether legend label should indicate the \code{wvar} objects are robust or classical
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of confidence interval. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines for wavelet variance and the edge of confidence interval, respectively. If not \code{NULL}, length of vector must equal to 2.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of point
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of point 
#' @param title A \code{string} that indicates the title of the graph
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param axis.x.label A \code{string} that indicates the label on x axis
#' @param axis.y.label A \code{string} that indicates the label on y axis
#' @param facet.label.size An \code{integer} that indicates the size of facet label
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend
#' @param legend.title.size An \code{integer} that indicates the size of title on legend
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param nrow An \code{integer} that indicates number of rows
#' @author JJB, Wenchao
#' @note 
#' Common error "Error in grid.Call(L_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : polygon edge not found"
#' Just run your code again or open R using the R GUI instead of RStudio.
#' @examples
#' \dontrun{
#' #1. Compare two objects
#' N1 = 1000
#' N2 = 2000
#' data.ar = gen.gts(AR1(phi = .32, sigma2=.01), N1)
#' data.arma = gen.gts(ARMA(ar=c(.8,.1), ma=c(.3), sigma2=1), N2)
#' wvar1 = wvar(data.ar)
#' wvar2 = wvar(data.arma, robust = T)
#' compare.wvar(wvar1, wvar2)
#' compare.wvar(wvar1, wvar2, split=F)
#' compare.wvar(wvar1, wvar2, CI = F)
#' compare.wvar(wvar1, wvar2, split=F, CI = F)
#' #2. Compare multiple objects
#' N1 = 1000
#' N2 = 2000
#' N3 = 4000
#' N4 = 3500
#' data1 = gen.gts(AR1(phi = .32, sigma2=.01), N1)
#' data2 = gen.gts(ARMA(ar=c(.8,.1), ma=c(.3), sigma2=1), N2)
#' data3 = gen.gts(AR1(phi = .32, sigma2=1), N3)
#' data4 = gen.gts(ARMA(ar=c(.8,.1), ma=c(.5), sigma2=1), N4)
#' wvar1 = wvar(data1)
#' #wvar1 = wvar(data1, robust = T)
#' wvar2 = wvar(data2)
#' wvar3 = wvar(data3)
#' wvar4 = wvar(data4)
#' compare.wvar(wvar1,wvar2,wvar3,wvar4, nrow = 2)
#' compare.wvar(wvar1,wvar2,wvar3,wvar4, split = F , CI = F)
#' #3. Change default setting
#' compare.wvar(wvar1, wvar2, wvar3,wvar4, CI.color = c('green','red','blue','black'))
#' compare.wvar(wvar1, wvar2, wvar3,wvar4, CI.color = c('green','red','blue','black'), 
#' facet.label.size = 9)
#' compare.wvar(wvar1, wvar2, wvar3,wvar4, CI.color = c('green','red','blue','black'), 
#' legend.label = c('1','2','3','4'))
#' compare.wvar(wvar1, wvar2, wvar3,wvar4, CI.color = c('green','red','blue','black'), 
#' legend.label = c('1','2','3','4'), split = F)
#' compare.wvar(wvar1, wvar2, wvar3,wvar4, CI.color = c('green','red','blue','black'), 
#' legend.label = c('1','2','3','4'), split = F, CI = F)
#' }
compare.wvar = function(..., background = 'white', split = TRUE, CI = TRUE, auto.label.wvar = T, transparence = 0.1, line.color = NULL, 
                        CI.color = NULL, line.type = NULL,  point.size = NULL, point.shape = NULL,
                        title = "Haar Wavelet Variance Representation", title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = expression(paste("Scale ", tau)),
                        axis.y.label = expression(paste("Wavelet Variance ", nu)),
                        facet.label.size = 13,
                        legend.label = NULL,
                        legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                        legend.text.size = 13, nrow = 1 ){
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  obj_list = list(...)
  numObj = length(obj_list)
  
  # check parameter
  params = c('line.type', 'line.color', 'CI.color', 'point.size', 'point.shape', 'legend.label')
  requireLength = c(2, numObj, numObj, numObj, numObj, numObj)
  default = list(c('solid','dotted'), NULL,  NULL, rep(5, numObj), rep(20, numObj), NULL)
  nullIsFine = c(rep(T,6))
  for (i in 1:length(params)){
    one_param = params[i]
    if( length(get(one_param))!=requireLength[i]){
      isNull = is.null(get(one_param))
      if(isNull && nullIsFine[i]){}else{
        warning(paste('Parameter', one_param, 'requires', requireLength[i],'elements,','but', length(get(one_param)),
                      'is supplied.','Default setting is used.'))
      }
      assign(one_param, default[[i]])
    }
  }
  
  if (numObj == 0){
    stop('At least one wvar object should be given')
  }
  else if (numObj == 1){
    ## just plot
    plot(...)
  }
  else  {
    
    if(is.null(legend.label)){
      legend.label = c()
      for (i in 1:numObj){
        legend.label[i] = paste('Dataset',i)
      }
    }
    
    if(auto.label.wvar){
      for (i in 1:numObj){
        legend.label[i] = paste(legend.label[i], if(obj_list[[i]]$robust) '(Robust)' else '(Classical)')
      }
    }
    
    total.len = 0
    each.len = numeric(numObj)
    for (i in 1:numObj){
      each.len[i] = length(obj_list[[i]]$variance)
      total.len = total.len + each.len[i]
    }
    #Initialize empty data frame with right number of rows
    obj = data.frame(WV = numeric(total.len),
                     scales = numeric(total.len),
                     low = numeric(total.len),
                     high = numeric(total.len),
                     dataset = 'XYZ', stringsAsFactors=FALSE)
    
    #put data into data frame
    t = 1
    for (i in 1:numObj){
      d = each.len[i]
      obj[t:(t+d-1),] = data.frame(WV = obj_list[[i]]$variance,
                                   scales = obj_list[[i]]$scales,
                                   low = obj_list[[i]]$ci_low,
                                   high = obj_list[[i]]$ci_high,
                                   dataset = legend.label[i], stringsAsFactors=FALSE)
      t = t +d
    }
    
    
    if (numObj == 2 ){
      if(is.null(line.color)){
        line.color = c("#003C7D","#F47F24")
      }
      if(is.null(CI.color)){
        CI.color = c("#003C7D","#F47F24")
      }
#       if(is.null(line.type)){
#         line.type = c('solid','dotted')
#       }
      
      autoplot.wvarComp(obj, split = split, CI = CI, background = background, transparence = transparence, line.color =line.color, 
                        CI.color = CI.color, line.type = line.type,  point.size = point.size, point.shape = point.shape,
                        title = title, title.size= title.size, 
                        axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                        axis.x.label = axis.x.label,
                        axis.y.label = axis.y.label,
                        facet.label.size = facet.label.size,
                        legend.label = legend.label,
                        legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                        legend.text.size = legend.text.size,
                        nrow = nrow)
    }
    else{
#       if(is.null(line.type)){
#         line.type = c('solid','dotted')
#       }
      autoplot.wvarComp(obj, split = split, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                        CI.color = CI.color, line.type = line.type,  point.size = point.size, point.shape = point.shape,
                        title = title, title.size= title.size, 
                        axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                        axis.x.label = axis.x.label,
                        axis.y.label = axis.y.label,
                        facet.label.size = facet.label.size,
                        legend.label = legend.label,
                        legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                        legend.text.size = legend.text.size,
                        nrow = nrow)
    }
    
  }
  
}