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

#' @title Create a GMWM TS Object based on data
#' @description Setups a time series oriented object that works well with graphing and summary utilities
#' @param data A one-column \code{matrix}, \code{data.frame}, or a numeric \code{vector}.
#' @param start A \code{numeric} that provides the time of the first observation.
#' @param end A \code{numeric} that provides the time of the last observation.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is \code{NULL}.
#' @return A \code{gts} object with the following attributes:
#' \describe{
#'   \item{start}{The time of the first observation}
#'   \item{end}{The time of the last observation}
#'   \item{freq}{Numeric representation of frequency}
#'   \item{unit}{String representation of the unit}
#'   \item{name}{Name of the dataset}
#' }
#' @author JJB, Wenchao
#' @examples
#' m = data.frame(rnorm(50))
#' x = gts(m, unit = 'sec', name = 'example')
#' x
#' plot(x)
#' 
#' x = gen.gts(WN(sigma2=1), 50)
#' x = gts(x, freq = 100, unit = 'sec')
#' plot(x)
gts = function(data, start = 0, end = NULL, freq = 1, unit = NULL, name = NULL){
  
  # 1. requirement for 'data'
  # Force data.frame to matrix  
  if (is.data.frame(data)){ 
    data = data.matrix(data)
  }
  
  # Check if the data is in matrix form
  if (is.matrix(data)) {
    # Check ncol
    ncolumn = ncol(data)
    if(ncolumn != 1){
      stop("'data' must have one column.")
    }
    
  } else {
    data = data.matrix(data) # convert vector to matrix
  }
  
  ndata = nrow(data)
  colnames(data) = if(is.null(name)) 'observed' else name 
  
  if(ndata == 0) {
    stop("Not a valid data object! Please supply a data set with one column that is in either a data.frame, matrix, or numeric object.")
  }
  
  # 2. requirement for 'freq'
  if(!is(freq,"numeric") || length(freq) != 1){ stop("'freq' must be one numeric number.") }
  if(freq <= 0) { stop("'freq' must be larger than 0.") }
  
  # 3. requirements for 'start' and 'end'
  if( is.numeric(start)==F && is.numeric(end)==F){
    stop("'start' or 'end' must be specified.")}
  
  if(is.null(start)==F && is.null(end)==F && (end-start)!= ((ndata-1)/freq) ){
    stop("end-start == (ndata-1)/freq must be TRUE.")
  }
  
  if ( is.null(end) ){
    end = start + (ndata - 1)/freq} # freq conversion (unit conversion is handled in graphical function)
  else if ( is.null(start) ){
    start = end - (ndata - 1)/freq}
  
  # 4. requirement for 'unit'
  if(!is.null(unit)){
    if(!unit %in% c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
    }
  }
  
  # x = 0:(ndata-1)
  # x = seq(from = 0, to = (ndata-1), length.out = ndata)
  # x = x/freq ###when generate the object, not deal with freq
  
  out = structure(data, 
                start = start, 
                end= end, # start and end will not be null now
                freq = freq,
                unit = unit,
                name = name, 
                class = c("gts","matrix"))
  
  out
}


#' @title Create a GMWM TS Object based on model
#' @description Create a \code{gts} object based on a supplied time series model.
#' @param model A \code{ts.model} or \code{gmwm} object containing one of the allowed models.
#' @param N An \code{interger} containing the amount of observations for the time series.
#' @param start A \code{numeric} that provides the time of the first observation.
#' @param end A \code{numeric} that provides the time of the last observation.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is \code{NULL}.
#' @return A \code{gts} object with the following attributes:
#' \describe{
#'   \item{start}{The time of the first observation}
#'   \item{end}{The time of the last observation}
#'   \item{freq}{Numeric representation of frequency}
#'   \item{unit}{String representation of the unit}
#'   \item{name}{Name of the dataset}
#' }
#' @details
#' This function accepts either a \code{ts.model} object (e.g. AR1(phi = .3, sigma2 =1) + WN(sigma2 = 1)) or a \code{gmwm} object.
#' @examples
#' # Set seed for reproducibility
#' set.seed(1336)
#' 
#' n = 1000
#' 
#' # AR1 + WN
#' model = AR1(phi = .5, sigma2 = .1) + WN(sigma2=1)
#' x = gen.gts(model, n)
#' x
#' plot(x)
#' 
#' set.seed(1336)
#' # GM + WN
#' # Convert from AR1 to GM values
#' m = ar1_to_gm(c(.5,.1),10)
#' # Beta = 6.9314718, Sigma2_gm = 0.1333333
#' model = GM(beta = m[1], sigma2_gm = m[2]) + WN(sigma2=1)
#' x2 = gen.gts(model, n, freq = 10, unit = 'sec')
#' x2
#' 
#' plot(x2, to.unit = 'min')
#' 
#' # Same time series
#' all.equal(x, x2, check.attributes = FALSE)
gen.gts = function(model, N = 1000, start = 0, end = NULL, freq = 1, unit = NULL, name = NULL){
  
  # 1. Do we have a valid model?
  if(!(is.ts.model(model) || is.gmwm(model))){
    stop("model must be created from a ts.model or gmwm object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  if(is.gmwm(model)){
    model = model$model.hat
  }
  
  # 2. requirement for 'freq'
  if(!is(freq,"numeric") || length(freq) != 1){ stop("'freq' must be one numeric number.") }
  if(freq <= 0) { stop("'freq' must be larger than 0.") }
  
  # 3. requirements for 'start' and 'end'
  if( is.numeric(start)==F && is.numeric(end)==F){
    stop("'start' or 'end' must be specified.")}
  
  if(is.null(start)==F && is.null(end)==F && (end-start)!= ((N-1)/freq) ){
    stop("end-start == (N-1)/freq must be TRUE.")
  }
  
  if ( is.null(end) ){
    end = start + (N - 1)/freq} # freq conversion (unit conversion is handled in graphical function)
  else if ( is.null(start) ){
    start = end - (N - 1)/freq}
  
  # 4. 'unit'
  if(!is.null(unit)){
    if(!unit %in% c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
    }
  }
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(!model$starting){
    
    theta = model$theta
    
    # Convert from AR1 to GM
    if(any(model$desc == "GM")){
      theta = conv.gm.to.ar1(theta, model$process.desc, freq)
    }
    
    out = .Call('gmwm_gen_model', PACKAGE = 'gmwm', N, theta, desc, obj)
  }else{
    stop("Need to supply initial values within the ts.model object.")
  }
  
  colnames(out) = if(is.null(name)) 'observed' else name 
  
  out = structure(.Data = out, 
                  start = start, 
                  end   = end, # start and end will not be null now
                  freq  = freq,
                  unit  = unit,
                  name  = name, 
                  class = c("gts","matrix"))
  
  out
  
}

#' @title Plot Time Series Data
#' @description This function is implemented with ggplot2.
#' @method plot gts
#' @export
#' @param x A \code{gts} object
#' @param to.unit A \code{string} indicating the unit which the data is converted to. The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size An \code{integer} that indicates the size of points on lines.
#' @param point.shape An \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of time series.
#' @author Wenchao
#' @examples
#' x = gen.gts(WN(sigma2=1), 50, unit = 'sec', freq = 1)
#' plot(x)
#' plot(x, to.unit = 'ns', line.type = 'dashed', line.color = 'black', point.size = 2)
plot.gts = function(x, to.unit = NULL, background = 'white',  
                    line.type = 'solid', line.color = '#003C7D',
                    point.size = 0, point.shape = 20,
                    title = NULL, title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, 
                    axis.x.label = NULL, axis.y.label = NULL, ... ){
  
  autoplot.gts(object = x, to.unit = to.unit, background = background,  
               line.type = line.type, line.color = line.color,
               point.size = point.size, point.shape = point.shape,
               title = title, title.size= title.size, 
               axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
               axis.x.label = axis.x.label, axis.y.label = axis.y.label )
  
}

#' @title Plot Time Series Data
#' @description This function is implemented with ggplot2.
#' @method autoplot gts
#' @export
#' @keywords internal
#' @param object A \code{gts} object
#' @param to.unit A \code{string} indicating the unit which the data is converted to. The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size An \code{integer} that indicates the size of points on lines.
#' @param point.shape An \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of time series.
#' @author Wenchao
#' @examples
#' x = gen.gts(WN(sigma2=1), 50, unit = 'sec', freq = 1)
#' autoplot(x)
#' autoplot(x, to.unit = 'ns', line.type = 'dashed', line.color = 'black', point.size = 2)
autoplot.gts = function(object, to.unit = NULL, background = 'white',  
                        line.type = 'solid', line.color = '#003C7D',
                        point.size = 0, point.shape = 20,
                        title = NULL, title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = NULL, axis.y.label = NULL, ... ){
  x1=y1=NULL
  from.unit = attr(object, 'unit')
  start =  attr(object, 'start')
  end = attr(object, 'end')
  ndata = nrow(object)
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  if(!is(object,"gts")){stop('object must be a gts object. Use function gts() or gen.gts() to create it.')}
  
  df = data.frame(y1 = object, x1 = seq( start, end, length.out = ndata) )
  colnames(df) = c('y1', 'x1')
  
  # NO unit conversion
  if( is.null(from.unit) && is.null(to.unit)==F ){
    warning('Unit of object is NULL. Unit conversion was not done.')
  }

  if (is.null(from.unit) == F){
    if (is.null(to.unit) == F){
      start_end = c(start, end)
      obj = unitConversion(start_end, from.unit, to.unit)
      
      if (obj$converted) {
        # YES unit conversion
        df$x1 = seq( obj$x[1], obj$x[2], length.out = ndata)
        message(paste0('Unit of object is converted from ', from.unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  p = ggplot(data = df, mapping = aes(x = x1,y = y1)) + geom_line(linetype = line.type, color = line.color) + 
    geom_point(color = line.color, size = point.size, shape = point.shape)
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  # axis label: default setting
  if(is.null(axis.x.label)){
    if(!is.null(from.unit) && !is.null(to.unit) && obj$converted){ 
      axis.x.label = paste('time(',to.unit,')', sep = '')
    }else if(is.null(from.unit) == F){
      axis.x.label =paste('time(',from.unit,')', sep = '')
    }else{
      axis.x.label = 'time'
    }
  }
  if(is.null(axis.y.label)){
    axis.y.label = 'observations'
  }
    
  p = p +  
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size))  
  p
}


#' @title Convert Unit of Time Series Data
#' @description Manipulate the units of time to different ones
#' @keywords internal
#' @param x A \code{vector} containing the values on x-axis.
#' @param from.unit A \code{string} indicating the unit which the data is converted from.
#' @param to.unit A \code{string} indicating the unit which the data is converted to.
#' @details
#' The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' Make sure \code{from.unit} and \code{to.unit} are not \code{NULL} before it is passed to this function.
#' @return A \code{list} with the following structure:
#' \describe{
#'  \item{x}{Data}
#'  \item{converted}{A \code{boolean} indicating whether conversion is made}
#' }
#' @examples
#' x = seq(60, 3600, 60)
#' unitConversion(x, 'sec', 'min')
#' y = 1:10
#' unitConversion(y, 'hour', 'sec')
unitConversion = function(x, from.unit, to.unit){
  
  #ns, ms, second, min, hour, day, month, year
  unit = c(ns = 1, ms = 2,se = 3, mi = 4, ho = 5, da = 6, mo = 7, ye = 8)
  
  #assume 1 month = 30 days
  ratio = c(1E6, 1E3, 60, 60, 24, 30, 12)
  from.unit.1 = substr(from.unit, 1, 2)
  to.unit.1 = substr(to.unit, 1, 2)
  
  #check unit:
  no.convert = F
  if(from.unit.1 == to.unit.1){no.convert = T}
  if(is.na(unit[from.unit.1]) ) {
    message = paste('No such unit: ', from.unit, '. Supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year". Conversion is terminated.', sep = '')
    warning(message); no.convert = T}
  if(is.na(unit[to.unit.1]) ) {
    message = paste('No such unit: ', to.unit, '. Supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year". Conversion is terminated.', sep = '')
    warning(message); no.convert = T}
  
  if(!no.convert){
    #print out warning when day is convert to month, or month is converted to day.
    conversionRange = unit[from.unit.1] : unit[to.unit.1]
    if(6 %in% conversionRange && 7 %in% conversionRange){
      warning('Unit conversion might be wrong because this function simply assumes 1 month = 30 days.')
    }
  }
  
  if(!no.convert){
    if(unit[from.unit.1] > unit[to.unit.1]){
      temp = ratio[unit[to.unit.1]: (unit[from.unit.1]-1)]
      multiplier = prod(temp)
      x = x*multiplier
    }else{
      temp = ratio[unit[from.unit.1]: (unit[to.unit.1]-1) ]
      multiplier = prod(temp)
      x = x/multiplier
    }
  }
  obj = list(x = x, converted = !no.convert)  
  return(obj)
}
