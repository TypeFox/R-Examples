#' @title Update Object Attribute
#' @description Update the attributes of \code{lts}, \code{gts} and \code{imu} object
#' @param object A \code{lts}, \code{gts} or \code{imu} object
#' @param type A \code{string} that contains the attribute to be updated
#' @param new The updated value for the attribute
#' @param keep.start A \code{boolean} value that indicates whether 'start' or 'end' should remain the same when 'freq' is updated
#' @param ... Further arguments passed to or from other methods.
#' @return An object with the updated attribute.
#' @export
#' @details
#' This function is able to update some attributes for \code{gts}, \code{lts} and \code{imu} objects. 
#' For \code{lts} object, the attributes that can be updated are 'start', 'end', 'freq', 'unit', 'name' and 'process'.
#' For \code{gts} object, the attributes that can be updated are 'start', 'end', 'freq', 'unit' and 'name'.
#' For \code{imu} object, the attributes that can be updated are 'axis', 'freq', 'unit' and 'name'.
#' 
#' If one between 'start' and 'end' is updated, the other one will also be updated, since \code{end-start == (N-1)/freq} must be TRUE, where \code{N}
#' is the number of observations in the object. 
#' 
#' If 'freq' is updated, by default 'start' will remain the same, and 'end' will be updated at the same time,
#' unless you set 'keep.start = F'.
#' 
#' If 'unit' is updated, the old unit will be replaced by the new one, and other attributes will remain the same.
#' It is different from the unit conversion feature.
#' 
#' @examples
#' gts1 = gts(rnorm(50), freq = 1, unit = 'sec', name = 'test1')
#' gts2 = update(gts1, 'unit', 'min')
#' attr(gts2, 'unit')
#' 
#' gts3 = update(gts1, 'name', 'test2')
#' attr(gts3, 'name')
update.lts = function(object, type, new, keep.start = T, ...){
  if(! (type%in%c('start', 'end', 'freq', 'unit', 'name', 'process')) ){
    stop("For lts object, you can only update 'start', 'end', 'freq', 'unit', 'name' and 'process'.")
  }
  update_obj(object, type, new, keep.start)
}

#' @rdname update.lts
#' @export
update.gts = function(object, type, new, keep.start = T, ...){
  if(! (type%in%c('start', 'end', 'freq', 'unit', 'name')) ){
    stop("For gts object, you can only update 'start', 'end', 'freq', 'unit' and 'name'.")
  }
  update_obj(object, type, new, keep.start)
}

#' @rdname update.lts
#' @export
update.imu = function(object, type, new, ...){
  if(! (type%in%c('axis', 'freq', 'unit', 'name')) ){
    stop("For imu object, you can only update 'axis', 'freq', 'unit', 'name'.")
  }
  update_obj(object, type, new)
}

#' @title Update the Attributes of Objects
#' @description  Internal Function to Update the Attributes of Objects
#' @param object A \code{lts}, \code{gts} or \code{imu} object
#' @param type A \code{string} that contains the attribute to be updated
#' @param new The updated value for the attribute
#' @param keep.start A \code{boolean} value that indicates whether 'start' or 'end' should remain the same when 'freq' is updated
#' @return An object with the updated attribute.
#' @keywords internal
update_obj = function(object, type, new, keep.start = T){
  
  if(type == 'start'){
    freq = attr(object, 'freq')
    N = nrow(object)

    attr(object, 'start') = new
    attr(object, 'end') = new + (N - 1)/freq
  
  }else if(type == 'end'){
    freq = attr(object, 'freq')
    N = nrow(object)
    
    attr(object, 'end') = new
    attr(object, 'start') = new - (N - 1)/freq
    
  }else if(type == 'freq'){
    if(!is(new,"numeric") || length(new) != 1){ stop("Frequency must be one numeric number.") }
    if(new <= 0) { stop("Frequency must be larger than 0.") }
    
    if(keep.start){
      start = attr(object, 'start')
      attr(object, 'freq') = new
      N = nrow(object)
      attr(object, 'end') = start + (N - 1)/new
    }else{
      end = attr(object, 'start')
      attr(object, 'freq') = new
      N = nrow(object)
      attr(object, 'start') = end - (N - 1)/new
    }
    
  }else if(type == 'unit'){
    if(!(new %in% c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year'))){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
    }
    
    attr(object, 'unit') = new
    
  } else if(type == 'name'){
    attr(object, 'name') = new
    
  }else if(type == 'process'){
    old.process = attr(object, 'process')
    if (length(old.process)!= length(new)){
      stop(paste0("'process' must be a vector that contains ", length(old.process), " elements.") )
    }
    attr(object, 'process') = new
    
  }else if(type == 'axis'){
    old.axis = attr(object, 'axis')
    if (length(old.axis)!= length(new)){
      stop(paste0("'axis' must be a vector that contains ", length(old.axis), " elements.") )
    }
    attr(object, 'axis') = new
    
  }else{
    stop("'type' is not recognized.")
  }
  
  return(object)
}
