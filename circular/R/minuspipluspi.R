#############################################################
#                                                           
#   MinusPiPlusPiRad function                                  
#   Author: Claudio Agostinelli                             
#   E-mail: claudio@unive.it                                
#   Date: October, 14, 2007                                  
#   Version: 0.1                                          
#                                                           
#   Copyright (C) 2007 Claudio Agostinelli                  
#                                                           
#############################################################

minusPiPlusPi <- function(x)
{
  if (is.circular(x)) {
    datacircularp <- circularp(x)
  } else {
    datacircularp <- list(type = "angles", units = "radians", template = "none",
                          modulo = "asis", zero = 0, rotation = "counter")
  }
  dc <- list()
  if (is.null(dc$type)) 
    dc$type <- datacircularp$type
  if (is.null(dc$units)) 
    dc$units <- datacircularp$units
  if (is.null(dc$template)) 
    dc$template <- datacircularp$template
  if (is.null(dc$modulo)) 
    dc$modulo <- datacircularp$modulo
  if (is.null(dc$zero)) 
    dc$zero <- datacircularp$zero
  if (is.null(dc$rotation)) 
    dc$rotation <- datacircularp$rotation
  
  	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  	attr(x, "class") <- attr(x, "circularp") <- NULL
	x <- MinusPiPlusPiRad(x)
	x <- conversion.circular(circular(x), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)

	return(x)
}


MinusPiPlusPiRad  <- function(x)
{
	x <- .C("MinusPiPlusPiRad",x=as.double(x),n=as.integer(length(x)))$x
 	return(x) 
} 

