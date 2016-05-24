#############################################################
#                                                           #
#   as.circular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 31, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-1                                           #
#############################################################

as.circular <- function (x, control.circular=list(), ...) {
   if (is.circular(x))
      return(x)
   else if (!is.null(xcircularp <- circularp(x)))
           circular(x, type=xcircularp$type, units=xcircularp$units, template=xcircularp$template, modulo=xcircularp$modulo, zero=xcircularp$zero, rotation=xcircularp$rotation)
         else {
           warntype <- warnunits <- warntemplate <- warnrotation <- warnmodulo <- warnzero <- ""
           printwarn <- FALSE
           dotc <- list(..., expand=TRUE)
           dc <- control.circular
           if (is.null(dc$type)) {
              if (!is.null(dotc$type)) 
                 dc$type <- dotc$type
              else {
                 dc$type <- "angles"
                 warntype <- "  type: 'angles'\n"
                 printwarn <- TRUE
              }
           }
           if (is.null(dc$units)) {
              if (!is.null(dotc$units))
                 dc$units <- dotc$units
              else {
                 dc$units <- "radians"
                 warnunits <- "  units: 'radians'\n"
                 printwarn <- TRUE
              }
           }
           if (is.null(dc$template)) {
              if (!is.null(dotc$template))
                 dc$template <- dotc$template
              else {
                 dc$template <- "none"
                 warntemplate <- "  template: 'none'\n"
                 printwarn <- TRUE
              }
           }
           if (is.null(dc$modulo)) {
              if (!is.null(dotc$modulo))
                 dc$modulo <- dotc$modulo
              else {
                 dc$modulo <- "asis"
                 warnmodulo <- "  modulo: 'asis'\n"
                 printwarn <- TRUE
              }
           }
           if (is.null(dc$zero)) {
              if (!is.null(dotc$zero))
                 dc$zero <- dotc$zero
              else {
                 dc$zero <- 0
                 warnzero <- "  zero: 0\n"
                 printwarn <- TRUE
              }
           }
           if (is.null(dc$rotation)) {
              if (!is.null(dotc$rotation))
                 dc$rotation <- dotc$rotation
              else {
                 dc$rotation <- "counter"
                 warnrotation <- "  rotation: 'counter'\n"
                 printwarn <- TRUE
              }
           }
           if (printwarn) {
               warn <- paste("an object is coerced to the class 'circular' using default value for the following components:\n", warntype, warnunits, warntemplate, warnmodulo, warnzero, warnrotation, sep="")
              warning(warn, sys.call(-1))
           }
           circular(x, type=dc$type, units=dc$units, template=dc$template, modulo=dc$modulo, zero=dc$zero, rotation=dc$rotation)
         }
}

