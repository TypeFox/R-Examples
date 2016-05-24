# generic pour toutes les classes Time*DataFrame
#-----------------------------------------------
# S4
setGeneric (name='nrow')
setGeneric (name='ncol')
setGeneric (name='lapply')

# fonctions génériques (time properties)
#---------------------------------------
# TODO : dispatch in the classe definition files ?
unit <- function(x, ...) UseMethod('unit')
'unit<-' <- function(object, value) UseMethod('unit<-')
setGeneric ('unit', function (x, ...) standardGeneric ('unit') )
setGeneric ('unit<-', function (object, value) standardGeneric ('unit<-') )

duration <- function(x, ...) UseMethod('duration')
setGeneric ('duration', function (x, ...) standardGeneric ('duration') )

timezone <- function (object) UseMethod('timezone')
'timezone<-' <- function(object, value) UseMethod('timezone<-')
setGeneric (name='timezone',
	    def=function (object) standardGeneric ('timezone') )
setGeneric (name='timezone<-',
	    def=function (object, value) standardGeneric ('timezone<-') )

of <- function(x, ...) UseMethod('of')

continuous <- function(x, ...) UseMethod('continuous')
setGeneric (name='continuous',
	    def=function(x, ...) standardGeneric('continuous'))
'continuous<-' <- function(x, value) UseMethod('continuous<-')
setGeneric (name='continuous<-',
	    def=function(x, value) standardGeneric('continuous<-'))

'homogeneous' <- function(x, ...) UseMethod('homogeneous')
setGeneric (name='homogeneous',
	    def=function(x, ...) standardGeneric('homogeneous'))

'overlapping' <- function(x, ...) UseMethod('overlapping')
setGeneric (name='overlapping',
	    def=function(x, ...) standardGeneric('overlapping'))

'period' <- function(x, ...) UseMethod('period')
setGeneric (name='period', def=function(x, ...) standardGeneric('period'))

'when' <- function(x, ...) UseMethod('when')
setGeneric (name='when', def=function(x, ...) standardGeneric('when'))

'regular' <- function(x, ...) UseMethod('regular')
setGeneric (name='regular', def=function(x, ...) standardGeneric('regular'))

'interval' <- function(x, ...) UseMethod('interval')
setGeneric (name='interval', def=function(x, ...) standardGeneric('interval'))

# generic pour interval (S3)
#---------------------------
'%intersect%' <- function(i1, i2) UseMethod ('%intersect%')

'%included%' <- function(i1, i2) UseMethod ('%included%')

# generic de coercition
#----------------------
as.POSIXcti <- function(from, ...) UseMethod ('as.POSIXcti')
as.POSIXctp <- function(from, ...) UseMethod ('as.POSIXctp')
as.TimeIntervalDataFrame <- function(from, ...) UseMethod ('as.TimeIntervalDataFrame')
as.TimeInstantDataFrame  <- function(from, ...) UseMethod ('as.TimeInstantDataFrame')
as.SubtimeDataFrame  <- function(x, unit, of, ...)
	UseMethod ('as.SubtimeDataFrame')

