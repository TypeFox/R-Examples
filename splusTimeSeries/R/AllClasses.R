setClass( "seriesVirtual" )

setClass( "series",
	  representation( 
			  data = "ANY",
			  positions = "positions",
			  start.position = "positions",
			  end.position = "positions",
			  future.positions = "positions",
			  units = "character",
			  title = "character",
			  documentation = "character",
			  attributes = "ANY" ),
         contains = "seriesVirtual",
         validity = function(object)
         {
           TRUE
         })


setClass( "timeSeries",
	  representation( "fiscal.year.start" = "numeric",
			  "type" = "character" ),
         contains = c("series", "seriesVirtual"),
         prototype = prototype(
           data = matrix(nrow=0, ncol=0),
           positions = timeDate(),
           start.position = timeDate(),
           end.position = timeDate(),
           future.positions = timeDate(),
           fiscal.year.start = 1,
           type = "",
           attributes = NULL )
         )

setClass( "signalSeries",
         representation( units.position = "character" ),
         contains = c("series", "seriesVirtual"),
         prototype = prototype(
           data = matrix(nrow=0, ncol=0),
           positions = numericSequence(),
           start.position = numericSequence(),
           end.position = numericSequence(),
           future.positions = numericSequence(),
           units.position = ""))

