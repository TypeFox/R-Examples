setClass("groupVec",
         representation(columns="list", names="character", classes="character"),
         prototype=prototype(columns=list(),
           names=character(),
           classes=character()),
         validity = function(object)
         {
           lens = c(length(object@columns), length(names(object)),
             length(object@classes))
           if(any(lens != lens[1]))
             return("slots columns, names and classes should have same length")
           else if(lens[1]>0){
             clens = sapply(object@columns, length)
             if(any(clens != clens[1]))
               return("each element in columns should have the same length")
           }
           TRUE
         })

setClass( "positionsCalendar", representation("VIRTUAL" ))

setClass( "timeZone", representation("VIRTUAL"))
setClass( "timeZoneC",
    representation( name = "character"),
         contains="timeZone" ,
    prototype=list( name = "utc" ),
    validity = function( object )
	  {
	    if( length( object@name ) != 1 )
	      return( "Valid timeZoneC objects contain a single name" )
	    if( !any( object@name == 
		     c( "st/newzealand", "newzealand", "st/caroline",
		        "st/eaustralia", "aust/nsw", "aust/tasmania",
		        "aust/victoria", "st/caustralia", "aust/south",
		        "st/japan", "st/china", "aust/western", "hongkong",
		        "singapore", "st/saigon", "st/kazakh", "st/pakistan",
		        "st/caspian", "st/moscow", "st/eeurope", 
		        "europe/east", "st/ceurope", "europe/central",
		        "utc", "britain", "europe/west", "st/azores",
		        "st/oscar", "st/wgreenland", "can/newfoundland",
		        "st/atlantic", "can/atlantic", "st/eastern",
		        "us/eastern", "can/eastern", "st/central",
		        "us/central", "can/central", "st/mountain",
		        "us/mountain", "can/mountain", "st/pacific",
		        "us/pacific", "can/pacific", "st/alaska",
		        "us/alaska", "can/yukon", "st/hawaii",
		        "us/hawaii", "st/samoa" )))
	      return( "unknown time zone" )
	    TRUE
	  })

setClass( "timeZoneR",
         representation( offset = "integer", rules = "ANY" ),
         contains="timeZone",
         prototype = list( offset = 0L, rules = data.frame() ))

setClass( "timeDate",
         representation(format = "character",
                        time.zone = "character" ),
         contains=c("groupVec", "positionsCalendar"),
         prototype=prototype(
           columns=list(numeric(), numeric()),
           names=c( "julian.day", "milliseconds" ),
           classes=c( "integer", "integer" ),
	   format = "%02m/%02d/%04Y %02H:%02M:%02S.%03N",
	   time.zone = "GMT"),
         validity = function(object)
         {
           lens = c(length(object@columns), length(names(object)),
             length(object@classes))
           if(any(lens != 2))
             return("slots columns, names and classes should have length 2")
           else {
             clens = sapply(object@columns, length)
             if(any(clens != clens[1]))
               return("each element in columns should have the same length")
           }
           TRUE
         })
         
           
setClass( "timeInterval")
setClass( "timeSpan",
         representation( format = "character"),
         contains = c( "groupVec", "timeInterval"),
         prototype=prototype(
           columns=list(numeric(), numeric()),
           names=c( "julian.day", "milliseconds" ),
           classes=c( "integer", "integer" ),
	   format = "%dd %Hh %Mm %Ss %NMS"),
         validity = function(object)
         {
           lens = c(length(object@columns), length(names(object)),
             length(object@classes))
           if(any(lens != 2))
             return("slots columns, names and classes should have length 2")
           else {
             clens = sapply(object@columns, length)
             if(any(clens != clens[1]))
               return("each element in columns should have the same length")
           }
           TRUE
         })
	

setClass( "timeRelative", 
         representation( Data="character", holidays="positionsCalendar" ),
         contains="timeInterval",
         prototype = prototype( list(Data = character(0),
           holidays = new("timeDate"))))

setClass( "timeEvent", 
         contains = "groupVec",
         prototype=prototype(
           names = c( "start", "end", "IDs" ),
           classes =  c( "positionsCalendar", "positionsCalendar", "ANY" ),
           columns = list( new("timeDate"), new("timeDate"), character())))

setClass( "timeSequence", 
	  representation(                     
			  from = "timeDate",
			  to = "timeDate",
			  by = "timeInterval",
			  length = "integer",
			  exceptions = "timeEvent",
			  additions = "positionsCalendar",
			  format = "character",
			  time.zone = "character"
			 ),
         contains="positionsCalendar",
         prototype = prototype( 
	      from = new("timeDate", columns=list(0,0)),
          to = new("timeDate", columns=list(0,0)),
	      by = new("timeSpan"),
	      length = 0L,
	      exceptions = new("timeEvent"),
	      additions = new("timeDate"),
	      format = "%02m/%02d/%04Y %02H:%02M:%02S.%03N",
	      time.zone = "GMT"))

setClass( "numericSequence", 
	 representation( 
			 from = "numeric",
			 to = "numeric", 
			 by = "numeric",
			 length = "integer" ),
         prototype = prototype(
           from = 0, to = 0, by = numeric( 0 ), length = 0L ))

setClassUnion("positionsNumeric", c("numericSequence", "numeric"))
setClassUnion( "positions" , c("positionsCalendar", "positionsNumeric"))

