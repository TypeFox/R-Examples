
timeZoneC <- function( name )
{
  if( missing( name ) || !length(name)) name <- "utc"
  if(length(name) > 1 )
    stop( "name can only contain one time zone name")

  new( "timeZoneC", name = name )
}

setAs("character", "timeZoneC", function(from) timeZoneC(from))

.time.zone.list <- list(
   Atlantic = timeZoneC( "can/atlantic" ),
   ADT = timeZoneC( "can/atlantic" ),
   AST = timeZoneC( "can/atlantic" ),
   Halifax = timeZoneC( "can/atlantic" ),
   PuertoRico = timeZoneC( "st/atlantic" ),

   Eastern = timeZoneC( "us/eastern" ),
   EST = timeZoneC( "us/eastern" ),
   EDT = timeZoneC( "us/eastern" ),
   EST5EDT = timeZoneC( "us/eastern" ),
   "EST/EDT" = timeZoneC( "us/eastern" ),
   Indiana = timeZoneC( "st/eastern" ),
   Toronto = timeZoneC( "can/eastern" ),

   Central = timeZoneC( "us/central" ),
   CST = timeZoneC( "us/central" ),
   CDT = timeZoneC( "us/central" ),
   CST6CDT = timeZoneC( "us/central" ),
   "CST/CDT" = timeZoneC( "us/central" ),
   Chicago = timeZoneC( "us/central" ),
   Winnipeg = timeZoneC( "can/central" ),
   
   Mountain = timeZoneC( "us/mountain" ),
   MST = timeZoneC( "us/mountain" ),
   MDT = timeZoneC( "us/mountain" ),
   MST7MDT = timeZoneC( "us/mountain" ),
   "MST/MDT" = timeZoneC( "us/mountain" ),
   Arizona = timeZoneC( "st/mountain" ),
   Denver = timeZoneC( "us/mountain" ),
   Edmonton = timeZoneC( "can/mountain" ),

   Pacific = timeZoneC( "us/pacific" ),
   PST = timeZoneC( "us/pacific" ),
   PDT = timeZoneC( "us/pacific" ),
   PST8PDT = timeZoneC( "us/pacific" ),
   "PST/PDT" = timeZoneC( "us/pacific" ),
   Vancouver = timeZoneC( "can/pacific" ),

   Alaska = timeZoneC( "us/alaska" ),
   AKST = timeZoneC( "us/alaska" ),
   AKDT = timeZoneC( "us/alaska" ),
   "AKST/AKDT" = timeZoneC( "us/alaska" ),

   Aleutian = timeZoneC( "us/hawaii" ),
   HST = timeZoneC( "st/hawaii" ),
   Hawaii = timeZoneC( "st/hawaii" ),

   Midway = timeZoneC( "st/samoa" ),
   Samoa = timeZoneC( "st/samoa" ),
   SST = timeZoneC( "st/samoa" ),

   Japan = timeZoneC( "st/japan" ),
   Tokyo = timeZoneC( "st/japan" ),
   JST = timeZoneC( "st/japan" ),

   China = timeZoneC( "st/china" ),
   HongKong = timeZoneC( "hongkong" ),
   Singapore = timeZoneC( "singapore" ),

   Sydney = timeZoneC( "aust/nsw" ),
   Hobart = timeZoneC( "aust/tasmania" ),
   Melbourne = timeZoneC( "aust/victoria" ),
   Adelaide = timeZoneC( "aust/south" ),
   Darwin = timeZoneC( "st/caustralia" ),
   Perth = timeZoneC( "aust/western" ),

   Marshall=timeZoneC("st/newzealand"),
   Wake=timeZoneC("st/newzealand"),
   IDLE=timeZoneC("st/newzealand"),

   Chamorro=timeZoneC("st/eaustralia"),
   ChST=timeZoneC("st/eaustralia"),
   Yap=timeZoneC("st/eaustralia"),
   YAPT=timeZoneC("st/eaustralia"),

   Caroline=timeZoneC("st/caroline"),
   LINT=timeZoneC("st/caroline"),

   Auckland = timeZoneC( "newzealand" ),
   NZST = timeZoneC( "newzealand" ),
   NZDT = timeZoneC( "newzealand" ),

   UTC = timeZoneC( "utc" ),
   GMT = timeZoneC( "utc" ),
   London = timeZoneC( "britain" ),
   GDT = timeZoneC( "britain" ),
   BST = timeZoneC( "britain" ),
   WET = timeZoneC( "europe/west" ),
   Wes = timeZoneC( "europe/west" ),
   WEST = timeZoneC( "europe/west" ),
   "WET/WEST" = timeZoneC( "europe/west" ),
   WED = timeZoneC( "europe/west" ),
   WEDT = timeZoneC( "europe/west" ),

   CET = timeZoneC( "europe/central" ),
   CEST = timeZoneC( "europe/central" ),
   MET = timeZoneC( "europe/central" ),
   MEST = timeZoneC( "europe/central" ),
   "MET/MEST" = timeZoneC( "europe/central" ),

   EET = timeZoneC( "europe/east" ),
   EEST = timeZoneC( "europe/east" ),
   "EET/EEST" = timeZoneC( "europe/east" )
)

timeZoneList <- function( ... )
{
  # Function for examining and changing the time zone list
  #
  # No arguments: return the current time zone list
  # List as argument: named components put into time zone list, and 
  #   previous value of the list is returned
  #   Note that each list component *must* be a time zone object.
  # One or more named arguments: works like list as argument
  #
  # The time zone list is stored as .time.zone.list in splusTimeDate
  # package

  current <- if(exists(".time.zone.list", envir=.splusTimeDateEnv)) {
    .splusTimeDateEnv$.time.zone.list
  } else {
    .time.zone.list
  }

  if( nargs() == 0 )
    return( current )
  
  # make a list out of the arguments
  arglist <- list( ... )

  # convert back if it was a list before
  if(( length(arglist) == 1 ) && !is( arglist[[1]], "timeZone" ))
    arglist <- arglist[[1]]

  # special case user passed in an empty list
  if( length( arglist ) == 0 )
    return(current)

  nam <- names( arglist )
  if( is.null( nam ))
    stop( "Time zones must have names" )
  if( !all( sapply( arglist, "is", "timeZone" )))
    stop( "All arguments or list components must be time zones" )


  oldzones <- current
  current[ nam ] <- arglist

  assign(".time.zone.list", current, envir=.splusTimeDateEnv)
  
  invisible(oldzones)
}


setMethod( "show", "timeZoneC",
function( object ) cat( "timeZoneC(\"", object@name, "\")\n",
		        sep = "") 
	  )

setMethod( "summary", "timeZoneC", 
	   function( object, ... ) 
	   {
	     ret <- object@name
	     names( ret ) <- "name"
	     ret <- matrix( ret, nrow = 1, dimnames = list("", "name" ))
	     oldClass( ret ) <- "table"
	     ret
	   })


timeZoneR <- function( offset=0, yearfrom=integer(0), yearto=integer(0), 
		       hasdaylight=logical(0), dsextra=integer(0),
		       monthstart=integer(0), codestart=integer(0), 
		       daystart=integer(0), xdaystart=integer(0), 
		       timestart=integer(0), monthend=integer(0), 
		       codeend=integer(0), dayend=integer(0), 
		       xdayend=integer(0), timeend=integer(0),
		       rules)
{

  if( length( offset ) != 1 )
    stop( "Offset must have length 1" )

  if( missing( rules ))
  {

    lenstend <- c(length( monthstart ), length( codestart ),
		  length( daystart ), length(xdaystart),
		  length( timestart ), 
		  length( monthend ), length( codeend ),
		  length( dayend ), length(xdayend),
		  length( timeend ), 
		  length( yearfrom ), length( yearto ))

    howmany <- length( hasdaylight )

    # see if the daylight savings part of rule is actually being used
    # and validate it if so
    if( howmany && any( hasdaylight ))
    {
      nodaylight <- FALSE

      if( any( lenstend != howmany ) || ( length( dsextra ) != howmany ))
	stop( "Rule inputs must have same length" )

      rules <- data.frame(yearfrom = as( yearfrom, "integer" ),
			  yearto = as( yearto, "integer" ),
			  hasdaylight = as( hasdaylight, "integer" ),
			  dsextra = as( dsextra, "integer" ),
			  monthstart = as( monthstart, "integer" ),
			  codestart = as( codestart, "integer" ),
			  daystart = as( daystart, "integer" ),
			  xdaystart = as( xdaystart, "integer" ),
			  timestart = as( timestart, "integer" ),
			  monthend = as( monthend, "integer" ),
			  codeend = as( codeend, "integer" ),
			  dayend = as( dayend, "integer" ),
			  xdayend = as( xdayend, "integer" ),
			  timeend = as( timeend, "integer" ))

    } else rules <- data.frame()
  } else
  {
    # rules supplied -- make sure it's a data frame
    if( !is( rules, "data.frame" ))
      rules <- as.data.frame( rules )
    if( !dim( rules )[[1]])
      rules <- data.frame()
    else if( !length( colnames( rules )))
      colnames(rules) <- c("yearfrom", "yearto", "hasdaylight",
			 "dsextra", "monthstart", "codestart",
			 "daystart", "xdaystart", "timestart",
			 "monthend", "codeend", "dayend", "xdayend",
			 "timeend")
  }

  # validate a bit more

  if( dim( rules )[[1]])
  {
    codes <- c( rules$codestart, rules$codeend )
    if( any(is.na( codes )) || max(codes) > 4 || min(codes) < 1 )
      stop( "Invalid codes" )

    if( any(( rules$yearto < rules$yearfrom ) & ( rules$yearto != -1 )))
      stop( "yearfrom years must be on or after yearto years" )
  }

  ret <- new("timeZoneR")
  ret@offset <- as( offset, "integer" )
  ret@rules <- rules
  ret  
}

setMethod( "show", "timeZoneR",
function( object )
{
  cat( "offset:", object@offset[1], "\n" )
  cat( "rules:\n" )
  show( object@rules )
})

setMethod( "summary", "timeZoneR", 
	   function( object, ... ) 
	   {
	     if( !is( object@rules, "data.frame" ))
	       stop( "invalid time zone object -- rules must be data frame")
	     ret <- c( object@offset, dim( object@rules )[[1]])
	     ret <- matrix( ret, nrow = 1, 
			   dimnames = list("", c( "offset", "rules")))
	     oldClass( ret ) <- "table"
	     ret
	   })

