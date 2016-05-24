#.mergeLists----------------------------2012-12-20
#  Taken from R Curl - merge x and y, if elements
#  match in both x and y, then elements from y are
#  stored (overwriting x's elements).
#----------------------------------------------ACB
.mergeLists <- function( x, y ) 
{
    if (length(x) == 0) 
        return(y)
    if (length(y) == 0) 
        return(x)
    i = match(names(y), names(x))
    i = is.na(i)
    if (any(i)) 
        x[names(y)[which(i)]] = y[which(i)]
    x
}
#--------------------------------------.mergeLists


#PBSoptions-class-----------------------2012-12-20
setClass ("PBSoptions", representation( instance = "list" ) )

setMethod( f="initialize", signature="PBSoptions",
definition=function(.Object, filename, initial.options = list(), gui.prefix = "option" )
{
	if( missing( filename ) ) stop( "class initializer requires a filename" )
	if( is.character( filename ) == FALSE ) stop( "filename must be a character vector (string)" )
	if( is.list( initial.options ) == FALSE ) stop( "initial.options must be a list" )
	if( is.character( gui.prefix ) == FALSE ) stop( "GUI prefix must be a character vector (string)" )

	options <- initial.options

	#create functions within this namespace - so that each time a new instance of this class is created
	#these functions will point to a new set of variables (i.e. just how non-static variables of a c++ class work)
	load <- function( fname, prompt = FALSE )
	{
		old_filename <- filename
		if( missing( fname ) == FALSE )
			setFileName( fname )
		if( prompt == TRUE ) {
			f <- selectFile( initialfile=filename, mode="open" )
			if( is.null( f ) == TRUE ) {
				setFileName( old_filename ) #restore old filename if user cancels
				return()
			}
			setFileName( f )
		}
		if( file.exists( filename ) == FALSE )
			return()
		tmp <- readList( filename )
		eval(parse(text="options <<- .mergeLists( options, tmp )"))
	}

	setFileName <- function( name )
	{
		eval(parse(text="filename <<- name"))
	}

	getFileName <- function()
	{
		return( filename )
	}

	saveAs <- function() #prompts user to select a file name
	{
		selected <- selectFile( mode="save")
		if( length( selected ) == 0 || selected == "" )
			return()

		#TODO alert if file already exists - do we really want to overwrite?

		eval(parse(text="filename <<- selected"))
		save()
	}
	
	#silent save
	save <- function( fname, prompt = FALSE )
	{
		old_filename <- filename
		if( missing( fname ) == FALSE )
			setFileName( fname )
		if( prompt == TRUE ) {
			f <- selectFile( initialfile=filename, mode="save" )
			if( is.null( f ) == TRUE ) {
				setFileName( old_filename ) #restore old filename if user cancels
				return()
			}
			setFileName( f )
		}
		writeList( options, fname = filename )
	}

	#get( some_key ) or get() for all key/values
	get <- function( key )
	{
		if( missing( key ) )
			return( options )
		return( options[[ key ]] )
	}

	#set( key = value )
	set <- function( ... )
	{
		v <- list( ... )
                if( length( v ) == 1 && is.list( v[[1]] ) && length( names( v ) ) == 0
                   && length( names( v[[1]] ) ) > 0)
                        # assume that a single, unnamed list that contains a named list of
                        # length greater than 0 is actually a list of options
                        v <- v[[1]]
		if( length( v ) == 0 ) return()
		if( is.null( names( v ) ) || any( names( v ) == "" ) )
			stop( "values must be named" )

		for( i in 1:length( v ) ) {
			eval(parse(text="options[ names(v)[ i ] ] <<- list( v[[ i ]] )"))
		}
		
	}

	#save GUI values to this option list (R)
	saveGUI <- function()
	{
		values <- getWinVal()
		widgets <- names( values )
		opts <- paste( gui.prefix, names( get() ), sep="" )
		m <- match( widgets, opts )
		to_update = widgets[ !is.na( m ) ] #names of widgets which correspond to a value

		for( w in to_update ) {
			k <- substring( w, nchar(gui.prefix)+1 ) #remove prefix to get option key
			eval(parse(text="options[[ k ]] <<- values[[ w ]]"))
		}
			
	}

	#load GUI values from this option list (R)
	loadGUI <- function()
	{
		values <- getWinVal()
		widgets <- names( values )
		opts <- paste( gui.prefix, names( get() ), sep="" )
		m <- match( widgets, opts )
		to_update = widgets[ !is.na( m ) ] #names of widgets which correspond to a value

		opts <- list()
		for( w in to_update ) {
			k <- substring( w, nchar(gui.prefix)+1 ) #remove prefix to get option key
			opts[[ w ]] <- get( k )
		}
		setWinVal( opts )
	}

	#get the gui prefix
	getPrefix <- function() { return( gui.prefix ) }

	#set the gui prefix
	setPrefix <- function( prefix ) { eval(parse(text="gui.prefix <<- prefix")) }

	#load all functions into the instance list to return
	instance <- list()
	items <- ls() 
	for( i in items ) {
		v <- base::get( i )
		if( is.function( v ) == FALSE ) next
		instance[[ i ]] <- v
	}

	load()
	.Object@instance <- instance
	return( .Object )
}
)
#---------------------------------PBSoptions-class


#.showOptions---------------------------2012-12-20
.showOptions <- function( object )
{
	cat( "filename:", object@instance$getFileName(), "\n" )
	cat( "gui.prefix:", object@instance$getPrefix(), "\n" )
	opts <- object@instance$get()
	if( length( opts ) == 0 )
		cat( "Options: None\n" )
	else {
		cat( "Options:\n" )
		str( opts, no.list=T)
	}
}
#-------------------------------------.showOptions


#PBSotions-methods----------------------2012-12-20
setMethod( "print", signature="PBSoptions",
definition=function( x, ... )
{
	.showOptions( x )
}
)
setMethod( "show", signature="PBSoptions",
definition=function( object )
{
	.showOptions( object )
}
)
#--------------------------------PBSotions-methods


#getOptions-----------------------------2012-12-20
getOptions <-function( option.object, key )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$get( key )
}
#---------------------------------------getOptions


#getOptionsFileName---------------------2012-12-20
getOptionsFileName <-function( option.object )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$getFileName()
}
#-------------------------------getOptionsFileName


#getOptionsPrefix-----------------------2012-12-20
getOptionsPrefix <-function( option.object )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$getPrefix()
}
#---------------------------------getOptionsPrefix


#loadOptions----------------------------2012-12-20
loadOptions <-function( option.object, fname, prompt = FALSE )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$load( fname, prompt )
}
#--------------------------------------loadOptions


#loadOptionsGUI-------------------------2012-12-20
loadOptionsGUI <-function( option.object )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$loadGUI()
}
#-----------------------------------loadOptionsGUI


#saveOptions----------------------------2012-12-20
saveOptions <-function( option.object, fname, prompt = FALSE )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$save( fname, prompt )
}
#--------------------------------------saveOptions


#saveOptionsGUI-------------------------2012-12-20
saveOptionsGUI <-function( option.object )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$saveGUI()
}
#-----------------------------------saveOptionsGUI


#setOptions-----------------------------2012-12-20
setOptions <-function( option.object, ... )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$set( ... )
}
#---------------------------------------setOptions


#setOptionsFileName---------------------2012-12-20
setOptionsFileName <-function( option.object, name )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$setFileName( name )
}
#-------------------------------setOptionsFileName


#setOptionsPrefix-----------------------2012-12-20
setOptionsPrefix <-function( option.object, prefix )
{
	if( is( option.object, "PBSoptions" ) == FALSE ) stop( "option.object must be a pbsmodelling option class" )
	option.object@instance$setPrefix( prefix )
}
#---------------------------------setOptionsPrefix


#===== THE END ===================================

