#=================================================
# Set classes used by presenttalk
# -------------------------------
setClass( "text", representation( text = "character", "break" = "logical" ) )
setClass( "file", representation( name = "character", filename = "character", "break" = "logical", button = "logical", col = "integer" ) )
setClass( "code", representation( show = "logical", print = "logical", code = "character", "break" = "character", eval = "logical" ) )
setClass( "break", representation( "NULL" ) ) #this prints a message - wtf R?!!! give me an empty class that isn't virtual

setClass( "section", representation( name = "character", items = "list", button = "logical", col = "integer", section_id = "integer" ) ) #items should be a list of the above 4 s4 classes
setClass( "talk", representation( name = "character", sections = "list", files = "list" ) )


setValidity ("code", function( object )
{
	valid_breaks <- c( "show", "print", "all", "none" )
	if( any( object@"break" == valid_breaks ) == FALSE )
		return( paste( "code's break arugment must be one of: ", paste( valid_breaks, collapse=" " ), sep="" ) )
	return( TRUE )
}
)
#=================================================


#.parseTalkFile-------------------------2012-12-20
#  Returns a talk s4 class - given a filename to an xml talk
#----------------------------------------------ACB
.parseTalkFile <- function( talk_fname )
{
	.processSection <- function( node )
	{
		x <- new( "section", 
				name = node$attributes[ "name" ],
				button = as.logical( xmlGetAttr( node, "button", FALSE ) ),
				col = as.integer( xmlGetAttr( node, "col", 2 ) )
				 )
		for( i in xmlChildren( node ) ) {
			
			#TODO refactor with a do.call
			if( xmlName( i ) == "text" ) item <- .processText( i )
			else if( xmlName( i ) == "file" ) item <- .processFile( i )
			else if( xmlName( i ) == "code" ) item <- .processCode( i )
			#else if( xmlName( i ) == "break" ) item <- processBreak( i )
			else if( xmlName( i ) == "comment" ) { next } #do nothing
			else stop( paste( "Error unhandled xml tag found:", xmlName(i) ) )
			
			x@items[[ length( x@items ) + 1 ]] <- item
		}

		if( length( x@items ) == 0 ) {
			stop( paste( "Error parsing ", talk_fname, ": section \"", x@name, "\" does not contain any children (text,file,code) items", sep="" ) )
		}
	
		#scan through items, and place breaks accordingly (it would be easier to just have a <break/>)
		items <- x@items
		x@items <- list()
		i <- 1
		for( item in items ) {
			if( inherits( item, "text" ) || inherits( item, "file" ) ) {
				#insert item
				x@items[[ i ]] <- item
				i <- i + 1
	
				#insert break
				if( item@"break" == TRUE ) {
					x@items[[ i ]] <- new( "break" )
					i <- i + 1
				}
			} else if( inherits( item, "code" ) ) {
				#just print code first
				item@eval = FALSE
				x@items[[ i ]] <- item
				i <- i + 1
	
				#break
				if( any( item@"break" == c( "show", "all" ) ) ) {
					x@items[[ i ]] <- new( "break" )
					i <- i + 1
				}
	
				#just eval code
				item@eval = TRUE
				x@items[[ i ]] <- item
				i <- i + 1
	
				#break
				if( any( item@"break" == c( "print", "all" ) ) ) {
					x@items[[ i ]] <- new( "break" )
					i <- i + 1
				}
	
			}
		}
	
		return( x )
	}
	
	.processText <- function( node )
	{
		return( new( "text", text = xmlValue( node ), "break" = as.logical( xmlGetAttr( node, "break", TRUE ) ) ) )
	}
	
	.processFile <- function( node )
	{
		return( new( "file", 
			name = xmlGetAttr( node, "name", "" ), 
			"break" = as.logical( xmlGetAttr( node, "break", TRUE ) ),
			filename = xmlValue( node ),
			button = as.logical( xmlGetAttr( node, "button", FALSE ) ),
			col = as.integer( xmlGetAttr( node, "col", 3 ) )
			 ) )
	}
	
	.processCode <- function( node )
	{
		return( new( "code", 
				show = as.logical( xmlGetAttr( node, "show", TRUE ) ), 
				print = as.logical( xmlGetAttr( node, "print", TRUE ) ), 
				code = xmlValue( node ), 
				"break" = tolower( xmlGetAttr( node, "break", "print" ) )
				) )
	}

	#start parsing below

	talk_xml <- xmlTreeParse( talk_fname )
	stopifnot( !is.null( talk_xml$doc$children$talk ) )

	#create root element
	talk_node <- talk_xml$doc$children$talk
	name <- xmlGetAttr( talk_node, "name" )
	stopifnot( !is.null( name ) )
	talk <- new( "talk", name = name )

	#parse child objects
	for( i in xmlChildren( talk_node ) ) {
		if( xmlName( i ) == "comment" ) {
			next
		} else if( xmlName( i ) == "file" ) {
			talk@files[[ length( talk@files ) + 1 ]] <- .processFile( i )
		} else if( xmlName( i ) == "section" ) {
			section <- .processSection( i )
			talk@sections[[ length( talk@sections ) + 1 ]] <- section
		} else {
			stop( paste( "unhandled xml tag:", xmlName( i ) ) )
		}
	}

	return( talk )
}
#-----------------------------------.parseTalkFile


#.getSectionNames-----------------------2012-12-20
#given a talk, return a vector of all section names
#----------------------------------------------ACB
.getSectionNames <- function( talk )
{
	stopifnot( any( class( talk ) == "talk" ) )
	section_names <- character( length( talk@sections ) )
	i <- 1
	for( s in talk@sections ) {
		section_names[ i ] <- s@name
		i <- i + 1
	}
	return( section_names )
}
#---------------------------------.getSectionNames


#.getTalkIndexes------------------------2012-12-20
#retuns a list of 2 element vectors (i,j) where i is the section index, and j is the items index
#each element of the list corresponds to a break point
#----------------------------------------------ACB
.getTalkIndexes <- function( talk )
{
	stopifnot( !missing( talk ) )
	stopifnot( inherits( talk, "talk" ) )
	i <- 1
	breaks <- list()
	for( section in talk@sections ) {
		j <- 1
		start <- TRUE
		for( item in section@items ) {
			is_break <- inherits( item, "break" )
			if( start && is_break == FALSE )
				breaks[[ length( breaks ) + 1 ]] <- c( i, j )
			start <- is_break
			j <- j + 1
		}
		i <- i + 1
	}
	return( breaks )
}
#----------------------------------.getTalkIndexes


#.getIndexForSection--------------------2012-12-20
#  Get the slide index which corresponds to the 
#  first slide in a given section.
#----------------------------------------------ACB
.getIndexForSection <- function( talk, section_id )
{
	indices <- .getTalkIndexes( talk )
	i <- 1
	for( index in indices ) {
		if( section_id == index[1] )
			break
		i <- i + 1
	}
	return( i )
}
#------------------------------.getIndexForSection


#.getButton-----------------------------2012-12-20
#  Returns win description for a button for a file or section.
#----------------------------------------------ACB
.getButton <- function( talk_name, obj )
{
	if( inherits( obj, "section" ) ) {
		b <- paste( "button text=\"", obj@name, "\" function=.setsection action=\"", talk_name, ":", obj@section_id,"\" width=10 padx=4 pady=4 fg=red3 bg=whitesmoke",sep="")
		return( b )
	}
	if( inherits( obj, "file" ) ) {
		b <- paste( "button text=\"", obj@name, "\" function=.presentTalkOpenFile action=\"", obj@filename,"\" width=10 padx=4 pady=4 fg=blue bg=whitesmoke",sep="")
		return( b )
	}
}
#---------------------------------------.getButton


#.getButtons----------------------------2012-12-20
#  Gets widget descriptions for file and section buttons.
#----------------------------------------------ACB
.getButtons <- function( talk )
{
	#create a list of buttons
	section_but <- list()
	file_but <- list()

	#process top level buttons (under talk)
	for( f in talk@files ) {
		i <- length( file_but ) + 1
		file_but[[ i ]] <- f
	}

	#process buttons under sections
	sect_id <- 1
	for( s in talk@sections ) {
		if( s@button == TRUE ) {
			i <- length( section_but ) + 1
			s@section_id <- as.integer( sect_id )
			section_but[[ i ]] <- s
		}
		for( item in s@items ) {
			if( inherits( item, "file" ) && item@button == TRUE ) {
				i <- length( file_but ) + 1
				file_but[[ i ]] <- item
			}
		}
		sect_id <- sect_id + 1
	}

	#create win desc corresponding to list of buttons

	n_files <- length( file_but ) #files on the left
	n_sections <- length( section_but ) #sections on the right
	cols <- 4 #number of cols to use
	min_rows <- 9999
	min_row_i <- 0
	for( i in 1:cols ) {
		#determine optimal column assignment
		rows <- ceiling( max( n_files / i, n_sections / ( cols - i ) ) )
		if( i == 0 || rows <= min_rows ) {
			min_rows <- rows
			min_row_i <- i
		}
	}
	
	stopifnot( min_rows > 0 )
	w <- paste( "grid ", min_rows, " ", cols + 1, " relief=groove sticky=NEW", sep="" )
	file_i <- 1
	section_i <- 1
	#iterate over each (i,j) position in button grid, insert file buttons on the left,
	#then a center null widget, then section buttons on the right
	for( j in 1:min_rows ) {
		for( i in 1:(cols+1) ) {
			if( i <= min_row_i ) {
				if( file_i > length( file_but ) ) {
					w <- append( w, "null padx=34" )
				} else {
					w <- append( w, .getButton( talk@name, file_but[[ file_i ]] ) )
					file_i <- file_i + 1
				}
			} else if( i == ( min_row_i + 1 ) ) {
				w <- append( w, "null padx=5" )
			} else {
				if( section_i > length( section_but ) ) {
					w <- append( w, "null padx=20" )
				} else {
					w <- append( w, .getButton( talk@name, section_but[[ section_i ]] ) )
					section_i <- section_i + 1
				}
			}
		}
	}
	return( w )
}
#--------------------------------------.getButtons


#.getMenus------------------------------2012-12-20
#  Get widget description for menus.
#----------------------------------------------ACB
.getMenus <- function( talk )
{
	stopifnot( inherits( talk, "talk" ) )
	sections <- c()
	files <- c()
	sect_id <- 1
	for( s in talk@sections ) {
		#save section under menu
		i <- length( sections ) + 1
		sections[ i ] <- paste( "menuitem label=\"", s@name, "\" function=.setsection action=\"", talk@name, ":", sect_id, "\"", sep="" )

		#look for files
		for( item in s@items ) {
			if( inherits( item, "file" ) ) {
				i <- length( files ) + 1
				files[ i ] <- paste( "menuitem label=\"", item@name, "\" function=.presentTalkOpenFile action=\"", item@filename, "\"", sep="" )
			}
		}
		sect_id <- sect_id + 1
	}

	if( length( sections ) > 0 ) {
		w <- paste( "menu nitems=", length( sections ), " label=Sections", sep="" )
		w <- append( w, sections )
	}

	if( length( files ) > 0 ) {
		w <- append( w, paste( "menu nitems=", length( files ), " label=Files", sep="" ) )
		w <- append( w, files )
	}

	return( w )
}
#----------------------------------------.getMenus


#.presentTalkOpenFile-------------------2012-12-20
#  Open files from the win act, and supports multiple files.
#----------------------------------------------ACB
.presentTalkOpenFile <- function()
{
	f <- getWinAct()[ 1 ]
	f <- strsplit( f, "\\s+" )
	print( f )
	sapply( f, openFile )
}
#-----------------------------.presentTalkOpenFile


#.updateSlide---------------------------2012-12-20
.updateSlide <- function( talk )
{
	tget(.PBSmod)
	index <- .PBSmod[[ ".presentTalk" ]][[ talk@name ]]$index 
	indicies <- .getTalkIndexes( talk )
	section_id <- indicies[[ index ]][ 1 ]
	item_id <- indicies[[ index ]][ 2 ]
	items <- talk@sections[[ section_id ]]@items
	num_sections <- length( talk@sections )

	#make sure the correct section is visible
	section_names <- .getSectionNames( talk )
	setWinVal( list( section = section_names[ section_id ] ) )

	#set slide label
	#setWinVal( list( slide_num = paste( "slide: ", index, "/", length( indicies ) ) ) )
	setWinVal( list( slides = index ) )

	if( index > length( indicies ) ) {
		cat( "end of talk" )
		return();
	}

	#setWidgetState( "next", ifelse( index >= length( indicies ), "disabled", "normal" ) );
	#setWidgetState( "prev", ifelse( index <= 1, "disabled", "normal" ) );
	#setWidgetState( "start", ifelse( indicies[[ index ]][ 2 ] == 1, "disabled", "normal" ) );

	setWidgetState( "prev", ifelse( section_id == 1, "disabled", "normal" ) );
	setWidgetState( "next", ifelse( section_id == num_sections, "disabled", "normal" ) );
	#setWidgetState( "start", ifelse( section_id == 1 && item_id == 1, "disabled", "normal" ) );

	setWidgetState( "go", ifelse( index >= length( indicies ), "disabled", "normal" ) );
	setWidgetState( "back", ifelse( index <= 1, "disabled", "normal" ) );


	
	#get next item to iterate to, if next item is in the next section, stop at last in this section
	if( index < length( indicies ) )
		m <- max( length( items ), indicies[[ index + 1 ]][ 2 ] )
	else
		m <- length( items ) #last slide case

	text <- c()
	while( item_id <= m ) {
		item <- items[[ item_id ]]
		item_id <- item_id + 1
		if( inherits( item, "text" ) ) {
			#text items
			cat( item@text, "\n" )
			text <- append( text, item@text )
		} else if( inherits( item, "file" ) ) {
			#file items
			openFile( item@filename )
		} else if( inherits( item, "code" ) ) {
			#code items
			code <- strsplit( item@code, "\n" )[[ 1 ]]
			if( item@print == TRUE && item@eval == FALSE ) {
				#print to consol
				code_cat <- paste( "> ", code, "\n", sep="" )
				cat( code_cat, sep="" )
			}
			if( item@eval == FALSE ) next
			res <- capture.output( eval( parse( text = code ), envir = globalenv() ) )
			#print results
			if( item@show == TRUE ) {
				res <- paste( c( res, ""), collapse = "\n" )
				cat( res )
			}
		} else if( inherits( item, "break" ) ) {
			break
		}
	}

	if( index < length( indicies ) )
		cat( "-----------------------------------< Press Go to continue >---------\n" )
	else
		cat( "--------------------------------------------< End of talk >---------\n" )
}
#-------------------------------------.updateSlide


#.startSlide----------------------------2012-12-20
.startSlide <- function( talk ) {
	tget(.PBSmod)
	#eval(parse(text=".PBSmod[[ \".presentTalk\" ]][[ talk@name ]]$index <<- 1"))
	.PBSmod[[ ".presentTalk" ]][[ talk@name ]]$index <- 1
	tput(.PBSmod)
	talk <- .PBSmod[[ ".presentTalk" ]][[ talk@name ]]$talk
	.updateSlide( talk )
}
#--------------------------------------.startSlide


#.prevSlide-----------------------------2012-12-20
.prevSlide <- function() {
	talk_name <- getWinAct()[1]
	tget(.PBSmod)
	index <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index
	talk <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$talk
	#eval(parse(text=".PBSmod[[ \".presentTalk\" ]][[ talk_name ]]$index <<- index - 1"))
	.PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index <- index - 1
	tput(.PBSmod)
	.updateSlide( talk )
}
#---------------------------------------.prevSlide


#.nextSlide-----------------------------2012-12-20
.nextSlide <- function() {
	talk_name <- getWinAct()[1]
	tget(.PBSmod)
	index <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index
	talk <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$talk
	#eval(parse(text=".PBSmod[[ \".presentTalk\" ]][[ talk_name ]]$index <<- index + 1"))
	.PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index <- index + 1
	tput(.PBSmod)
	.updateSlide( talk )
}
#---------------------------------------.nextSlide


#.slidedrop-----------------------------2012-12-20
.slidedrop <- function() {
	#get talk
	talk_name <- getWinAct()[1]
	tget(.PBSmod)
	talk <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$talk
	index <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index
	section_id <- .getTalkIndexes( talk )[[ index ]][ 1 ]
	new_index <- getWinVal()$slides.id
	#do nothing if current section is re-selected
	if( index == new_index )
		return()
	#eval(parse(text=".PBSmod[[ \".presentTalk\" ]][[ talk_name ]]$index <<- new_index"))
	.PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index <- new_index
	tput(.PBSmod)
	.updateSlide( talk )
}
#---------------------------------------.slidedrop


#.sectiondrop---------------------------2012-12-20
.sectiondrop <- function() {
	#get talk
	talk_name <- getWinAct()[1]
	tget(.PBSmod)
	talk <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$talk
	index <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index
	section_id <- .getTalkIndexes( talk )[[ index ]][ 1 ]
	new_sect_id <- getWinVal()$section.id
	#do nothing if current section is re-selected
	if( section_id == new_sect_id )
		return()
	#eval(parse(text=".PBSmod[[ \".presentTalk\" ]][[ talk_name ]]$index <<- .getIndexForSection( talk, new_sect_id )"))
	.PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index <- .getIndexForSection( talk, new_sect_id )
	tput(.PBSmod)
	.updateSlide( talk )
}
#-------------------------------------.sectiondrop


#.setsection----------------------------2012-12-20
.setsection <- function()
{
	act = getWinAct()[ 1 ]
	act = unlist( strsplit( act, ":" ) )
	talk_name <- act[ 1 ]
	act <- act[ 2 ]
	tget(.PBSmod)
	talk <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$talk
	index <- .PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index
	indicies <- .getTalkIndexes( talk )
	section_id <- indicies[[ index ]][ 1 ]
	if( act == "+1" )
		index <- .getIndexForSection( talk, section_id + 1 )
	else if( act == "-1" )
		index <- .getIndexForSection( talk, section_id - 1 )
	else if( act == "0" )
		index <- .getIndexForSection( talk, section_id )
	else
		index <- .getIndexForSection( talk, as.integer( act ) )
	#eval(parse(text=".PBSmod[[ \".presentTalk\" ]][[ talk_name ]]$index <<- index"))
	.PBSmod[[ ".presentTalk" ]][[ talk_name ]]$index <- index
	tput(.PBSmod)
	.updateSlide( talk )
}
#--------------------------------------.setsection


#presentTalk----------------------------2015-01-20
# Present an interactive talk within R.
#-------------------------------------------ACB/RH
presentTalk <- function( talk )
{
	#if (!requireNamespace("XML", quietly=TRUE))
	#	stop( "The XML R package is required to use `presentTalk' - please install it first. Linux users might have to install `libxml-dev' (via apt) before installing the R XML package." )
	.initPBSoptions()
	tget(.PBSmod)
	#setup .PBSmod$.talk (should be seperate package)
	if( !is.null( .PBSmod[[ ".presentTalk" ]] ) ) {
		#eval(parse(text=".PBSmod[[ \".presentTalk\" ]] <<- list()"))
		.PBSmod[[ ".presentTalk" ]] <- list()
		tput(.PBSmod) }
	#parse XML into a DOM
	talk <- .parseTalkFile( talk )
	indicies <- .getTalkIndexes( talk )
	vals <- 1:(length( indicies ) )
	nvals = length(vals)

	#save parsed talk
	name <- talk@name
	#eval(parse(text=".PBSmod[[ \".presentTalk\" ]][[ name ]] <<- list( index = 0, talk = talk )"))
	.PBSmod[[ ".presentTalk" ]][[ name ]] <- list( index = 0, talk = talk )
	tput(.PBSmod)

	#create a GUI for it
	win_desc <- c(
	paste( "window name=presentwin onclose=.win.restoreCWD title=\"", .addslashes( name ), "\"", sep="" ),
	.getMenus( talk ),
	"grid 1 2 pady=\"0 5\"",
		"grid 1 1 relief=groove",
		"grid 2 1 padx=3 pady=3",
		"grid 1 2 sticky=E",
		"label \"section:\" font=8 sticky=W",
		paste( "droplist name=section values=\"\" function=.sectiondrop bg=skyblue width=15 action=\"", name, "\"", sep="" ),
		"grid 1 3",
		paste( "button name=prev text=\"< Prev\" bg=skyblue sticky=S function=.setsection action=\"",name,":-1\" width=7", sep="" ),
		paste( "button name=curr text=\"Restart\" bg=skyblue sticky=S function=.setsection action=\"",name,":0\" width=7", sep="" ),
		paste( "button name=next text=\"Next >\" bg=skyblue sticky=S function=.setsection action=\"",name,":+1\" width=7", sep="" ),
	
		"grid 1 1 relief=groove padx=\"5 0\"",
		"grid 2 1 padx=3 pady=3",
		"grid 1 3 sticky=E",
		"label \"slide:\" font=\"8\" sticky=W",
		paste( "droplist name=slides values=\"\" function=.slidedrop bg=greenyellow width=9 action=\"",name,"\"", sep="" ),
		paste("label \"/ ",nvals,"\" name=slidecount font=\"8\" sticky=W",sep=""),
		"grid 1 2",
		paste( "button name=back text=\"< Back\" bg=greenyellow sticky=S function=.prevSlide action=\"",name,"\" width=10", sep="" ),
		paste( "button name=go text=\"Go >\" bg=greenyellow sticky=S function=.nextSlide action=\"",name,"\" width=10", sep="" ),
	
	.getButtons( talk ),
	""
	)
	createWin( win_desc, astext = TRUE )
	
	#initialize section droplist
	section_names <- .getSectionNames( talk )
	setWinVal( list( section.values = section_names ) )
	setWinVal( list( section = section_names[1] ) )
	
	#initialize slide droplist
	indicies <- .getTalkIndexes( talk )
	vals <- 1:(length( indicies ) )
	setWinVal( list( slides.values = vals ) )
	setWinVal( list( slides = vals[1] ) )
	#setWinVal( list( slidecount = paste( "/", length( indicies ) ) ) )  # cannot uptdate a label widget

	#move to first slide
	.startSlide( talk )
}
#--------------------------------------presentTalk


#===== THE END ===================================

