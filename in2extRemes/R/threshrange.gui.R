threshrange.gui <- function( base.txt) {
    
    #
    # Function to provide guis for the 'gpd.fitrange' fcn of Stuart Coles.
    #
    
    # Initialize tcl variables.
    
    umin.value <- tclVar("")
    umax.value <- tclVar("")
    type.value <- tclVar("GP")
    nint.value <- tclVar("")
    na.value <- tclVar("na.fail")
    
    # Internal functions...
    
    refresh <- function() {
    #
    # function to determine which variables exist given data object chosen.
    
    	# Remove recent contents of 'var.listbox'.
    	tkdelete( var.listbox, 0.0, "end")
    
    	# Obtain recent data object chosen.
    	data.select <- as.numeric( tkcurselection( data.listbox))+1
    	dd <- get( full.list[ data.select])
    
    	# put column names into 'var.listbox'.
    	tmp <- colnames( dd$data)
    	for( i in 1:length( tmp))
    		tkinsert( var.listbox, "end", paste( tmp[i]))
    
    	invisible()
    	} # end of refresh fcn
    
    submit <- function() {
        
        # Function that is executed when the "OK" button is hit.  Actually calls
        # the 'gpd.fitrange' fcn.
        
        # Obtain argument values entered into gui.
        	umin.val <- as.numeric( tclvalue( umin.value))
        	umax.val <- as.numeric( tclvalue( umax.value))
        	nint.val <- as.numeric( tclvalue( nint.value))
        
        # Obtain data to use in 'gpd.fitrange' fcn.
        	if( !is.nothing) {
        		data.select <- as.numeric( tkcurselection( data.listbox))+1
        		dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")
        		} else  dd.cmd <- "dd <- in2extRemesData"
        	eval( parse( text=dd.cmd))
        	write( dd.cmd, file="in2extRemes.log", append=TRUE)
        
        	var.select <- as.numeric( tkcurselection( var.listbox))+1
        	var.val.cmd <- paste( "var.val <- dd[[\"data\"]][, ", var.select, "]", sep="")
        	eval( parse( text=var.val.cmd))
        	write( var.val.cmd, file="in2extRemes.log", append=TRUE)
        
        	fitrange.cmd <- paste("threshrange.plot(x = var.val, r = c(", umin.val, ", ", umax.val,
                                        "), type = \"", tclvalue(type.value), 
    				    "\", nint = ", nint.val,
    				    ", na.action = ", tclvalue(na.value), 
    				    ")", sep="")
    
            eval( parse(text=fitrange.cmd))
            write( fitrange.cmd, file="in2extRemes.log", append=TRUE)
    
            tkdestroy( base)

    	    invisible()

    	} # end of internal 'submit' fcn
    
    gpdfitrangehelp <- function() {

    	cat("\n", "Invokes the \'extRemes\' function: \'threshrange.plot\'\n")
    	help(threshrange.plot)

    	invisible()

    } # end of gpdfitrangehelp fcn
    
    endprog <- function() {
    	tkdestroy( base)
    	} # end of endprog fcn
    
    #######################
    # Frame/button setup
    #######################
    
    base <- tktoplevel()
    tkwm.title( base, "Fit a POT model for a range of thresholds")
    
    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")
    
    # Top frame to select data object.
    
    data.listbox <- tklistbox( top.frm,
    			yscrollcommand=function(...) tkset( data.scroll, ...),
    			selectmode="single",
    			width=20,
    			height=5,
    			exportselection=0)
    
    data.scroll <- tkscrollbar( top.frm, orient="vert",
    			command=function(...) tkyview( data.listbox, ...))
    
    temp <- ls( all.names=TRUE, name=".GlobalEnv")
    full.list <- character(0)
    is.nothing <- TRUE
    for( i in 1:length( temp)) {
    	if( is.null( class( get( temp[i])))) next
    	if( (class( get( temp[i]))[1] == "in2extRemesDataObject")) {
    		tkinsert( data.listbox, "end", paste( temp[i]))
    		full.list <- c( full.list, temp[i])
    		is.nothing <- FALSE
    			}
    	} # end of for i loop
    
    tkpack( tklabel( top.frm, text="Data Object", padx=4), side="left")
    tkpack( data.listbox, side="left")
    tkpack( data.scroll, side="right", fill="y")
    tkpack( top.frm)
    
    # place binding on 'data.listbox' to reflect the chosen data from list.
    tkbind( data.listbox, "<Button-1>", "")
    tkbind( data.listbox, "<ButtonRelease-1>", refresh)
    
    # Middle frame tochoose which columns (vars) of data to use and other args
    # to 'gpd.fitrange' fcn.
    
    midleft <- tkframe( mid.frm, borderwidth=2, relief="groove")
    
    var.listbox <- tklistbox( midleft,
    			yscrollcommand=function(...) tkset( var.scroll, ...),
    			selectmode="single",
    			width=15,
    			height=6,
    			exportselection=0)
    
    var.scroll <- tkscrollbar( midleft, orient="vert",
    			command=function(...) tkyview( var.listbox, ...))
    
    # Insert column names of 'in2extRemesData' if no other data objects exist.
    # Otherwise, begin with no variables.
    if( is.nothing) {
    	for( i in 1:length( colnames( in2extRemesData$data)))
    		tkinsert( var.listbox, "end",
    				paste( colnames( in2extRemesData$data)[i]))
    } else tkinsert( var.listbox, "end", " ")
    
    tkpack( tklabel( midleft, text="Select Variable", padx=4), side="top")
    tkpack( var.listbox, side="left")
    tkpack( var.scroll, side="right")
    
    # Other args to 'gpd.fitrange' fcn.
    
    midright <- tkframe( mid.frm, borderwidth=2, relief="groove")
    
    umin.frm <- tkframe( midright, borderwidth=2, relief="flat")
    umin.entry <- tkentry( umin.frm, textvariable=umin.value, width=5)
    
    umax.frm <- tkframe( midright, borderwidth=2, relief="flat")
    umax.entry <- tkentry( umax.frm, textvariable=umax.value, width=5)
  
    type.frm <- tkframe(midright, borderwidth = 2, relief = "flat")
    type.entry <- tkentry(type.frm, textvariable = type.value, width = 2)
  
    nint.frm <- tkframe( midright, borderwidth=2, relief="flat")
    nint.entry <- tkentry( nint.frm, textvariable=nint.value, width=5)

    na.frm <- tkframe(midright, borderwidth = 2, relief = "flat")
    na.entry <- tkentry(na.frm, textvariable = na.value, width = 20)
    
    tkpack( tklabel( umin.frm, text="Minimum Threshold", padx=4), umin.entry, side="left", fill="y", anchor="w")

    tkpack( tklabel( umax.frm, text="Maximum Threshold", padx=4), umax.entry, side="left", fill="y", anchor="w")

    tkpack(tklabel(type.frm, text = "Type of POT model", padx = 4), type.entry, side = "left", fill = "y", anchor = "w") 

    tkpack( tklabel( nint.frm, text="Number of thresholds", padx=4), nint.entry, side="left", fill="y", anchor="w")

    tkpack(tklabel(na.frm, text = "NA Action", padx = 4), na.entry, side = "left", fill = "y", anchor = "w")
    
    tkpack( umin.frm)
    tkpack( umax.frm)
    tkpack(type.frm)
    tkpack( nint.frm)
    tkpack(na.frm)
    
    tkpack( midleft, side="left")
    tkpack( midright, side="right")
    
    # Bottom frame for execution or cancellation.
    
    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton( bot.frm, text="Help", command=gpdfitrangehelp)
    
    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")
    
    # Place bindings on "OK" and "Cancel" so that Return key will execute them.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", gpdfitrangehelp)
    
    tkpack( top.frm, side="top")
    tkpack( mid.frm, fill="x")
    tkpack( bot.frm, side="bottom")

} # end of 'threshrange.gui' fcn
