atdf.gui <- function( base.txt) {
    
    #
    # Function to provide guis for the 'gpd.fitrange' fcn of Stuart Coles.
    #
    
    # Initialize tcl variables.
    
    threshold.value <- tclVar("0.8")
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
        	threshold.val <- as.numeric( tclvalue( threshold.value))
        
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
       
	atdf.cmd <- paste("atdf(x = var.val, u = ", threshold.val, ", na.action = ", tclvalue(na.value), ")", sep = "") 
    
        eval( parse(text=atdf.cmd))
        write( atdf.cmd, file="in2extRemes.log", append=TRUE)
    
        tkdestroy( base)

    	invisible()

    } # end of internal 'submit' fcn
    
    atdfhelp <- function() {

    	cat("\n", "Invokes the \'extRemes\' function: \'atdf\'\n")
    	help(threshrange.plot)

    	invisible()

    } # end of atdfhelp fcn
    
    endprog <- function() {
    	tkdestroy( base)
    	} # end of endprog fcn
    
    #######################
    # Frame/button setup
    #######################
    
    base <- tktoplevel()
    tkwm.title( base, "Auto Tail-Dependence Function")
    
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
    
    threshold.frm <- tkframe( midright, borderwidth=2, relief="flat")
    threshold.entry <- tkentry( threshold.frm, textvariable=threshold.value, width=5)
    
    na.frm <- tkframe(midright, borderwidth = 2, relief = "flat")
    na.entry <- tkentry(na.frm, textvariable = na.value, width = 20)
    
    tkpack( tklabel( threshold.frm, text="Quantile Threshold", padx=4), threshold.entry, side="left", fill="y", anchor="w")

    tkpack(tklabel(na.frm, text = "NA Action", padx = 4), na.entry, side = "left", fill = "y", anchor = "w")
    
    tkpack( threshold.frm)
    tkpack(na.frm)
    
    tkpack( midleft, side="left")
    tkpack( midright, side="right")
    
    # Bottom frame for execution or cancellation.
    
    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton( bot.frm, text="Help", command=atdfhelp)
    
    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")
    
    # Place bindings on "OK" and "Cancel" so that Return key will execute them.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", atdfhelp)
    
    tkpack( top.frm, side="top")
    tkpack( mid.frm, fill="x")
    tkpack( bot.frm, side="bottom")

} # end of 'atdf.gui' fcn
