scrubber.gui <- function( base.txt) {
    
    ##################################
    # Internal functions
    ##################################

    # Refresh fcn 
    refresh <- function() {
    	tkdelete( fit.listbox, 0.0, "end")
    	tkdelete( resp.listbox, 0.0, "end")
    	if( !is.nothing) {
                    data.select <- as.numeric( tkcurselection( data.listbox))+1
                    dd <- get( full.list[ data.select])
                    } else stop("fitdiag.gui: Must load a data object!")
    	for( i in 1:ncol(dd$data)) tkinsert( resp.listbox, "end", paste( colnames( dd$data)[i]))
    	models.fit <- names( dd$models)
    	for( i in 1:length( models.fit))
    		tkinsert( fit.listbox, "end", paste( models.fit[i]))
    	invisible()
    } # end of refresh fcn.
    
    submit <- function() {
    	if( !is.nothing) {
                    data.select <- as.numeric( tkcurselection( data.listbox))+1
    		data.name <- full.list[ data.select]
                    dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")
                    } else stop("scrubber.gui: Must load a data object!")
    	eval( parse( text=dd.cmd))
    	write( dd.cmd, file="in2extRemes.log", append=TRUE)
    
    	# Scrub any selected data columns.
    	resp.select <- as.numeric( tkcurselection( resp.listbox))+1
    	if( tclvalue( tkcurselection( resp.listbox)) != "") {
    		NameCMD <- "tmp.names <- colnames( dd[[\"data\"]])"
    		eval( parse( text=NameCMD))
    		write( NameCMD, file="in2extRemes.log", append=TRUE)
    		if( length( resp.select) == ncol( dd$data)) {
    			warning("scrubber.gui: Removing entire data set!")
    			dd.cmd <- "dd[[\"data\"]] <- NULL"
    			eval( parse( text=dd.cmd))
    			write( dd.cmd, file="in2extRemes.log", append=TRUE)
    		} else {
    			if( length( resp.select) > 1) {
    				respos <- character(0)
    				for( i in 1:(length( resp.select)-1)) respos <- paste( respos, resp.select[i], ", ", sep="")
    				respos <- paste( "c( ", respos, resp.select[ length( resp.select)], ")", sep="")
    			} else respos <- resp.select
    			# tmp <- dd$data[,-resp.select]
    			tmpCMD <- paste( "tmp <- dd[[\"data\"]][, -", respos, "]", sep="")
    			eval( parse( text=tmpCMD))
    			write( tmpCMD, file="in2extRemes.log", append=TRUE)
    			colnamesCMD <- paste( "colnames( tmp) <- tmp.names[-", respos, "]", sep="")
    			eval( parse( text=colnamesCMD))
                            write( colnamesCMD, file="in2extRemes.log", append=TRUE)
    			dd.cmd <- "dd[[\"data\"]] <- tmp"
    			eval( parse( text=dd.cmd))
                            write( dd.cmd, file="in2extRemes.log", append=TRUE)
    			} # end of if else remove all data stmt.
    		} # end of if a data column selected or not stmt.
    
    	# Scrub any selected model fits.
    	fit.select <- as.numeric( tkcurselection( fit.listbox))+1
    	if( tclvalue( tkcurselection( fit.listbox)) != "" ) {
    		if( length( fit.select) == length( dd$models)) {
    			warning("scrubber.gui: Removing all model fits!")
    			dd.cmd <- "dd[[\"models\"]] <- NULL"
    			eval( parse( text=dd.cmd))
    			write( dd.cmd, file="in2extRemes.log", append=TRUE)
    		} else {
    		cat("\n", "Sorry, this feature no longer available.  May remove all fits only.\n")
    		do.it <- FALSE
    		if( do.it) {
    		if( length( fit.select) > 1) {
    		   fitsos <- character(0)
    		   for( i in 1:(length( fit.select)-1)) fitsos <- paste( fitsos, fit.select[i], ", ", sep="")
    		   fitsos <- paste( "c( ", fitsos, fit.select[ length( fit.select)], ")", sep="")
    		} else fitsos <- fit.select
    
    		tmpCMD <- paste( "tmp.names <- names( dd[[\"models\"]])[ -", fitsos, "]", sep="")
    		eval( parse( text=tmpCMD))
    		write( tmpCMD, file="in2extRemes.log", append=TRUE)
    
    		tmpCMD <- "tmp <- list()"
    		eval( parse( text=tmpCMD))
    		write( tmpCMD, file="in2extRemes.log", append=TRUE)
    
    		cmd <- "tmp2 <- dd[[\"models\"]]"
    		eval( parse( text=cmd))
                    write( cmd, file="in2extRemes.log", append=TRUE)
    
    		for( i in 1:length( tmp.names)) {
    		   tmpCMD <- paste( "tmp[[", i, "]] <- tmp2[[\"", tmp.names[i], "\"]]", sep="")
    		   eval( parse( text=tmpCMD))
    		   write( tmpCMD, file="in2extRemes.log", append=TRUE)
    		} # end of for 'i' loop.
    
    		tmpCMD <- "names( tmp) <- tmp.names"
    		eval( parse( text=tmpCMD))
    		write( tmpCMD, file="in2extRemes.log", append=TRUE)
    
    		cmd <- "dd[[\"models\"]] <- NULL"
    		eval( parse( text=cmd))
    		write( cmd, file="in2extRemes.log", append=TRUE)
    
    		cmd <- "dd[[\"models\"]] <- list()"
    		eval( parse( text=cmd))
                    write( cmd, file="in2extRemes.log", append=TRUE)
    
    		dd.cmd <- "dd[[\"models\"]] <- tmp"
    		eval( parse( text=dd.cmd))
    		write( dd.cmd, file="in2extRemes.log", append=TRUE)
    		} # end of if 'do.it' stmts.
    	   } # end of if else remove all models stmt.  
    	} # end of if any models selected or not stmt.
    	assignCMD <- paste( "assign( \"", data.name, "\", dd, pos=\".GlobalEnv\")", sep="")
    	eval( parse( text=assignCMD))
    	write( assignCMD, file="in2extRemes.log", append=TRUE)
    	tkdestroy(base)
        # tkconfigure(base.txt,state="disabled")
    } # end of submit fcn.
    
    scrubberhelp <- function() {
    	# tkconfigure( base.txt, state="normal")
    	nl1 <- paste(" ", "---------------------------------", " ", sep="\n")
            help.msg1 <- paste( " ", "Scrubber Info:", " ", "Select data columns and/or model fits for deletion.", " ", sep="\n")
    	help.msg2 <- paste(" ",
    		"To remove an entire data object use the \'rm\' function: (use \'help( rm)\' for more help)", " ",
    				" ", "from the R command prompt.", sep="\n")
    	# tkinsert( base.txt, "end", nl1)
    	cat( nl1)
            # tkinsert( base.txt, "end", help.msg1)
    	cat( help.msg1)
    	# tkinsert( base.txt, "end", help.msg2)
    	cat( help.msg2)
    	# tkinsert( base.txt, "end", nl1)
    	cat( nl1)
            # tkconfigure( base.txt, state="disabled")
            invisible()
    	} # end of scrubberhelp fcn
    
    endprog <- function() {
    	tkdestroy( base)
    	}
    # Function to plot diagnostic plots for various fits.
    
    #####################
    # Frame/button setup.
    #####################
    
    base <- tktoplevel()
    tkwm.title( base, "Scrubber")
    
    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    top.l <- tkframe( base, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")
    
    ## Data Objects.
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
    
    tkpack( tklabel( top.frm, text="Data Object:  ", padx=4), side="left")
    tkpack( data.listbox, data.scroll,  side="left", fill="y")
    
    # Place bindings on data listbox to update fit listbox.
    tkbind( data.listbox, "<Button-1>", "")
    tkbind( data.listbox, "<ButtonRelease-1>", refresh)
    
    ## Data Columns (to be scrubbed).
    resp.listbox <-
            tklistbox(top.l,yscrollcommand=function(...)tkset(resp.scroll,...),
                            selectmode="multiple",width=35,height=4,exportselection=0)
    resp.scroll <- tkscrollbar(top.l,orient="vert",
                            command=function(...)tkyview(resp.listbox,...))
    
    if( is.nothing) {
    for( i in 1:ncol(in2extRemesData$data))
            tkinsert( resp.listbox, "end", paste(colnames(in2extRemesData$data)[i]))
    # end of for i loop
            } else tkinsert( resp.listbox, "end", "")
    
    tkpack(tklabel(top.l,text="Data Columns:",padx=4), side="left")
    tkpack(resp.listbox,side="left")
    tkpack(resp.scroll,side="right",fill="y")
    
    ## Middle frame for choosing which fit to plot.
    fit.listbox <- tklistbox( mid.frm,
    			yscrollcommand=function(...) tkset( fit.scroll, ...),
    			selectmode="multiple",
    			width=20,
    			height=5,
    			exportselection=0)
    
    fit.scroll <- tkscrollbar( mid.frm, orient="vert",
    			command=function(...) tkyview( fit.listbox, ...))
    tkinsert( fit.listbox, "end", "")
    
    tkpack( tklabel( mid.frm, text="Model Fits: ", padx=4), side="left")
    tkpack( fit.listbox, fit.scroll, side="left", fill="y")
    
    # Bottom frame for execution and cancellation.
    
    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton( bot.frm, text="Help", command=scrubberhelp)
    
    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")
    
    # place bindings on buttons.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", scrubberhelp)
    
    tkpack( top.frm, top.l, mid.frm, bot.frm, side="top", fill="x")
    invisible()
    
} # end of 'scrubber.gui' fcn
