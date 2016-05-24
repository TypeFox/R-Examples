indicatorTransform.gui <- function(base.txt) {

    # This function provides a gui interface for taking an indicator transformation 
    # of data, 'x'.
    # That is, functions creates a column in the 'in2extRemesDataObject' list object's data
    # matrix with the value 'z=I_{x > u}', where 'I' is an indicator that is 1 if
    # 'x > u' and 0 otherwise.
    # The gui will list all objects of class "in2extRemesDataObject" and the column names of the
    # user selected object.
    #
    # After taking the indicator transformation of the selected data, will return a 
    # new column with the same column name(s), but with a ".indU" extension; where
    # 'U' is the value used for 'u'
    
    # Set the tcl variables
    
    u <- tclVar("0")
    # Internal functions
    
    refresh <- function() {
    
        # When data is selected, this function fills the lists for the columns of
        # the data set so that the user can select which column(s) to transform.
        
        tkdelete( col.listbox, 0.0, "end")
        
        if( !is.nothing) {

            data.select <- as.numeric( tkcurselection( data.listbox))+1
            dd <- get( full.list[ data.select])

        } else dd <- in2extRemesData
        
            for( i in 1:ncol(dd$data)) tkinsert(col.listbox, "end", paste( colnames( dd$data)[i]))

        invisible()

    } # end of internal 'refresh' fcn
    
    indicator.trans <- function() {

        # Function invoked when the "ok" button is pressed.  Actually takes indicator
        # transformation of the data.
    
        if(!is.nothing) {

    	    data.select <- as.numeric(tkcurselection(data.listbox)) + 1
    	    dd.cmd <- paste("dd <- get(\"", full.list[ data.select], "\")", sep="")

        } else dd.cmd <- "dd <- in2extRemesData"

    	eval( parse( text=dd.cmd))
    	write( dd.cmd, file="in2extRemes.log", append=TRUE)
    
        cols.selected.cmd <- "cols.selected <- numeric(0)"
        eval( parse( text=cols.selected.cmd))
        write( cols.selected.cmd, file="in2extRemes.log", append=TRUE)
    
        # Gather column names to add indicator transform(s) to end.
        cnames.cmd <- "cnames <- colnames( dd[[\"data\"]])"
        eval( parse( text=cnames.cmd))
        write( cnames.cmd, file="in2extRemes.log", append=TRUE)
        
        temp <- as.numeric( tkcurselection( col.listbox)) + 1
        
        # Make sure a column has been selected.
        if( is.na( temp)) return()
        
        # cols.selected <- cnames[temp]
        for( i in 1:length( temp)) {

            cols.selected.cmd <- paste( "cols.selected <- c( cols.selected, ", temp, ")", sep="") 
            eval( parse( text=cols.selected.cmd))
            write( cols.selected.cmd, file="in2extRemes.log", append=TRUE)

        }

        cols.selected.cmd <- "cols.selected <- cnames[cols.selected]"
        eval( parse( text=cols.selected.cmd))
        write( cols.selected.cmd, file="in2extRemes.log", append=TRUE)
        
        u.value <- as.numeric(tclvalue(u))
        z.cmd <- "z <- dd[[\"data\"]][, cols.selected]"
        eval( parse( text=z.cmd))
        write( z.cmd, file="in2extRemes.log", append=TRUE)

        ind.cmd <- paste( "ind <- z > ", u.value, sep="")
        eval( parse( text=ind.cmd))
        write( ind.cmd, file="in2extRemes.log", append=TRUE)

        z.cmd <- "z[ind] <- 1"
        eval( parse( text=z.cmd))
        write( z.cmd, file="in2extRemes.log", append=TRUE)

        z.cmd <- "z[!ind] <- 0"
        eval( parse( text=z.cmd))
        write( z.cmd, file="in2extRemes.log", append=TRUE)
        
        newnames.cmd <- "newnames <- character(0)"
        eval( parse( text=newnames.cmd))
        write( newnames.cmd, file="in2extRemes.log", append=TRUE)

        for( i in 1:length( temp)) {

            newnames.cmd <- paste( "newnames <- c( newnames, \"", paste( cnames[temp[i]], ".ind", u.value, sep=""),
        					"\")", sep="")
            eval( parse( text=newnames.cmd))
            write( newnames.cmd, file="in2extRemes.log", append=TRUE)

        } # end of for 'i' loop.

        # cnames <- c( cnames, newnames)
        cnames.cmd <- "cnames <- c( cnames, newnames)"
        eval( parse( text=cnames.cmd))
        write( cnames.cmd, file="in2extRemes.log", append=TRUE)

        dd.cmd <- "dd[[\"data\"]] <- cbind( dd[[\"data\"]], z)"
        eval( parse( text=dd.cmd))
        write( dd.cmd, file="in2extRemes.log", append=TRUE)

        colnames.cmd <- "colnames( dd[[\"data\"]]) <- cnames"
        eval( parse( text=colnames.cmd))
        write( colnames.cmd, file="in2extRemes.log", append=TRUE)

        assignCMD <- paste( "assign( \"", full.list[ data.select], "\", dd, pos=\".GlobalEnv\")", sep="")
        eval( parse( text=assignCMD))
        write( assignCMD, file="in2extRemes.log", append=TRUE)
        
        # tkconfigure( base.txt, state="normal")
        msg <- paste( "\n", "Indicator transformation z = I_{x > ", u.value, "}",
        		" taken for", "\n",
        		" x= ", cols.selected, " and assigned to ", newnames, "\n", sep="")
        cat( msg)
        tkdestroy( base)

        invisible()

    } # end of indicator.trans fcn
    
    indicatorhelp <- function() {

    	help.msg1 <- paste( " ",
    	    "This is a simple function that takes a data set X and returns the",
    	    "indicator transformation of X.  That is", "  ", sep="\n")

    	help.msg2 <- paste(" ", "z=I_{X > u}", " ",
    	    "where \'I\' is 1 if X > u and 0 otherwise.", sep="\n")

        help.msg3 <- paste( " ", "Returns object of class \"in2extRemesDataObject\" with an extra column in the data",
    	    "indicated by a .indU extension where \'U\' is the value of u.", " ", sep="\n")

    	cat( help.msg1)
    	cat( help.msg2)
    	cat( help.msg3)

    	invisible()

    } # end of internal 'indicatorhelp' fcn
    
    endprog <- function() {

    	tkdestroy( base)

    }
    
    #####################
    # Frame/button setup
    #####################
    
    base <- tktoplevel()
    tkwm.title( base, "Indicator Transformation")
    
    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")
    
    # Top frame for data sets...
    
    data.listbox <- tklistbox( top.frm,
    			yscrollcommand=function(...) tkset(data.scroll, ...),
    			selectmode="single",
    			width=20,
    			height=5,
    			exportselection=0)
    
    data.scroll <- tkscrollbar( top.frm, orient="vert",
    			command=function(...) tkyview( data.listbox, ...))
    
    temp <- ls(all.names=TRUE, name=".GlobalEnv")
    full.list <- character(0)
    is.nothing <- TRUE
    for( i in 1:length( temp)) {
            if( is.null( class( get( temp[i])))) next
            if( (class(get( temp[i]))[1] == "in2extRemesDataObject")) {
                    tkinsert( data.listbox, "end", paste( temp[i]))
            	full.list <- c( full.list, temp[i])
    		is.nothing <- FALSE
    		}
    } # end of for i loop
    
    tkpack( tklabel( top.frm, text="Data Object", padx=4), side="top")
    tkpack( data.listbox, data.scroll, side="left", fill="y")
    # tkpack( data.scroll, side="right", fill="y")
    tkpack( top.frm)
    tkbind( data.listbox, "<Button-1>", "")
    tkbind( data.listbox, "<ButtonRelease-1>", refresh)
    
    # Middle frame for columns of chosen dataset...
    
    var.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
    col.listbox <- tklistbox( var.frm,
    			yscrollcommand=function(...) tkset( col.scroll, ...),
    			selectmode="multiple",
    			width=20,
    			height=5,
    			exportselection=0)
    
    col.scroll <- tkscrollbar( var.frm, orient="vert",
    			command=function(...) tkyview( col.listbox, ...))
    
    if( is.nothing) {
    for( i in 1:ncol( in2extRemesData$data))
    	tkinsert( col.listbox, "end", paste( colnames( dd$data)[i]))
    # end of for i loop
    	} else tkinsert( col.listbox, "end", "")
    
    tkpack( tklabel( var.frm, text="Variables to Transform", padx=4),
    		side="top")
    tkpack( col.listbox, side="left")
    tkpack( col.scroll, side="right")
    tkpack( var.frm, side="top")
    
    u.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    u.entry <- tkentry( u.frm, textvariable=u, width=4)
    tkpack( tklabel( u.frm, text="threshold (u)", padx=4), u.entry, side="left")
    tkpack( u.frm)
    
    ok.but <- tkbutton( bot.frm, text="OK", command=indicator.trans)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    
    help.but <- tkbutton( bot.frm, text="Help", command=indicatorhelp)
    
    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")
    
    # place bindings for return key.
    tkbind( ok.but, "<Return>", indicator.trans)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", indicatorhelp)
    
    tkpack( top.frm, fill="x")
    tkpack( mid.frm, fill="x")
    tkpack( bot.frm, side="bottom")

} # end of 'indicator.gui' fcn
