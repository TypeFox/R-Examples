DataSummaryGUI <- function(base.txt) {

    # This function provides a gui interface for simply printing a summary of all
    # of the data contained in the 'data' component of an 'in2extRemesDataObject'.

    # Internal functions

    submit <- function() {

        # Function invoked when the "ok" button is pressed.
        # Invokes the summary command.

        if( !is.nothing) {

	    data.select <- as.numeric( tkcurselection( data.listbox))+1
	    dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")

        } else dd.cmd <- "dd <- in2extRemesData"

        eval( parse( text=dd.cmd))
        write( dd.cmd, file="in2extRemes.log", append=TRUE)

        summaryCMD <- "print( summary( dd[[\"data\"]]))"
        eval( parse( text=summaryCMD))
        write( summaryCMD, file="in2extRemes.log", append=TRUE)

        tkdestroy( base)
        invisible()
    } # end of neg.trans fcn

    DataSummaryHelp <- function() {

	# tkconfigure( base.txt, state="normal")

	help.msg1 <- paste( " ",
            "This is a simple function that summarizes the entire selected data set and prints the summary to the console", "  ",
			sep="\n")

	# tkinsert( base.txt, "end", help.msg1)

	cat( help.msg1)
	# tkconfigure( base.txt, state="disabled")

	invisible()

    } # end of internal 'DataSummaryHelp' fcn

    endprog <- function() {

	tkdestroy( base)

    } # end of internal 'endprog' function.

    #####################
    # Frame/button setup
    #####################

    base <- tktoplevel()
    tkwm.title( base, "Data Summary")

    top.frm <- tkframe( base, borderwidth=2, relief="groove")
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
    tkpack( top.frm)

    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)

    help.but <- tkbutton( bot.frm, text="Help", command=DataSummaryHelp)

    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")

    # place bindings for return key.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", DataSummaryHelp)

    tkpack( top.frm, fill="x")
    tkpack( bot.frm, side="bottom")

    } # end of 'DataSummaryGUI' fcn
