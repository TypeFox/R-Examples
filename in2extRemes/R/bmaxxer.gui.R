bmaxxer.gui <- function(base.txt) {

    # This function provides a gui interface for taking the block maxima a dataset.
    # The gui will list all objects of class
    # "in2extRemesDataObject" and the column names of the user selected object.
    # After taking the block maxima of the selected column, will return a new data frame,
    # which will keep the same columns and names as the original data frame.

    # Set the tcl variables
    save.as.value <- tclVar("")

    # Internal functions

    refresh <- function() {

        # When data is selected, this function fills the lists for the columns of
        # the data set so that the user can select which column(s) to decluster.

	tkdelete(col.listbox, 0.0, "end")
	tkdelete(groupby.listbox, 0.0, "end")

	if( !is.nothing) {

	    data.select <- as.numeric(tkcurselection(data.listbox))+1
	    dd <- get(full.list[data.select])

	} else dd <- in2extRemesData

	for(i in 1:ncol(dd$data)) {

	    tkinsert(col.listbox, "end", paste( colnames( dd$data)[i])) 
	    tkinsert(groupby.listbox, "end", paste( colnames( dd$data)[i]))

	} # end of for 'i' loop.

	invisible()

    } # end of internal 'refresh' fcn

    submit <- function() {

        # Function invoked when the "ok" button is pressed.  Actually declusters
        # the data.

        if( !is.nothing) {

    	    data.select <- as.numeric( tkcurselection( data.listbox))+1
	    dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")
	    # M <- dim( dd$data)[1]

	} else dd.cmd <- "dd <- in2extRemesData"

        eval(parse(text=dd.cmd))
        write(dd.cmd, file="in2extRemes.log", append=TRUE)

	cnames.cmd <- "cnames <- colnames( dd[[\"data\"]])"
	eval(parse(text=cnames.cmd))
	write(cnames.cmd, file="in2extRemes.log", append=TRUE)

	temp <- as.numeric(tkcurselection(col.listbox)) + 1

	# Make sure a column has been selected.
	if(is.na(temp)) return()

	# Grab column to degroup.by...
	temp2 <- as.numeric(tkcurselection(groupby.listbox))+1

	group.by.col <- cnames[temp2]

	print( paste( "Finding block maxima ...", sep=""))

	out.cmd <- paste("out <- blockmaxxer(x = dd[[\"data\"]], blocks = dd[[\"data\"]][[\"", group.by.col, "\"]], which = \"", cnames[temp], "\")", sep="")
	eval(parse(text = out.cmd))
	write(out.cmd, file = "in2extRemes.log", append = TRUE)

	## Save the new data frame.
	save.as <- tclvalue(save.as.value)

        if(save.as != "") {

	    cmd <- "dd2 <- as.in2extRemesDataObject(out)"
	    eval(parse(text = cmd))
	    write(cmd, file = "in2extRemes.log", append = TRUE)

            assignCMD <- paste( "assign( \"", save.as, "\", dd2, pos=\".GlobalEnv\")", sep="")
            eval( parse( text=assignCMD))
            write( assignCMD, file="in2extRemes.log", append=TRUE)

            cat("\n", "Saving workspace (may take a few moments) ...\n")
            saveCMD <- "save.image()"
            eval( parse( text=saveCMD))
            write( saveCMD, file="in2extRemes.log", append=TRUE)

	} else cat("\n", "Block maxima taken, but not saved!\n")

 	endprog()

        invisible()

    } # end of internal 'submit' fcn

    bmhelp <- function() {

	help(bmaxxer)

	invisible()

    } # end of internal 'lthelp' function.

    endprog <- function() {

	tkdestroy(base)

    }

    #####################
    # Frame/button setup
    #####################

    base <- tktoplevel()
    tkwm.title( base, "Block Maxima")

    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    data.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
    groupby.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")
    save.frm <- tkframe(bot.frm, borderwidth = 2, relief = "groove")

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

    for(i in 1:length( temp)) {

        if(is.null( class( get( temp[i])))) next

        if((class(get(temp[i]))[1] == "in2extRemesDataObject")) {

                tkinsert(data.listbox, "end", paste( temp[i]))
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

    # leftmid.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
    col.listbox <- tklistbox( data.frm,
			yscrollcommand=function(...) tkset( col.scroll, ...),
			selectmode="single",
			width=20,
			height=5,
			exportselection=0)

    col.scroll <- tkscrollbar( data.frm, orient="vert",
			command=function(...) tkyview( col.listbox, ...))

    if( is.nothing) {

        for(i in 1:ncol(in2extRemesData$data)) tkinsert( col.listbox, "end", paste( colnames( dd$data)[i]))

    } else tkinsert(col.listbox, "end", "")

    tkpack( tklabel( data.frm, text="Find block maxima of ...", padx=4), side="top")
    # tkpack( col.listbox, side="left")
    # tkpack( col.scroll, side="right")
    tkpack( col.listbox, col.scroll, side="left", fill="y")

    # group.by frame...
    groupby.listbox <- tklistbox( groupby.frm,
                        yscrollcommand=function(...) tkset(groupby.scroll, ...),
                        selectmode="single",
                        width=20,
                        height=5,
                        exportselection=0)
    groupby.scroll <- tkscrollbar( groupby.frm, orient="vert",
                        command=function(...) tkyview( groupby.listbox, ...))
    if(is.nothing) {

	for( i in 1:ncol( in2extRemesData$data)) tkinsert( groupby.listbox, "end", paste( colnames( dd$data)[i]))

    } else tkinsert( col.listbox, "end", "")

    tkpack(tklabel(groupby.frm, text="Blocks", padx=4), side="top")
    tkpack(groupby.listbox, groupby.scroll, side="left", fill="y")

    tkpack(data.frm, groupby.frm, side="left", fill="x")

    save.as.entry <- tkentry( save.frm, textvariable=save.as.value, width=20) 
    tkpack(tklabel(save.frm, text = "Save As (in R)", padx = 4), save.as.entry, side = "left")
    tkpack(save.frm, side = "top")

    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)

    help.but <- tkbutton( bot.frm, text="Help", command=bmhelp)

    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")

    # place bindings for return key.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", bmhelp)

    tkpack( top.frm, fill="x")
    tkpack( mid.frm, fill="x")
    tkpack( bot.frm, side="bottom")

} # end of bmaxxer.gui fcn
