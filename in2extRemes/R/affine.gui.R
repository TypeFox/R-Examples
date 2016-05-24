affine.gui <- function(base.txt) {

    # This function provides a gui interface for taking an affine transformation of data, 'x'.
    # That is, (x-c)/b, where 'c' and 'b' are constants.  The gui will list all objects of class
    # "in2extRemesDataObject" and the column names of the user selected object.
    #
    # After taking the affine transformation of the selected data, will return a new column with
    # the same column name(s), but with a ".cNbM" extension where 'N' is the value used for 'c' and
    # 'M' the value used for 'b'.

    # Set the tcl variables

    c <- tclVar("0")
    b <- tclVar("1") # Internal functions

    refresh <- function() {

        # When data is selected, this function fills the lists for the columns of
        # the data set so that the user can select which column(s) to transform.

	tkdelete(col.listbox, 0.0, "end")

	if( !is.nothing) {

	    data.select <- as.numeric(tkcurselection(data.listbox))+1
	    dd <- get(full.list[data.select])

	} else dd <- in2extRemesData

	for(i in 1:ncol( dd$data)) tkinsert( col.listbox, "end", paste( colnames( dd$data)[i]))

	invisible()

    } # end of refresh fcn

    affine.trans <- function() {

        # Function invoked when the "ok" button is pressed.  Actually takes affine
        # transformation of the data.

        if(!is.nothing) {

	    data.select <- as.numeric( tkcurselection( data.listbox))+1
	    dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")

	} else dd.cmd <- "dd <- in2extRemesData"

        eval(parse(text=dd.cmd))
        write(dd.cmd, file="in2extRemes.log", append=TRUE)

        # Gather column names to add affine transform(s) to end.
        cnames.cmd <- "cnames <- colnames( dd[[\"data\"]])"
        eval(parse( text=cnames.cmd))
        write(cnames.cmd, file="in2extRemes.log", append=TRUE)

        temp <- as.numeric(tkcurselection( col.listbox)) + 1

        # Make sure a column has been selected.
        if(is.na(temp)) return()

	if(length(temp) > 1) {

	    tempos <- character(0)
	    for(i in 1:(length( temp)-1)) tempos <- paste( tempos, temp, ", ", sep="")
	    tempos <- paste( "c( ", tempos, temp[length( temp)], ")", sep="")

	} else tempos <- temp

	cols.selected.cmd <- paste( "cols.selected <- cnames[", tempos, "]", sep="")
	eval(parse(text=cols.selected.cmd))
	write(cols.selected.cmd, file="in2extRemes.log", append=TRUE)

	c.value <- as.numeric(tclvalue(c))
	b.value <- as.numeric(tclvalue(b))
	# cb <- (dd$data[ ,cols.selected]-c.value)/b.value
	cb.cmd <- paste( "cb <- (dd[[\"data\"]][, cols.selected] - ", c.value, ")/", b.value, sep="")
	# cb.cmd <- paste( "tmp <- (", full.list[ data.select], "$data[, cols.selected]-", c.value, ")/", b.value, sep="")
	eval(parse(text=cb.cmd))
	write(cb.cmd, file="in2extRemes.log", append=TRUE)

	newnames <- paste( cnames[temp], ".c", c.value, "b", b.value, sep="")

	for( i in 1:length( newnames)) {

	    cnames.cmd <- paste("cnames <- c( cnames, \"", newnames[i], "\")", sep="")
	    eval(parse(text=cnames.cmd))
	    write(cnames.cmd, file="in2extRemes.log", append=TRUE)

	} # end of for 'i' loop.

	# cnames <- c(cnames, newnames)
	dd.cmd <- "dd[[\"data\"]] <- cbind( dd[[\"data\"]], cb)"
	eval(parse( text=dd.cmd))
	write(dd.cmd, file="in2extRemes.log", append=TRUE)
	colnames.cmd <- "colnames( dd[[\"data\"]]) <- cnames"
	eval(parse( text=colnames.cmd))
	write(colnames.cmd, file="in2extRemes.log", append=TRUE)
	assign.cmd <- paste("assign( \"", full.list[ data.select], "\", dd, pos=\".GlobalEnv\")", sep="")
	# assign.cmd <- paste(full.list[ data.select], "$data <- cbind( ", full.list[ data.select], "$data, tmp)", sep="")
	eval(parse( text=assign.cmd))
	write(assign.cmd, file="in2extRemes.log", append=TRUE)

	# tkconfigure( base.txt, state="normal")
	print(paste( "Affine transformation"))
	cat("\n")

	print( paste( "(x - ", c.value, ")/", b.value, sep=""))
	cat("\n")
	print( paste( "taken for: "))
	for( i in 1:length( cols.selected)) print( paste( cols.selected[i]))
	print( paste( " and assigned to "))
	for( i in 1:length( cols.selected)) print( paste( newnames[i]))
	# tkinsert( base.txt, "end", msg)
	tkdestroy( base)
	# tkconfigure( base.txt, state="disabled")
	invisible()

    } # end of internal 'affine.trans' fcn

    affinehelp <- function() {

	# tkconfigure(base.txt, state="normal")

	print(paste("This is a simple function that takes a data set X and returns the"))
	print(paste("affine transformation of X.  That is"))

	cat("\n")

	print(paste("(X - c)/b"))

	cat("\n")

	print(paste("Returns object of class \"in2extRemesDataObject\" with an extra column in the data"))
	print(paste("indicated by a .cNbM extension where N and M are the values of b and c respectively."))

	# tkinsert( base.txt, "end", help.msg1)
	# tkinsert( base.txt, "end", help.msg2)
	# tkinsert( base.txt, "end", help.msg3)
	# tkconfigure( base.txt, state="disabled")

	invisible()

    } # end of internal 'affinehelp' function.

    endprog <- function() {

	tkdestroy( base)

    } # end of internal 'endprog' function.

    #####################
    # Frame/button setup
    #####################

    base <- tktoplevel()
    tkwm.title(base, "Affine Transformation")

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

    for(i in 1:length(temp)) {

        if(is.null(class(get(temp[i])))) next

        if((class(get( temp[i]))[1] == "in2extRemesDataObject")) {

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

        for( i in 1:ncol( in2extRemesData$data)) tkinsert( col.listbox, "end", paste( colnames( dd$data)[i]))

    } else tkinsert( col.listbox, "end", "")

    tkpack( tklabel( var.frm, text="Variables to Transform", padx=4),
		side="top")
    tkpack( col.listbox, side="left")
    tkpack( col.scroll, side="right")
    tkpack( var.frm, side="top")

    c.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    c.entry <- tkentry( c.frm, textvariable=c, width=4)
    tkpack( tklabel( c.frm, text="(X-c", padx=4), c.entry, side="left")

    b.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    b.entry <- tkentry( b.frm, textvariable=b, width=4)
    tkpack( tklabel( b.frm, text=")/b", padx=4), b.entry, side="left")
    tkpack( c.frm, b.frm, side="left")

    ok.but <- tkbutton( bot.frm, text="OK", command=affine.trans)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)

    help.but <- tkbutton( bot.frm, text="Help", command=affinehelp)

    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")

    # place bindings for return key.
    tkbind( ok.but, "<Return>", affine.trans)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", affinehelp)

    tkpack( top.frm, fill="x")
    tkpack( mid.frm, fill="x")
    tkpack( bot.frm, side="bottom")

} # end of affine.gui fcn
