declusterGUI <- function(base.txt) {

    # This function provides a gui interface for declustering a dataset.
    # The gui will list all objects of class
    # "in2extRemesDataObject" and the column names of the user selected object.
    # After declustering the selected data, will return a new column
    # to the data.  
    # Declustered data will keep the same column name as the original
    # clustered column, but with a ".dc" extension.

    # Set the tcl variables

    threshold.value <- tclVar("")
    r.value <- tclVar("1")
    add.plot.value <- tclVar(0)

    # Internal functions

    refresh <- function() {

        # When data is selected, this function fills the lists for the columns of
        # the data set so that the user can select which column(s) to decluster.

	tkdelete(col.listbox, 0.0, "end")
	tkdelete(clustby.listbox, 0.0, "end")

	if( !is.nothing) {

	    data.select <- as.numeric(tkcurselection(data.listbox))+1
	    dd <- get(full.list[data.select])

	} else dd <- in2extRemesData

	for(i in 1:ncol(dd$data)) {

	    tkinsert(col.listbox, "end", paste( colnames( dd$data)[i])) 
	    tkinsert(clustby.listbox, "end", paste( colnames( dd$data)[i]))

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

        M.cmd <- "M <- dim( dd[[\"data\"]])[1]"
        eval(parse(text=M.cmd))
        write(M.cmd, file="in2extRemes.log", append=TRUE)

        # Grab data to decluster.
	# cols.selected.cmd <- "cols.selected <- character(0)"
	# eval(parse(text=cols.selected.cmd))
	# write(cols.selected.cmd, file="in2extRemes.log", append=TRUE)

	# Gather column names to add decluster to end.
	cnames <- colnames(dd[["data"]])

	temp <- as.numeric( tkcurselection( col.listbox)) + 1
	# cols.selected.cmd <- paste( "cols.selected <- c( cols.selected, \"", cnames[temp], "\")", sep="")
	# eval(parse( text=cols.selected.cmd))
	# write(cols.selected.cmd, file="in2extRemes.log", append=TRUE)

	# Make sure a column has been selected.
	if(is.na(temp)) return()

	# Grab column to decluster by...
	temp2 <- as.numeric(tkcurselection(clustby.listbox)) + 1
	if(length(temp2) > 0) cluster.by.val <- cnames[temp2]
	else cluster.by.val <- NULL

	if(!is.null(cluster.by.val)) {

	    dcol <- cnames[temp]

	    cmd <- paste("datf <- dd[[\"data\"]]", sep = "")
	    eval(parse(text=cmd))
            write(cmd, file="in2extRemes.log", append=TRUE)

	}
	    
	tmp.data.cmd <- paste("tmp.data <- dd[[\"data\"]][[\"", cnames[temp], "\"]]", sep = "")
	eval(parse(text=tmp.data.cmd))
	write(tmp.data.cmd, file="in2extRemes.log", append=TRUE)

	r.val <- as.numeric(tclvalue(r.value))

        # Obtain the threshold values (is it a number or a name?).
	if(!is.na(as.numeric(tclvalue(threshold.value)))) {

	    threshold.val <- as.numeric( tclvalue( threshold.value))
	    thresh.name <- as.character( threshold.val)
	    if(is.null(cluster.by.val)) newnames <- paste( cnames[temp], ".u",thresh.name, "r", r.val, "dc", sep="")
	    else newnames <- paste( cnames[temp], ".u",thresh.name, "r", r.val, "dc", "by", cnames[temp2], sep="")

	} else { 

	    threshold.val <- tclvalue(threshold.value)
	    thresh.name <- tclvalue( threshold.value)
	    if(is.null(cluster.by.val)) newnames <- paste( cnames[temp], ".u",thresh.name, "r", r.val, "dc", sep="")
	    else newnames <- paste(cnames[temp], ".u",thresh.name, "r", r.val, "dc", "by", cnames[temp2], sep="")

	} # end of if else 'threshold.value' is a number or a name stmts.


	## Make sure this new decluster doesn't delete an existing column.
	## The assumption is that if the column exists, then the procedure
	## has already been performed, but this may not be the case.

	for(i in 1:length(cnames)) {

            if(newnames == cnames[i]) {

		warn.msg <- paste(" ", "***************", "Warning: ",
				"This declustering procedure appears to have already been performed.",
				"No declustering applied.  Maybe change column name.", "***************", " ", sep="\n")
                stop(warn.msg)

            } # end of if stmt

        } # end of for i loop

	print( paste( "Declustering ...", sep=""))

	if(!is.null(cluster.by.val)) out.cmd <- paste("out <- decluster(datf, threshold = ", threshold.val,
			", which.cols = c(\"", dcol, "\", \"", cluster.by.val, "\"), r = ", r.val, ")", sep="")
	else out.cmd <- paste("out <- decluster(x = tmp.data, threshold = ", threshold.val, ", r = ", r.val, ")", sep="")
	eval(parse(text = out.cmd))
	write(out.cmd, file = "in2extRemes.log", append = TRUE)

	tmp.dc.cmd <- "tmp.dc <- c(out)"
	eval(parse(text=tmp.dc.cmd))
	write(tmp.dc.cmd, file="in2extRemes.log", append = TRUE)

	if(tclvalue(add.plot.value)==1) {

	    plotCMD <- "plot(out)"
	    eval(parse(text=plotCMD))
	    write(plotCMD, file="in2extRemes.log", append=TRUE)

	}

	eiCMD <- paste("print(extremalindex(tmp.data, threshold = ", threshold.val, ", run.length = ", r.val, ", method = \"runs\"))", sep="")
	eval(parse(text=eiCMD))
	write(eiCMD, file = "in2extRemes.log", append = TRUE)

	eiCMD <- paste("print(extremalindex(tmp.data, threshold = ", threshold.val, "))", sep = "")
	eval(parse(text=eiCMD))
	write(eiCMD, file = "in2extRemes.log", append = TRUE)

        cnames.cmd <- paste( "cnames <- c( cnames, \"", newnames, "\")", sep="")
        eval(parse(text=cnames.cmd))
        write(cnames.cmd, file="in2extRemes.log", append=TRUE)

        dd.cmd <- "dd[[\"data\"]] <- cbind( dd[[\"data\"]], tmp.dc)"
        eval(parse( text=dd.cmd))
        write(dd.cmd, file="in2extRemes.log", append=TRUE)

        colnames.cmd <- "colnames( dd[[\"data\"]]) <- cnames"
        eval(parse(text=colnames.cmd))
        write(colnames.cmd, file="in2extRemes.log", append=TRUE)

        assignCMD <- paste("assign( \"", full.list[data.select], "\", dd, pos=\".GlobalEnv\")", sep="")
        eval(parse(text=assignCMD))
        write(assignCMD, file="in2extRemes.log", append=TRUE)

        msg1 <- paste("declustering performed for:", sep="")
        msg2 <- paste(cnames[temp], " and assigned to ", newnames, sep="")

        print(msg1)
        print(msg2)

        invisible()

    } # end of internal 'submit' fcn

    dchelp <- function() {

	# tkconfigure(base.txt, state="normal")
	nl1 <- paste("*******************", sep="")

	help.msg <- paste( " ", "Runs declustering", "",
		"This is a simple function that takes a data set X and declusters",
		"X based on a given threshold and the number of subsequent obs that",
		"fall below the threshold.  That is, once an observation exceeds a", 
		"threshold (say u) a cluster begins.  After r observations fall below u,",
		"that cluster terminates and a new one begins.  The",
		"maximum of each of these clusters is then returned to the data object",
		"with a new data column terminating with a .u[u]r[r]dc extension.  The rest of",
		"the new column is filled with the lesser of the minimum of the original data",
		"or one less than the chosen threshold for dimensional compatibility.",
		"","Also displays estimates for the extremal index (for both original and",
		"declustered data) based on both n_c/N, where N is the number of exceedances",
		"of the threshold and n_c are the number of clusters found from runs declustering,",
		"and on the intervals estimator from Ferro and Segers (2003) (see the help file for",
		"extremalindex).", "",
		"See Coles (2001) for more information on runs declustering.",
		"See the help file for \'decluster\' (i.e., \'help(decluster)\') as this is the primary ",
		"function utilized.", " ", sep="\n")

		cat(nl1)
		cat(help.msg)

		invisible()

	} # end of internal 'lthelp' function.

    endprog <- function() {

	tkdestroy(base)

    }

    #####################
    # Frame/button setup
    #####################

    base <- tktoplevel()
    tkwm.title( base, "Declustering")

    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    data.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
    clustby.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
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

    tkpack( tklabel( data.frm, text="Variable to Decluster", padx=4), side="top")
    # tkpack( col.listbox, side="left")
    # tkpack( col.scroll, side="right")
    tkpack( col.listbox, col.scroll, side="left", fill="y")

    # cluster by frame...
    clustby.listbox <- tklistbox( clustby.frm,
                        yscrollcommand=function(...) tkset(clustby.scroll, ...),
                        selectmode="single",
                        width=20,
                        height=5,
                        exportselection=0)
    clustby.scroll <- tkscrollbar( clustby.frm, orient="vert",
                        command=function(...) tkyview( clustby.listbox, ...))
    if(is.nothing) {

	for( i in 1:ncol( in2extRemesData$data)) tkinsert( clustby.listbox, "end", paste( colnames( dd$data)[i]))

    } else tkinsert( col.listbox, "end", "")

    tkpack(tklabel(clustby.frm, text="Decluster by", padx=4), side="top")
    tkpack(clustby.listbox, clustby.scroll, side="left", fill="y")

    tkpack(data.frm, clustby.frm, side="left", fill="x")

    # Select a threshold.
    rightmid.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")

    threshold.frm <- tkframe( rightmid.frm, borderwidth=2, relief="flat")
    threshold.entry <- tkentry( threshold.frm, textvariable=threshold.value,
			width=10)
    tkpack( threshold.entry, tklabel( threshold.frm, text="Threshold(s)", padx=4), side="right")

    # Select r.
    r.frm <- tkframe( rightmid.frm, borderwidth=2, relief="groove")
    r.entry <- tkentry( r.frm, textvariable=r.value, width=4)
    tkpack(r.entry, tklabel(r.frm, text="Run Length", padx=4), side="right")

    # Add plot checkbutton...
    add.plot.checkbutton <- tkcheckbutton( rightmid.frm, text="Plot data",
				variable=add.plot.value)
    tkpack(threshold.frm, r.frm, add.plot.checkbutton, side="top")
    tkpack(rightmid.frm, fill="x")

    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)

    help.but <- tkbutton( bot.frm, text="Help", command=dchelp)

    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")

    # place bindings for return key.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", dchelp)

    tkpack( top.frm, fill="x")
    tkpack( mid.frm, fill="x")
    tkpack( bot.frm, side="bottom")

} # end of declusterGUI fcn
