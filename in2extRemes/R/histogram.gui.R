histogram.gui <- function( base.txt) {

    breaks.val <- tclVar("Sturges")

    # Refresh fcn 
    refresh <- function() {

	tkdelete(fit.listbox, 0.0, "end")

	if(!is.nothing) {

                data.select <- as.numeric( tkcurselection( data.listbox))+1
                dd <- get( full.list[ data.select])

        } else stop("histogram.gui: Must load a data object!")

	models.fit <- names( dd$models)
	for( i in 1:length( models.fit))
		tkinsert( fit.listbox, "end", paste( models.fit[i]))

	invisible()

    } # end of internal 'refresh' fcn.

    submit <- function() {

	if(!is.nothing) {

                data.select <- as.numeric( tkcurselection( data.listbox))+1
                dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")

        } else stop("fitdiag.gui: Must load a data object!")

	eval(parse(text = dd.cmd))
	write(dd.cmd, file="in2extRemes.log", append=TRUE)

	fit.select <- as.numeric(tkcurselection(fit.listbox))+1
	breaks.choices <- c("Sturges", "Scott", "FD")
	breaks.select <- as.numeric(tkcurselection(breaks.listbox)) + 1
	breaks.value <- breaks.choices[ breaks.select]

	histCMD <- paste("plot( dd[[\"models\"]][[ ", fit.select, "]], type = \"hist\"",
			", hist.args = list(breaks = \"", breaks.value, "\", freq = FALSE))", sep="")
	eval(parse(text=histCMD))
	write(histCMD, file="in2extRemes.log", append=TRUE)

	invisible()

    } # end of internal 'submit' fcn.

    histhelp <- function(base.txt) {

	cat("\n", "Histogram with densities for fitted objects from the \'fevd\' function:\n")
	cat("See the help files for \'plot.fevd.mle\', and also for the R function \'hist\'.\n")

	help(hist)

    } # end of internal 'histhelp' fcn

    endprog <- function() {

	tkdestroy(base)

    }
    # Function to plot diagnostic plots for various fits.

    #####################
    # Frame/button setup.
    #####################

    base <- tktoplevel()
    tkwm.title( base, "Histogram Plots for Fitted Objects")

    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")

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

    # Middle frame for choosing which fit to plot and which breaks method to use.
    fit.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    fit.listbox <- tklistbox( fit.frm,
			yscrollcommand=function(...) tkset( fit.scroll, ...),
			selectmode="single",
			width=20,
			height=5,
			exportselection=0)

    fit.scroll <- tkscrollbar( fit.frm, orient="vert", command=function(...) tkyview( fit.listbox, ...))
    tkinsert( fit.listbox, "end", "")

    tkpack( tklabel( fit.frm, text="Select a fit: ", padx=4), side="left")
    tkpack( fit.listbox, fit.scroll, side="left", fill="y")

    breaks.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    breaks.listbox <- tklistbox( breaks.frm,
                        yscrollcommand=function(...) tkset( breaks.scroll, ...),
                        selectmode="single",
                        width=20,
                        height=1,
                        exportselection=0)

    breaks.scroll <- tkscrollbar( breaks.frm, orient="vert",
                        command=function(...) tkyview( breaks.listbox, ...))

    tkinsert( breaks.listbox, "end", "Sturges")
    tkinsert( breaks.listbox, "end", "Scott")
    tkinsert( breaks.listbox, "end", "Friedman-Diaconis")

    tkpack( tklabel( breaks.frm, text="Breaks Algorithm", padx=4), side="left")
    tkpack( breaks.listbox, breaks.scroll, side="left", fill="y")

    tkpack( fit.frm, breaks.frm, side="top")

    # Bottom frame for execution and cancellation.

    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton( bot.frm, text="Help", command=histhelp)

    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")

    # place bindings on buttons.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", histhelp)

    tkpack( top.frm, mid.frm, bot.frm, side="top", fill="x")

    invisible() 

} # end of 'histogram.gui' fcn
