fitsummary.gui <- function(base.txt) {

    closevar <- tclVar(1)

    # Refresh fcn 
    refresh <- function() {

	tkdelete( fit.listbox, 0.0, "end")

	if(!is.nothing) {

                data.select <- as.numeric(tkcurselection(data.listbox))+1
                dd <- get(full.list[ data.select])

        } else stop("fitdiag.gui: Must load a data object!")

	models.fit <- names(dd$models)

	for(i in 1:length(models.fit)) tkinsert(fit.listbox, "end", paste(models.fit[i]))

	invisible()

    } # end of internal 'refresh' fcn.

    submit <- function() {

	if(!is.nothing) {

                data.select <- as.numeric( tkcurselection( data.listbox))+1
                dd.cmd <- paste("dd <- get( \"", full.list[ data.select], "\")", sep="")
		eval(parse(text = dd.cmd))
		write(dd.cmd, file = "in2extRemes.log", append = TRUE)

        } else stop("fitsummary.gui: Must load a data object!")

	fit.select <- as.numeric(tkcurselection(fit.listbox))+1
	z.cmd <- paste("z <- dd[[\"models\"]][[ ", fit.select, "]]", sep="")
	eval(parse(text = z.cmd))
	write(z.cmd, file="in2extRemes.log", append=TRUE)
	summaryCMD <- "print(summary(z))"
	eval(parse(text = summaryCMD))
	write(summaryCMD, file = "in2extRemes.log", append=TRUE)

	if(as.numeric(tclvalue(closevar))==1) endprog()

	invisible()

    } # end of internal 'submit' fcn.

    summaryhelp <- function( base.txt) {

	cat("\n", "prints summary information about an object to the console.", "\n")
	cat("Use \'help(summary.fevd)\' for more information.\n")
	help(summary.fevd)

    } # end of summaryhelp fcn

    endprog <- function() {

	tkdestroy(base)

    }

    # Function to give fit summaries for various fits.

    #####################
    # Frame/button setup.
    #####################

    base <- tktoplevel()
    tkwm.title(base, "Fitted Objects Summaries")

    top.frm <- tkframe(base, borderwidth=2, relief="groove")
    mid.frm <- tkframe(base, borderwidth=2, relief="groove")
    lowermidfrm <- tkframe(base, borderwidth=2, relief="groove")
    bot.frm <- tkframe(base, borderwidth=2, relief="groove")

    data.listbox <- tklistbox(top.frm,
                        yscrollcommand=function(...) tkset( data.scroll, ...),
                        selectmode="single",
                        width=20,
                        height=5,
                        exportselection=0)

    data.scroll <- tkscrollbar(top.frm, orient="vert", command=function(...) tkyview( data.listbox, ...))

    temp <- ls(all.names=TRUE, name=".GlobalEnv")
    full.list <- character(0)
    is.nothing <- TRUE

    for(i in 1:length(temp)) {

        if(is.null(class(get(temp[i])))) next

        if((class(get(temp[i]))[1] == "in2extRemesDataObject")) {

                tkinsert(data.listbox, "end", paste( temp[i]))
                full.list <- c(full.list, temp[i])
                is.nothing <- FALSE

        }

    } # end of for i loop

    tkpack(tklabel(top.frm, text="Data Object:  ", padx=4), side="left")
    tkpack(data.listbox, data.scroll,  side="left", fill="y")

    # Place bindings on data listbox to update fit listbox.
    tkbind(data.listbox, "<Button-1>", "")
    tkbind(data.listbox, "<ButtonRelease-1>", refresh)

    # Middle frame for choosing which fit to plot.
    fit.listbox <- tklistbox( mid.frm,
			yscrollcommand=function(...) tkset( fit.scroll, ...),
			selectmode="single",
			width=20,
			height=5,
			exportselection=0)

    fit.scroll <- tkscrollbar( mid.frm, orient="vert",
			command=function(...) tkyview( fit.listbox, ...))
    tkinsert(fit.listbox, "end", "")

    tkpack(tklabel( mid.frm, text="Select a fit: ", padx=4), side="left")
    tkpack(fit.listbox, fit.scroll, side="left", fill="y")

    # Close on execute?
    close.but <- tkcheckbutton(lowermidfrm, text="Close window on execution", variable=closevar)
    tkpack(close.but, fill="x")

    # Bottom frame for execution and cancellation.

    ok.but <- tkbutton(bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton(bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton(bot.frm, text="Help", command=summaryhelp)

    tkpack(ok.but, cancel.but, side="left")
    tkpack(help.but, side="right")

    # place bindings on buttons.
    tkbind(ok.but, "<Return>", submit)
    tkbind(cancel.but, "<Return>", endprog)
    tkbind(help.but, "<Return>", summaryhelp)

    tkpack(top.frm, mid.frm, lowermidfrm, bot.frm, side="top", fill="x")

    invisible()


} # end of 'fitsummary.gui' fcn
