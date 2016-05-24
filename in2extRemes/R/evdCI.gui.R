evdCI.gui <- function( base.txt) {

    # Initialize tcl variables.

    m.value <- tclVar("100")
    rl.lowval.value <- tclVar("")
    rl.upval.value <- tclVar("")
    conf.value <- tclVar("0.95")
    nint.value <- tclVar("20")
    makeplot <- tclVar(1)
    do.rl <- tclVar(1)
    do.par <- tclVar(1)
    R.value <- tclVar(502)

    # Internal functions.

    # Refresh fcn
    refresh <- function() {

        tkdelete(fit.listbox, 0.0, "end")

        if(!is.nothing) {

            data.select <- as.numeric(tkcurselection(data.listbox))+1
            dd <- get(full.list[data.select])

        } else stop("fitdiag.gui: Must load a data object!")

        models.fit <- names( dd$models)
        for(i in 1:length(models.fit)) tkinsert(fit.listbox, "end", paste(models.fit[i]))

        invisible()

    } # end of internal 'refresh' fcn.

#     refresh2 <- function() {
# 
# 	tkdelete(wpar.listbox, 0.0, "end")
# 
# 	if(!is.nothing) {
# 
# 	    wpar.select <- as.numeric(tkcurselection(fit.listbox)) + 1
# 	    fittedobj <- get(full.list[ data.select])
# 	    parnames <- names(fittedobj$results$par)
# 	    for(i in 1:length(parnames)) tkinsert(wpar.listbox, "end", paste(parnames[i]))
# 	}
# 
# 	invisible()
# 
#     } # end of 'refresh2' internal function.

    submit <- function() {

	do.rl.value <- ifelse(as.numeric(tclvalue(do.rl)) == 1, TRUE, FALSE)
        do.par.value <- ifelse(as.numeric(tclvalue(do.par)) == 1, TRUE, FALSE)

	if(!do.rl.value && !do.par.value) stop("You must select at least one of return level or parameter.")

	citype <- ifelse(do.rl.value, "return.level", "parameter")

	# Grab the data object and make sure it has a 'gev.fit' component.
	if( !is.nothing) {

	    data.select <- as.numeric( tkcurselection( data.listbox))+1
	    dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")

	} else dd.cmd <- "dd <- in2extRemesData"

	eval(parse(text=dd.cmd))
	write(dd.cmd, file="in2extRemes.log", append=TRUE)

	fit.cmd <- paste("fit <- dd[[\"models\"]][[", as.numeric(tkcurselection(fit.listbox))+1, "]]", sep="")
	eval(parse(text=fit.cmd))
	write(fit.cmd, file="in2extRemes.log", append=TRUE)

	# Collect inputs for fcn args.
	m.val <- as.numeric(tclvalue(m.value))

	conf.val <- as.numeric(tclvalue(conf.value))
	R.val <- as.numeric(tclvalue(R.value))
	nint.val <- as.numeric(tclvalue(nint.value))

	methval <- as.numeric(tkcurselection(method.listbox)) + 1
	method.val <- c("normal", "boot", "proflik")[methval]

	makeplot2 <- ifelse(as.numeric(tclvalue(makeplot)) == 1, TRUE, FALSE)

	estRLup <- tclvalue(rl.upval.value)
	if(estRLup == "") estRLup <- NULL
	else estRLup <- as.numeric(estRLup)

	estRLdn <- tclvalue(rl.lowval.value)
	if(estRLdn == "") estRLdn <- NULL
	else estRLdn <- as.numeric(estRLdn)

	if(method.val == "proflik") {

            estRLdn <- ifelse(!is.null(estRLdn), estRLdn, "NULL")
            estRLup <- ifelse(!is.null(estRLup), estRLup, "NULL")
	    if(estRLdn != "NULL" && estRLup != "NULL") xrange.val <- paste("c(", estRLdn, ", ", estRLup, ")", sep = "")
	    else xrange.val <- "NULL"

	} else {

	    xrange.val <- "NULL"

	}


	if(do.par.value) {

	    if(method.val != "proflik") {

		if(method.val == "boot") {

	            ci.cmd <- paste("ciout <- ci(fit, alpha = ", 1 - conf.val,
			", type = \"parameter\", R = ", R.val,
			", method = \"", method.val, "\")", sep = "")

		} else {

		    ci.cmd <- paste("ciout <- ci(fit, alpha = ", 1 - conf.val,
			", type = \"parameter\", method = \"", method.val, "\")", sep = "")

		}

	        eval(parse(text=ci.cmd))
                write( ci.cmd, file="in2extRemes.log", append=TRUE)

                printCMD <- "print(ciout)"
                eval(parse(text = printCMD))
	        write(printCMD, file = "in2extRemes.log", append = TRUE)

	    } else {

		np <- length(fit$results$par)

		if(makeplot2) {

		    if(do.rl.value) Np <- np + 1
		    else Np <- np

		    if(Np == 2) par(mfrow = c(2, 1))
		    else if(Np == 3) par(mfrow = c(3,1))
		    else if(Np == 4) par(mfrow = c(2,2))
		    else if(Np < 6) par(mfrow = c(2,3))
		    else if(Np < 9) par(mfrow = c(2,4))
		    else if(Np < 10) par(mfro = c(3, 3))

		    else {

			cat("\n", "May be too many paramters to plot on one device.\n")
			cat("You may want to use pdf/jpeg/png/postscript and dev.off functions before and after plotting\n")
			cat("In order to print them to a pdf, etc. file for subsequent viewing.\n")

		    }
		}

		for(i in 1:np) {

		    ci.cmd <- paste("ciout <- ci(fit, alpha = ", 1 - conf.val,
                        ", type = \"parameter\"",
                        ", which.par = ", i,
                        ", method = \"", method.val, "\", nint = ", nint.val,
                        ", verbose = ", makeplot2, ")", sep = "")

                    eval(parse(text=ci.cmd))
                    write( ci.cmd, file="in2extRemes.log", append=TRUE)

                    printCMD <- "print(ciout)"
                    eval(parse(text = printCMD))
                    write(printCMD, file = "in2extRemes.log", append = TRUE)

		} # end of for 'i' loop.
	    }

	}

	if(do.rl.value) {

            ci.cmd <- paste("ciout <- ci(fit, alpha = ", 1 - conf.val,
                        ", return.period = ", m.val,
                        ", R = ", R.val,
                        ", method = \"", method.val,
                        "\", xrange = ", xrange.val,
                        ", nint = ", nint.val,
                        ", verbose = ", makeplot2, ")", sep = "")

            eval(parse(text=ci.cmd))
            write( ci.cmd, file="in2extRemes.log", append=TRUE)

            printCMD <- "print(ciout)"
            eval(parse(text = printCMD))
            write(printCMD, file = "in2extRemes.log", append = TRUE)

        }


	invisible()

    } # end of internal 'submit' fcn

    endprog <- function() {

	tkdestroy(base)

    }

    cihelp <- function() {

	msg1 <- paste("Estimates confidence intervals for return levels or parameters for fits to EVD distributions.",
			sep="")

	cat(msg1)

	help(ci.fevd.mle)

    } # end of internal 'cihelp' function.

    #####################
    # Frame/button setup.
    #####################

    base <- tktoplevel()
    tkwm.title( base, "Estimated confidence limits for EVD fit parameters and return levels")

    top.frm <- tkframe(base, borderwidth=2, relief="groove")
    data.frm <- tkframe(top.frm, borderwidth=2, relief="groove")
    fit.frm <- tkframe(top.frm, borderwidth=2, relief="groove")
    mid.frm <- tkframe(base, borderwidth=2, relief="groove")
    bot.frm <- tkframe(base, borderwidth=2, relief="groove")

    # Top frame to select data object.

    data.listbox <- tklistbox( data.frm,
			yscrollcommand=function(...) tkset( data.scroll, ...),
			selectmode="single",
			width=20,
			height=5,
			exportselection=0)

    data.scroll <- tkscrollbar(data.frm, orient="vert", command=function(...) tkyview( data.listbox, ...))

    temp <- ls( all.names=TRUE, name=".GlobalEnv")
    full.list <- character(0)
    is.nothing <- TRUE

    for( i in 1:length( temp)) {

	if(is.null(class(get(temp[i])))) next

	if((class(get(temp[i]))[1] == "in2extRemesDataObject")) {

	    tkinsert(data.listbox, "end", paste(temp[i]))
	    full.list <- c(full.list, temp[i])
	    is.nothing <- FALSE

	}

    } # end of for i loop

    tkpack(tklabel(data.frm, text="Data Object", padx=4), side="top")
    tkpack(data.listbox, side="left")
    tkpack(data.scroll, side="right", fill="y")

    # Place bindings on data listbox to update fit listbox.
    tkbind(data.listbox, "<Button-1>", "")
    tkbind(data.listbox, "<ButtonRelease-1>", refresh)
    tkpack(data.frm, fit.frm, side="left")

    # Frame for choosing which fit to plot.
    fit.listbox <- tklistbox(fit.frm,
                        yscrollcommand=function(...) tkset(fit.scroll, ...),
                        selectmode="single",
                        width=20,
                        height=5,
                        exportselection=0)

    fit.scroll <- tkscrollbar(fit.frm, orient="vert", command=function(...) tkyview( fit.listbox, ...))
    tkinsert(fit.listbox, "end", "")

    tkpack(tklabel(fit.frm, text="Select a fit: ", padx=4), side="top")
    tkpack(fit.listbox, fit.scroll, side="left", fill="y")
    tkpack(data.frm, fit.frm, fill="y", side="left")

    # Middle frame to enter arguments for 'ci' fcn.

    # Frame for 'lowval' and 'upval' args.
    options.frm <- tkframe(mid.frm, borderwidth=2, relief="groove")
    method.frm <- tkframe(options.frm, borderwidth=2, relief="flat")
    proflik.frm <- tkframe(mid.frm, borderwidth=2, relief="groove")
    valrange.frm <- tkframe(proflik.frm, borderwidth=2, relief="groove")

    # Frame for m-year return level.
    m.frm <- tkframe(options.frm, borderwidth=2, relief="groove")
    m.entry <- tkentry(m.frm, textvariable=m.value, width=5)

    rl.lowval.frm <- tkframe( valrange.frm, borderwidth=2, relief="flat")
    rl.lowval.entry <- tkentry( rl.lowval.frm, textvariable=rl.lowval.value, width=5)

    rl.upval.frm <- tkframe(valrange.frm, borderwidth=2, relief="flat")
    rl.upval.entry <- tkentry(rl.upval.frm, textvariable=rl.upval.value, width=5)

    conf.frm <- tkframe(options.frm, borderwidth=2, relief="flat")
    conf.entry <- tkentry(conf.frm, textvariable=conf.value, width=5)

    # CI method listbox.

    tkpack(tklabel(options.frm, text = "Options", padx = 4), side = "top")

    method.listbox <- tklistbox(method.frm, yscrollcommand = function(...) tkset(method.scroll, ...),
			selectmode = "single", width = 20, height = 3, exportselection = 0)
    method.scroll <- tkscrollbar(method.frm, orient = "vert", command = function(...) tkyview(method.listbox, ...))

    # tkinsert(method.listbox, "end", "")
    tkpack(tklabel(method.frm, text="Method ", padx=4), side="top")
    tkinsert(method.listbox, "end", "normal approximation")
    tkinsert(method.listbox, "end", "parametric bootstrap")
    tkinsert(method.listbox, "end", "profile likelihood")

    tkpack(method.listbox, method.scroll, side = "left", fill = "y")
    tkpack(method.frm, side = "left")
    # tkbind(method.listbox, "<Button-1>", "")
    # tkbind(method.listbox, "<ButtonRelease-1>", submit)
   
    R.frm <- tkframe(mid.frm, borderwidth = 2, relief = "flat")
    R.entry <- tkentry(R.frm, textvariable = R.value, width = 5) 

    nint.frm <- tkframe(proflik.frm, borderwidth=2, relief="flat")
    nint.entry <- tkentry(nint.frm, textvariable=nint.value, width=5)

    tkpack(tklabel(m.frm, text="return period", padx=4), m.entry, side="left", fill="x", anchor="w")
    tkpack(m.frm)

    # Buttons to indicate which parameters to find CI's for.
    do.frm <- tkframe(options.frm, borderwidth=2, relief="flat")
    rl.button <- tkcheckbutton(do.frm, text="Return Level", variable=do.rl)
    tkpack(rl.button, side="left")

    par.button <- tkcheckbutton(do.frm, text="Parameter", variable=do.par)
    tkpack(par.button, side="left")

    # wpar.listbox <- tklistbox(do.frm, yscrollcommand = function(...) tkset(wpar.scroll, ...),
# 			selectmode = "single", width = 10, height = 5, exportselection = 0)
 #    wpar.scroll <- tkscrollbar(do.frm, orient = "vert", command = function(...) tkyview(wpar.listbox, ...))
  #   tkinsert(wpar.listbox, "end", "")

   #  tkpack(wpar.listbox, wpar.scroll, side = "left", fill = "y")
    # tkbind(fit.listbox, "<Button-1>", "")
    # tkbind(fit.listbox, "<ButtonRelease-1>", refresh)
    tkpack(do.frm)

    tkpack(tklabel(rl.lowval.frm, text="Lower limit", padx=4), rl.lowval.entry, side="left")
    tkpack(tklabel(rl.upval.frm, text="Upper limit", padx=4), rl.upval.entry, side="left")
    tkpack(rl.upval.frm, side="bottom")
    tkpack(rl.lowval.frm, side="bottom")

    # Pack together the search ranges for both return level and shape parameter.
    tkpack(tklabel(proflik.frm, text = "Profile Likelihood Options", padx = 4), side = "top")
    tkpack(tklabel(valrange.frm, text="Profile Search Range", padx=4),
    tklabel(valrange.frm, text="(leave blank to try to find automatically)", padx=4), side="top")
    tkpack(valrange.frm, side="left")

    tkpack(tklabel(conf.frm, text="Confidence Value (* 100%)", padx=4), conf.entry, side="left", fill="x", anchor="w")
    tkpack(conf.frm)

    tkpack(tklabel(nint.frm, text="Number of points at which to calculate profile likelihood", padx=4),
	nint.entry, side="left", fill="x")
    tkpack(nint.frm)

    # create check button for plotting profile likelihoods (or not)

    makeplot.but <- tkcheckbutton(proflik.frm,text="Plot profile likelihoods",variable=makeplot)
    tkpack(makeplot.but, side="left")

    tkpack(options.frm, proflik.frm, side = "top")

    tkpack(tklabel(R.frm, text = "Parametric Bootstrap Options", padx = 4), side = "top")
    tkpack(tklabel(R.frm, text = "Number of Bootstrap Replicates", padx = 4), R.entry, side = "left", fill="x", anchor="w")
    tkpack(R.frm)

    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton( bot.frm, text="Help", command=cihelp)

    tkpack(ok.but, cancel.but, side="left")
    tkpack(help.but, side="right")

    tkbind(ok.but, "<Return>", submit)
    tkbind(cancel.but, "<Return>", endprog)
    tkbind(help.but, "<Return>", cihelp)

    tkpack(top.frm, side="top")
    tkpack(mid.frm, fill="x")
    tkpack(bot.frm, side="bottom")

} # end of 'evdCI.gui' fcn
