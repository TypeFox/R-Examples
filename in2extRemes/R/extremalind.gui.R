extremalind.gui <- function( base.txt) {

    # Initialize tcl values.

    u.value <- tclVar("")
    conf.value <- tclVar("0.95")
    method.value <- tclVar("intervals")
    r.value <- tclVar("1")

    # Internal functions.

    refresh <- function() {

	tkdelete( var.listbox, 0.0, "end")

	data.select <- as.numeric( tkcurselection( data.listbox))+1
	dd <- get( full.list[ data.select])

	tmp <- colnames(dd$data)
	for(i in 1:length(tmp)) tkinsert(var.listbox, "end", paste( tmp[i]))

	invisible()

    } # end of internal 'refresh' fcn

    submit <- function() {

	if(tclvalue(u.value) == "") u.val <- NULL
	else if(!is.na(as.numeric(tclvalue(u.value)))) u.val <- as.numeric(tclvalue(u.value)) 
	else u.val <- tclvalue(u.value)

	conf.val <- as.numeric(tclvalue(conf.value))

	data.select <- as.numeric(tkcurselection( data.listbox))+1
	dd.cmd <- paste( "dd <- get(\"", full.list[ data.select], "\")", sep="")
	eval(parse(text=dd.cmd))
	write(dd.cmd, file="in2extRemes.log", append=TRUE)

	var.select <- as.numeric(tkcurselection(var.listbox))+1
	xdat.cmd <- paste("xdat <- dd[[\"data\"]][, ", var.select, "]", sep="")
	eval(parse(text = xdat.cmd))
	write( xdat.cmd, file="in2extRemes.log", append=TRUE)

	if(tclvalue(method.value) == "runs") {

	    cmd <- paste("print(ci(extremalindex(x = xdat, threshold = ", u.val, ", method = \"runs\", run.length = ",
			as.numeric(tclvalue(r.value)), ")))", sep="")

	} else {

	    cmd <- paste("print(ci(extremalindex(x = xdat, threshold = ", u.val, ", method =\"intervals\")))", sep = "")

	}

	eval(parse(text=cmd))
	write(cmd, file = "in2extRemes.log", append = TRUE)

	# ci.cmd <- paste("eiCI <- ci(ei, alpha = ", conf.val, ")", sep = "")
	# eval(parse(text = ci.cmd))
	# write(ci.cmd, file = "in2extRemes.log", append = TRUE)

	# printCMD <- "print(eiCI)"
	# eval(parse(text = printCMD))
	# write(printCMD, file = "in2extRemes.log", append = TRUE)

	tkdestroy(base)
	invisible()

    } # end of internal 'submit' fcn

    extindhelp <- function() {

	cat( "\n", paste("Invokes the function: \'extremalindex\'  ",
				"Use \'help(extremalindex)\' for more help.", " ", sep="\n"))
	help(extremalindex)

	invisible()

    }

    endprog <- function() {

	tkdestroy(base)

    }

    #####################
    # Frame/button setup
    #####################

    base <- tktoplevel()
    tkwm.title(base, "Extremal index estimation")

    top.frm <- tkframe(base, borderwidth=2, relief="groove")
    mid.frm <- tkframe(base, borderwidth=2, relief="groove")
    bot.frm <- tkframe(base, borderwidth=2, relief="groove")

    # Top frame to select data object.

    data.listbox <- tklistbox(top.frm,
			yscrollcommand=function(...) tkset( data.scroll,...),
			selectmode="single",
                        width=20,
                        height=5,
                        exportselection=0)

    data.scroll <- tkscrollbar(top.frm, orient="vert",
                        command=function(...) tkyview( data.listbox, ...))

    temp <- ls(all.names=TRUE, name=".GlobalEnv")
    full.list <- character(0)

    for(i in 1:length( temp)) {

        if(is.null(class(get(temp[i])))) next
        if((class(get(temp[i]))[1] == "in2extRemesDataObject")) {

                tkinsert( data.listbox, "end", paste( temp[i]))
        	full.list <- c( full.list, temp[i])

		}

    } # end of for i loop

    tkpack(tklabel(top.frm, text="Data Object", padx=4), side="left")
    tkpack(data.listbox, side="left")
    tkpack(data.scroll, side="right", fill="y")
    tkpack(top.frm)

    # place binding on data.listbox to reflect the chosen data from the list.
    tkbind(data.listbox, "<Button-1>", "")
    tkbind(data.listbox, "<ButtonRelease-1>", refresh)

    # Middle frame to choose which column of data to use and other args.

    midleft <- tkframe(mid.frm, borderwidth=2, relief="groove")

    var.listbox <- tklistbox( midleft,
			yscrollcommand=function(...) tkset( var.scroll, ...),
			selectmode="single",
			width=15,
			height=6,
			exportselection=0)

    var.scroll <- tkscrollbar( midleft, orient="vert",
                        command=function(...) tkyview( var.listbox, ...))

    tkinsert(var.listbox, "end", " ")
    tkpack(tklabel(midleft, text="Select Variable", padx=4), side="top")
    tkpack(var.listbox, side="left")
    tkpack(var.scroll, side="right")

    # Frame for threshold.
    midright <- tkframe(mid.frm, borderwidth=2, relief="groove")
    u.frm <- tkframe(midright, borderwidth=2, relief="flat")
    confFrm <- tkframe(midright, borderwidth=2, relief="flat")
    methodFrm <- tkframe(midright, borderwidth = 2, relief = "flat")
    rFrm <- tkframe(midright, borderwidth = 2, relief = "flat")

    u.entry <- tkentry(u.frm, textvariable=u.value, width=5)
    tkpack(tklabel(u.frm, text="Threshold", padx=4), u.entry, side="left")

    # method.entry <- tkentry(methodFrm, textvariable = method.value, width = 9)
    # tkpack(tklabel(methodFrm, text = "Method", padx = 4), method.entry, side = "left")

    tkpack(tkradiobutton(methodFrm, text = "intervals", value = "intervals", variable = method.value), anchor = "w")
    tkpack(tkradiobutton(methodFrm, text = "runs", value = "runs", variable = method.value), anchor = "w")
    
    r.entry <- tkentry(rFrm, textvariable = r.value, width = 2)
    tkpack(tklabel(rFrm, text = "Run Length", padx = 4), r.entry, side = "left")

    conf.entry <- tkentry(confFrm, textvariable=conf.value, width=5)
    tkpack(tklabel(confFrm, text="Confidence", padx=4), conf.entry, side="left")

    tkpack(u.frm, confFrm, methodFrm, rFrm, side="top")

    tkpack(midleft, midright, side="left")

    # Bottom frame for execution or cancellation.

    ok.but <- tkbutton(bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton(bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton(bot.frm, text="Help", command=extindhelp)

    tkpack(ok.but, cancel.but, side="left")
    tkpack(help.but, side="right")

    tkpack(top.frm, side="top")
    tkpack(mid.frm, fill="x")
    tkpack(bot.frm, side="bottom")

} # end of fcn
