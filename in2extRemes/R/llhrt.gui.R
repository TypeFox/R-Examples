llhrt.gui <- function( base.txt) {
    
    alpha.val <- tclVar("0.05")
    
    # Refresh fcn 
    refresh <- function() {

    	tkdelete( fit1.listbox, 0.0, "end")
    	tkdelete( fit2.listbox, 0.0, "end")

    	if( !is.nothing) {

            data.select <- as.numeric( tkcurselection( data.listbox))+1
            dd <- get( full.list[ data.select])

        } else stop("llhrt.gui: Must load a data object!")

    	models.fit <- names( dd$models)

    	for( i in 1:length( models.fit)) {

    	    tkinsert( fit1.listbox, "end", paste( models.fit[i]))
    	    tkinsert( fit2.listbox, "end", paste( models.fit[i]))

    	}

    	invisible()

    } # end of internal 'refresh' fcn.
    
    submit <- function() {

    	alpha <- as.numeric( tclvalue( alpha.val))

    	if( !is.nothing) {

            data.select <- as.numeric( tkcurselection( data.listbox))+1
            dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")

        } else stop("llhrt.gui: Must load a data object!")
    	eval( parse( text=dd.cmd))
    	write( dd.cmd, file="in2extRemes.log", append=TRUE)
    
    	fit1.select <- as.numeric( tkcurselection( fit1.listbox))+1
    	fit2.select <- as.numeric( tkcurselection( fit2.listbox))+1

    	m0.cmd <- paste( "m0 <- dd[[\"models\"]][[ ", fit1.select, "]]", sep="")
    	eval( parse( text=m0.cmd))
    	write( m0.cmd, file="in2extRemes.log", append=TRUE)

    	m1.cmd <- paste( "m1 <- dd[[\"models\"]][[ ", fit2.select, "]]", sep="")
    	eval( parse( text=m1.cmd))
    	write( m1.cmd, file="in2extRemes.log", append=TRUE)

    	out.cmd <- paste( "out <- lr.test(x = m0, y = m1, alpha = ", alpha, ")", sep="")
    	eval( parse( text=out.cmd))
    	write( out.cmd, file="in2extRemes.log", append=TRUE)

	printCMD <- "print(out)"
	eval(parse(text = printCMD))
	write(printCMD, file = "in2extRemes.log", append = TRUE)

    	tkyview.moveto( base.txt, 1.0)

    	invisible()

    } # end of submit fcn.
    
    devhelp <- function() {

    	nl1 <- paste(" ", "**************", " ", sep="\n")
            nl2 <- paste(" ", " ", sep="\n")
    	cat( nl1)
    	h1 <- paste( "Simply computes the likelihood-ratio test", "D=2*log(M1/M0)",
    			"where M0 contained in M1 are the likelihood functions.",
    			"Also computed is the 1-alpha quantile of",
    			"a Chi-squared distribution with degrees of freedom",
    			"equal to the difference in the number of parameters of M0 and M1.", sep="\n")
    	cat( h1)
    	cat( nl2)
    	cat( nl1)

    	tkyview.moveto( base.txt, 1.0)

    	invisible()

    } # end of internal 'devhelp' fcn
    
    endprog <- function() {

    	tkdestroy(base)

    }
    # Function to compare two fits M0 and M1, where M0 is contained in M1.
    
    #####################
    # Frame/button setup.
    #####################
    
    base <- tktoplevel()
    tkwm.title( base, "Likelihood-ratio test for comparing nested fits")
    
    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    d.frm <- tkframe( top.frm, borderwidth=2, relief="groove")
    alpha.frm <- tkframe( top.frm, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")
    
    data.listbox <- tklistbox( d.frm,
                            yscrollcommand=function(...) tkset( data.scroll, ...),
                            selectmode="single",
                            width=20,
                            height=5,
                            exportselection=0)
    
    data.scroll <- tkscrollbar( d.frm, orient="vert",
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
    
    tkpack( tklabel( d.frm, text="Data Object:  ", padx=4), side="top")
    tkpack( data.listbox, data.scroll,  side="left", fill="y")
    
    # Place bindings on data listbox to update fit listbox.
    tkbind( data.listbox, "<Button-1>", "")
    tkbind( data.listbox, "<ButtonRelease-1>", refresh)
    
    # alpha frame
    alpha.entry <- tkentry( alpha.frm, textvariable=alpha.val, width=5)
    tkpack( tklabel( alpha.frm, text="significance level (alpha)", padx=4), alpha.entry, side="left")
    tkpack( d.frm, alpha.frm, side="left")
    
    # Middle frame for choosing which fits to compare.
    midleft.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    fit1.listbox <- tklistbox( midleft.frm,
    			yscrollcommand=function(...) tkset( fit1.scroll, ...),
    			selectmode="single",
    			width=20,
    			height=5,
    			exportselection=0)
    
    fit1.scroll <- tkscrollbar( midleft.frm, orient="vert",
    			command=function(...) tkyview( fit1.listbox, ...))
    tkinsert( fit1.listbox, "end", "")
    
    tkpack( tklabel( midleft.frm, text="Select base fit (M0): ", padx=4), side="top")
    tkpack( fit1.listbox, fit1.scroll, side="left", fill="y")
    
    midright.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    fit2.listbox <- tklistbox( midright.frm,
                            yscrollcommand=function(...) tkset( fit2.scroll, ...),
                            selectmode="single",
                            width=20,
                            height=5,
                            exportselection=0)
    
    fit2.scroll <- tkscrollbar( midright.frm, orient="vert",
                            command=function(...) tkyview( fit2.listbox, ...))
    tkinsert( fit2.listbox, "end", "")
    
    tkpack( tklabel( midright.frm, text="Select comparison fit (M1): ", padx=4), side="top")
    tkpack( fit2.listbox, fit2.scroll, side="left", fill="y")
    
    tkpack( midleft.frm, midright.frm, side="left", fill="x")
    
    # Bottom frame for execution and cancellation.
    
    ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton( bot.frm, text="Help", command=devhelp)
    
    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")
    
    # place bindings on buttons.
    tkbind( ok.but, "<Return>", submit)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", devhelp)
    
    tkpack( top.frm, mid.frm, bot.frm, side="top", fill="x")
    invisible()
    
    
} # end of fcn
