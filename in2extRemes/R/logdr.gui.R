logdr.gui <- function( base.txt) {
    
    # This function provides a gui interface for finding the log daily returns
    # transformation of a dataset.  The gui will list all objects of class
    # "in2extRemesDataObject" and the column names of the user selected object.
    # After taking the log daily returns of the selected data, will return a new
    # column with the same column name(s), but with a ".ldr" extension.
    
    # Set the tcl variables
    
    ebase <- tclVar(1)
    bb <- tclVar( 10)
    
    # Internal functions
    
    refresh <- function() {
    
        # When data is selected, this function fills the lists for the columns of
        # the data set so that the user can select which column(s) to transform.
    
    	tkdelete( col.listbox, 0.0, "end")
    
    	if( !is.nothing) {

    	    dd.select <- as.numeric( tkcurselection( data.listbox))+1
    	    dd <- get( full.list[ dd.select])

    	} else dd <- in2extRemesData
    
    	for( i in 1:ncol( dd$data)) tkinsert( col.listbox, "end", paste( colnames( dd$data)[i]))

    	invisible()

    	} # end of internal 'refresh' fcn
    
    ldr.trans <- function() {
        
            # Function invoked when the "ok" button is pressed.  Actually takes log
            # transformation of the data.
            
            if( !is.nothing) {

            	data.select <- as.numeric( tkcurselection( data.listbox))+1
            	dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")
            	ddname <- full.list[ data.select]

            } else dd.cmd <- "dd <- in2extRemesData"
            	eval( parse( text=dd.cmd))
            	write( dd.cmd, file="in2extRemes.log", append=TRUE)
            
            if( tclvalue(ebase) == 1 ) log.b <- exp(1)
            else {
            	log.b <- as.numeric( tclvalue( bb))
            	bbext <- deparse( log.b)
            	}
            
            cols.selected.cmd <- "cols.selected <- character(0)"
            eval( parse( text=cols.selected.cmd))
            write( cols.selected.cmd, file="in2extRemes.log", append=TRUE)
            
            temp <- as.numeric( tkcurselection( col.listbox)) + 1
            
            # Make sure a column has been selected.
            if( is.na( temp)) return()
            
            # Keep column names for writing new ones.
            cnames.cmd <- "cnames <- colnames( dd[[\"data\"]])"
            eval( parse( text=cnames.cmd))
            write( cnames.cmd, file="in2extRemes.log", append=TRUE)
            
            for( i in 1:length( temp)) {

            	cols.selected.cmd <- paste(" cols.selected <- c( cols.selected, \"", cnames[ temp[i]], "\")", sep="")
            	eval( parse( text=cols.selected.cmd))
            	write( cols.selected.cmd, file="in2extRemes.log", append=TRUE)

            } # end of for 'i' loop.

            N.cmd <- "N <- dim( dd[[\"data\"]])[1]"
            eval( parse( text=N.cmd))
            write( N.cmd, file="in2extRemes.log", append=TRUE)
            
            tmp.ldr.cmd <- "tmp.ldr <- matrix( NA, nrow=N, ncol=length( cols.selected) )"
            eval( parse( text=tmp.ldr.cmd))
            write( tmp.ldr.cmd, file="in2extRemes.log", append=TRUE)
            
            tmp.ldr.cmd <- paste( "tmp.ldr[ 1:(N-1),] <- log( dd[[\"data\"]][2:N, cols.selected], base=", log.b, ") -",
            			" log( dd[[\"data\"]][1:(N-1), cols.selected], base=", log.b, ")", sep="")
            eval( parse( text=tmp.ldr.cmd))
            write( tmp.ldr.cmd, file="in2extRemes.log", append=TRUE)
            
            # Fill in last place (and anywhere else that there is an NA
            # with the average value so that there is not a missing value.
            CMD <- "for( i in 1: dim( tmp.ldr)[2]) tmp.ldr[is.na( tmp.ldr[,i]),i] <- mean(tmp.ldr[,i],na.rm=TRUE)"
            eval( parse( text=CMD))
            write( CMD, file="in2extRemes.log", append=TRUE)
            
            if( tclvalue(ebase) == 1) newnames <- paste( cnames[temp], ".ldr", sep="")
            else newnames <- paste( cnames[temp], ".", bbext, "ldr", sep="")

            for( i in 1:length( newnames)) {

            	# cnames <- c( cnames, newnames)
            	cnames.cmd <- paste( "cnames <- c( cnames, \"", newnames[i], "\")", sep="")
            	eval( parse( text=cnames.cmd))
            	write( cnames.cmd, file="in2extRemes.log", append=TRUE)

            } # end of for 'i' loop.
            
            newdata.cmd <- "dd[[\"data\"]] <- cbind( dd[[\"data\"]], tmp.ldr)"
            eval( parse( text=newdata.cmd))
            write( newdata.cmd, file="in2extRemes.log", append=TRUE)

            colnamesCMD <- "colnames( dd[[\"data\"]]) <- cnames"
            eval( parse( text=colnamesCMD))
            write( colnamesCMD, file="in2extRemes.log", append=TRUE)

            assignCMD <- paste( "assign( \"", full.list[ data.select], "\", dd, pos=\".GlobalEnv\")", sep="")
            eval( parse( text=assignCMD))
            write( assignCMD, file="in2extRemes.log", append=TRUE)
            
            # tkconfigure( base.txt, state="normal")
            msg <- paste( "\n", "log daily return (base= ", round( log.b, digits=6), 
            		") taken for", "\n",
            		" ", cols.selected, " and assigned to ", newnames, sep="")
            # tkinsert( base.txt, "end", msg)
            cat( msg)
            tkdestroy( base)
            # tkconfigure( base.txt, state="disabled")

            invisible()

    	} # end of internal 'ldr.trans' fcn
    
    ldrhelp <- function() {
    	# tkconfigure( base.txt, state="normal")
    	help.msg1 <- paste( " ",
    "This is a simple function that takes an n X 1 data set Y and computes the ",
    	"log (base bb) daily return.  That is:", " ", sep="\n")
    	help.msg2 <- paste("  ",
    "log( Y[2:n], base=bb) - log( Y[1:(n-1)], base=bb)", " ", sep="\n")
    	help.msg3 <- paste(" ",
    "Returns an object of class \"in2extRemesDataObject\" ", " ", sep="\n")
    	# tkinsert( base.txt, "end", help.msg1)
    	cat( help.msg1)
    	# tkinsert( base.txt, "end", help.msg2)
    	cat( help.msg2)
    	# tkinsert( base.txt, "end", help.msg3)
    	cat( help.msg3)
    	# tkconfigure( base.txt, state="disabled")
    	invisible()
    	} # end of ldrhelp fcn
    
    endprog <- function() {
    	tkdestroy( base)
    	}
    
    #####################
    # Frame/button setup
    #####################
    
    base <- tktoplevel()
    tkwm.title( base, "Log Daily Returns Transformation")
    
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
    tkpack( top.frm)
    tkbind( data.listbox, "<Button-1>", "")
    tkbind( data.listbox, "<ButtonRelease-1>", refresh)
    
    # Middle frame for columns of chosen dataset...
    
    leftmid.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
    
    col.listbox <- tklistbox( leftmid.frm,
    			yscrollcommand=function(...) tkset( col.scroll, ...),
    			selectmode="multiple",
    			width=20,
    			height=5,
    			exportselection=0)
    
    col.scroll <- tkscrollbar( leftmid.frm, orient="vert",
    			command=function(...) tkyview( col.listbox, ...))
    
    if( is.nothing) {
    for( i in 1:ncol( in2extRemesData$data))
    	tkinsert( col.listbox, "end", paste( colnames( dd$data)[i]))
    # end of for i loop
    	} else tkinsert( col.listbox, "end", "")
    
    tkpack( tklabel( leftmid.frm, text="Variables to Transform", padx=4),
    		side="top")
    tkpack( col.listbox, side="left")
    tkpack( col.scroll, side="right")
    
    # Choose which base if other than 'exp'...
    
    rightmid.frm <- tkframe( mid.frm, borderwidth=2, relief="groove")
    
    exp.base <- tkcheckbutton( rightmid.frm, text="Use Exponential Base",
    					variable=ebase)
    
    other.base <- tkentry( rightmid.frm, textvariable=bb, width=4)
    tkpack( exp.base, side="top")
    
    tkpack( other.base,
    	tklabel( rightmid.frm, text="Use other base", padx=4), side="bottom")
    tkpack( leftmid.frm, rightmid.frm, side="left")
    
    ok.but <- tkbutton( bot.frm, text="OK", command=ldr.trans)
    cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    
    help.but <- tkbutton( bot.frm, text="Help", command=ldrhelp)
    
    tkpack( ok.but, cancel.but, side="left")
    tkpack( help.but, side="right")
    
    # place bindings on buttons so that Return key executes them.
    tkbind( ok.but, "<Return>", ldr.trans)
    tkbind( cancel.but, "<Return>", endprog)
    tkbind( help.but, "<Return>", ldrhelp)
    
    tkpack( top.frm, fill="x")
    tkpack( mid.frm, fill="x")
    tkpack( bot.frm, side="bottom")

} # end of logtrans.gui fcn
