load.data <- function(base.txt) {

    #
    # This function loads the data set using tkfilefind to locate it.
    #

    ############################
    # Internal functions
    ############################
    
    readit <- function() {

    	# This is the name for the new data object in R.
    	save.name <- tclvalue( sname)
    	in2extRemesData.cmd <- "in2extRemesData <- list()"
    	eval( parse( text=in2extRemesData.cmd))
    	write( in2extRemesData.cmd, file="in2extRemes.log", append=TRUE)
    
    	# Data is actually read into R here.
    	hh <- as.logical(as.numeric(tclvalue(head))) 
    	if( tclvalue( file.type)=="common") {

    	    in2extRemesData.cmd <- paste( "read.table( \"", file.name, "\", header=", hh,
    					", sep=\"", tclvalue(delimiter), "\")", sep="")
    	    in2extRemesData$data <- eval( parse( text=in2extRemesData.cmd))
    	    in2extRemesData.cmd <- paste( "in2extRemesData$data <- ", in2extRemesData.cmd, sep="")
    	    write( in2extRemesData.cmd, file="in2extRemes.log", append=TRUE)

    	} else {
    		in2extRemesData.cmd <- paste( "in2extRemesData[[\"data\"]] <- source( \"", file.name, "\")$value", sep="")
    		# in2extRemesData$data <- source( file.name)$value
    		eval( parse( text=in2extRemesData.cmd))
                    write( in2extRemesData.cmd, file="in2extRemes.log", append=TRUE)
    	}
    
    	if( is.null( dim( in2extRemesData$data))) {
    		nl <- length( in2extRemesData$data)
    		# in2extRemesData$data <- cbind( 1:nl, in2extRemesData$data)
    		in2extRemesData.cmd <- paste( "in2extRemesData[[\"data\"]] <- cbind( 1:", nl, ", in2extRemesData[[\"data\"]])",
    						sep="")
    		eval( parse( text=in2extRemesData.cmd))
                    write( in2extRemesData.cmd, file="in2extRemes.log", append=TRUE)
    		if( is.null( colnames( in2extRemesData$data))) {
    			colnames.cmd <- "colnames( in2extRemesData[[\"data\"]]) <- c(\"obs\", \"value\")"
    			# colnames( in2extRemesData$data) <- c("obs", "value")
    			eval( parse( text=colnames.cmd))
    			write( colnames.cmd, file="in2extRemes.log", append=TRUE)
    		} else {
    			colnames.cmd <- "colnames( in2extRemesData[[\"data\"]])[1] <- c(\"obs\", \"value\")"
    			# colnames( in2extRemesData$data)[1] <- "obs"
    			eval( parse( text=colnames.cmd))
                            write( colnames.cmd, file="in2extRemes.log", append=TRUE)
    		}
    	} else {
    		nc <- dim( in2extRemesData$data)[2]
    		if( is.null( colnames( in2extRemesData$data))) {
    			colnames.cmd <- paste( "colnames( in2extRemesData[[\"data\"]]) <- paste( \"V\", 1:", nc, ", sep=\"\")",
    						sep="")
    			# colnames( in2extRemesData$data) <- as.character(1:nc)
    			eval( parse( text=colnames.cmd))
                            write( colnames.cmd, file="in2extRemes.log", append=TRUE)
    			}
    		}
    	class.cmd <- "class( in2extRemesData) <- \"in2extRemesDataObject\""
    	eval( parse( text=class.cmd))
    	write( class.cmd, file="in2extRemes.log", append=TRUE)
    	# class( in2extRemesData) <- "in2extRemesDataObject"
    
    	in2extRemesData.cmd <- paste( "in2extRemesData$name <- strsplit( \"", file.name, "\", \"/\")[[1]][ ",
    					length(strsplit(file.name,"/")[[1]]), "]", sep="")
    	# in2extRemesData$name <- strsplit(file.name,"/")[[1]][length(strsplit(file.name,"/")[[1]])] 
    	eval( parse( text=in2extRemesData.cmd))
    	write( in2extRemesData.cmd, file="in2extRemes.log", append=TRUE)
    
    	in2extRemesData.cmd <- paste( "in2extRemesData$file.path <- \"", file.name, "\"", sep="")
    	# in2extRemesData$file.path <- file.name
    	eval( parse( text=in2extRemesData.cmd))
            write( in2extRemesData.cmd, file="in2extRemes.log", append=TRUE)
    	if( save.name == "") {
    		save.name <- "in2extRemesData"
    		saveit <- FALSE
    	} else saveit <- TRUE
    
    	# Assign data object the value of 'save.name' in R--default is "in2extRemesData".
    	assignCMD <- paste( "assign( \"", save.name, "\", in2extRemesData, pos=\".GlobalEnv\")", sep="")
    	eval( parse( text=assignCMD))
    	write( assignCMD, file="in2extRemes.log", append=TRUE)
    	in2extRemesData$default.ldata <- save.name
    	print( paste( "Successfully opened file: ", in2extRemesData$name, sep=""))
    	msg.cmd <- "print( summary( in2extRemesData[[\"data\"]]))"
    	eval( parse( text=msg.cmd))
    	write( msg.cmd, file="in2extRemes.log", append=TRUE)
    	if( saveit) {
                    cat("\n", "Saving workspace (may take a few moments for large workspaces) ...\n")
                    saveCMD <- "save.image()"
                    eval( parse( text=saveCMD))
                    write( saveCMD, file="in2extRemes.log", append=TRUE)
    		cat( "\n", "Workspace saved.\n")
            }

    	tkdestroy( tt)

    } # end of internal 'readit' fcn
    
    endprog <-function() {

        tkdestroy(tt)

    }
    
    
    ########################################### 
    # get the filename
    # file.name <-tkfilefind()
    file.name <- tclvalue( tkgetOpenFile()) 
     
    # make sure that a file was selected
    if (!is.null(file.name)) {
    	tt<-tktoplevel()
    	tkwm.title(tt,"Read File")
    	delimiter <- tclVar("")
    	head <- tclVar("0")
    	sname <- tclVar("")
    	file.type <- tclVar("common") # Other choice is 'Rsrc' (R source).
    
    ################################
    #  Frame/button setup
    ################################
    
    
    	spec.frm <- tkframe(tt, borderwidth=2)
    	
    	ftype.frm <- tkframe( spec.frm, relief="groove", borderwidth=2)
    	tkpack( tklabel( ftype.frm, text="File Type", padx=4), side="top")
    	types <- c("Common", "R source")
    	types2 <- c("common", "Rsrc")
    	for( i in 1:2) {
    		tmp <- tkradiobutton( ftype.frm, text=types[i],
    				value=types2[i], variable=file.type)
    		tkpack( tmp, anchor="w")
    		}
    
    	left.frm <- tkframe(spec.frm,relief="groove",borderwidth=2)
    	sep.entry <- tkentry(left.frm,textvariable=delimiter, width=1)
    	tkpack(tklabel(left.frm, text="Delimiter:", padx=0), sep.entry, 
    			anchor="w")
    	
    	right.frm <- tkframe(spec.frm, relief="groove", borderwidth=2)
    	header.cbut <- tkcheckbutton(right.frm, text="Header", variable=head) 
    	tkpack(header.cbut, anchor="w")
    
    	tkpack( ftype.frm, side="left", fill="y")
    	tkpack(left.frm,side="left",fill="y")
    	tkpack(right.frm, side="left",fill="y")
    	
    	txt.frm<-tkframe(tt,relief="groove",borderwidth=2)
    	txt<-tktext(txt.frm,bg="white",font="courier")
    	scr.txt<-tkscrollbar(txt.frm,command=function(...)tkyview(txt,...))
    	tkconfigure(txt,yscrollcommand=function(...)tkset(scr.txt,...))
    	tkpack(txt,side="left",fill="both",expand=TRUE)
    	tkpack(scr.txt,side="right",fill="y")
    # File is read in here preliminarily.  This took a tremendous amount of time.
    # Thus it has been removed.
    # data.file<-tkcmd("open",file.name)
    #	tkinsert(txt,"end",tkcmd("read",data.file))
    # 	tkcmd("close",data.file)
    	tkconfigure(txt,state="disabled")
       
    # left.frm <- tkframe(spec.frm,relief="groove",borderwidth=2)
    #        sep.entry <- tkentry(left.frm,textvariable=delimiter, width=1)
    #        tkpack(tklabel(left.frm, text="Separator:", padx=0), sep.entry,
    #                        anchor="w")
    
    	save.frm <- tkframe(spec.frm, relief="groove", borderwidth=2)
    	save.entry <- tkentry(save.frm, textvariable=sname, width=20) 
    	tkpack( tklabel( save.frm, text="Save As (in R)", padx=0), save.entry,
    			anchor="w")
    	tkpack( save.frm, side="left", fill="y")
    
    	sub.but <- tkbutton(spec.frm,text="OK",command=readit)
    	# Place binding on 'sub.but' so that user may simply hit return key.
    	# However, only works if user "TABS" over to the 'sub.but' first.
    	tkbind( sub.but, "<Return>", readit)
    
    	quit.but <- tkbutton(spec.frm,text="Cancel", command=endprog)
    	tkbind( quit.but, "<Return>", endprog)
    	tkpack(save.entry, sub.but, quit.but,fill="x", anchor="n")         
    	tkpack(spec.frm)
      
    	tkpack(txt.frm)
    	tkwait.window(tt) 
    
    	} # end of if file actually read stmt
}
