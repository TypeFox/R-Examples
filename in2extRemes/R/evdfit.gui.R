evdfit.gui <- function( base.txt) {

    #  Set the tcl variables
    plot.diags <- tclVar(0)
    threshold.value <- tclVar("")
    mu.fun <- tclVar("1")
    sig.fun <- tclVar("1")
    gam.fun <- tclVar("1")
    threshold.type <- tclVar("constant")
    tiidunits <- tclVar("days")
    mu.link <- tclVar("identity")
    sig.link <- tclVar("identity")
    gam.link <- tclVar("identity")
    log.value <- tclVar(0)
    save.as.value <- tclVar("")
    units.value <- tclVar("")
    na.value <- tclVar("na.fail")

    #########################################
    # internal functions
    #########################################

    refresh <- function() {

        # when a data object is chosen, this function will reset the various lists
        # to have the correct covariates, etc... (at least in theory)

        if(!is.nothing) {

    	    data.select <- as.numeric(tkcurselection(data.listbox))+1
	    dd <- get(full.list[data.select])

	} else dd <- in2extRemesData

	tkdelete(resp.listbox, 0.0, "end")

	for(i in 1:ncol(dd$data)) tkinsert(resp.listbox, "end", paste(colnames(dd$data)[i])) 

    } # end of internal 'refresh' fcn

    redolists<-function() {
        #
        # When a response variable is selected, this function eliminates it
        #  as a covariate option from the other lists.

        if(!is.nothing) {

	    data.select <- as.numeric(tkcurselection( data.listbox))+1
	    dd <- get(full.list[ data.select])

	} else dd <- in2extRemesData

        dd2 <- dd$data[, as.numeric( tkcurselection( resp.listbox))+1]

        resp.name <- colnames(dd$data)[as.numeric(tkcurselection(resp.listbox))+1] 


	for(i in colnames(dd$data)) {

		if(i != resp.name) {


      		} # end of if i != resp.name stmt

	} # end of for i loop    

    } # end of internal 'redolists' fcn

    submit <- function() {
        #
        # The meat of this program.  Actually fits the EVD df/model.
        #

        # names of the covariates used (if any) 
        cov.names.cmd <- "cov.names <- character(0)"
	eval(parse(text=cov.names.cmd))
	write(cov.names.cmd, file="in2extRemes.log", append=TRUE)

	if(!is.nothing) {

		data.select <- as.numeric(tkcurselection(data.listbox))+1
		dd.cmd <- paste( "dd <- get( \"", full.list[data.select], "\")", sep="")

	} else dd.cmd <- "dd <- in2extRemesData"
	eval(parse(text = dd.cmd))
	write(dd.cmd, file = "in2extRemes.log", append = TRUE)

        resp.select <- as.numeric(tkcurselection(resp.listbox))+1

	xdat.cmd <- "xdat <- dd[[\"data\"]]"
	eval(parse(text = xdat.cmd))
	write(xdat.cmd, file = "in2extRemes.log", append = TRUE)
 
        # make sure that a response was selected
        if(is.na(resp.select)) return()

        type.list <- c("GEV", "GP", "PP", "Gumbel")
        type.select <- as.numeric(tkcurselection(type.listbox))+1
        if(length(type.select) == 0) {
                type.value <- "GEV"
        } else type.value <- type.list[type.select]

        units.val <- as.numeric( tclvalue( units.value))
        utype <- tclvalue( threshold.type)
        if( utype == "constant") {

	    threshold.val.cmd <- paste("threshold.val <- ", as.numeric(tclvalue(threshold.value)), sep="")
	    u.fun <- "~1"

	} else if( utype == "vector") {
	
	    threshold.val.cmd <- paste("threshold.val <- ", tclvalue(threshold.value), sep="") 
	    u.fun <- "~1"

	} else {

	    threshold.val.cmd <- "threshold.val <- \"NULL\""
	    u.fun <- paste("~ ", tclvalue(threshold.value), sep = "")

	} # end of if else utype stmts

	eval(parse(text=threshold.val.cmd))
	write(threshold.val.cmd, file="in2extRemes.log", append=TRUE)

        loc.fun <- paste("~ ", tclvalue(mu.fun), sep = "")
	sc.fun <- paste("~ ", tclvalue(sig.fun), sep = "")
	sh.fun <- paste("~ ", tclvalue(gam.fun), sep = "")

	if(tclvalue(log.value) == 0) log.val <- FALSE
	else log.val <- TRUE

	number.of.models <- length(dd$models)
	names.of.models <- names(dd$models)
	if(is.null(names.of.models)) names.of.models <- character(0)

	names.of.models <- c(names.of.models, paste("fit", length(dd$models) + 1, sep=""))

	# fit the df/model
	response.name <- colnames(dd[["data"]])[resp.select]

	if(type.value == "GEV") {

	    cmd <- paste( "dd[[\"models\"]][[\"fit", length(dd$models) + 1, "\"]] <- fit <- fevd( ", response.name,
			", data = xdat",
			", location.fun = ", loc.fun,
                        ", scale.fun = ", sc.fun,
                        ", shape.fun = ", sh.fun,
                        ", use.phi = ", log.val,
                        ", units = \"", tclvalue(units.value), "\"",
			", type = \"", type.value, "\"",
                        ", na.action = ", tclvalue(na.value), ")", sep="")
	} else {

	    cmd <- paste( "dd[[\"models\"]][[\"fit", length(dd$models) + 1, "\"]] <- fit <- fevd( ", response.name, 
			", data = xdat", 
			", threshold = threshold.val, threshold.fun = ", u.fun,
			", location.fun = ", loc.fun,
			", scale.fun = ", sc.fun,
			", shape.fun = ", sh.fun,
			", use.phi = ", log.val,
			", units = \"", tclvalue(units.value), "\"", 
			", time.units = \"", tclvalue(tiidunits), "\"",
			", type = \"", type.value, "\"", 
			", na.action = ", tclvalue(na.value), ")", sep="")

	}
	eval(parse(text=cmd))
	write(cmd, file="in2extRemes.log", append=TRUE)

	if(is.null(dd$models[[number.of.models+1]])) {

	    # failure to fit

	    print(paste("Fit failed."))

        } else {

	    printCMD <- "print(fit)"
	    eval(parse(text = printCMD))
	    write(printCMD, file = "in2extRemes.log", append = TRUE)

	    # plot diagnostics if requested

            if(tclvalue(plot.diags)==1) {

		plotCMD <- paste( "plot(fit)", sep="")
		eval(parse(text=plotCMD))
                write(plotCMD, file="in2extRemes.log", append=TRUE)
            }

	    if(is.nothing) assignCMD <- "assign( \"in2extRemesData\", dd, pos=\".GlobalEnv\")"
            else assignCMD <- paste("assign( \"", full.list[ data.select], "\", dd, pos=\".GlobalEnv\")", sep="")
	    eval( parse( text=assignCMD))
	    write( assignCMD, file="in2extRemes.log", append=TRUE)

	    print(paste("Model name: ", names.of.models[number.of.models+1], sep=""))
            tkyview.moveto( base.txt, 1.0)

	} # end of failure to fit stmt
 
	tkdestroy( base)

    } # end of internal 'submit' fcn

    fevdhelp <- function() {

	cat("\n", "Invokes \'fevd\'", "\n")
	cat( "Use \'help(fevd)\' for more information.\n")
	# help(fevd)

	invisible()

    } # end of internal 'fevdhelp' fcn

    endprog<-function() {

	tkdestroy(base)

    } # end of internal 'endprog' fcn

    #################################
    # Frame/button setup
    #################################


    base<-tktoplevel()
    tkwm.title(base,"Fit EV Model")

    data.frm <- tkframe( base, borderwidth=2, relief="groove")
    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")
    args.frm <- tkframe( base, borderwidth=2, relief="groove")
    threshold.frm <- tkframe( args.frm, borderwidth=2, relief="groove")
    units.frm <- tkframe( args.frm, borderwidth=2, relief="groove")
    units2.frm <- tkframe( units.frm, borderwidth=2, relief="flat")
    tiidunits.frm <- tkframe( units.frm, borderwidth=2, relief="flat")
    optim.frm <- tkframe( base, borderwidth=2, relief="groove")
    type.frm <- tkframe( optim.frm, borderwidth=2, relief="flat")
    plot.frm <- tkframe( optim.frm, borderwidth=2, relief="flat")
    na.frm <- tkframe(optim.frm, borderwidth = 2, relief = "flat")

    # Choose which data object to use (set the listbox to contain all objects of
    # class "in2extRemesDataObject").

    data.listbox <- tklistbox(data.frm,
			yscrollcommand=function(...) tkset(data.scroll,...),
			selectmode="single",
			width=20,
			height=5,
			exportselection=0)

    data.scroll <- tkscrollbar( data.frm, orient="vert",
			command=function(...)tkyview(data.listbox,...))

    # initialize variables for data list.
    # 'temp' is list of everything in global environment.
    # 'full.list' will be list of all objects in '.GlobalEnv' of class "in2extRemesDataObject".
    temp <- ls(all.names=TRUE, name=".GlobalEnv")
    is.nothing <- TRUE
    full.list <- character(0)

    for( i in 1:length( temp)) {

	if( is.null( class( get( temp[i])))) next

	if( (class( get( temp[i]))[1] == "in2extRemesDataObject")) {

		tkinsert( data.listbox, "end", paste( temp[i]))
		full.list <- c(full.list, temp[i])
		is.nothing <- FALSE

	}

    } # end of for i loop

    tkpack(tklabel(data.frm, text="Data Object", padx=4), side="left")
    tkpack(data.listbox, side="left")
    tkpack(data.scroll, side="right", fill="y")
    tkpack(data.frm)

    # place binding on data.listbox to reflect the chosen data from the list.
    tkbind(data.listbox, "<Button-1>", "")
    tkbind(data.listbox, "<ButtonRelease-1>", refresh)

    # top frame for response variable

    top.r <- tkframe(top.frm,borderwidth=2)
    top.l <- tkframe(top.frm,borderwidth=2)
    param.frm <- tkframe(top.r, borderwidth=2, relief="groove")
    mu.frm <- tkframe(param.frm, borderwidth=2, relief="flat")
    sig.frm <- tkframe(param.frm, borderwidth=2, relief="flat")
    gam.frm <- tkframe(param.frm, borderwidth=2, relief="flat")

    resp.listbox <- tklistbox(top.l,yscrollcommand=function(...)tkset(resp.scroll,...),
			selectmode="single",width=35,height=4,exportselection=0)
    resp.scroll <- tkscrollbar(top.l,orient="vert", command=function(...)tkyview(resp.listbox,...))

    if(is.nothing) {

	for(i in 1:ncol(in2extRemesData$data)) tkinsert( resp.listbox, "end", paste(colnames( in2extRemesData$data)[i]))  

    } else tkinsert(resp.listbox, "end", "")

    tkpack(tklabel(top.l,text="Response:",padx=4), side="top")
    tkpack(resp.listbox,side="left")
    tkpack(resp.scroll,side="right",fill="y")

    # place binding on resp.listbox to eliminate the response from the lists of covs.
    tkbind(resp.listbox,"<ButtonRelease-1>",redolists)

    # threshold frame

    tkpack(tklabel(threshold.frm, text="Threshold", padx=4), side="top")
    for(i in c("constant", "vector")) {

	tmp <- tkradiobutton(threshold.frm, text=i, value=i, variable=threshold.type)
	tkpack(tmp, anchor="w")

    } # end of for i loop

    threshold.entry <- tkentry(threshold.frm, textvariable=threshold.value, width=5)

    tkpack(tklabel(threshold.frm, text="Value(s)", padx=4), threshold.entry, side="left")

    # tiidunits frame

    tiidunits.entry <- tkentry(tiidunits.frm, textvariable=tiidunits, width=6)
    tkpack(tklabel(tiidunits.frm, text="Time Units", padx=4), tiidunits.entry, side="left")
    units.entry <- tkentry(units2.frm, textvariable=units.value, width=6)
    tkpack(tklabel(units2.frm, text="Response Units", padx=4), units.entry, side="left")
    tkpack(units2.frm, tiidunits.frm, side = "top")
    tkpack(threshold.frm, units.frm, side="left", fill="both")

    # mu frame

    mu.entry <- tkentry(mu.frm, textvariable = mu.fun, width = 25)
    tkpack(tklabel(mu.frm, text = "Location parameter function ~", padx = 4), mu.entry, side = "left")

    # sigma frame

    sig.entry <- tkentry(sig.frm, textvariable = sig.fun, width = 25)
    tkpack(tklabel(sig.frm, text = "Scale parameter function ~", padx = 4), sig.entry, side = "left")
 
    # gamma frame

    gam.entry <- tkentry(gam.frm, textvariable = gam.fun, width = 25)
    tkpack(tklabel(gam.frm, text = "Shape parameter function ~", padx = 4), gam.entry, side = "left")

    log.but <- tkcheckbutton(param.frm, text = "log-link for scale parameter", variable = log.value)
    tkpack(log.but, side = "left")

    tkpack(param.frm, side="left")
    tkpack(top.l,top.r,side="top")
 
    # type frame

    type.listbox <- tklistbox( type.frm,
                        yscrollcommand=function(...)tkset(type.scroll,...),
                        selectmode="single",
                        width=50,
                        height=3,
                        exportselection=0)

    type.scroll <- tkscrollbar( type.frm, orient="vert",
                command=function(...)tkyview(type.listbox, ...))

    tkinsert(type.listbox, "end", "Generalized Extreme Value (GEV)")
    tkinsert(type.listbox, "end", "Generalized Pareto (GP)")
    tkinsert(type.listbox, "end", "Point Process (PP)")
    tkinsert(type.listbox, "end", "Gumbel")

    tkpack(tklabel(type.frm, text="Model Type", padx=4), side="left")
    tkpack(type.listbox, type.scroll, side="left")

    plot.but<- tkcheckbutton(plot.frm,text="Plot diagnostics",variable=plot.diags)
    tkpack(plot.but, side = "left")

    na.entry <- tkentry(na.frm, textvariable = na.value, width = 10)
    tkpack(tklabel(na.frm, text = "NA Action", padx = 4), na.entry, side = "left")

    # bottom frame
    ok.but <- tkbutton( bot.frm, text="OK", command=submit)  
    quit.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    help.but <- tkbutton( bot.frm, text="Help", command=fevdhelp)

    tkpack( ok.but, quit.but, side="left")
    tkpack( help.but, side="right")

    # place bindings on "OK", "Cancel" and "Help" buttons so that user can hit the return key to execute them.
    tkbind(ok.but, "<Return>", submit)
    tkbind(quit.but, "<Return>", endprog)
    tkbind(help.but, "<Return>", fevdhelp)

    tkpack(top.frm, side="top", fill="x")
    tkpack(type.frm, na.frm, plot.frm, fill="x")
    tkpack(tklabel(optim.frm, text="Options", padx=4), side="top", before=type.frm)

    tkpack(optim.frm, args.frm, side="top", fill="both")
    tkpack(mu.frm, sig.frm, gam.frm, side="top")
    tkpack(bot.frm, side="top", fill="x")

} # end of ppfit.gui fcn
