if(getRversion() >= "2.15.1") utils::globalVariables(c("ablineCMD", "assignCMD", "assign.cmd", "atdf.cmd", "bmaxxer",
    "cb.cmd", "ci.cmd", "class.cmd", "cmd", "cnames", "cnamesCMD", "cnames.cmd", "colnamesCMD", "colnames.cmd", "cols.selected",
    "cols.selected.cmd", "cosTcmd", "cov.names.cmd",
    "dataCMD", "dd", "dd.cmd", "eiCMD", "fit", "fit.cmd", "fitrange.cmd",
    "gev.sim.cmd", "gp.sim.cmd", "histCMD", "in2extRemesData", "ind.cmd", "lt.cmd",
    "m0.cmd", "m1.cmd", "M.cmd", "mod.fit", "mod.fit.cmd", "mrlplotCMD", "msg.cmd",
    "NameCMD", "N.cmd", "negCMD", "newdata.cmd", "newnames",
    "out.cmd", "p.cmd", "plotCMD", "printCMD", "putfitCMD", "response.cmd", "rlplotCMD",
    "saveCMD", "sinTcmd", "summaryCMD",
    "threshold.val.cmd", "tmpCMD", "tmp.data.cmd", "tmp.ldr.cmd", "tmp.names", "trendCMD", "trendNameCMD", 
    "var.name.cmd", "var.val.cmd", "Xcmd", "xdat.cmd", "z.cmd"))

in2extRemes <- function () {

    ev.dataexists <- function() {

        tmp1 <- ls(all.names=TRUE, pos=".GlobalEnv")
        n <- length(tmp1)

        # If a data object of class "in2extRemesDataObject" exists, will return 1.
        # Otherwise, will return 0.

        data.exists <- FALSE

	for( i in 1:n) {

	    if(is.null(class(get(tmp1[i])))) next

	    if(class(get(tmp1[i]))[1] == "in2extRemesDataObject") {

	        data.exists <- TRUE
		break

	    }

	} # end of for i loop

        if(!data.exists) {

	    cat("\n", "Must load or simulate a data set.\n")
    
    	    return(0)

        } else return(1)

    } # end of internal 'ev.dataexists' fcn

    savework <- function() {

	cat( "\n", "Saving R workspace ...\n")
	cat( "This may take a few moments if the workspace is large.\n") 
	saveCMD <- "save.image()"

	eval( parse( text=saveCMD))
	write( saveCMD, file="extRemes.log", append=TRUE)

	cat( "\n", "Workspace saved.\n")

	invisible()

    } # end of internal 'savework' function.

    endprog <- function() {

        tkdestroy(base)

    } # end of 'endprog' internal function.

    readdata <- function() {

        load.data(txt)

    } # end of 'readdata' internal fcn

    simgevdata <- function() {

        simgev.gui(txt)
        tkyview.moveto(txt, 1)

    } # end of internal 'simgevdata' fcn

    simgpdata <- function() {

	simgp.gui(txt)
	tkyview.moveto(txt, 1)

    } # end of internal 'simgpdata' fcn.

    bmaxxer <- function() {

        if( ev.dataexists()) {
            bmaxxer.gui(txt)
            tkyview.moveto(txt, 1)
        }

    } # end of internal 'bmaxxer' fcn

    decluster <- function() {

	if( ev.dataexists()) {
            declusterGUI(txt)
            tkyview.moveto(txt, 1)
        }

    } # end of internal 'decluster' fcn
	
    negtrans <- function() {

	if( ev.dataexists()) {

            negtrans.gui(txt)
            tkyview.moveto(txt, 1)

        }

    } # end of internal 'negtrans' fcn

    logtrans <- function() {

	if( ev.dataexists()) {

	    logtrans.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'logtrans' function

    logdr <- function() {

	if( ev.dataexists()) {

	    logdr.gui( txt)
	    tkyview.moveto( txt, 1)

	}

    } # end of internal 'logdr' function
	
    affinetrans <- function() {

	if( ev.dataexists()) {

            affine.gui( txt)
            tkyview.moveto( txt, 1)

	}

    } # end of internal 'affinetrans' fcn
	
    indicatortrans <- function() {

        if( ev.dataexists()) {

            indicatorTransform.gui( txt)
            tkyview.moveto( txt, 1)

        }

    } # end of internal 'indicatortrans' fcn

    trigtrans <- function() {

        if( ev.dataexists()) {

	    trigtrans.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'trigtrans' fcn

    DataSummary <- function() {

	if( ev.dataexists()) {

	    DataSummaryGUI(txt)
	    tkyview.moveto( txt, 1)

	}

    } # end of internal 'DataSummary' fcn.

    scrubber <- function() {

	if( ev.dataexists()) {

            scrubber.gui(txt)
            tkyview.moveto(txt, 1)

        }

    } # end of internal 'scrubber' function

    clearlogfile <- function() {

		clearlog(txt)

    } # end of internal 'clearlogfile' function

    evdfit <- function() {

        if (ev.dataexists()) {

            evdfit.gui(txt)
            tkyview.moveto(txt, 1)

        }

    } # end of internal 'evdfit' fcn

    poissonfit <- function() {

        if (ev.dataexists()) {

            poisson.gui(txt)
            tkyview.moveto(txt, 1)

        }

    } # end of internal 'poissonfit' fcn

    extremalind <- function() {

	if( ev.dataexists()) {

	    extremalind.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'extremalind' fcn

    deviancecomparison <- function() {

	if(ev.dataexists()) {

	    llhrt.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'deviancecomparison' fcn

    plotdata <- function() {

	if( ev.dataexists()) {

	    scatterplot.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'plotdata' fcn

    mrlplot <- function() {

	if( ev.dataexists()) {

	    mrlplot.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'mrlplot' fcn

    threshrange <- function() {

	if( ev.dataexists()) {

	    threshrange.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'threshrange' fcn

    atdffun <- function() {

        if( ev.dataexists()) {

            atdf.gui(txt)
            tkyview.moveto(txt, 1)

        }

    } # end of internal 'atdffun' fcn 

    fitdiag <- function() {

	if(ev.dataexists()) {

            fitdiag.gui(txt)
            tkyview.moveto(txt, 1)

        }

    } # end of internal 'fitdiag' fcn
	
    # histplot <- function() {

# 	if(ev.dataexists()) {

# 	    histogram.gui(txt)
# 	    tkyview.moveto(txt, 1)

# 	}

 #    } # end of internal 'histplot' fcn

    rlplot <- function() {

        if(ev.dataexists()) {

            rlplot.gui(txt)
            tkyview.moveto(txt, 1)

        }

    } # end of internal 'rlplot' fcn

    evdCI <- function() {

	if(ev.dataexists()) {

	    evdCI.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'evdCI' fcn

    fitsummary <- function() {

	if(ev.dataexists()) {

	    fitsummary.gui(txt)
	    tkyview.moveto(txt, 1)

	}

    } # end of internal 'fitsummary' fcn

    base <- tktoplevel()
    tkwm.title(base, "Into the extRemes Package")
    top.frm <- tkframe(base, borderwidth = 2)
    bottom.frm <- tkframe(base, borderwidth = 2)

    ##
    ## Menu buttons...  
    ##

    # File menu...

    fmenu.but <- tkmenubutton(top.frm, text = "File", relief = "raised", borderwidth = 2)

    file.menu <- tkmenu(fmenu.but)
    tkconfigure(fmenu.but, menu = file.menu)
    tkadd(file.menu, "command", label = "Read Data", command = readdata)

    SimMenu <- tkmenu(file.menu, tearoff=FALSE)
    tkadd(SimMenu, "command", label = "Generalized Extreme Value (GEV)", command=simgevdata)
    tkadd(SimMenu, "command", label = "Generalized Pareto (GP)", command = simgpdata)
    tkadd(file.menu, "cascade", label="Simulate Data", menu=SimMenu)

    tkadd(file.menu, "separator")
    tkadd(file.menu, "command", label="Block Maxima", command=bmaxxer)
    tkadd(file.menu, "command", label="Decluster", command=decluster)
    TransMenu <- tkmenu(file.menu, tearoff=FALSE)

    tkadd(TransMenu, "command", label="Negative", command=negtrans)
    tkadd(TransMenu, "command", label="Logarithm", command=logtrans)
    tkadd(TransMenu, "command", label="Log Daily Returns", command=logdr)
    tkadd(TransMenu, "command", label="Affine Transformation", command=affinetrans)
    tkadd(TransMenu, "command", label="Indicator Transformation", command=indicatortrans)
    tkadd(TransMenu, "command", label="Trigonometric Transformation", command=trigtrans)
    tkadd(file.menu, "cascade", label="Transform Data", menu=TransMenu)

    tkadd(file.menu, "separator")
    tkadd(file.menu, "command", label="Data Summary", command=DataSummary)
    tkadd(file.menu, "separator")
    tkadd(file.menu, "command", label="Scrubber", command=scrubber)
    tkadd(file.menu, "command", label="Clear log file", command=clearlogfile)
    tkadd(file.menu, "separator")
    tkadd(file.menu, "command", label="Save", command=savework)
    tkadd(file.menu, "separator")
    tkadd(file.menu, "command", label = "Exit", command = endprog)

    # Plot menu...

    pmenu.but <- tkmenubutton( top.frm, text="Plot", relief="raised", borderwidth=2)
    plot.menu <- tkmenu( pmenu.but)
    tkconfigure( pmenu.but, menu=plot.menu)
    tkadd(plot.menu, "command", label = "Scatter Plot", command = plotdata)
    tkadd(plot.menu, "separator")
    tkadd(plot.menu, "command", label="Mean Residual Life Plot", command=mrlplot)
    tkadd(plot.menu, "command", label="Fit POT model to a range of thresholds", command=threshrange)
    tkadd(plot.menu, "command", label="Auto tail-dependence function", command=atdffun)
    tkadd(plot.menu, "separator")
    tkadd(plot.menu, "command", label="Fit Diagnostics", command=fitdiag)
    # tkadd(plot.menu, "command", label="Fitted Model Density with Histogram", command=histplot)
    tkadd(plot.menu, "command", label="Return Level Plot", command=rlplot)

    # Analyze menu...

    amenu.but <- tkmenubutton(top.frm, text = "Analyze", relief = "raised", borderwidth = 2)
    ana.menu <- tkmenu(amenu.but)
    tkconfigure(amenu.but, menu = ana.menu)

    tkadd(ana.menu, "command", label="Extreme Value Distributions", command=evdfit)
    tkadd( ana.menu, "command", label="Poisson Distribution", command=poissonfit)

    tkadd( ana.menu, "separator")

    tkadd( ana.menu, "command", label="Parameter Confidence Intervals", command = evdCI)

    tkadd( ana.menu, "command", label="Likelihood-ratio test", command=deviancecomparison)
    tkadd( ana.menu, "separator")
    tkadd( ana.menu, "command", label="Fit Summary", command=fitsummary)
    tkadd( ana.menu, "separator")
    tkadd( ana.menu, "command", label="Extremal Index", command=extremalind)
    tkpack(top.frm, side = "top", fill = "x")
    tkpack(fmenu.but, pmenu.but, amenu.but, side = "left")

# Base frame...  (for text messages)
	txt <- tktext(bottom.frm, bg = "darkorange", font = "courier")
	# txt <- tktext(bottom.frm, bg = "white", font = "courier")
	# txt <- tktext(bottom.frm, bg = "gray", font = "courier")
	# scr <- tkscrollbar(bottom.frm, repeatinterval = 5,
	# 			command = function(...) tkyview(txt, ...))
	# tkconfigure(txt, yscrollcommand = function(...) tkset(scr, ...))
	tkconfigure(txt, state = "disabled")
	tkpack(txt, side = "left", fill = "both")
	# tkpack(scr, side = "right", fill = "y")
	tkpack(bottom.frm, side = "bottom")
	write.in2extRemesMainMessage(txt=txt)

} # end of 'in2extRemes' function
