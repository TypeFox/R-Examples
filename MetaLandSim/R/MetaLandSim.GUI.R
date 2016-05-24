MetaLandSim.GUI <- function() {
	ls.label <- NULL
	output <- NULL
    mainW <- tktoplevel(width = 700, height = 20)
    loading <- tclVar()
    tcl("image", "create", "photo", loading, file = paste(system.file(package = "MetaLandSim"), 
        "/logo/logo.gif", sep = ""))
    background.logo <- tklabel(mainW, image = loading)
    tkpack(background.logo)
    tktitle(mainW) <- "MetaLandSim - Landscape and Metapopulation Simulation"
    topMenu <- tkmenu(mainW)
    tkconfigure(mainW, menu = topMenu)
    fileM <- tkmenu(topMenu, tearoff = FALSE)
    paraM <- tkmenu(topMenu, tearoff = FALSE)
    simM <- tkmenu(topMenu, tearoff = FALSE)
    rangeM <- tkmenu(topMenu, tearoff = FALSE)
    viewM <- tkmenu(topMenu, tearoff = FALSE)
    aboutM <- tkmenu(topMenu, tearoff = FALSE)
    Newsession<- function(gisBase, override) {
      initGRASS(gisBase = gisBase, home = tempdir(), 
	  override = as.logical(override))
    }
tkadd(rangeM, "command", label = "Start GRASS session",
  command = function() guiv(Newsession, exec = "Execute",
  argText = c(gisBase = "Path to GRASS binaries and libraries:",
  override = "override GRASS session (TRUE or FALSE):"),
  helpsFunc = "initGRASS"))
	choose.df <- function(df.entry, dfnr.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (is.data.frame(xobj)) {
                tkinsert(tlb, "end", x1)
                cbind(nrow(xobj))
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(df.entry, 0, "end")
        tkinsert(df.entry, "end", choix)
        tkconfigure(dfnr.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    convert.graph.out <- function(data1, data2, data3, outfile) {
        exp.eval <- paste(outfile, "<<- convert.graph(dframe = ", data1, ", mapsize = ", 
            data2, ", dispersal = ", data3, ")", sep = "")
        eval(parse(text = exp.eval))
    }
    convert.graph.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "Import data")
        IFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        dfvar <- tclVar()
        df.entry <- tkentry(IFrame, textvariable = dfvar)
        dfnr.label <- tklabel(IFrame, width = 5)
        choosevect.but <- tkbutton(IFrame, text = "Set", command = function() choose.df(df.entry, 
            dfnr.label))
        tkgrid(tklabel(IFrame, text = "Patch occupancy data frame"), df.entry, choosevect.but, 
            dfnr.label, sticky = "w")
        sizevar <- tclVar()
        size.entry <- tkentry(IFrame, textvariable = sizevar)
        tkgrid(tklabel(IFrame, text = "Landscape mosaic side length (meters)"), size.entry, 
            sticky = "w")
        dispvar <- tclVar()
        disp.entry <- tkentry(IFrame, textvariable = dispvar)
        tkgrid(tklabel(IFrame, text = "Species dispersal ability (meters)"), disp.entry, 
            sticky = "w")
        tkgrid(IFrame)
        OFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for output: "), out.entry, sticky = "w")
        tkgrid(OFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            df1 <- tclvalue(dfvar)
            mapsize <- tclvalue(sizevar)
            dispersal <- tclvalue(dispvar)
            outname <- tclvalue(outvar)
            substitute(convert.graph.out(data1 = df1, data2 = mapsize, data3 = dispersal, 
                outfile = outname))
        }
        "reset" <- function() {
            tclvalue(dfvar) <- ""
            tkconfigure(dfnr.label, text = "")
            tclvalue(sizevar) <- ""
            tclvalue(dispvar) <- ""
            tclvalue(outvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(tt, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(tt))
        tkgrid(submit.but, reset.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tt, "<KeyPress-Return>", function() execcomp())
        tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(fileM, "command", label = "Import data", command = function() convert.graph.gui())
    import.shape.gui <- function(filename, species.col, ID.col, area.col, dispersal, 
        output) {
        result1 <- import.shape(filename = filename, species.col = species.col, ID.col = ID.col, 
            area.col = area.col, dispersal = dispersal)
        exp.eval <- paste(output, "<<- result1", sep = "")
        eval(parse(text = exp.eval))
    }
    importguiCallback <- function(arg) {
        if (arg == "filename") {
            columnames <- names(readShapePoly(guiGetValue("filename"))@data)
            print(columnames)
            guiSet("columnames", columnames)
            setListElements("ID.col", columnames)
            setListElements("area.col", columnames)
            setListElements("species.col", columnames)
        }
    }
    tkadd(fileM, "command", label = "Import shapefile", command = function() guiv(import.shape.gui, 
        exec = "Execute", title = "Import shapefile", argFilename = list(filename = NULL), 
        argList = list(ID.col = NULL, species.col = NULL, area.col = NULL), callback = importguiCallback, 
        argText = c(filename = "Chose your file in ESRI *.shp format", species.col = "Name of the column with species data", 
            ID.col = "Name of the column with the ID (if available)", area.col = "Name of column with the area", 
            dispersal = "Species dispersal ability", output = "Output name"), helpsFunc = "import.shape"))
    tkadd(fileM, "command", label = "Quit", command = function() tkdestroy(mainW))
    tkadd(topMenu, "cascade", label = "File", menu = fileM)
    choose.metapop <- function(df.entry, dfnr.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (class(xobj) == "metapopulation") {
                tkinsert(tlb, "end", x1)
                cbind(xobj$number.patches)
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(df.entry, 0, "end")
        tkinsert(df.entry, "end", choix)
        tkconfigure(dfnr.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    parameter.estimate.out <- function(data1, data2, data3, data4, outfile) {
        exp.eval <- paste(outfile, "<<- parameter.estimate(sp = ", data1, ", method = data2, alpha = ", 
            data3, ", nsnap = ", data4, ")", sep = "")
        eval(parse(text = exp.eval))
    }
    parameter.estimate.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "SPOM parameter estimation")
        IFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        dfvar <- tclVar()
        df.entry <- tkentry(IFrame, textvariable = dfvar)
        dfnr.label <- tklabel(IFrame, width = 5)
        choosevect.but <- tkbutton(IFrame, text = "Set", command = function() choose.metapop(df.entry, 
            dfnr.label))
        tkgrid(tklabel(IFrame, text = "Patch occupancy data frame"), df.entry, choosevect.but, 
            dfnr.label, sticky = "w")
        methodvar <- tclVar()
        method.entry <- tkentry(IFrame, textvariable = methodvar)
        tkgrid(tklabel(IFrame, text = "Method for estimation (see help for available methods)"), 
            method.entry, sticky = "w")
        alphavar <- tclVar()
        alpha.entry <- tkentry(IFrame, textvariable = alphavar)
        tkgrid(tklabel(IFrame, text = "Alpha parameter initial value"), alpha.entry, 
            sticky = "w")
        snapvar <- tclVar()
        snap.entry <- tkentry(IFrame, textvariable = snapvar)
        tkgrid(tklabel(IFrame, text = "number of snapshots (available for method='Rsnap_x')"), 
            snap.entry, sticky = "w")
        tkgrid(IFrame)
        OFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for output: "), out.entry, sticky = "w")
        tkgrid(OFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            df1 <- tclvalue(dfvar)
            method <- as.character(tclvalue(methodvar))
            if (tclvalue(alphavar) != "") {
                alpha <- tclvalue(alphavar)
            }
            else alpha <- NULL
            nsnap <- tclvalue(snapvar)
            outname <- tclvalue(outvar)
            substitute(parameter.estimate.out(data1 = df1, data2 = method, data3 = alpha, 
                data4 = nsnap, outfile = outname))
        }
        "reset" <- function() {
            tclvalue(dfvar) <- ""
            tkconfigure(dfnr.label, text = "")
            tclvalue(sizevar) <- ""
            tclvalue(dispvar) <- ""
            tclvalue(outvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(tt, relief = "groove")
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(tt))
        submit.but <- tkbutton(RCSFrame, text = "  Run  ", default = "active", command = function() execcomp())
        tkgrid(cancel.but, submit.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tt, "<KeyPress-Return>", function() execcomp())
        tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(paraM, "command", label = "SPOM parameter estimation", command = function() parameter.estimate.gui())
    create.parameter.df.gui <- function(alpha, x, y, e, output) {
        result1 <- create.parameter.df(alpha = alpha, x = x, y = y, e = e)
        exp.eval <- paste(output, "<<- result1", sep = "")
        eval(parse(text = exp.eval))
    }
    tkadd(paraM, "command", label = "Create parameter data frame", command = function() guiv(create.parameter.df.gui, 
        exec = "Execute", title = "Create parameter data frame", argText = c(alpha = "alpha", 
            x = "x", y = "y", e = "e", output = "Output name"), helpsFunc = "create.parameter.df"))
    tkadd(topMenu, "cascade", label = "Parameter Estimation", menu = paraM)
    rland.graph.gui <- function(mapsize, dist_m, areaM, areaSD, Npatch, disp, plotG = TRUE, 
        output) {
        result1 <- rland.graph(mapsize = mapsize, dist_m = dist_m, areaM = areaM, 
            areaSD = areaSD, Npatch = Npatch, disp = disp, plotG)
        exp.eval <- paste(output, "<<- result1", sep = "")
        eval(parse(text = exp.eval))
    }
    tkadd(simM, "command", label = "Random landscape generation", command = function() guiv(rland.graph.gui, 
        exec = "Execute", title = "Random landscape generation", argText = c(mapsize = "Landscape mosaic side length (meters)", 
            dist_m = "Minimum distance allowed between patches", areaM = "Mean area (hectares)", 
            areaSD = "Standard deviation for the area of patches", Npatch = "Number of patches. WARNING: limited by the minimum distance parameter configuration", 
            disp = "Mean dispersal ability of focal species(meters)", plotG = "Do plot", 
            output = "Output name"), argOption = list(plotG = c("TRUE", "FALSE")), 
        helpsFunc = "rland.graph"))
    choose.landscape <- function(df.entry, dfnr.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (class(xobj) == "landscape") {
                tkinsert(tlb, "end", x1)
                cbind(xobj$number.patches)
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(df.entry, 0, "end")
        tkinsert(df.entry, "end", choix)
        tkconfigure(dfnr.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    species.graph.out <- function(data1, data2, data3, data4, data6, outfile) {
        exp.eval <- paste(outfile, "<<- species.graph(rl = ", data1, ", method = data2, parm = ", 
            data3, ", nsew = data4, plotG = ", data6, ")", sep = "")
        eval(parse(text = exp.eval))
    }
    species.graph.gui <- function() {
        sg <- tktoplevel()
        tkwm.title(sg, "Generate species occupation")
        IFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        dfvar <- tclVar()
        df.entry <- tkentry(IFrame, textvariable = dfvar)
        dfnr.label <- tklabel(IFrame, width = 5)
        choosevect.but <- tkbutton(IFrame, text = "Set", command = function() choose.landscape(df.entry, 
            dfnr.label))
        tkgrid(tklabel(IFrame, text = "Object of class landscape"), df.entry, choosevect.but, 
            dfnr.label, sticky = "w")
        methodvar <- tclVar()
        method.entry <- tkentry(IFrame, textvariable = methodvar)
        tkgrid(tklabel(IFrame, text = "Species occupancy method ('percentage' by default)"), 
            method.entry, sticky = "w")
        paramvar <- tclVar()
        param.entry <- tkentry(IFrame, textvariable = paramvar)
        tkgrid(tklabel(IFrame, text = "Parameter to specify the species occupancy"), 
            param.entry, sticky = "w")
        pointvar <- tclVar()
        point.entry <- tkentry(IFrame, textvariable = pointvar)
        tkgrid(tklabel(IFrame, text = "Point of entrance for the species. Default 'none'"), 
            point.entry, sticky = "w")
        tkgrid(IFrame)
        LFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(LFrame, text = "Do plot"))
        loadvar <- tclVar(1)
        tkgrid(tkradiobutton(LFrame, text = "TRUE", value = 1, variable = loadvar), 
            columnspan = 3)
        tkgrid(tkradiobutton(LFrame, text = "FALSE", value = 2, variable = loadvar), 
            columnspan = 3)
        tkgrid(LFrame)
        OFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for output: "), out.entry, sticky = "w")
        tkgrid(OFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            df1 <- tclvalue(dfvar)
            if (tclvalue(methodvar) != "") {
                method <- tclvalue(methodvar)
            }
            else method <- "percentage"
            param <- tclvalue(paramvar)
            if (tclvalue(pointvar) != "") {
                point <- tclvalue(pointvar)
            }
            else point <- "none"
            available <- c(TRUE, FALSE)
            plotting <- available[as.numeric(tclvalue(loadvar))]
            outname <- tclvalue(outvar)
            substitute(species.graph.out(data1 = df1, data2 = method, data3 = param, 
                data4 = point, data6 = plotting, outfile = outname))
        }
        "reset" <- function() {
            tclvalue(dfvar) <- ""
            tkconfigure(dfnr.label, text = "")
            tclvalue(methodvar) <- ""
            tclvalue(paramvar) <- ""
            tclvalue(pointvar) <- ""
            tclvalue(minavar) <- ""
            tclvalue(outvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(sg, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(sg))
        tkgrid(submit.but, reset.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(sg, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(sg, "<KeyPress-Return>", function() execcomp())
        tkbind(sg, "<KeyPress-Escape>", function() tkdestroy(sg))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(simM, "command", label = "Generate species occupation", command = function() species.graph.gui())
    iterate.graph.gui <- function(iter, mapsize, dist_m, areaM, areaSD, Npatch, disp, 
        span, par1, par2, par3, par4, par5, method, parm, nsew, succ, param_df, 
        kern, conn, colnz, ext, beta1, b, c1, c2, z, R, graph, outname) {
        result1 <- iterate.graph(iter, mapsize, dist_m, areaM, areaSD, Npatch, disp, 
            span, par1, par2, par3, par4, par5, method, parm, nsew, succ, param_df, 
            kern, conn, colnz, ext, beta1, b, c1, c2, z, R, graph)
        exp.eval <- paste(outname, "<<- result1", sep = "")
        eval(parse(text = exp.eval))
    }
    choose.df <- function(df.entry, dfnr.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (is.data.frame(xobj)) {
                tkinsert(tlb, "end", x1)
                cbind(nrow(xobj))
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(df.entry, 0, "end")
        tkinsert(df.entry, "end", choix)
        tkconfigure(dfnr.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    iterate.GUI <- function() {
        it1 <- tktoplevel()
        tkwm.title(it1, "Persistence in dynamic landscapes")
        IOFrame <- tkframe(it1, relief = "groove", borderwidth = 2)
        IOFrame2 <- tkframe(it1, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IOFrame2, foreground = "blue"), columnspan = 5)
        var1 <- tclVar()
        var2 <- tclVar()
        var3 <- tclVar()
        var4 <- tclVar()
        var5 <- tclVar()
        var6 <- tclVar()
        var7 <- tclVar()
        var8 <- tclVar()
        var9 <- tclVar(init = "none")
        var10 <- tclVar(init = "NULL")
        var11 <- tclVar(init = "NULL")
        var12 <- tclVar(init = "NULL")
        var13 <- tclVar(init = "NULL")
        var14 <- tclVar(init = "percentage")
        var15 <- tclVar()
        var16 <- tclVar(init = "none")
        var17 <- tclVar(init = "none")
        var18 <- tclVar()
        var19 <- tclVar(init = "op1")
        var20 <- tclVar(init = "op1")
        var21 <- tclVar(init = "op1")
        var22 <- tclVar(init = "op1")
        var23 <- tclVar(init = "NULL")
        var24 <- tclVar(init = 1)
        var25 <- tclVar(init = "NULL")
        var26 <- tclVar(init = "NULL")
        var27 <- tclVar(init = "NULL")
        var28 <- tclVar(init = "NULL")
        it.entry1 <- tkentry(IOFrame, textvariable = var1)
        it.entry2 <- tkentry(IOFrame, textvariable = var2)
        it.entry3 <- tkentry(IOFrame, textvariable = var3)
        it.entry4 <- tkentry(IOFrame, textvariable = var4)
        it.entry5 <- tkentry(IOFrame, textvariable = var5)
        it.entry6 <- tkentry(IOFrame, textvariable = var6)
        it.entry7 <- tkentry(IOFrame, textvariable = var7)
        it.entry8 <- tkentry(IOFrame, textvariable = var8)
        it.entry9 <- tkentry(IOFrame, textvariable = var9)
        it.entry10 <- tkentry(IOFrame, textvariable = var10)
        it.entry11 <- tkentry(IOFrame, textvariable = var11)
        it.entry12 <- tkentry(IOFrame, textvariable = var12)
        it.entry13 <- tkentry(IOFrame, textvariable = var13)
        it.entry14 <- tkentry(IOFrame, textvariable = var14)
        it.entry15 <- tkentry(IOFrame2, textvariable = var15)
        it.entry16 <- tkentry(IOFrame2, textvariable = var16)
        it.entry17 <- tkentry(IOFrame2, textvariable = var17)
        df.entry <- tkentry(IOFrame2, textvariable = var18)
        dfnr.label <- tklabel(IOFrame2, width = 5)
        choosevect.but <- tkbutton(IOFrame2, text = "Set", command = function() choose.df(df.entry, 
            dfnr.label))
        it.entry19 <- tkentry(IOFrame2, textvariable = var19)
        it.entry20 <- tkentry(IOFrame2, textvariable = var20)
        it.entry21 <- tkentry(IOFrame2, textvariable = var21)
        it.entry22 <- tkentry(IOFrame2, textvariable = var22)
        it.entry23 <- tkentry(IOFrame2, textvariable = var23)
        it.entry24 <- tkentry(IOFrame2, textvariable = var24)
        it.entry25 <- tkentry(IOFrame2, textvariable = var25)
        it.entry26 <- tkentry(IOFrame2, textvariable = var26)
        it.entry27 <- tkentry(IOFrame2, textvariable = var27)
        it.entry28 <- tkentry(IOFrame2, textvariable = var28)
        but1 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Number of repetitions of the simulation.", title = "iter", 
            type = "ok", parent = IOFrame))
        but2 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Landscape mosaic side length, in meters. To be internally passed to rland.graph.", 
            title = "mapsize", type = "ok", parent = IOFrame))
        but3 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Minimum distance between patches (centroid). To be internally passed to rland.graph.", 
            title = "dist_m", type = "ok", parent = IOFrame))
        but4 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Mean area (in hectares). To be internally passed to rland.graph.", 
            title = "areaM", type = "ok", parent = IOFrame))
        but5 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "SD of the area of patches, in order to give variability to the patches area. To be internally passed to rland.graph.", 
            title = "areaSD", type = "ok", parent = IOFrame))
        but6 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Number of patches (might be impaired by the dist_m, see the 'Note' section). To be internally passed to rland.graph.", 
            title = "Npatch", type = "ok", parent = IOFrame))
        but7 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Species mean dispersal ability, in meters. To be internally passed to rland.graph.", 
            title = "disp", type = "ok", parent = IOFrame))
        but8 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Number of time steps (e.g. years) to simulate. To be internally passed to span.graph", 
            title = "span", type = "ok", parent = IOFrame))
        but9 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Method to produce habitat disturbance.", title = "par1", type = "ok", 
            parent = IOFrame))
        but10 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter specifying details for the options in par1: percentage of patches do delete (if par1 = 'hab'); distance, in meters (if par1 = 'dincr'); percentage of increase/decrease of the mean area of patches (if par1 = 'area'); percentage of new patches (if par1 = 'stoc'); 'northerndness' of created patches (if par1 = 'ncsd'); percentage of destroyed patches (if par1 = 'aggr'). To be internally passed to span.graph. Default NULL.", 
            title = "par2", type = "ok", parent = IOFrame))
        but11 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Additional parameter specifying details for the options in par1: percentage of destroyed patches (if par1 = 'stoc'); 'southerndness' of destroyed patches (if par1 = 'ncsd'); aggregation of destruction (if par1 = 'aggr'). Minimum area for patch deletion, in hectares (if par1 = 'darea'). To be internally passed to span.graph. Default NULL.", 
            title = "par3", type = "ok", parent = IOFrame))
        but12 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Percentage of created patches (if par1 = 'ncsd'). To be internally passed to span.graph. Default NULL.", 
            title = "par4", type = "ok", parent = IOFrame))
        but13 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Percentage of destroyed patches (if par1 = 'ncsd'). To be internally passed to span.graph. Default NULL.", 
            title = "par5", type = "ok", parent = IOFrame))
        but14 <- tkbutton(IOFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "One of the following (default 'percentage'): click - individually select the patches with occurrence of the species by clicking on the map. Use only for individual landscape simulations. However, this option should not be used with iterate.graph. percentage - percentage of the patches to by occupied by the species. number - number of patches to be occupied by the species. To be internally passed to species.graph.", 
            title = "method", type = "ok", parent = IOFrame))
        but15 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "parameter to specify the species occurrence - either percentage of occupied patches or number of occupied patches, depending on the method chosen. To be internally passed to species.graph.", 
            title = "parm", type = "ok", parent = IOFrame2))
        but16 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "'N', 'S', 'E', 'W' or none - point of entry of the species in the landscape. By default set to 'none'. To be internally passed to species.graph.", 
            title = "nsew", type = "ok", parent = IOFrame2))
        but17 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Sucessional preference of the species (Default:'none'). Options: 'early', 'mid' and 'late'.", 
            title = "succ", type = "ok", parent = IOFrame2))
        but18 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter data frame delivered by parameter.estimate.", title = "param_df", 
            type = "ok", parent = IOFrame2))
        but19 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "'op1' or 'op2'. Dispersal kernel. See details in the spom function. To be internally passed to spom.", 
            title = "kern", type = "ok", parent = IOFrame2))
        but20 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "'op1' or 'op2'. Connectivity function. See details in the spom function. To be internally passed to spom.", 
            title = "conn", type = "ok", parent = IOFrame2))
        but21 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "'op1', 'op2' or 'op3'. Colonization function. See details in the spom function. To be internally passed to spom.", 
            title = "colnz", type = "ok", parent = IOFrame2))
        but22 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "'op1', 'op2' or 'op3'. Extinction function. See details in the spom function. To be internally passed to spom.", 
            title = "ext", type = "ok", parent = IOFrame2))
        but23 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter affecting long distance dispersal probability (if the Kern = 'op2'). To be internally passed to spom.", 
            title = "beta1", type = "ok", parent = IOFrame2))
        but24 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter scaling emigration with patch area (if conn = 'op1' or 'op2'). To be internally passed to spom. By default set to 1.", 
            title = "b", type = "ok", parent = IOFrame2))
        but25 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter scaling immigration with the focal patch area (if conn = 'op2'). To be internally passed to spom.", 
            title = "c1", type = "ok", parent = IOFrame2))
        but26 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter c in the option 3 of the colonization probability (if colnz = 'op3'). To be internally passed to spom.", 
            title = "c2", type = "ok", parent = IOFrame2))
        but27 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter giving the strength of the Allee effect (if colnz = 'op3'). To be internally passed to spom.", 
            title = "z", type = "ok", parent = IOFrame2))
        but28 <- tkbutton(IOFrame2, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Parameter giving the strength of the Rescue effect (if ext = 'op3'). To be internally passed to spom.", 
            title = "R", type = "ok", parent = IOFrame2))
        tkgrid(tklabel(IOFrame, text = "Number of simulations"), it.entry1, but1, 
            sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Landscape mosaic side length (meters)"), 
            it.entry2, but2, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Minimum distance between patches (centroid)"), 
            it.entry3, but3, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Mean area (hectares)"), it.entry4, but4, 
            sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Standard deviation of the area of patches"), 
            it.entry5, but5, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Number of patches. WARNING: limited by the minimum distance parameter configuration"), 
            it.entry6, but6, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Mean dispersal ability of focal species(meters)"), 
            it.entry7, but7, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Number of time steps (e.g. years) to simulate"), 
            it.entry8, but8, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Type of habitat disturbance"), it.entry9, 
            but9, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Parameter specifying details for habitat disturbance (see manual)"), 
            it.entry10, but10, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Additional parameter specifying details for habitat disturbance (see manual)"), 
            it.entry11, but11, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Percentage of created patches (if par1 = 'ncsd')"), 
            it.entry12, but12, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Percentage of destroyed patches (if par1 = 'ncsd')"), 
            it.entry13, but13, sticky = "w")
        tkgrid(tklabel(IOFrame, text = "Species occupancy method"), it.entry14, but14, 
            sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Parameter to specify the species occupancy"), 
            it.entry15, but15, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Point of entrance for the species. Default 'none'"), 
            it.entry16, but16, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Succesional preference of the species"), 
            it.entry17, but17, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Parameter data frame delivered by parameter.estimate: "), 
            df.entry, choosevect.but, dfnr.label, but18, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Dispersal Kernel"), it.entry19, but19, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Connectivity"), it.entry20, but20, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Colonization function"), it.entry21, but21, 
            sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Extinction function"), it.entry22, but22, 
            sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "beta - parameter affecting long distance dispersal (if kern='op2')"), 
            it.entry23, but23, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Parameter scaling emigration with patch area (b)"), 
            it.entry24, but24, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Parameter scaling immigration with focal patch area (if conn = op2, parameter c1)"), 
            it.entry25, but25, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Parameter c2 (if colnz='op3')"), it.entry26, 
            but26, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Strength of Allee effect (if colnz='op4', parameter z)"), 
            it.entry27, but27, sticky = "w")
        tkgrid(tklabel(IOFrame2, text = "Strength of Rescue effect (if ext='op3', parameter R)"), 
            it.entry28, but28, sticky = "w")
        tkgrid(IOFrame, IOFrame2)
        LFrame <- tkframe(it1, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(LFrame, text = "Show graphic: "))
        loadvar <- tclVar(1)
        but_plot <- tkbutton(LFrame, text = "?", command = function() tkmessageBox(icon = "info", 
            message = "Show graphic output?", title = "Do plot", type = "Ok", parent = IOFrame))
        tkgrid(tkradiobutton(LFrame, text = "yes", value = 1, variable = loadvar), 
            columnspan = 3)
        tkgrid(tkradiobutton(LFrame, text = "No", value = 2, variable = loadvar), 
            columnspan = 3)
        tkgrid(LFrame)
        OFrame <- tkframe(it1, relief = "groove", borderwidth = 2)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Output name: "), out.entry, sticky = "w")
        tkgrid(OFrame)
        vnr = NULL
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- as.numeric(tclvalue(var1))
            vect2 <- as.numeric(tclvalue(var2))
            vect3 <- as.numeric(tclvalue(var3))
            vect4 <- as.numeric(tclvalue(var4))
            vect5 <- as.numeric(tclvalue(var5))
            vect6 <- as.numeric(tclvalue(var6))
            vect7 <- as.numeric(tclvalue(var7))
            vect8 <- as.numeric(tclvalue(var8))
            vect9 <- tclvalue(var9)
            vect10 <- as.numeric(tclvalue(var10))
            vect11 <- as.numeric(tclvalue(var11))
            vect12 <- as.numeric(tclvalue(var12))
            vect13 <- as.numeric(tclvalue(var13))
            vect14 <- tclvalue(var14)
            vect15 <- as.numeric(tclvalue(var15))
            vect16 <- tclvalue(var16)
            vect17 <- tclvalue(var17)##########
            vect18 <- parse(text = tclvalue(var18))[[1]]
            vect19 <- tclvalue(var19)
            vect20 <- tclvalue(var20)
            vect21 <- tclvalue(var21)
            vect22 <- tclvalue(var22)
            vect23 <- as.numeric(tclvalue(var23))
            vect24 <- as.numeric(tclvalue(var24))
            vect25 <- as.numeric(tclvalue(var25))
            vect26 <- as.numeric(tclvalue(var26))
            vect27 <- as.numeric(tclvalue(var27))
            vect28 <- as.numeric(tclvalue(var28))
            available <- c(TRUE, FALSE)
            choice.plot <- available[as.numeric(tclvalue(loadvar))]
            outname <- tclvalue(outvar)
            substitute(iterate.graph.gui(vect1, vect2, vect3, vect4, vect5, vect6, 
                vect7, vect8, vect9, vect10, vect11, vect12, vect13, vect14, vect15, 
                vect16, vect17, vect18, vect19, vect20, vect21, vect22, vect23, vect24, 
                vect25, vect26, vect27, vect28, graph = choice.plot, outname))
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(it1, relief = "groove")
        cancel.but <- tkbutton(RCSFrame, text = "Cancel", command = function() tkdestroy(it1))
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        tkgrid(submit.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(it1, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(it1, "<KeyPress-Return>", function() execcomp())
        tkbind(it1, "<KeyPress-Escape>", function() tkdestroy(it1))
        if (tclvalue(done) == "2") 
            return()
    }
   tkadd(simM, "command", label = "Persistence in dynamic landscapes", command = function() iterate.GUI())
    manage.simulations.out <- function(data1, data2, outfile) {
        exp.eval <- paste(output, " <<- manage.simulations(par_df = ", data1, ", parameters_spom = ", 
            data2, ")", sep = "")
        eval(parse(text = exp.eval))
    }
    manage.simulations.gui <- function() {
        sg <- tktoplevel()
        tkwm.title(sg, "Batch simulation")
        IFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        dfpar <- tclVar()
        dfpar.entry <- tkentry(IFrame, textvariable = dfpar)
        dfparnr.label <- tklabel(IFrame, width = 5)
        choosepar.but <- tkbutton(IFrame, text = "Set", command = function() choose.df(dfpar.entry, 
            dfparnr.label))
        tkgrid(tklabel(IFrame, text = "Arguments data frame to be parsed by iterate.graph"), 
            dfpar.entry, choosepar.but, dfparnr.label, sticky = "w")
        dfvar <- tclVar()
        df.entry <- tkentry(IFrame, textvariable = dfvar)
        dfnr.label <- tklabel(IFrame, width = 5)
        choosevect.but <- tkbutton(IFrame, text = "Set", command = function() choose.df(df.entry, 
            dfnr.label))
        tkgrid(tklabel(IFrame, text = "Parameter data frame (see manual)"), df.entry, 
            choosevect.but, dfnr.label, sticky = "w")
        tkgrid(IFrame)
        OFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for output: "), out.entry, sticky = "w")
        tkgrid(OFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            df1 <- tclvalue(dfpar)
            df2 <- tclvalue(dfvar)
            outname <- tclvalue(outvar)
            substitute(manage.simulations.out(data1 = df1, data2 = df2, outfile = outname))
        }
        "reset" <- function() {
            tclvalue(dfpar) <- ""
            tkconfigure(dfparnr.label, text = "")
            tclvalue(dfvar) <- ""
            tkconfigure(dfnr.label, text = "")
            tclvalue(outvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(sg, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(sg))
        tkgrid(submit.but, reset.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(sg, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(sg, "<KeyPress-Return>", function() execcomp())
        tkbind(sg, "<KeyPress-Escape>", function() tkdestroy(sg))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(simM, "command", label = "Batch simulation", command = function() manage.simulations.gui())
    tkadd(topMenu, "cascade", label = "Landscape simulation", menu = simM)
    range.expansion.out <- function(data1, data2, data4, data5, data6, data7, 
        outfile) {
        exp.eval <- paste(outfile, " <<- range_expansion(rl = ", data1, ", percI = ", 
            data2, ", param = ", data4, ", b = ", data5, ", tsteps = ", 
            data6, ", iter = ", data7, ")", sep = "")
        eval(parse(text = exp.eval))
    }
    range.expansion.gui <- function() {
        sg <- tktoplevel()
        tkwm.title(sg, "Generate dispersal model")
        IFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        dfvar <- tclVar()
        df.entry <- tkentry(IFrame, textvariable = dfvar)
        dfnr.label <- tklabel(IFrame, width = 5)
        choosevect.but <- tkbutton(IFrame, text = "Set", command = function() choose.landscape(df.entry, 
            dfnr.label))
        tkgrid(tklabel(IFrame, text = "Starting landscape"), df.entry, choosevect.but, 
            dfnr.label, sticky = "w")
        percvar <- tclVar()
        perc.entry <- tkentry(IFrame, textvariable = percvar)
        tkgrid(tklabel(IFrame, text = "Percentage of occupation for starting landscape"), 
            perc.entry, sticky = "w")
        paramvar <- tclVar()
        param.entry <- tkentry(IFrame, textvariable = paramvar)
        param.label <- tklabel(IFrame, width = 5)
        chooseparam.but <- tkbutton(IFrame, text = "Set", command = function() choose.df(param.entry, 
            param.label))
        tkgrid(tklabel(IFrame, text = "Parameter data frame (see manual)"), param.entry, 
            chooseparam.but, param.label, sticky = "w")
        parambvar <- tclVar()
        paramb.entry <- tkentry(IFrame, textvariable = parambvar)
        tkgrid(tklabel(IFrame, text = "Parameter scaling emigration with patch area (b)"), 
            paramb.entry, sticky = "w")
        stepsvar <- tclVar()
        steps.entry <- tkentry(IFrame, textvariable = stepsvar)
        tkgrid(tklabel(IFrame, text = "Number of time steps"), steps.entry, sticky = "w")
        itervar <- tclVar()
        iter.entry <- tkentry(IFrame, textvariable = itervar)
        tkgrid(tklabel(IFrame, text = "Number of iterations"), iter.entry, sticky = "w")
        tkgrid(IFrame)
        OFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for output: "), out.entry, sticky = "w")
        tkgrid(OFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            dfvar <- tclvalue(dfvar)
            perc <- tclvalue(percvar)
            param <- tclvalue(paramvar)
            paramb <- tclvalue(parambvar)
            steps <- tclvalue(stepsvar)
            iter <- tclvalue(itervar)
            outname <- tclvalue(outvar)
            substitute(range.expansion.out(data1 = dfvar, data2 = perc,  
                data4 = param, data5 = paramb, data6 = steps, data7 = iter, outfile = outname))
        }
        "reset" <- function() {
            tclvalue(dfvar) <- ""
            tkconfigure(dfnr.label, text = "")
            tclvalue(percvar) <- ""
            tclvalue(paramvar) <- ""
            tkconfigure(param.label, text = "")
            tclvalue(parambvar) <- ""
            tclvalue(stepsvar) <- ""
            tclvalue(itervar) <- ""
            tclvalue(outvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(sg, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(sg))
        tkgrid(cancel.but, submit.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(sg, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(sg, "<KeyPress-Return>", function() execcomp())
        tkbind(sg, "<KeyPress-Escape>", function() tkdestroy(sg))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(rangeM, "command", label = "Generate dispersal model", command = function() range.expansion.gui())
    range.to.raster.gui <- function(presences.map, re.out, mask.map, output) {
        result1 <- range_raster(presences.map = presences.map, re.out = re.out, mask.map = mask.map, 
            plot.directions = TRUE)
        exp.eval <- paste(output, "<<- result1", sep = "")
        eval(parse(text = exp.eval))
    }
    tkadd(rangeM, "command", label = "Compute raster from dispersal models", command = function() guiv(range.to.raster.gui, 
        exec = "Execute", title = "Estimate raster from dispersal models", argFilename = list(presences.map = NULL, 
            mask.map = NULL), argText = c(presences.map = "Raster file with occurrences", 
            re.out = "Object of class 'expansion'", mask.map = "Mask map", output = "Output name"), 
        helpsFunc = "range_raster"))
    tkadd(topMenu, "cascade", label = "Range expansion simulation", menu = rangeM)
    save.graphic.gui <- function() {
        tf <- tktoplevel()
        tkwm.title(tf, "Save graphic as ...")
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        frame2 <- tkframe(tf, relief = "groove", borderwidth = 2)
        frame3 <- tkframe(tf, relief = "groove", borderwidth = 2)
        devframe <- tkframe(frame2, relief = "groove", borderwidth = 2)
        done <- tclVar(0)
        formatvar <- tclVar(1)
        widthvar <- tclVar(32)
        heightvar <- tclVar(32)
        "savefig" <- function(formatvar, widthvar, heightvar) {
            outform <- tclvalue(formatvar)
            width <- as.numeric(tclvalue(widthvar))
            height <- as.numeric(tclvalue(heightvar))
            if (outform == 1) {
                filename <- tclvalue(tkgetSaveFile(initialfile = "MetaLandSimplot.tif", 
                  defaultextension = ".tif", title = "Save graphic as ...", filetypes = "{tiff {.tif .tiff}} {{All Files} {*.*}}"))
                if (filename != "") {
                  dev2bitmap(file = filename, type = "tiff24nc", width = width, height = height, 
                    units = "cm")
                }
            }
            else if (outform == 2) {
                filename <- tclvalue(tkgetSaveFile(initialfile = "MetaLandSimplot.png", 
                  defaultextension = ".png", title = "Save graphic as ...", filetypes = "{PNG {.png}} {{All Files} {*.*}}"))
                if (filename != "") {
                  dev2bitmap(file = filename, type = "png16m", width = width, height = height, 
                    units = "cm")
                }
            }
            else if (outform == 3) {
                filename <- tclvalue(tkgetSaveFile(initialfile = "MetaLandSimplot.jpeg", 
                  defaultextension = ".jpeg", title = "Save graphic as ...", filetypes = "{JPEG {.jpeg .jpg}} {{All Files} {*.*}}"))
                if (filename != "") {
                  dev2bitmap(file = filename, type = "jpeg", width = width, height = height, 
                    units = "cm")
                }
            }
            tkdestroy(tf)
        }
        tkgrid(tklabel(tf, text = "Save current graphic as ..."), columnspan = 2)
        tkgrid(tklabel(frame2, text = "Output format: "), sticky = "n")
        tkgrid(tkradiobutton(frame2, text = "tiff", value = 1, variable = formatvar), 
            sticky = "w")
        tkgrid(tkradiobutton(frame2, text = "png", value = 2, variable = formatvar), 
            sticky = "w")
        tkgrid(tkradiobutton(frame2, text = "jpeg", value = 3, variable = formatvar), 
            sticky = "w")
        tkgrid(frame2, rowspan = 2, sticky = "n")
        tkgrid(tklabel(frame3, text = "Output size: "))
        width.entry <- tkentry(frame3, textvariable = widthvar, width = 10)
        height.entry <- tkentry(frame3, textvariable = heightvar, width = 10)
        tkgrid(tklabel(frame3, text = "Width: "), width.entry)
        tkgrid(tklabel(frame3, text = "Height: "), height.entry)
        tkgrid(frame3, column = 1, row = 1, sticky = "n")
        save.but <- tkbutton(frame1, text = "Save as ...", command = function() savefig(formatvar, 
            widthvar, heightvar))
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        tkgrid(save.but, cancel.but)
        tkgrid(frame1, column = 1, row = 2, sticky = "n")
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() savefig(formatvar, widthvar, heightvar))
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        tkdestroy(tf)
    }
    tkadd(viewM, "command", label = "Save graphic as ...", command = function() save.graphic.gui())
    choose.multi <- function(multi.entry, multi.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (class(xobj) == "metapopulation" | class(xobj) == "landscape") {
                tkinsert(tlb, "end", x1)
                cbind(length(xobj))
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(multi.entry, 0, "end")
        tkinsert(multi.entry, "end", choix)
        tkconfigure(multi.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    plot_graph.gui <- function() {
        sg <- tktoplevel()
        tkwm.title(sg, "Plot single landscape")
        IFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        multivar <- tclVar()
        multi.entry <- tkentry(IFrame, textvariable = multivar)
        multi.label <- tklabel(IFrame, width = 5)
        choosemulti.but <- tkbutton(IFrame, text = "Set", command = function() choose.multi(multi.entry, 
            multi.label))
        tkgrid(tklabel(IFrame, text = "Object of class 'landscape' or 'metapopulation'"), 
            multi.entry, choosemulti.but, multi.label, sticky = "w")
        tkgrid(IFrame)
        LFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(LFrame, text = "Show links between patches"))
        loadvar <- tclVar(1)
        tkgrid(tkradiobutton(LFrame, text = "TRUE (not for 'metapopulation' class)", 
            value = 1, variable = loadvar), columnspan = 3)
        tkgrid(tkradiobutton(LFrame, text = "FALSE", value = 2, variable = loadvar), 
            columnspan = 3)
        tkgrid(LFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            data1 <- parse(text = tclvalue(multivar))[[1]]
            test.class <- get(tclvalue(multivar), envir = globalenv())
            if (class(test.class) == "metapopulation") {
                data2 <- TRUE
            }
            else data2 <- FALSE
            links <- c(TRUE, FALSE)
            data3 <- links[as.numeric(tclvalue(loadvar))]
            substitute(plot_graph(rl = data1, species = data2, links = data3))
        }
        "reset" <- function() {
            tclvalue(lsvar) <- ""
            tkconfigure(ls.label, text = "")
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(sg, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(sg))
        tkgrid(submit.but, reset.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(sg, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(sg, "<KeyPress-Return>", function() execcomp())
        tkbind(sg, "<KeyPress-Escape>", function() tkdestroy(sg))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(viewM, "command", label = "Plot single landscape", command = function() plot_graph.gui())
    choose.list <- function(list.entry, list.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (is.list(xobj)) {
                tkinsert(tlb, "end", x1)
                cbind(length(xobj))
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(list.entry, 0, "end")
        tkinsert(list.entry, "end", choix)
        tkconfigure(list.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    plotL.graph.gui <- function() {
        sg <- tktoplevel()
        tkwm.title(sg, "Plot landscape from list")
        IFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        lsvar <- tclVar()
        df.entry <- tkentry(IFrame, textvariable = lsvar)
        dfnr.label <- tklabel(IFrame, width = 5)
        choosels.but <- tkbutton(IFrame, text = "Set", command = function() choose.landscape(df.entry, 
            dfnr.label))
        tkgrid(tklabel(IFrame, text = "Object of class landscape"), df.entry, choosels.but, 
            dfnr.label, sticky = "w")
        listvar <- tclVar()
        list.entry <- tkentry(IFrame, textvariable = listvar)
        list.label <- tklabel(IFrame, width = 5)
        chooselist.but <- tkbutton(IFrame, text = "Set", command = function() choose.list(list.entry, 
            list.label))
        tkgrid(tklabel(IFrame, text = "List returned by span.graph"), list.entry, 
            chooselist.but, list.label, sticky = "w")
        nrvar <- tclVar()
        nr.entry <- tkentry(IFrame, textvariable = nrvar)
        tkgrid(tklabel(IFrame, text = "Time step position (numeric)"), nr.entry, 
            sticky = "w")
        tkgrid(IFrame)
        LFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(LFrame, text = "Show links between patches"))
        loadvar <- tclVar(1)
        tkgrid(tkradiobutton(LFrame, text = "TRUE", value = 1, variable = loadvar), 
            columnspan = 3)
        tkgrid(tkradiobutton(LFrame, text = "FALSE", value = 2, variable = loadvar), 
            columnspan = 3)
        tkgrid(LFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            data1 <- parse(text = tclvalue(lsvar))[[1]]
            data2 <- parse(text = tclvalue(listvar))[[1]]
            data3 <- as.numeric(tclvalue(nrvar))
            links <- c(TRUE, FALSE)
            data4 <- links[as.numeric(tclvalue(loadvar))]
            substitute(plotL.graph(rl = data1, rlist = data2, nr = data3, species = FALSE, 
                links = data4))
        }
        "reset" <- function() {
            tclvalue(lsvar) <- ""
            tkconfigure(ls.label, text = "")
            tclvalue(listvar) <- ""
            tkconfigure(list.label, text = "")
            tclvalue(nrvar) <- ""
            tclvalue(loadvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(sg, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(sg))
        tkgrid(submit.but, reset.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(sg, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(sg, "<KeyPress-Return>", function() execcomp())
        tkbind(sg, "<KeyPress-Escape>", function() tkdestroy(sg))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(viewM, "command", label = "Plot landscape from list", command = function() plotL.graph.gui())
    choose.exp <- function(exp.entry, exp.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (class(xobj) == "expansion") {
                tkinsert(tlb, "end", x1)
                cbind(length(xobj))
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(exp.entry, 0, "end")
        tkinsert(exp.entry, "end", choix)
        tkconfigure(exp.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    plot_expansion.gui <- function() {
        sg <- tktoplevel()
        tkwm.title(sg, "Plot single landscape settings")
        IFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        expvar <- tclVar()
        exp.entry <- tkentry(IFrame, textvariable = expvar)
        exp.label <- tklabel(IFrame, width = 5)
        chooseexp.but <- tkbutton(IFrame, text = "Set", command = function() choose.exp(exp.entry, 
            exp.label))
        tkgrid(tklabel(IFrame, text = "Object of class 'expansion'"), exp.entry, 
            chooseexp.but, exp.label, sticky = "w")
        tkgrid(IFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            data1 <- parse(text = tclvalue(expvar))[[1]]
            substitute(plot_expansion(exp = data1))
        }
        "reset" <- function() {
            tclvalue(expvar) <- ""
            tkconfigure(exp.label, text = "")
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(sg, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(sg))
        tkgrid(submit.but, reset.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(sg, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(sg, "<KeyPress-Return>", function() execcomp())
        tkbind(sg, "<KeyPress-Escape>", function() tkdestroy(sg))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(viewM, "command", label = "Expansion graphics", command = function() plot_expansion.gui())
    choose.raster <- function(raster.entry, raster.label) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
        vnr <- NULL
        vnc <- NULL
        numi <- 1
        tlb <- tklistbox(tf)
        scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb, 
            ...))
        tkconfigure(tlb, yscrollcommand = function(...) tkset(scr, ...))
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        submit.but <- tkbutton(frame1, text = "Choose", default = "active", command = function() tclvalue(done) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1, side = "bottom")
        tkpack(tlb, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        obj <- ls(globalenv())
        flb <- function(x1) {
            xobj <- get(x1, envir = globalenv())
            if (class(xobj) == "RasterLayer") {
                tkinsert(tlb, "end", x1)
                cbind(length(xobj))
            }
        }
        v <- unlist(lapply(obj, flb))
        if (length(v) > 0) {
            vnr <- v[seq(from = 1, to = length(v), by = 1)]
        }
        tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(raster.entry, 0, "end")
        tkinsert(raster.entry, "end", choix)
        tkconfigure(raster.label, text = as.character(vnr[numi]))
        tkdestroy(tf)
    }
    plot.map.expansion <- function(map, name) {
        eval(parse(text = paste("plot(", map, ", main = name, xlab = 'longitude', ylab = 'latitude')", 
            sep = "")))
        eval(parse(text = paste("contour(", map, ", add = TRUE)", sep = "")))
    }
    plot.raster.gui <- function() {
        sg <- tktoplevel()
        tkwm.title(sg, "Plot expansion map settings")
        IFrame <- tkframe(sg, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        rastervar <- tclVar()
        raster.entry <- tkentry(IFrame, textvariable = rastervar)
        raster.label <- tklabel(IFrame, width = 5)
        chooseraster.but <- tkbutton(IFrame, text = "Set", command = function() choose.raster(raster.entry, 
            raster.label))
        tkgrid(tklabel(IFrame, text = "Expansion map"), raster.entry, chooseraster.but, 
            raster.label, sticky = "w")
        namevar <- tclVar()
        name.entry <- tkentry(IFrame, textvariable = namevar)
        tkgrid(tklabel(IFrame, text = "Title for map"), name.entry, sticky = "w")
        tkgrid(IFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            data1 <- tclvalue(rastervar)
            data2 <- tclvalue(namevar)
            substitute(plot.map.expansion(data1, data2))
        }
        "reset" <- function() {
            tclvalue(rastervar) <- ""
            tkconfigure(raster.label, text = "")
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(sg, relief = "groove")
        submit.but <- tkbutton(RCSFrame, text = "Execute", default = "active", command = function() execcomp())
        reset.but <- tkbutton(RCSFrame, text = " Reset ", default = "active", command = function() reset())
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(sg))
        tkgrid(submit.but, reset.but, cancel.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(sg, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(sg, "<KeyPress-Return>", function() execcomp())
        tkbind(sg, "<KeyPress-Escape>", function() tkdestroy(sg))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(viewM, "command", label = "Plot expansion map", command = function() plot.raster.gui())
    tkadd(topMenu, "cascade", label = "View tools", menu = viewM)
    about.message <- function() {
        AboutW <- tktoplevel(width = 700, height = 200)
        tktitle(AboutW) <- "About the package MetaLandSim"
        file.about <- paste(system.file(package = "MetaLandSim"), "/logo/logo2.gif", 
            sep = "")
        tcl("image", "create", "photo", "imageID", file = file.about)
        l <- ttklabel(AboutW, image = "imageID", compound = "image")
        tkpack(l)
        msg <- "MetaLandSim was developed as part of a PhD by Frederico Mestre entitled:\n 'Synergistic effects of climate change and habitat fragmentation on species range shifts and metapopulation persistence'.\n\n This work was funded by FEDER through COMPETE and by national funds through FCT, Fundacao para a Ciencia e a Tecnologia,\n under the project NETPERSIST (PTDC/AAG-MAA/3227/2012) and a PhD Fellowship (Ref. SFRH/BD/73768/2010).\n\n It was developed by:\n\n Mestre, F.*; Canovas, F.**; Pita, R.*; Mira, A.* and P. Beja***\n\n *Research Center in Biodiversity and Genetic Resources-Univ. Evora - CIBIO/InBio-UE\n\n **Center of Marine Sciences-Univ. Algarve - CCMAR\n\n ***Research Center in Biodiversity and Genetic Resources-Univ. Porto - CIBIO/InBio-UP\n\n PORTUGAL"
        tkpack(tklabel(AboutW, text = msg))
        Dismiss.but <- tkbutton(AboutW, text = "OK", command = function() tkdestroy(AboutW))
        tkpack(Dismiss.but)
    }
    tkadd(aboutM, "command", label = "About the package", command = function() about.message())
    tkadd(aboutM, "command", label = "Key references", command = function() tkmessageBox(title = "Key references", 
        message = "\n\tBogaert, J., Ceulemans, R., & Salvador-Van Eysenrode, D. (2004). Decision tree algorithm for detection of spatial \n\t\tprocesses in landscape transformation. Environmental Management, 33(1): 6273.\n\tHanski, I. (1994). A practical model of metapopulation dynamics. Journal of Animal Ecology, 63: 151-162. \n\tHanski, I. (1999). Metapopulation Ecology. Oxford University Press. 313 pp.\n\tHanski, I. & Gaggiotti, O.E. (Eds.) (2004) Ecology, Genetics, and Evolution of Metapopulations. Elsevier Academic\n\t\tPress. 696 pp.\n\tHanski, I., & Ovaskainen, O. (2000). The metapopulation capacity of a fragmented landscape. Nature, 404: 755-758.\n\tMinor, E. S., and Urban, D. L. (2007). Graph theory as a proxy for spatially explicit population models in \n\t\tconservation planning. Ecological Applications, 17(6): 1771-1782.\n\tMinor, E. S., and Urban, D. L. (2008). A Graph-Theory Framework for Evaluating Landscape Connectivity and \n\t\tConservation Planning. Conservation Biology, 22(2): 297-307.\n\tMoilanen, A., & Nieminen, M. (2002). Simple connectivity measures in spatial ecology. Ecology, 83(4): 1131-1145.\n\tMoilanen, A. (2004). SPOMSIM: software for stochastic patch occupancy models of metapopulation dynamics. Ecological\n\t\tModelling, 179(4): 533-550.\n\tOvaskainen, O. & Hanski, I. (2001). Spatially structured metapopulation models: global and local assessment of \n\t\tmetapopulation capacity. Theoretical Population Biology, 60(4), 281-302.\n\tOvaskainen, O., & Hanski, I. (2002). Transient dynamics in metapopulation response to perturbation. Theoretical \n\t\tPopulation Biology, 61(3): 285-295.\n\tPascual-Hortal, L., and Saura, S. (2006). Comparison and development of new graph-based landscape connectivity indices: \n\t\ttowards the prioritization of habitat patches and corridors for conservation. Landscape Ecology, \n\t\t21(7): 959-967.\n\tSaura, S., and Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity in landscape \n\t\tconservation planning: comparison with existing indices and application to a case study. Landscape and Urban \n\t\tPlanning, 83(2): 91-103.\n\tShaw, M.W., (1994). Simulation of population expansion and spatial pattern when individual dispersal distributions do \n\t\tnot decline exponentially with distance. Proc. R. Soc. London B: 259, 243-248.\n\tUrban, D., and Keitt, T. (2001). Landscape connectivity: a graph-theoretic perspective. Ecology, 82(5): 1205-1218.\n    ", 
        icon = "info", type = "ok"))
    tkadd(aboutM, "command", label = "Manual", command = function() openPDF(file = paste(system.file(package = "MetaLandSim"), 
        "/doc/MetaLandSim.pdf", sep = "")))
    tkadd(aboutM, "command", label = "Vignette on Model Parameter Estimation", command = function() openPDF(file = paste(system.file(package = "MetaLandSim"), 
        "/doc/parametrization.pdf", sep = "")))
    tkadd(aboutM, "command", label = "Vignette on Landscape Occupation Simulation in Dynamic Landscapes", 
        command = function() openPDF(file = paste(system.file(package = "MetaLandSim"), 
            "/doc/landscape_simulation.pdf", sep = "")))
    tkadd(aboutM, "command", label = "Vignette on Range Expansion Simulation", command = function() openPDF(file = paste(system.file(package = "MetaLandSim"), 
        "/doc/range_expansion.pdf", sep = "")))
	tkadd(aboutM, "command", label = "Reference Paper", command = function() tkmessageBox(title = "Reference Paper", message = "Mestre,F.; Canovas,F.; Pita, R.; Mira, A. and Beja, P. (2016). 'An R Package for simulating metapopulation dynamics and range expansion under environmental change', Environmental Modelling and Software, 81:40-44.", icon = "info", type = "ok"))
    tkadd(topMenu, "cascade", label = "About", menu = aboutM)
}