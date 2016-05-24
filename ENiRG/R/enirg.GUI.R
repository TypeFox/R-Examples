enirg.GUI <-
function() {
    mainW <- tktoplevel(width = 700, height = 20)
    loading <- tclVar()
    file.back <- paste(system.file(package='ENiRG'), '/img/ENiRG.gif', sep='')
    tcl("image", "create", "photo", loading, file = file.back)
    background.logo <- tklabel(mainW, image = loading)
    tkpack(background.logo)
    tktitle(mainW) <- "ENiRG - Ecological Niche in R-GRASS"
    topMenu <- tkmenu(mainW)
    tkconfigure(mainW, menu = topMenu)
    DataM <- tkmenu(topMenu, tearoff = FALSE)
    ImportS <- tkmenu(DataM, tearoff = FALSE)
    ExampleS <- tkmenu(ImportS, tearoff = FALSE)
    AnalysisM <- tkmenu(topMenu, tearoff = FALSE)
    ViewM <- tkmenu(topMenu, tearoff = FALSE)
    AboutM <- tkmenu(topMenu, tearoff = FALSE)
    Newsession<- function(gisBase, new) {
      initGRASS(gisBase = gisBase, home = tempdir(), override = as.logical(new))
    }
    tkadd(DataM, "command", label = "New GRASS session ...", command = function() guiv(Newsession, 
          exec = "Execute",
          argText = c(gisBase = "Path to GRASS binaries and libraries:", 
            new = "override GRASS session (TRUE or FALSE):"), 
          helpsFunc = "initGRASS"))
    tkadd(DataM, "cascade", label = "Import from ...", menu = ImportS)
    import.text.gui <- function(filename, long.col, lat.col, presence.col, output)
      {
         temp.table <- read.table(filename, sep = "\t", header = TRUE)
         data.table <- data.frame(longitude = temp.table[, long.col], latitude = temp.table[, lat.col],
                  presences = temp.table[, presence.col])
         exp.eval <- paste(output," <<- data.table", sep="")
         eval(parse(text = exp.eval))
      }
    importtxtCallback <- function(arg)
      {
        if( arg == "filename" )
          {
            columnames <- names(read.table(guiGetValue("filename"), sep = "\t", header= TRUE))
            print(columnames)
            guiSet("columnames", columnames)
            setListElements("long.col", columnames)
            setListElements("lat.col", columnames)
            setListElements("presence.col", columnames)            
          }
      }
    tkadd(ImportS, "command", label = "tab separated text columns",command = function() guiv(import.text.gui,
          exec="Execute",
		  title="Import from tab separated text columns",
          argFilename=list(filename=NULL),
          argList=list(long.col=NULL, lat.col=NULL, presence.col=NULL),
          callback=importtxtCallback,
          argText=c(filename="Chose text file:",
          long.col="Name of the column with longitude/x:",
          lat.col="Name of the column with latitude/y:",
          presence.col="Name of column with presence/abundance:",
          output="Output name")))
    import.excel.gui <- function(filename, long.col, lat.col, presence.col, output)
      {
         temp.table <- read.xls(filename)
         data.table <- data.frame(longitude = temp.table[, long.col], latitude = temp.table[, lat.col],
                  presences = temp.table[, presence.col])
         exp.eval <- paste(output," <<- data.table", sep="")
         eval(parse(text = exp.eval))
      }
    importxlsCallback <- function(arg)
      {
        if( arg == "filename" )
          {
            columnames <- names(read.xls(guiGetValue("filename")))
            print(columnames)
            guiSet("columnames", columnames)
            setListElements("long.col", columnames)
            setListElements("lat.col", columnames)
            setListElements("presence.col", columnames)            
          }
      }
    tkadd(ImportS, "command", label = "excel",command = function() guiv(import.excel.gui,
          exec="Execute",
		  title="Import from excel",
          argFilename=list(filename=NULL),
          argList=list(long.col=NULL, lat.col=NULL, presence.col=NULL),
          callback=importxlsCallback,
          argText=c(filename="Chose excel file:",
          long.col="Name of the column with longitude/x:",
          lat.col="Name of the column with latitude/y:",
          presence.col="Name of column with presence/abundance:",
          output="Output name")))
    tkadd(DataM, "command", label = "Import egvs", command = function() guiv(import.egvs, 
        exec = "Execute", argFilename = list(filenames = NULL),
        argText = c(filenames = "Chose your raster file (check GDAL library for supported formats):", 
        output.names = "Output names for imported maps separated by a comma:"), 
        helpsFunc = "import.egvs"))
    tkadd(DataM, "command", label = "Map's list", command = function() gui(list.maps, 
        exec = "Get", argText = c(prefix = "Filtering pattern (regular expressions):")
        , helpsFunc = "list.maps"))
    load.grass <- function() {
        grassenv <- tktoplevel()
        tktitle(grassenv) <- "Load GRASS mapset info"
        available.maps <- execGRASS("g.list", type = "raster", intern = TRUE)
        tl <- tklistbox(grassenv, height = length(available.maps), selectmode = "single", 
            background = "white")
        tkgrid(tklabel(grassenv, text = "Choose a map from the list:"))
        tkgrid(tl)
        for (i in (1:length(available.maps))) {
            tkinsert(tl, "end", available.maps[i])
        }
        OnOK <- function() {
            mapChoice <- available.maps[as.numeric(tkcurselection(tl)) + 1]
            tkdestroy(grassenv)
            metadata <- map.info(mapChoice)
            msg <- paste(attr(metadata, "resOut"), collapse = "\n")
            metadataW <- tktoplevel(width = 700, height = 300)
            tktitle(metadataW) <- "Metadata"
            tkgrid(tklabel(metadataW, text = msg, font = tkfont.create(family = "courier", 
                size = 12, weight = "bold")))
            Dismiss.but <- tkbutton(metadataW, text = "Dismiss", command = function() tkdestroy(metadataW))
            tkgrid(Dismiss.but)
        }
        OK.but <- tkbutton(grassenv, text = "  Get  ", command = OnOK)
        Cancel.but <- tkbutton(grassenv, text = "Dismiss", command = function() tkdestroy(grassenv))
        tkgrid(OK.but, Cancel.but)
    }
    tkadd(DataM, "command", label = "Map info", command = function() load.grass())
    stdz.maps.gui <- function() {
        grassenv <- tktoplevel()
        tktitle(grassenv) <- "Load GRASS mapset info"
        available.maps <- execGRASS("g.list", type = "raster", intern = TRUE)
        tl <- tklistbox(grassenv, height = length(available.maps), selectmode = "multiple", 
            background = "white")
        tkgrid(tklabel(grassenv, text = "Choose map/s from the list:"))
        tkgrid(tl)
        for (i in (1:length(available.maps))) {
            tkinsert(tl, "end", available.maps[i])
        }
        OnOK <- function() {
            mapChoice <- available.maps[as.numeric(tkcurselection(tl)) + 1]
            tkdestroy(grassenv)
            stdz.maps(mapChoice)
            tkmessageBox(title = "info", message = "Standardization completed successfuly!", 
                icon = "info", type = "ok")
        }
        OK.but <- tkbutton(grassenv, text = "Execute", command = OnOK)
        Cancel.but <- tkbutton(grassenv, text = "Dismiss", command = function() tkdestroy(grassenv))
        tkgrid(OK.but, Cancel.but)
    }
    tkadd(DataM, "command", label = "Standardization", command = function() stdz.maps.gui())
    tkadd(DataM, "command", label = "Quit", command = function() tkdestroy(mainW))
    tkadd(topMenu, "cascade", label = "Data", menu = DataM)
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
    choose.maps <- function(qt.entry, qt.label) {
        qtmap <- tktoplevel()
        tkwm.title(qtmap, "Choose:")
        done.map <- tclVar(0)
        available.maps <- execGRASS("g.list", type = "raster", intern = TRUE)
        vnc.map <- NULL
        numi.map <- 1
        tlb.map <- tklistbox(qtmap, selectmode = "multiple")
        scr.map <- tkscrollbar(qtmap, repeatinterval = 5, command = function(...) tkyview(tlb.map, 
            ...))
        tkconfigure(tlb.map, yscrollcommand = function(...) tkset(scr.map, ...))
        frame1.qt <- tkframe(qtmap, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1.qt, text = "Dismiss", command = function() tkdestroy(qtmap))
        submit.but <- tkbutton(frame1.qt, text = "Choose", default = "active", command = function() tclvalue(done.map) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1.qt, side = "bottom")
        tkpack(tlb.map, side = "left", fill = "both", expand = TRUE)
        tkpack(scr.map, side = "right", fill = "y")
        for (i in (1:length(available.maps))) {
            tkinsert(tlb.map, "end", available.maps[i])
        }
        tkbind(tlb.map, "<Double-ButtonPress-1>", function() tclvalue(done.map) <- 1)
        tkbind(qtmap, "<Destroy>", function() tclvalue(done.map) <- 2)
        tkbind(qtmap, "<KeyPress-Return>", function() tclvalue(done.map) <- 1)
        tkbind(qtmap, "<KeyPress-Escape>", function() tkdestroy(qtmap))
        tkwait.variable(done.map)
        if (tclvalue(done.map) == "2") 
            return(0)
        numc.map <- tclvalue(tkcurselection(tlb.map))
        if (numc.map == "") {
            tkdestroy(qtmap)
            return(0)
        }
        choix.maps <- available.maps[as.numeric(strsplit(numc.map, " ")[[1]]) + 1]
        tkdelete(qt.entry, 0, "end")
        tkinsert(qt.entry, "end", choix.maps)
        tkconfigure(qt.label, text = length(choix.maps))
        tkdestroy(qtmap)
    }
    choose.ql.maps <- function(ql.entry, ql.label) {
        qlmap <- tktoplevel()
        tkwm.title(qlmap, "Choose:")
        done.ql.map <- tclVar(0)
        available.maps <- execGRASS("g.list", type = "raster", intern = TRUE)
        vnc.ql.map <- NULL
        numi.ql.map <- 1
        tlb.ql.map <- tklistbox(qlmap, selectmode = "multiple")
        scr.ql.map <- tkscrollbar(qlmap, repeatinterval = 5, command = function(...) tkyview(tlb.ql.map, 
            ...))
        tkconfigure(tlb.ql.map, yscrollcommand = function(...) tkset(scr.ql.map, 
            ...))
        frame1.ql <- tkframe(qlmap, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1.ql, text = "Dismiss", command = function() tkdestroy(qlmap))
        submit.but <- tkbutton(frame1.ql, text = "Choose", default = "active", command = function() tclvalue(done.ql.map) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1.ql, side = "bottom")
        tkpack(tlb.ql.map, side = "left", fill = "both", expand = TRUE)
        tkpack(scr.ql.map, side = "right", fill = "y")
        for (i in (1:length(available.maps))) {
            tkinsert(tlb.ql.map, "end", available.maps[i])
        }
        tkbind(tlb.ql.map, "<Double-ButtonPress-1>", function() tclvalue(done.ql.map) <- 1)
        tkbind(qlmap, "<Destroy>", function() tclvalue(done.ql.map) <- 2)
        tkbind(qlmap, "<KeyPress-Return>", function() tclvalue(done.ql.map) <- 1)
        tkbind(qlmap, "<KeyPress-Escape>", function() tkdestroy(qlmap))
        tkwait.variable(done.ql.map)
        if (tclvalue(done.ql.map) == "2") 
            return(0)
        numc.ql.map <- tclvalue(tkcurselection(tlb.ql.map))
        if (numc.ql.map == "") {
            tkdestroy(qlmap)
            return(0)
        }
        choix.ql.maps <- available.maps[as.numeric(strsplit(numc.ql.map, " ")[[1]]) + 
            1]
        tkdelete(ql.entry, 0, "end")
        tkinsert(ql.entry, "end", choix.ql.maps)
        tkconfigure(ql.label, text = length(choix.ql.maps))
        tkdestroy(qlmap)
    }
    enirg.gui.out <- function(data1, data2, data3, data4, outfile) {
        exp.eval <- paste(outfile,
            " <<- enirg(presences.table = ", data1,
            ", qtegv.maps = data2",
            ", qlegv.maps = data3",
            ", col.w = NULL, scannf = TRUE, load.maps = TRUE",
            ", species.name = data4, method = 'normal')",
            sep="")
        eval(parse(text = exp.eval))
    }
    enirg.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "enirg settings")
        IFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- input -", foreground = "blue"), columnspan = 5)
        dfvar <- tclVar()
        df.entry <- tkentry(IFrame, textvariable = dfvar)
        dfnr.label <- tklabel(IFrame, width = 5)
        choosevect.but <- tkbutton(IFrame, text = "Set", command = function() choose.df(df.entry, 
            dfnr.label))
        tkgrid(tklabel(IFrame, text = "Observations' data frame: "), df.entry, choosevect.but, 
            dfnr.label, sticky = "w")
        tkgrid(IFrame)
        OFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for enirg output: "), out.entry, sticky = "w")
        speciesvar <- tclVar()
        species.entry <- tkentry(OFrame, textvariable = speciesvar)
        tkgrid(tklabel(OFrame, text = "Species name for map names: "), species.entry, 
            sticky = "w")
        tkgrid(OFrame)
        QTFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(QTFrame, text = "Quantitative map/s", foreground = "blue"), 
            columnspan = 5)
        qtvar <- tclVar()
        qt.entry <- tkentry(QTFrame, textvariable = qtvar)
        qt.label <- tklabel(QTFrame, width = 5)
        chooseqt.but <- tkbutton(QTFrame, text = "Set", command = function() choose.maps(qt.entry, 
            qt.label))
        tkgrid(tklabel(QTFrame, text = "Number of map/s: "), qt.label, qt.entry, 
            chooseqt.but, sticky = "w")
        tkgrid(QTFrame)
        QLFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(QLFrame, text = "Qualitative map/s", foreground = "blue"), 
            columnspan = 5)
        qlvar <- tclVar()
        ql.entry <- tkentry(QLFrame, textvariable = qlvar)
        ql.label <- tklabel(QLFrame, width = 5)
        chooseql.but <- tkbutton(QLFrame, text = "Set", command = function() choose.ql.maps(ql.entry, 
            ql.label))
        tkgrid(tklabel(QLFrame, text = "Number of map/s: "), ql.label, ql.entry, 
            chooseql.but, sticky = "w")
        tkgrid(QLFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- tclvalue(dfvar)
            if (tclvalue(qtvar) != "") {
                vect2 <- strsplit(tclvalue(qtvar), " ")[[1]]
            }
            else vect2 <- NULL
            if (tclvalue(qlvar) != "") {
                vect3 <- strsplit(tclvalue(qlvar), " ")[[1]]
            }
            else vect3 <- NULL
            if (tclvalue(speciesvar) != "") {
                spn <- tclvalue(speciesvar)
            }
            else spn <- "Species"
            if (tclvalue(outvar) != "") {
                outname <- tclvalue(outvar)
            }
            else outname <- "enirg.out"
            substitute(enirg.gui.out(data1 = vect1, data2 = vect2, data3 = vect3, 
                data4 = spn, outfile = outname))
        }
        "reset" <- function() {
            tclvalue(dfvar) <- ""
            tkconfigure(dfnr.label, text = "")
            tclvalue(speciesvar) <- ""
            tclvalue(outvar) <- ""
            tclvalue(qlvar) <- ""
            tkconfigure(ql.label, text = "")
            tclvalue(qtvar) <- ""
            tkconfigure(qt.label, text = "")
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
    tkadd(AnalysisM, "command", label = "enirg", command = function() enirg.gui())
    enirg.predict.out <- function(data1, data2, data3, data4, outfile) {
        exp.eval <- paste(outfile,
            " <<- enirg.predict(enirg.results = ", data1,
            ", qtegv.maps = data2",
            ", qlegv.maps = data3",
            ", prediction.name = data4",
            ", method = 'normal', load.map = TRUE)",
            sep="")
        eval(parse(text = exp.eval))
    }
    enirg.predict.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "Prediction settings")
        IFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- Input -", foreground = "blue"), columnspan = 5)
        enirgvar <- tclVar()
        enirg.entry <- tkentry(IFrame, textvariable = enirgvar)
        choose.enirg.but <- tkbutton(IFrame, text = "Set", command = function() choose.enirg(enirg.entry))
        tkgrid(tklabel(IFrame, text = "enirg object: "), enirg.entry, choose.enirg.but, 
            sticky = "w")
        tkgrid(IFrame)
        OFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- Output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for prediction output: "), out.entry, 
            sticky = "w")
        outnamevar <- tclVar()
        outname.entry <- tkentry(OFrame, textvariable = outnamevar)
        tkgrid(tklabel(OFrame, text = "Name for prediction maps: "), outname.entry, 
            sticky = "w")
        tkgrid(OFrame)
        QTFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(QTFrame, text = "Quantitative map/s", foreground = "blue"), 
            columnspan = 5)
        qtvar <- tclVar()
        qt.entry <- tkentry(QTFrame, textvariable = qtvar)
        qt.label <- tklabel(QTFrame, width = 5)
        chooseqt.but <- tkbutton(QTFrame, text = "Set", command = function() choose.maps(qt.entry, 
            qt.label))
        tkgrid(tklabel(QTFrame, text = "Number of map/s: "), qt.label, qt.entry, 
            chooseqt.but, sticky = "w")
        tkgrid(QTFrame)
        QLFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(QLFrame, text = "Qualitative map/s", foreground = "blue"), 
            columnspan = 5)
        qlvar <- tclVar()
        ql.entry <- tkentry(QLFrame, textvariable = qlvar)
        ql.label <- tklabel(QLFrame, width = 5)
        chooseql.but <- tkbutton(QLFrame, text = "Set", command = function() choose.ql.maps(ql.entry, 
            ql.label))
        tkgrid(tklabel(QLFrame, text = "Number of map/s: "), ql.label, ql.entry, 
            chooseql.but, sticky = "w")
        tkgrid(QLFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- tclvalue(enirgvar)
            if (tclvalue(qtvar) != "") {
                vect2 <- strsplit(tclvalue(qtvar), " ")[[1]]
            }
            else vect2 <- NULL
            if (tclvalue(qlvar) != "") {
                vect3 <- strsplit(tclvalue(qlvar), " ")[[1]]
            }
            else vect3 <- NULL
            if (tclvalue(outvar) != "") {
                outname <- tclvalue(outvar)
            }
            else outname <- "enirg.pred"
            if (tclvalue(outnamevar) == "") {
                outnamemap <- NULL
            }
            else outnamemap <- tclvalue(outnamevar)
            substitute(enirg.predict.out(data1 = vect1, data2 = vect2, data3 = vect3, 
                data4 = outnamemap, outfile = outname))
        }
        "reset" <- function() {
            tclvalue(enirgvar) <- ""
            tclvalue(outvar) <- ""
            tclvalue(outnamevar) <- ""
            tclvalue(qlvar) <- ""
            tkconfigure(ql.label, text = "")
            tclvalue(qtvar) <- ""
            tkconfigure(qt.label, text = "")
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
    tkadd(AnalysisM, "command", label = "enirg prediction", command = function() enirg.predict.gui())
    choose.hsm <- function(hsm.entry) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
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
            if (class(xobj) == "hsm") {
                tkinsert(tlb, "end", x1)
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
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(hsm.entry, 0, "end")
        tkinsert(hsm.entry, "end", choix)
        tkdestroy(tf)
    }
    pre.boyce <- function(results, outcat) {
        dataobs <- results[["predictions"]]
        dataobs <- dataobs$predicted
        datavect <- results[[4]]@data@values
        datavect <- datavect[which(datavect > 0)]
        boyce(prediction = dataobs, prediction.map = datavect, outcat = outcat)
    }
    boyce.GUI <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "Classification settings")
        hsmvar <- tclVar()
        IOFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IOFrame, text = "- Input -", foreground = "blue"), columnspan = 5)
        hsm.entry <- tkentry(IOFrame, textvariable = hsmvar)
        choosehsm.but <- tkbutton(IOFrame, text = "Set", command = function() choose.hsm(hsm.entry))
        tkgrid(tklabel(IOFrame, text = "Prediction object: "), hsm.entry, choosehsm.but, sticky = "w")
        outvar <- tclVar()
        cat.entry <- tkentry(IOFrame, textvariable = outvar)
        tkgrid(tklabel(IOFrame, text = "Name to store classification info: "), cat.entry, 
            sticky = "w")
        tkgrid(IOFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- parse(text = tclvalue(hsmvar))[[1]]
            outpn <- tclvalue(outvar)
            substitute(pre.boyce(results = vect1, outcat = outpn))
        }
        "reset" <- function() {
            tclvalue(hsmvar) <- ""
            tclvalue(outvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(tt, relief = "groove")
        reset.but <- tkbutton(RCSFrame, text = " Reset ", command = reset)
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(tt))
        submit.but <- tkbutton(RCSFrame, text = "  Set  ", default = "active", command = function() execcomp())
        tkgrid(cancel.but, submit.but, reset.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tt, "<KeyPress-Return>", function() execcomp())
        tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(AnalysisM, "command", label = "Boyce classification", command = function() boyce.GUI())
    choose.map <- function(map.entry, map.label) {
        qtmap <- tktoplevel()
        tkwm.title(qtmap, "Choose:")
        done.map <- tclVar(0)
        available.maps <- execGRASS("g.list", type = "raster", intern = TRUE)
        vnc.map <- NULL
        numi.map <- 1
        tlb.map <- tklistbox(qtmap, selectmode = "multiple")
        scr.map <- tkscrollbar(qtmap, repeatinterval = 5, command = function(...) tkyview(tlb.map, 
            ...))
        tkconfigure(tlb.map, yscrollcommand = function(...) tkset(scr.map, ...))
        frame1.qt <- tkframe(qtmap, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(frame1.qt, text = "Dismiss", command = function() tkdestroy(qtmap))
        submit.but <- tkbutton(frame1.qt, text = "Choose", default = "active", command = function() tclvalue(done.map) <- 1)
        tkpack(cancel.but, submit.but, side = "left")
        tkpack(frame1.qt, side = "bottom")
        tkpack(tlb.map, side = "left", fill = "both", expand = TRUE)
        tkpack(scr.map, side = "right", fill = "y")
        for (i in (1:length(available.maps))) {
            tkinsert(tlb.map, "end", available.maps[i])
        }
        tkbind(tlb.map, "<Double-ButtonPress-1>", function() tclvalue(done.map) <- 1)
        tkbind(qtmap, "<Destroy>", function() tclvalue(done.map) <- 2)
        tkbind(qtmap, "<KeyPress-Return>", function() tclvalue(done.map) <- 1)
        tkbind(qtmap, "<KeyPress-Escape>", function() tkdestroy(qtmap))
        tkwait.variable(done.map)
        if (tclvalue(done.map) == "2") 
            return(0)
        numc.map <- tclvalue(tkcurselection(tlb.map))
        if (numc.map == "") {
            tkdestroy(qtmap)
            return(0)
        }
        choix.maps <- available.maps[as.numeric(strsplit(numc.map, " ")[[1]]) + 1]
        tkdelete(map.entry, 0, "end")
        tkinsert(map.entry, "end", choix.maps)
        tkconfigure(map.label, text = length(choix.maps))
        tkdestroy(qtmap)
    }
    choose.cbi <- function(cbi.entry) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
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
            if (class(xobj) == "CBI") {
                tkinsert(tlb, "end", x1)
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
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(cbi.entry, 0, "end")
        tkinsert(cbi.entry, "end", choix)
        tkdestroy(tf)
    }
    classify.map.out <- function(data1, data2, data3, outfile) {
        exp.eval <- paste(outfile,
            " <<- classify.map(map = data1",
            ", suit.classes = data2",
            ", output.name = data3",
            ", load.map = TRUE)",
            sep="")
        eval(parse(text = exp.eval))
    }
    classify.map.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "Classification settings")
        IFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- Input -", foreground = "blue"), columnspan = 5)
        mapvar <- tclVar()
        map.label <- tklabel(IFrame, width = 2)
        map.entry <- tkentry(IFrame, textvariable = mapvar)
        choose.map.but <- tkbutton(IFrame, text = "Set", command = function() choose.map(map.entry, map.label))
        tkgrid(tklabel(IFrame, text = "Suitability map: "), map.entry, choose.map.but, 
            sticky = "w")
        tkgrid(IFrame)
        OFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(OFrame, text = "- Output -", foreground = "blue"), columnspan = 5)
        outvar <- tclVar()
        out.entry <- tkentry(OFrame, textvariable = outvar)
        tkgrid(tklabel(OFrame, text = "Name for classified map: "), out.entry, 
            sticky = "w")
        outfilevar <- tclVar()
        outfile.entry <- tkentry(OFrame, textvariable = outfilevar)
        tkgrid(tklabel(OFrame, text = "Name for output: "), outfile.entry, sticky = "w")
        tkgrid(OFrame)
        QTFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(QTFrame, text = "Suitability intervals: ", foreground = "blue"), 
            columnspan = 5)
        cbivar <- tclVar()
        cbi.entry <- tkentry(QTFrame, textvariable = cbivar)
        choose.cbi.but <- tkbutton(QTFrame, text = "Set", command = function() choose.cbi(cbi.entry))
        tkgrid(tklabel(QTFrame, text = "CBI object: "), cbi.entry, choose.cbi.but, 
            sticky = "w")
        tkgrid(QTFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- tclvalue(mapvar)
            vect2 <- parse(text = tclvalue(cbivar))[[1]]
            if (tclvalue(outvar) != "") 
                vect3 <- tclvalue(outvar)
            else vect3 <- NULL
            outname <- tclvalue(outfilevar)
            substitute(classify.map.out(data1 = vect1, data2 = vect2, data3 = vect3, 
                outfile = outname))
        }
        "reset" <- function() {
            tclvalue(mapvar) <- ""
            tclvalue(outvar) <- ""
            tclvalue(outfilevar) <- ""
            tclvalue(cbivar) <- ""
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
    tkadd(AnalysisM, "command", label = "Classify map", command = function() classify.map.gui())
    tkadd(topMenu, "cascade", label = "Analysis", menu = AnalysisM)
    choose.enirg <- function(enirg.entry) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
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
            if (class(xobj) == "enirg") {
                tkinsert(tlb, "end", x1)
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
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(enirg.entry, 0, "end")
        tkinsert(enirg.entry, "end", choix)
        tkdestroy(tf)
    }
    enirg.plot.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "plot settings")
        IOFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IOFrame, text = "- Input -", foreground = "blue"), columnspan = 5)
        enirgvar <- tclVar()
        enirg.entry <- tkentry(IOFrame, textvariable = enirgvar)
        choose.enirg.but <- tkbutton(IOFrame, text = "Set", command = function() choose.enirg(enirg.entry))
        tkgrid(tklabel(IOFrame, text = "enirg object: "), enirg.entry, choose.enirg.but, 
            sticky = "w")
        tkgrid(IOFrame)
        LFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(LFrame, text = "Method: "))
        loadvar <- tclVar(1)
        tkgrid(tkradiobutton(LFrame, text = "simplified", value = 1, variable = loadvar), 
            columnspan = 3)
        tkgrid(tkradiobutton(LFrame, text = "extended", value = 2, variable = loadvar), 
            columnspan = 3)
        tkgrid(LFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- parse(text = tclvalue(enirgvar))[[1]]
            print(vect1)
            available <- c("simplified", "extended")
            Choices <- available[as.numeric(tclvalue(loadvar))]
            print(Choices)
            substitute(enirg.plot(enirg.results = vect1, mar.col = "grey", spe.col = "black", 
                method = Choices))
        }
        "reset" <- function() {
            tclvalue(enirgvar) <- ""
            tclvalue(loadvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(tt, relief = "groove")
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(tt))
        submit.but <- tkbutton(RCSFrame, text = "  Plot  ", default = "active", command = function() execcomp())
        tkgrid(cancel.but, submit.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tt, "<KeyPress-Return>", function() execcomp())
        tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(ViewM, "command", label = "Plot enirg", command = function() enirg.plot.gui())
    choose.hsm.map <- function(hsm.entry) {
        tf <- tktoplevel()
        tkwm.title(tf, "Choose:")
        done <- tclVar(0)
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
            if (class(xobj) == "hsm") {
                tkinsert(tlb, "end", x1)
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
        numc <- tclvalue(tkcurselection(tlb))
        numi <- as.integer(numc) + 1
        if (numc == "") {
            tkdestroy(tf)
            return(0)
        }
        choix <- tclvalue(tkget(tlb, numc))
        tkdelete(hsm.entry, 0, "end")
        tkinsert(hsm.entry, "end", choix)
        tkdestroy(tf)
    }
    plot.hsm.map<- function(object.hsm) {
        if(length(class(object.hsm)) == 2) {
            plot(object.hsm[[1]], legend = FALSE, col = c("red", "yellow", "blue", "green"))
            tkmessageBox(title = "Locate legend.",
                message = "Please, click the left bottom of the mouse to choose the position of the legend, and then the right buttom to place it.", icon = "info", type = "ok")
            position <- locator()
            legend(position, legend = c("unsuitable", "marginal", "suitable", "optimal"), fill = c("red", "yellow", "blue", "green"))
        }
        else {
            plot(object.hsm[[1]], main = names(object.hsm)[1])
            contour(object.hsm[[1]], add=T)
        }
    }
    plot.hsm.map.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "plot settings")
        IOFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IOFrame, text = "- Input -", foreground = "blue"), columnspan = 5)
        hsmvar <- tclVar()
        hsm.entry <- tkentry(IOFrame, textvariable = hsmvar)
        choose.hsm.but <- tkbutton(IOFrame, text = "Set", command = function() choose.hsm.map(hsm.entry))
        tkgrid(tklabel(IOFrame, text = "hsm object: "), hsm.entry, choose.hsm.but, 
            sticky = "w")
        tkgrid(IOFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- parse(text = tclvalue(hsmvar))[[1]]
            substitute(plot.hsm.map(object.hsm = vect1))
        }
        "reset" <- function() {
            tclvalue(hsmvar) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(tt, relief = "groove")
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(tt))
        submit.but <- tkbutton(RCSFrame, text = "  Plot  ", default = "active", command = function() execcomp())
        tkgrid(cancel.but, submit.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tt, "<KeyPress-Return>", function() execcomp())
        tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(ViewM, "command", label = "Plot HSM map", command = function() plot.hsm.map.gui())
    plot.map<- function(mapname) {
        eval(parse(text = "plot(raster(readRAST(mapname)), main = mapname)"))
    }
    plot.map.gui <- function() {
        tt <- tktoplevel()
        tkwm.title(tt, "plot settings")
        IFrame <- tkframe(tt, relief = "groove", borderwidth = 2)
        tkgrid(tklabel(IFrame, text = "- Input -", foreground = "blue"), columnspan = 5)
        mapname <- tclVar()
        map.entry <- tkentry(IFrame, textvariable = mapname)
        map.label <- tklabel(IFrame, width = 2)
        choose.map.but <- tkbutton(IFrame, text = "Set", command = function() choose.map(map.entry, map.label))
        tkgrid(tklabel(IFrame, text = "Choose map: "), map.entry, choose.map.but, 
            sticky = "w")
        tkgrid(IFrame)
        vnr = NULL
        vnc = NULL
        numi = 1
        done <- tclVar(0)
        "build" <- function() {
            vect1 <- tclvalue(mapname)
            substitute(plot.map(vect1))
        }
        "reset" <- function() {
            tclvalue(mapname) <- ""
        }
        "execcomp" <- function() {
            cmd <- build()
            eval.parent(cmd)
        }
        RCSFrame <- tkframe(tt, relief = "groove")
        cancel.but <- tkbutton(RCSFrame, text = "Dismiss", command = function() tkdestroy(tt))
        submit.but <- tkbutton(RCSFrame, text = "  Plot  ", default = "active", command = function() execcomp())
        tkgrid(cancel.but, submit.but, ipadx = 20)
        tkgrid(RCSFrame)
        tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tt, "<KeyPress-Return>", function() execcomp())
        tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
        if (tclvalue(done) == "2") 
            return()
    }
    tkadd(ViewM, "command", label = "Plot map", command = function() plot.map.gui())
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
                filename <- tclvalue(tkgetSaveFile(initialfile = "ENiRGplot.tif", 
                  defaultextension = ".tif", title = "Save graphic as ...", 
                  filetypes = "{tiff {.tif .tiff}} {{All Files} {*.*}}"))
                if (filename != "") {
                  dev2bitmap(file = filename, type = "tiff24nc", width = width, height = height, units="cm")
                }
            }
            else if (outform == 2) {
                filename <- tclvalue(tkgetSaveFile(initialfile = "ENiRGplot.png", 
                  defaultextension = ".png", title = "Save graphic as ...", 
                  filetypes = "{PNG {.png}} {{All Files} {*.*}}"))
                if (filename != "") {
                  dev2bitmap(file = filename, type = "png16m", width = width, height = height, units="cm")
                }
            }
            else if (outform == 3) {
                filename <- tclvalue(tkgetSaveFile(initialfile = "ENiRGplot.jpeg", 
                  defaultextension = ".jpeg", title = "Save graphic as ...", 
                  filetypes = "{JPEG {.jpeg .jpg}} {{All Files} {*.*}}"))
                if (filename != "") {
                  dev2bitmap(file = filename, type = "jpeg", width = width, height = height, units="cm")
                }
            }
            tkdestroy(tf)
        }
        tkgrid(tklabel(tf, text = "Save current graphic as ..."), columnspan = 2)
        tkgrid(tklabel(frame2, text = "Output format: "), sticky = "n")
        tkgrid(tkradiobutton(frame2, text = "tiff", value = 1, 
            variable = formatvar), sticky = "w")
        tkgrid(tkradiobutton(frame2, text = "png", value = 2, 
            variable = formatvar), sticky = "w")
        tkgrid(tkradiobutton(frame2, text = "jpeg", value = 3, 
            variable = formatvar), sticky = "w")
        tkgrid(frame2, rowspan = 2, sticky = "n")
        tkgrid(tklabel(frame3, text = "Output size: "))
        width.entry <- tkentry(frame3, textvariable = widthvar, 
            width = 10)
        height.entry <- tkentry(frame3, textvariable = heightvar, 
            width = 10)
        tkgrid(tklabel(frame3, text = "Width: "), width.entry)
        tkgrid(tklabel(frame3, text = "Height: "), height.entry)
        tkgrid(frame3, column = 1, row = 1, sticky = "n")
        save.but <- tkbutton(frame1, text = "Save as ...", command = function() savefig(formatvar, 
            widthvar, heightvar))
        cancel.but <- tkbutton(frame1, text = "Dismiss", command = function() tkdestroy(tf))
        tkgrid(save.but, cancel.but)
        tkgrid(frame1, column = 1, row = 2, sticky = "n")
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkbind(tf, "<KeyPress-Return>", function() savefig(formatvar, 
            widthvar, heightvar))
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkwait.variable(done)
        if (tclvalue(done) == "2") 
            return(0)
        tkdestroy(tf)
    }
    tkadd(ViewM, "command", label = "Save graphic as ...", command = function() save.graphic.gui())
    tkadd(topMenu, "cascade", label = "Visualization", menu = ViewM)
    about.message <- function() {
        AboutW <- tktoplevel(width = 701, height = 232)
        tktitle(AboutW) <- "About"
        file.about <- paste(system.file(package='ENiRG'), '/img/ENiRG_1.gif', sep='')
        tcl("image","create","photo", "imageID", file=file.about)
        l <- ttklabel(AboutW, image="imageID", compound="image")
        tkpack(l)
        msg <- "Ecological Niche in R and GRASS (ENiRG)\nVersion: 1.0-1\nDate: 2016-05-03\n\nThis package was developed by:\n\nF. Canovas *\nC. Magliozzi *\nF. Mestre **\nJ.A. Palazon-Ferrando ***\nM. Gonzalez-Wanguemert *\n\n* Centro de Ciencias do Mar do Algarve (CCMAR, Portugal)\n** Centro de Investigacao em Biodiversidade e Recursos Geneticos,\nUniv. Evora (CIBIO/InBio-Evora, Portugal)\n*** Universidad de Murcia (Spain)\n\nMaintainer: F. Canovas <fcgarcia@ualg.pt>"
        tkpack(tklabel(AboutW, text = msg))
        Dismiss.but <- tkbutton(AboutW, text = "Dismiss", command = function() tkdestroy(AboutW))
        tkpack(Dismiss.but)
    }
    tkadd(AboutM, "command", label = "About the package", command = function() about.message())
    references.message <- function() {
        ReferencesW <- tktoplevel(width = 700, height = 235)
        tktitle(ReferencesW) <- "Main references"
        file.references <- paste(system.file(package='ENiRG'), '/img/ENiRG_2.gif', sep='')
        tcl("image","create","photo", "imageRef", file=file.references)
        ref <- ttklabel(ReferencesW, image="imageRef", compound="image")
        tkpack(ref)
        references <- "\n\nBasille, M., Calenge, C., Marboutin, E., Andersen, R. and Gaillard, J.M. (2008) Assessing habitat\nselection using multivariate statistics: some refinements of the ecological-niche factor\nanalysis. Ecological Modelling, 211, 233-240.\n\nBoyce, M.S.,Vernier, P.R.,Nielsen,S.E.,Schmiegelow, F.K.A. (2002) Evaluating resource selection\nfunctions Ecological Modelling 157, 281-300.\n\nCanovas, F., Magliozzi, C., Palazon-Ferrando, J.A., Mestre, F. and Gonzalez-Wanguemert, M. (2015)\nENiRG: R-GRASS interface to efficiently characterize the ecological niche of species and predict\nhabitat suitability. Ecography. DOI: 10.1111/ecog.01426.\n\nHirzel, A.H., Hausser, J., Chessel, D. and Perrin, N. (2002) Ecological-niche factor analysis: How\nto compute habitat-suitability maps without absence data? Ecology, 83, 2027-2036.\n\nHutchinson, G.E. (1957) Concluding Remarks. Cold Spring Harbor Symposium on Quantitative\nBiology, 22: 415-427. master's thesis.\n\n"
        tkpack(tklabel(ReferencesW, text = references))
        Dismiss.but <- tkbutton(ReferencesW, text = "Dismiss", command = function() tkdestroy(ReferencesW))
        tkpack(Dismiss.but)
    }
    tkadd(AboutM, "command", label = "References", command = function() references.message())
    tkadd(topMenu, "cascade", label = "About", menu = AboutM)
}
