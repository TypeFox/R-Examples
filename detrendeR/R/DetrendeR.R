detrender = function () DetrendeR()



DetrendeR = function () {

	detrendeRversion <- "detrendeR 1.0.4"
    try(tkdestroy(.detrendeRGUI), silent = TRUE)
	pos <- match("detrenderEnv", search())
	if (!is.na(pos)){
	if (exists (".detrendeRGUI", envir=as.environment(pos))) try(rm(.detrendeRGUI, pos=pos), silent = TRUE)
    }
    if (is.na(pos)) 


        .optionsDefault()

    listDataSets = function(envir = .GlobalEnv, ...) {
        Vars <- ls(envir = envir, all.names = TRUE)
        if (length(Vars) == 0) 
            return(Vars)
        out = names(which(sapply(Vars, function(.x) is.data.frame(get(.x, 
            envir = envir)) || is.matrix(get(.x, envir = envir)))))
        out
    }
    

    
    
   # .assign(".RemoveTrend", .RemoveTrend)
   # .assign(".as.logic", .as.logic)
   # .assign(".SPLINE", .SPLINE)
   # .assign(".NegExp", .NegExp)
   # .assign(".SPLINE", .SPLINE)
   # .assign(".GetDetrendMethod", .GetDetrendMethod)

    
    .assign("fname", "")
    .assign("TempDataBase", "")
	.assign(".detrendeRGUI",.detrendeRGUI <- tktoplevel())
#tkwm.resizable(.detrendeRGUI, 1, 0)
    tkwm.title(.detrendeRGUI, detrendeRversion)
    tkwm.geometry(.detrendeRGUI, "+0+0")
    frame0 <- tkframe(.detrendeRGUI, relief = "groove", borderwidth = 2)
    topMenuFile <- tkmenubutton(frame0, text = "File")
    fileMenu <- tkmenu(topMenuFile, tearoff = FALSE)
    tkconfigure(topMenuFile, menu = fileMenu)
    tkadd(fileMenu, "command", label = "Read compact", command = function() {
        try(eval(parse(text = "readCompact()")), silent = T)
    })
    tkadd(fileMenu, "command", label = "Read file", command = function() {
        readTable()
    })
    tkadd(fileMenu, "command", label = "Read rwl", command = function() {
        readRwlFile()
    })
    tkadd(fileMenu, "command", label = "Read crn", command = function() {
        readCrnFile()
    })
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Save compact", command = function() {
        if (DataBaseChoice != "<No active dataset>") {
            fname = tclvalue(tkgetSaveFile(initialfile = DataBaseChoice, 
                defaultextension = ".cpt", filetypes = " {{All Files} {*.*}} {CPT {.cpt}}"))
            if (fname != "") 
                try(eval(parse(text = paste("dplR:::write.rwl(", 
                  DataBaseChoice, ", fname=fname, format=\"compact\")"))), 
                  silent = T)
        }
    })
    tkadd(fileMenu, "command", label = "Save csv", command = function() {
        if (DataBaseChoice != "<No active dataset>") {
            fname = tclvalue(tkgetSaveFile(initialfile = DataBaseChoice, 
                defaultextension = ".csv", filetypes = " {{All Files} {*.*}} {CSV {.csv}}"))
            if (fname != "") {
                try(eval(parse(text = paste("YEAR<-row.names(", 
                  DataBaseChoice, ")"))), silent = T)
                try(eval(parse(text = paste("TempDataBase=data.frame(YEAR,", 
                  DataBaseChoice, ")"))), silent = T)
                try(write.table(TempDataBase, file = fname, quote = F, 
                  sep = ";", na = "", row.names = FALSE))
            }
        }
    })
    saveRwlMenu <- tkmenu(fileMenu, tearoff = FALSE)
    tkadd(saveRwlMenu, "command", label = "[0.01]", command = function() {
        if (DataBaseChoice != "<No active dataset>") {
            fname = tclvalue(tkgetSaveFile(initialfile = DataBaseChoice, 
                defaultextension = ".rwl", filetypes = " {{All Files} {*.*}} {RWL {.rwl}}"))
            if (fname != "") 
                try(eval(parse(text = paste("write.rwl(", DataBaseChoice, 
                  ", fname=fname, long.names=TRUE)"))), silent = T)
        }
    })
    tkadd(saveRwlMenu, "command", label = "[0.001]", command = function() {
        if (DataBaseChoice != "<No active dataset>") {
            fname = tclvalue(tkgetSaveFile(initialfile = DataBaseChoice, 
                defaultextension = ".rwl", filetypes = " {{All Files} {*.*}} {RWL {.rwl}}"))
            if (fname != "") 
                try(eval(parse(text = paste("write.rwl(", DataBaseChoice, 
                  ", fname=fname, prec=0.001,long.names=TRUE)"))), 
                  silent = T)
        }
    })
    tkadd(fileMenu, "cascade", label = "Save rwl", menu = saveRwlMenu)
    tkadd(fileMenu, "command", label = "Save crn", command = function() {
        if (DataBaseChoice != "<No active dataset>") {
            .assign("fname", tclvalue(tkgetSaveFile(initialfile = DataBaseChoice, 
                defaultextension = ".crn", filetypes = " {{All Files} {*.*}} {CRN {.crn}}")))
            if (fname != "") {
                try(eval(parse(text = paste("dplR:::write.crn(", 
                  DataBaseChoice, "[,1:2], fname = fname)"))), 
                  silent = T)
            }
        }
    })
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Change dir...", command = function() {
        try(setwd(toString(tkchooseDirectory())), silent = TRUE)
    })
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Save Workspace...", command = function() {
        try(save.image(tclvalue(tkgetSaveFile(title = "Save image...", 
            defaultextension = ".RData", filetypes = " {{All Files} {*.*}} {R_Image {.RData}}"))), 
            silent = T)
    })
    exitR = function() {
        response <- tclvalue(tkmessageBox(message = "Save workspace image?", 
            type = "yesnocancel", title = "Question", icon = "question"))
        if (response == "cancel") 
            return(invisible(response))
        if (response == "yes") {
            try(save.image(tclvalue(tkgetSaveFile(title = "Save image...", 
                defaultextension = ".RData", filetypes = " {{All Files} {*.*}} {R_Image {.RData}}"))), 
                silent = T)
            return(invisible(q(save = "no")))
        }
        if (response == "no") 
            return(invisible(q(save = "no")))
    }
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Quit DetrendeR", command = function() tkdestroy(.detrendeRGUI))
    tkadd(fileMenu, "command", label = "Quit R", command = function() exitR())
    call.detrender = function() {
        flag <- detrenderGUI()
        try(if (flag) 
            ARSTAN(), silent = TRUE)
    }
    topMenuTools <- tkmenubutton(frame0, text = "Tools")
    WinTools <- tkmenu(topMenuTools, tearoff = FALSE)
    tkconfigure(topMenuTools, menu = WinTools)
    tkadd(WinTools, "command", label = "Define settings", command = function() detrenderGUI())
    tkadd(WinTools, "command", label = "Batch mode", command = function() call.detrender())
    tkpack(topMenuFile, topMenuTools, side = "left")
    tkpack(frame0, fill = "x")
    frame1.0 <- tkframe(.detrendeRGUI, relief = "groove", borderwidth = 2)
    DataDetTcl <- tclVar("<Please select a dataset>")
    .assign("DataBaseChoice", c("<No active dataset>"))
    foo <- function() {
        DataBases <- listDataSets()
        if (length(DataBases) == 0) {
            .assign("DataBaseChoice", c("<No active dataset>"))
            tkconfigure(DataBaseSelectedBt, textvariable = tclVar(DataBaseChoice), 
                bg = "red", foreground = "black")
            tkfocus(.detrendeRGUI)
            return()
        }
        .assign("DataBaseChoice", tk_select.list(sort(listDataSets()), 
            title = "Select one"))
        if (DataBaseChoice != "") 
            tkconfigure(DataBaseSelectedBt, textvariable = tclVar(DataBaseChoice), 
                background = "grey80", foreground = "blue")
        if (DataBaseChoice == "") {
            .assign("DataBaseChoice", c("<No active dataset>"))
            tkconfigure(DataBaseSelectedBt, textvariable = tclVar(DataBaseChoice), 
                bg = "red", foreground = "black")
        }
    }
    DataBaseSelectedBt <- tkbutton(frame1.0, text = as.character(tclvalue(DataDetTcl)), 
        command = foo, width = 25, height = 1, background = "grey90", 
        foreground = "black")
    DELETE_DATASET = function() {
        if (DataBaseChoice != "<No active dataset>") {
            DatasetsToDelete <- tk_select.list(listDataSets(), 
                preselect = DataBaseChoice, multiple = TRUE, 
                title = "Select datasets to delete")
            if (length(DatasetsToDelete) > 0) {
                for (i in 1:length(DatasetsToDelete)) a <- try(eval(parse(text = paste("rm(", 
                  DatasetsToDelete[i], ",envir = .GlobalEnv)"))), 
                  silent = TRUE)
                if (any(DatasetsToDelete == DataBaseChoice)) {
                  .assign("DataBaseChoice", c("<No active dataset>"))
                  tkconfigure(DataBaseSelectedBt, textvariable = tclVar(DataBaseChoice), 
                    bg = "red", foreground = "black")
                }
            }
            tkfocus(.detrendeRGUI)
        }
    }
    DeleteDatasetBt <- tkbutton(frame1.0, text = "  Delete  ", 
        command = DELETE_DATASET)
    tkgrid(DataBaseSelectedBt, DeleteDatasetBt, sticky = "w")
    tkpack(frame1.0, fill = "x")
    frame2 <- tkframe(.detrendeRGUI, relief = "groove", borderwidth = 2)
    SERIES_INFORMATION = function() {
        if (DataBaseChoice != "<No active dataset>") {
            cat("\n[", DataBaseChoice, "]\n", sep = "")
            a <- try(eval(parse(text = sprintf("seriesInfo(%s)", 
                DataBaseChoice))), silent = TRUE)
        }
    }
    Information.but <- tkbutton(frame2, text = " Information ", 
        command = SERIES_INFORMATION)
    TREE_IDS = function() {
        if (DataBaseChoice != "<No active dataset>") {
            cat("\n[", DataBaseChoice, "]", sep = "")
            a <- try(eval(parse(text = paste("TreeIds(", DataBaseChoice, 
                ", stc=c(", stc[1], ",", stc[2], ",", stc[3], 
                "))"))), silent = TRUE)
        }
    }
    TreeIds.but <- tkbutton(frame2, text = "   TreeIds   ", command = TREE_IDS)
    MISSING_RINGS = function() {
        if (DataBaseChoice != "<No active dataset>") {
            cat("\n[", DataBaseChoice, "]", sep = "")
            a <- try(eval(parse(text = paste("TrwLessThan(", 
                DataBaseChoice, ",TRW=0)"))), silent = TRUE)
        }
    }
    MISSING_RINGS.but <- tkbutton(frame2, text = "Missing rings ", 
        command = MISSING_RINGS)
    RWL_INFO = function() {
        if (DataBaseChoice != "<No active dataset>") {
            cat("\n[", DataBaseChoice, "]\n", sep = "")
            a <- try(eval(parse(text = paste("RwlInfo(", DataBaseChoice, 
                ")"))), silent = TRUE)
        }
    }
    RWL_INFO.but <- tkbutton(frame2, text = "   RwlInfo ", command = RWL_INFO)
    SEG_PLOT = function() {
        if (DataBaseChoice != "<No active dataset>") {
            a <- try(eval(parse(text = paste("seg.plot(", DataBaseChoice, 
                ", main=\"", paste(DataBaseChoice), "\")"))), 
                silent = TRUE)
        }
    }
    SEG_PLOT.but <- tkbutton(frame2, text = " Segment plot ", 
        command = SEG_PLOT)
    RWL_PLOT = function() {
        if (DataBaseChoice != "<No active dataset>") {
            a <- try(eval(parse(text = paste("plotRwl(", DataBaseChoice, 
                ", main=\"", paste(DataBaseChoice), "\", save.csv=F)"))), 
                silent = TRUE)
        }
    }
    RWL_PLOT.but <- tkbutton(frame2, text = " Rwl plot ", command = RWL_PLOT)
    tkgrid(Information.but, TreeIds.but, MISSING_RINGS.but, RWL_INFO.but, 
        SEG_PLOT.but, RWL_PLOT.but)
    tkpack(frame2, fill = "x")
    frame3 <- tkframe(.detrendeRGUI, relief = "groove", borderwidth = 2)
    frame3.1 <- tkframe(frame3, relief = "groove", borderwidth = 1)
    frame3.2 <- tkframe(frame3, relief = "groove", borderwidth = 0)
    DETRENDING = function(TwoSteps = T, input = "", ...) {
        if (length(listDataSets()) == 0) 
            return()
        detrending(TwoSteps = TwoSteps, input = input)
    }
    topMenuDetrending <- tkmenubutton(frame3.1, text = "Detrending  ")
    detrendingMenu <- tkmenu(topMenuDetrending, tearoff = FALSE)
    tkconfigure(topMenuDetrending, menu = detrendingMenu)
    tkadd(detrendingMenu, "command", label = "1 step", command = function(...) {
        DETRENDING(TwoSteps = F, input = DataBaseChoice)
    })
    tkadd(detrendingMenu, "command", label = "2 steps", command = function(...) {
        DETRENDING(TwoSteps = T, input = DataBaseChoice)
    })
    AR.MODEL = function() {
        if (length(listDataSets()) == 0) 
            return()
        arMODEL(input = DataBaseChoice)
    }
    AR.but <- tkbutton(frame3.2, text = "  AR model ", command = AR.MODEL)
    makeCRONO = function() {
        if (length(listDataSets()) == 0) 
            return()
        CRONO(input = DataBaseChoice)
        tkfocus(.detrendeRGUI)
    }
    CRONObut <- tkbutton(frame3.2, text = "     Chrono     ", 
        command = makeCRONO)
    EPS = function(...) {
        if (length(listDataSets()) == 0) 
            return()
        interactiveEPS(input = DataBaseChoice)
    }
    EPSbut <- tkbutton(frame3.2, text = "     EPS     ", command = EPS)
    tkpack(topMenuDetrending)
    tkpack(frame3.1, side = "left")
    tkpack(AR.but, CRONObut, EPSbut, side = "left")
    tkpack(frame3.2, side = "left")
    tkpack(frame3, fill = "x")
    tkfocus(.detrendeRGUI)
    tkwm.resizable(.detrendeRGUI, 0, 0)
    #Sys.sleep(0.015)
    #.Geometry <- toString(tkwm.geometry(.detrendeRGUI))
    #.geometry <- strsplit(toString(tkwm.geometry(.detrendeRGUI)), "+")[[1]]
    #.n <- which(.geometry == "x")
    .assign(".width", 0)#as.integer(strtrim(.Geometry, .n - 1)))
    .assign(".heigth", 0)#as.integer(substr(.Geometry, .n + 1, which(.geometry == "+")[1] - 1)) + 26)
    return(invisible())
}

