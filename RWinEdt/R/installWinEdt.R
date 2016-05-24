delIni <- function(filename, ApplData){
    localfile <- file.path(ApplData, "ConfigEx", filename)
    if(file.exists(localfile))
        inifile <- readLines(localfile)
    else
        return(FALSE)
    if(filename == "Modes.ini")
        inifile <- gsub("*.rnw;*.snw;*.rd;", "", inifile, fixed=TRUE)
    Ws <- "//RWinEdt-start"
    We <- "//RWinEdt-end"
    pos1 <- grep(Ws, inifile)[1]
    pos2 <- grep(We, inifile)[1]
    if(length(pos1) && length(pos2))
        writeLines(inifile[c(1:(pos1-1), (pos2+1):length(inifile))], localfile)
}

updateIni <- function(filename, InstallRoot, ApplData, rwloc, force = FALSE){
    Ws <- "//RWinEdt-start"
    We <- "//RWinEdt-end"
    RWinEdtini <- c(Ws, readLines(file.path(rwloc, filename)), We, "[END]")
    filename <- gsub("[[:digit:]]", "", filename)
    localfile <- file.path(ApplData, "ConfigEx", filename)
    if(file.exists(localfile))
        inifile <- readLines(localfile)
    else {
        inifile <- readLines(file.path(InstallRoot, "ConfigEx", filename))
        dir.create(file.path(ApplData, "ConfigEx"), showWarnings=FALSE)
    }
    end <- "^\\[END\\]$"
    if(length(grep(Ws, inifile))){
        if(!force) return(FALSE)
        delIni(filename, ApplData)
        inifile <- readLines(localfile)
    }
    pos <- grep(end, inifile)[1]
    if(pos < 2) {
        inifile <- RWinEdtini
    } else {
        inifile <- c(inifile[1:(pos-1)], RWinEdtini, if(length(inifile) > pos) inifile[(pos+1):length(inifile)])
    }
    if(filename == "Modes.ini")
        inifile <- gsub("*.tex;*.ltx;*.texi;", "*.tex;*.ltx;*.texi;*.rnw;*.snw;*.rd;", inifile, fixed=TRUE)
    writeLines(inifile, localfile)
}


"installWinEdt" <-
function(InstallRoot, ApplData, WinEdtVersion, force = FALSE){
    if(is.null(force)){
        force <- "YES" == winDialog(type = "yesno",
            paste("User customized settings of R-WinEdt might be lost after resetting!\nReset?"))
    }
    rwloc <- normalizePath(file.path(system.file(package="RWinEdt"), "PlugIn"))
    rwloc6 <- normalizePath(file.path(system.file(package="RWinEdt"), "PlugIn6"))

    MR <- paste(InstallRoot, if(WinEdtVersion < 6) "\\macros\\R" else "\\Contrib\\R", sep = "")
    dirCreated <- file.exists(MR)
    if(!dirCreated) dirCreated <- dir.create(MR, recursive = TRUE)
    if((!(dirCreated || force)) ||  (file.access(MR, 2) == -1)){
        text <- paste("You need Administrator privileges the first time you run RWinEdt.",
                   "On Windows versions later than XP, please restart R by right click",
                   "and select 'Run as administrator'.", sep="\n")
        cat(text)
        stop(text)
    }
    if(WinEdtVersion < 6){
        for(i in dir(rwloc, "\\.edt$"))
            file.copy(file.path(rwloc, i, fsep = "\\"),
                file.path(InstallRoot, "macros", "R", i, fsep = "\\"), overwrite = force)
        for(i in dir(rwloc, "^R\\."))
            file.copy(file.path(rwloc, i, fsep = "\\"),
                file.path(InstallRoot, i, fsep = "\\"), overwrite = force)
        file.copy(file.path(rwloc, "send2R.edt", fsep = "\\"),
            file.path(InstallRoot, "send2R.edt", fsep = "\\"), overwrite = force)
        for(i in dir(rwloc, "\\.bmp$"))
            file.copy(file.path(rwloc, i, fsep = "\\"),
                file.path(InstallRoot, "bitmaps", "buttons", i, fsep = "\\"), overwrite = force)
        for(i in dir(rwloc, "^R_"))
            file.copy(file.path(rwloc, i, fsep = "\\"),
                file.path(InstallRoot, "bitmaps", "images", i, fsep = "\\"), overwrite = force)
        if(!missing(ApplData)){
            for(i in dir(rwloc, "^R\\."))
                file.copy(file.path(rwloc, i, fsep = "\\"),
                    file.path(ApplData, i, fsep = "\\"), overwrite = force)
            file.copy(file.path(rwloc, "send2R.edt", fsep = "\\"),
                file.path(ApplData, "send2R.edt", fsep = "\\"), overwrite = force)
        }
    } else {
        for(i in dir(rwloc))
            file.copy(file.path(rwloc, i, fsep = "\\"),
                file.path(InstallRoot, "Contrib", "R", i, fsep = "\\"), overwrite = force)
        for(i in dir(rwloc6))
            file.copy(file.path(rwloc6, i, fsep = "\\"),
                file.path(InstallRoot, "Contrib", "R", i, fsep = "\\"), overwrite = TRUE)
        MR <- file.path(ApplData, "Bitmaps", "Images")
        if(!file.exists(MR))
            dirCreated <- dir.create(MR, recursive = TRUE)
        for(i in dir(rwloc, "\\.bmp$"))
            file.copy(file.path(rwloc, i, fsep = "\\"),
                file.path(ApplData, "Bitmaps", "Images", i, fsep = "\\"), overwrite = force)
        for(i in dir(rwloc6, "\\.bmp$"))
            file.copy(file.path(rwloc6, i, fsep = "\\"),
                file.path(ApplData, "Bitmaps", "Images", i, fsep = "\\"), overwrite = TRUE)
        for(i in dir(rwloc6, "\\.png$"))
            file.copy(file.path(rwloc6, i, fsep = "\\"),
                file.path(ApplData, "Bitmaps", "Images", i, fsep = "\\"), overwrite = TRUE)
        file.copy(file.path(rwloc, "send2R.edt", fsep = "\\"),
                file.path(ApplData, "send2R.edt", fsep = "\\"), overwrite = force)
        file.copy(file.path(rwloc, "R.ver", fsep = "\\"),
                file.path(ApplData, "R.ver", fsep = "\\"), overwrite = force)
        inis <- c("MainMenu.ini", if(WinEdtVersion > 6) "Toolbar2.ini" else "Toolbar.ini", "FilterSets.ini",
                  "Keywords.ini", "Switches.ini", "Modes.ini")
        for(ini in inis)
            updateIni(ini, InstallRoot, ApplData, rwloc, force = force)
        startWinEdt(InstallRoot, ApplData, WinEdtVersion,
            args= paste(" [", paste("LoadConfig('.\\ConfigEx\\", inis, "');", sep="", collapse=""), "]", sep=""))
    }

    icon <- character(0)
    ## Statmenu icon?
    if(WinEdtVersion >= 6){
        if("YES" == winDialog(type = "yesno", "Create file type associations for .R, .Rd and .Rnw?"))
            startWinEdt(InstallRoot, ApplData, WinEdtVersion,
                args= " [Exe('%B\\Contrib\\R\\Filetypes.edt');]")
        winDialog("ok", "Starting from WinEdt 6, RWinEdt is integrated\nin the regular WinEdt, hence it is advisable\nnot to create any icons during the next questions,\nbut we offer this for completeness.")
    }
    icon <- if("YES" == winDialog(type = "yesno", "Create a StartMenu icon?"))
                "StartMenu"
    ## Desktop icon?
    icon <- if("YES" == winDialog(type = "yesno", "Create a Desktop icon?"))
                c(icon, "Desktop")
    for(i in icon){
        createIcon(InstallRoot, ApplData, i, WinEdtVersion)
    }
}
