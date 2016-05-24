
# The R Commander and command logger

# last modified 2016-03-20 by John Fox

# contributions by Milan Bouchet-Valat, Richard Heiberger, Duncan Murdoch, Erich Neuwirth, Brian Ripley

Commander <- function(){
    library(Rcmdr, quietly=TRUE)

    # set up RcmdrEnv
    RcmdrEnv.on.path <- getOption("Rcmdr")[["RcmdrEnv.on.path"]]
    if (is.null(RcmdrEnv.on.path)) RcmdrEnv.on.path <- FALSE
    if (RcmdrEnv.on.path){
        RcmdrEnv <- function() {
            pos <-  match("RcmdrEnv", search())
            if (is.na(pos)) { # Must create it
                RcmdrAttach <- base::attach
                RcmdrEnv <- list()
                RcmdrAttach(RcmdrEnv, pos = length(search()) - 1)
                rm(RcmdrEnv)
                pos <- match("RcmdrEnv", search())
            }
            return(pos.to.env(pos))
        }
        
        # the following two lines to be commented-out for debugging:
        assignInMyNamespace("RcmdrEnv", RcmdrEnv)        
        assignInMyNamespace(".RcmdrEnv", NULL)
        
    }
    RStudioP <- function() nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")) # to detect RStudio
    DESCRIPTION <- readLines(file.path(find.package("Rcmdr"), "DESCRIPTION")[1])
    RcmdrVersion <- trim.blanks(sub("^Version:", "",
        grep("^Version:", DESCRIPTION, value=TRUE)))
    putRcmdr("quotes", options(useFancyQuotes=FALSE))
    putRcmdr("messageNumber", 0)
    if (exists(".RcmdrEnv") && is.environment(RcmdrEnv()) &&
            exists("commanderWindow", RcmdrEnv()) &&
            !is.null(get("commanderWindow", RcmdrEnv()))) {
        return(invisible(NULL))
    }
    
    # check for auxiliary software
    putRcmdr("capabilities", RcmdrCapabilities())
    
    # the following function used to apply Rcmdr options with specified defaults
    #   if global == TRUE, store option
    setOption <- function(option, default, global=TRUE) {
        opt <- if (is.null(current[option][[1]])) default else current[option][[1]]
        if (global) putRcmdr(option, opt)
        opt
    }
    current <- getOption("Rcmdr")
    
    # define icons
    setOption("suppress.icon.images", FALSE)
    icon.images <- !getRcmdr("suppress.icon.images")
    tkimage.create("photo", "::image::RlogoIcon", file = system.file("etc", "R-logo.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::okIcon", 
        file = if (icon.images) system.file("etc", "ok.gif", package="Rcmdr") else system.file("etc", "blank.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::cancelIcon", file = if (icon.images) system.file("etc", "cancel.gif", package="Rcmdr") 
        else system.file("etc", "blank.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::helpIcon", file = if (icon.images) system.file("etc", "help.gif", package="Rcmdr")
        else system.file("etc", "blank.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::resetIcon", file = if (icon.images) system.file("etc", "reset.gif", package="Rcmdr")
        else system.file("etc", "blank.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::applyIcon", file = if (icon.images) system.file("etc", "apply.gif", package="Rcmdr")
        else system.file("etc", "blank.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::submitIcon", file = system.file("etc", "submit.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::editIcon", file = system.file("etc", "edit.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::viewIcon", file = system.file("etc", "view.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::dataIcon", file = system.file("etc", "data.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::modelIcon", file = system.file("etc", "model.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::removeIcon", file = system.file("etc", "remove.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::copyIcon", file = system.file("etc", "copy.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::cutIcon", file = system.file("etc", "cut.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::deleteIcon", file = system.file("etc", "delete.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::findIcon", file = system.file("etc", "find.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::pasteIcon", file = system.file("etc", "paste.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::redoIcon", file = system.file("etc", "redo.gif", package="Rcmdr"))
    tkimage.create("photo", "::image::undoIcon", file = system.file("etc", "undo.gif", package="Rcmdr"))
    
    # locate Rcmdr etc directory and directory for menus (usually the same)
    etc <- setOption("etc", system.file("etc", package="Rcmdr"))
    etcMenus <- setOption("etcMenus", etc)
    putRcmdr("etcMenus", etcMenus)
    
    # standard edit actions
    onCopy <- function(){
        focused <- tkfocus()
        if ((tclvalue(focused) != LogWindow()$ID) && (tclvalue(focused) != OutputWindow()$ID) && 
                (tclvalue(focused) != MessagesWindow()$ID) && (tclvalue(focused) != RmdWindow()$ID) && (tclvalue(focused) != RnwWindow()$ID))
            focused <- LogWindow()
        selection <- strsplit(tclvalue(tktag.ranges(focused, "sel")), " ")[[1]]
        if (is.na(selection[1])) return()
        text <- tclvalue(tkget(focused, selection[1], selection[2]))
        tkclipboard.clear()
        tkclipboard.append(text)
    }
    onDelete <- function(){
        focused <- tkfocus()
        if ((tclvalue(focused) != LogWindow()$ID) && (tclvalue(focused) != OutputWindow()$ID) && 
                (tclvalue(focused) != MessagesWindow()$ID) && (tclvalue(focused) != RmdWindow()$ID) && (tclvalue(focused) != RnwWindow()$ID))
            focused <- LogWindow()
        selection <- strsplit(tclvalue(tktag.ranges(focused, "sel")), " ")[[1]]
        if (is.na(selection[1])) return()
        tkdelete(focused, selection[1], selection[2])
    }
    onCut <- function(){
        onCopy()
        onDelete()
    }
    onPaste <- function(){
        onDelete()
        focused <- tkfocus()
        if ((tclvalue(focused) != LogWindow()$ID) && (tclvalue(focused) != OutputWindow()$ID)  && 
                (tclvalue(focused) != MessagesWindow()$ID) && (tclvalue(focused) != RmdWindow()$ID) && (tclvalue(focused) != RnwWindow()$ID))
            focused <- LogWindow()
        text <- tclvalue(.Tcl("selection get -selection CLIPBOARD"))
        if (length(text) == 0) return()
        tkinsert(focused, "insert", text)
    }
    onFind <- function(){
        focused <- tkfocus()
        if ((tclvalue(focused) != LogWindow()$ID) && (tclvalue(focused) != OutputWindow()$ID)  && 
                (tclvalue(focused) != MessagesWindow()$ID) && (tclvalue(focused) != RmdWindow()$ID) && (tclvalue(focused) != RnwWindow()$ID))
            focused <- LogWindow()
        initializeDialog(title=gettextRcmdr("Find"))
        textFrame <- tkframe(top)
        textVar <- tclVar(getRcmdr("last.search"))
        textEntry <- ttkentry(textFrame, width="20", textvariable=textVar)
        checkBoxes(frame="optionsFrame", boxes=c("regexpr", "case"), initialValues=c("0", "1"),
            labels=gettextRcmdr(c("Regular-expression search", "Case sensitive")))
        radioButtons(name="direction", buttons=c("foward", "backward"), labels=gettextRcmdr(c("Forward", "Backward")),
            values=c("-forward", "-backward"), title=gettextRcmdr("Search Direction"))
        onOK <- function(){
            text <- tclvalue(textVar)
            putRcmdr("last.search", text)
            if (text == ""){
                errorCondition(recall=onFind, message=gettextRcmdr("No search text specified."))
                return()
            }
            type <- if (tclvalue(regexprVariable) == 1) "-regexp" else "-exact"
            case <- tclvalue(caseVariable) == 1
            direction <- tclvalue(directionVariable)
            stop <- if (direction == "-forward") "end" else "1.0"
            where.txt <- if (case) tksearch(focused, type, direction, "--", text, "insert", stop)
            else tksearch(focused, type, direction, "-nocase", "--", text, "insert", stop)
            where.txt <- tclvalue(where.txt)
            if (where.txt == "") {
                Message(message=gettextRcmdr("Text not found."),
                    type="note")
                if (GrabFocus()) tkgrab.release(top)
                tkdestroy(top)
                tkfocus(CommanderWindow())
                return()
            }
            if (GrabFocus()) tkgrab.release(top)
            tkfocus(focused)
            tkmark.set(focused, "insert", where.txt)
            tksee(focused, where.txt)
            tkdestroy(top)
        }
        .exit <- function(){
            text <- tclvalue(textVar)
            putRcmdr("last.search", text)
            return("")
        }
        OKCancelHelp()
        tkgrid(labelRcmdr(textFrame, text=gettextRcmdr("Search for:")), textEntry, sticky="w")
        tkgrid(textFrame, sticky="w")
        tkgrid(optionsFrame, sticky="w")
        tkgrid(directionFrame, sticky="w")
        tkgrid(buttonsFrame, sticky="w")
        dialogSuffix(focus=textEntry)
    }
    onSelectAll <- function() {
        focused <- tkfocus()
        if ((tclvalue(focused) != LogWindow()$ID) && (tclvalue(focused) != OutputWindow()$ID) 
            && (tclvalue(focused) != MessagesWindow()$ID) && (tclvalue(focused) != RmdWindow()$ID) && (tclvalue(focused) != RnwWindow()$ID))
            focused <- LogWindow()
        tktag.add(focused, "sel", "1.0", "end")
        tkfocus(focused)
    }
    onClear <- function(){
        onSelectAll()
        onDelete()
    }
    onUndo <- function(){
        focused <- tkfocus()
        if ((tclvalue(focused) != LogWindow()$ID) && (tclvalue(focused) != OutputWindow()$ID) && 
                (tclvalue(focused) != MessagesWindow()$ID) && (tclvalue(focused) != RmdWindow()$ID) && (tclvalue(focused) != RnwWindow()$ID))
            focused <- LogWindow()
        tcl(focused, "edit", "undo")
    }
    onRedo <- function(){
        focused <- tkfocus()
        if ((tclvalue(focused) != LogWindow()$ID) && (tclvalue(focused) != OutputWindow()$ID) && 
                (tclvalue(focused) != MessagesWindow()$ID) && (tclvalue(focused) != RmdWindow()$ID) && (tclvalue(focused) != RnwWindow()$ID))
            focused <- LogWindow()
        tcl(focused, "edit", "redo")
    }
    
    # various initializations
    messageTag(reset=TRUE)
    putRcmdr("installed.packages", installed.packages())
    putRcmdr("RcmdrVersion", RcmdrVersion)
    putRcmdr(".activeDataSet", NULL)
    putRcmdr(".activeModel", NULL)
    putRcmdr("logFileName", NULL)
    putRcmdr("RmdFileName", "RcmdrMarkdown.Rmd")
    putRcmdr("RnwFileName", "RcmdrKnitr.Rnw")
    putRcmdr("outputFileName", NULL)
    putRcmdr("saveFileName", NULL)
    putRcmdr("modelNumber", 0)
    putRcmdr("reset.model", FALSE)
    putRcmdr("rgl", FALSE)
    putRcmdr("rgl.command", FALSE)
    putRcmdr("Identify3d", NULL)
    putRcmdr("open.dialog.here", NULL)
    putRcmdr("restoreTab", FALSE)
    putRcmdr("cancelDialogReopen", FALSE)
    putRcmdr("last.search", "")
    
    setOption("use.rgl", TRUE)

    # set up Rcmdr default and text (log) fonts, Tk scaling factor
    default.font.family.val <- tclvalue(.Tcl("font actual TkDefaultFont -family"))
    default.font.family.val <- gsub("\\{", "", gsub("\\}", "", default.font.family.val))
    default.font.family <- setOption("default.font.family", default.font.family.val)
    if (!("RcmdrDefaultFont" %in% as.character(.Tcl("font names")))){
        .Tcl(paste("font create RcmdrDefaultFont", tclvalue(tkfont.actual("TkDefaultFont"))))
        .Tcl("option add *font RcmdrDefaultFont")
    }
    
    .Tcl(paste("font configure RcmdrDefaultFont -family {", default.font.family, "}", sep=""))
    
    if (!("RcmdrTitleFont" %in% as.character(.Tcl("font names")))){
        .Tcl(paste("font create RcmdrTitleFont", tclvalue(tkfont.actual("TkDefaultFont"))))
    }
    .Tcl(paste("font configure RcmdrTitleFont -family {", default.font.family, "}", sep=""))
    if (!("RcmdrOutputMessagesFont" %in% as.character(.Tcl("font names")))){
        .Tcl(paste("font create RcmdrOutputMessagesFont", tclvalue(tkfont.actual("RcmdrTitleFont"))))
    }
    .Tcl(paste("font configure RcmdrTitleFont -family {", default.font.family, "}", sep=""))
    .Tcl(paste("font configure RcmdrOutputMessagesFont -family {", default.font.family, "}", sep=""))
    
    .Tcl(paste("font configure TkDefaultFont -family {",  default.font.family, "}", sep=""))
    log.font.family.val <- tclvalue(.Tcl("font actual TkFixedFont -family"))
    log.font.family.val <- gsub("\\{", "", gsub("\\}", "", log.font.family.val))
    log.font.family <- setOption("log.font.family", log.font.family.val)
    if (!("RcmdrLogFont" %in% as.character(.Tcl("font names")))){
        .Tcl(paste("font create RcmdrLogFont", tclvalue(tkfont.actual("TkFixedFont"))))
    }
    .Tcl(paste("font configure RcmdrLogFont -family {", log.font.family, "}", sep=""))
    .Tcl(paste("font configure TkFixedFont -family {",  log.font.family, "}", sep=""))
    putRcmdr("logFont", "RcmdrLogFont")    
    scale.factor <- current$scale.factor

    if (!is.null(scale.factor)) .Tcl(paste("tk scaling ", scale.factor, sep=""))
    # set various font sizes 
    if (WindowsP()){
      default.font.size.val <- abs(as.numeric(.Tcl("font actual TkDefaultFont -size")))
      if (is.na(default.font.size.val)) default.font.size.val <- 10
    }
    else default.font.size.val <- 10
    default.font.size <- setOption("default.font.size", default.font.size.val)
    tkfont.configure("RcmdrDefaultFont", size=default.font.size)
    tkfont.configure("RcmdrTitleFont", size=default.font.size)
    tkfont.configure("RcmdrOutputMessagesFont", size=default.font.size)
    tkfont.configure("TkDefaultFont", size=default.font.size)
    tkfont.configure("TkTextFont", size=default.font.size)
    tkfont.configure("TkCaptionFont", size=default.font.size)
    log.font.size <- setOption("log.font.size", 10)
    tkfont.configure("RcmdrLogFont", size=log.font.size)
    tkfont.configure("TkFixedFont", size=log.font.size)    
    
    .Tcl("ttk::style configure TButton -font RcmdrDefaultFont")
    .Tcl("ttk::style configure TLabel -font RcmdrDefaultFont")
    .Tcl("ttk::style configure TCheckbutton -font RcmdrDefaultFont")
    .Tcl("ttk::style configure TRadiobutton -font RcmdrDefaultFont")
    
    # set various options
    setOption("default.contrasts", c("contr.Treatment", "contr.poly"))
    standard.title.color <- as.character(.Tcl("ttk::style lookup TLabelframe.Label -foreground"))
    title.color <- setOption("title.color", standard.title.color) 
    if (tolower(title.color) == "black" || title.color == "#000000"){
        tkfont.configure("RcmdrTitleFont", weight="bold")
    }
    else tkfont.configure("RcmdrTitleFont", weight="normal")
    .Tcl("ttk::style configure TNotebook.Tab -font RcmdrDefaultFont")
    .Tcl(paste("ttk::style configure TNotebook.Tab -foreground", title.color))
    setOption("number.messages", TRUE)
    setOption("log.commands", TRUE)
    setOption("use.knitr", FALSE)
    setOption("use.markdown", !getRcmdr("use.knitr"))
    if ((!packageAvailable("markdown") && !packageAvailable("rmarkdown")) || (!packageAvailable("knitr"))) 
        putRcmdr("use.markdown", FALSE)
    if (!packageAvailable("knitr") || !getRcmdr("capabilities")$pdflatex) putRcmdr("use.knitr", FALSE)
    setOption("rmd.output.format", "html")
    putRcmdr("startNewCommandBlock", TRUE)
    putRcmdr("startNewKnitrCommandBlock", TRUE)
    putRcmdr("rmd.generated", FALSE)
    putRcmdr("rnw.generated", FALSE)
    setOption("RStudio", RStudioP())
    setOption("console.output", getRcmdr("RStudio"))
    setOption("retain.selections", TRUE)
    putRcmdr("dialog.values", list())
    putRcmdr("dialog.values.noreset", list())
    putRcmdr("savedTable", NULL)
    log.height <- as.character(setOption("log.height", if (!getRcmdr("log.commands")) 0 else 10))
    log.width <- as.character(setOption("log.width", 80))
    output.height <- as.character(setOption("output.height",
        if (getRcmdr("console.output")) 0
        else if ((as.numeric(log.height) != 0) || (!getRcmdr("log.commands"))) 2*as.numeric(log.height)
        else 20))
    messages.height <- as.character(setOption("messages.height", 3))
    putRcmdr("saveOptions", options(warn=1, contrasts=getRcmdr("default.contrasts"), width=as.numeric(log.width),
        na.action="na.exclude", graphics.record=TRUE))
    setOption("ask.to.exit", TRUE)
    setOption("ask.on.exit", TRUE)
    setOption("double.click", FALSE)
    setOption("sort.names", TRUE)
    setOption("grab.focus", TRUE)
    setOption("attach.data.set", FALSE)
    setOption("log.text.color", "black")
    setOption("command.text.color", "darkred")
    setOption("output.text.color", "darkblue")
    setOption("error.text.color", "red")
    setOption("warning.text.color", "darkgreen")
    setOption("prefixes", c("Rcmdr> ", "Rcmdr+ ", "RcmdrMsg: ", "RcmdrMsg+ "))
    setOption("multiple.select.mode", "extended")
    setOption("suppress.X11.warnings",
        interactive() && .Platform$GUI == "X11") # to address problem in X11 (Linux or Mac OS X)
    setOption("showData.threshold", 100)
    setOption("editDataset.threshold", 10000)
    setOption("retain.messages", TRUE)
    setOption("crisp.dialogs",  TRUE)
    setOption("length.output.stack", 10)
    setOption("length.command.stack", 10)
    setOption("quit.R.on.close", FALSE)
    putRcmdr("outputStack", as.list(rep(NA, getRcmdr("length.output.stack"))))
    putRcmdr("commandStack", as.list(rep(NA, getRcmdr("length.command.stack"))))
    setOption("variable.list.height", 6)
    setOption("variable.list.width", c(20, Inf))
    all.themes <- tk2theme.list()
    current.theme <- tk2theme()
    all.themes <- union(all.themes, current.theme)
    setOption("theme", current.theme)
    theme <- (getRcmdr("theme"))
    if (!(theme %in% all.themes)){
        warning(gettextRcmdr("non-existent theme"), ', "', theme,  '"\n  ', 
            gettextRcmdr("theme set to"), ' "', current.theme, '"')
        theme <- current.theme
    }
    putRcmdr("theme", theme)
    tk2theme(theme)
    placement <- setOption("placement", "", global=FALSE)
    
    putRcmdr("open.showData.windows", list())
    
    # platform-specific issues
    if (getRcmdr("suppress.X11.warnings")) {
        putRcmdr("messages.connection", file(open = "w+"))
        sink(getRcmdr("messages.connection"), type="message")
    }
    if (!(WindowsP())) {
        putRcmdr("oldPager", options(pager=RcmdrPager))
    }
    putRcmdr("restore.help_type", getOption("help_type"))
    setOption("help_type", "html")
    options(help_type=getRcmdr("help_type"))
#    putRcmdr("restore.use.external.help", FALSE)
    putRcmdr("restore.device", getOption("device"))
    if (RStudioP()){
        if (WindowsP()) options(device="windows")
        else if (MacOSXP()) options(device="quartz")
        else options(device="x11")
    }
    setOption("tkwait.dialog", FALSE)
    if (getRcmdr("tkwait.dialog")) putRcmdr("editDataset.threshold", 0)
    if (MacOSXP()){
        #       PATH <- system2("/usr/libexec/path_helper", "-s", stdout=TRUE)
        #       PATH <- sub("\"; export PATH;$", "", sub("^PATH=\\\"", "", PATH))
        #       Sys.setenv(PATH=PATH)
        PATH <- Sys.getenv("PATH")
        PATH <- unlist(strsplit(PATH, .Platform$path.sep, fixed=TRUE))
        if (MacOSXP("15.0.0")){
            if (length(grep("^/Library/TeX/texbin$", PATH)) == 0) {
                PATH[length(PATH) + 1] <- "/Library/TeX/texbin"
                Sys.setenv(PATH=paste(PATH, collapse=.Platform$path.sep))
            }
        }
        else{
            if (length(grep("^/usr/texbin$", PATH)) == 0) {
                PATH[length(PATH) + 1] <- "/usr/texbin"
                Sys.setenv(PATH=paste(PATH, collapse=.Platform$path.sep))
            }
        }
    }
    
    # source additional .R files, plug-ins preferred
    source.files <- list.files(etc, pattern="\\.[Rr]$")
    for (file in source.files) {
        source(file.path(etc, file))
        cat(paste(gettextRcmdr("Sourced:"), file, "\n"))
    }
    
    # collect plug-ins to be used
    Plugins <- options()$Rcmdr$plugins
    allPlugins <- listPlugins(loaded=TRUE)
    for (plugin in Plugins){
        if (!require(plugin, character.only=TRUE)){
            putRcmdr("commanderWindow", NULL)
            stop(sprintf(gettextRcmdr("the plug-in package %s is missing"), plugin))
        }
        if (!is.element(plugin, allPlugins)){
            putRcmdr("commanderWindow", NULL)
            stop(sprintf(gettextRcmdr("the package %s is not an Rcmdr plug-in"), plugin))
        }
    }
    
    # build Rcmdr menus
    Menus <- read.table(file.path(etcMenus, "Rcmdr-menus.txt"), colClasses = "character")
    addMenus <- function(Menus){
        removeMenus <- function(what){
            children <- Menus[Menus[,3] == what, 2]
            which <- what == Menus[,2] |  what == Menus[,5]
            Menus <<- Menus[!which,]
            for (child in children) removeMenus(child)
        }
        nms <- c("type", "menuOrItem", "operationOrParent", "label",
            "commandOrMenu", "activation", "install")
        names(Menus) <- nms
        for (plugin in Plugins) {
            MenusToAdd <- read.table(file.path(path.package(package=plugin)[1], "etc/menus.txt"),
                colClasses = "character")
            names(MenusToAdd) <- nms
            for (i in 1:nrow(MenusToAdd)){
                line <- MenusToAdd[i,]
                line[, "label"] <- gettext(line[,"label"], domain=paste("R=", plugin, sep=""))
                if (line[1, "type"] == "remove"){
                    removeMenus(line[1, "menuOrItem"])
                    next
                }
                if (line[1, "type"] == "menu"){
                    where <- if (line[1, "operationOrParent"] == "topMenu") 0
                    else max(which((Menus[, "type"] == "menu") &
                            (Menus[, "menuOrItem"] == line[1, "operationOrParent"])))
                }
                else if (line[1, "type"] == "item"){
                    if ((line[1, "operationOrParent"] == "command") || (line[1, "operationOrParent"] == "separator")){
                        which <- which(((Menus[, "operationOrParent"] == "command") | 
                                (Menus[, "operationOrParent"] == "separator")) &
                                (Menus[, "menuOrItem"] == line[1, "menuOrItem"]))
                        where <- if (length(which) == 0)
                            which((Menus[, "type"] == "menu")
                                & (Menus[, "menuOrItem"] == line[1, "menuOrItem"]))
                        else max(which)
                    }
                    else if (line[1, "operationOrParent"] == "cascade"){
                        where <- if (line[1, "menuOrItem"] != "topMenu")
                            max(which((Menus[, "operationOrParent"] == "cascade") &
                                    (Menus[, "menuOrItem"] == line[1, "menuOrItem"]) | (Menus[, "commandOrMenu"] == line[1, "menuOrItem"])))
                        else {
                            max(which((Menus[, "operationOrParent"] == "cascade") &
                                    (Menus[, "menuOrItem"] == "topMenu") &
                                    (Menus[, "commandOrMenu"] != "toolsMenu") &
                                    (Menus[, "commandOrMenu"] != "helpMenu")))
                        }
                    }
                    else stop(sprintf(gettextRcmdr('unrecognized operation, "%s", in plugin menu line %i'),
                        line[1, "operation"], i))
                }
                else stop(sprintf(gettextRcmdr('unrecognized type, "%s", in plugin menu line %i'),
                    line[1, "type"], i))
                Menus <- insertRows(Menus, line, where)
            }
        }
        Menus
    }
    Menus <- addMenus(Menus)
    menuNames <- Menus[Menus[,1] == "menu",]
    duplicateMenus <- duplicated(menuNames)
    if (any(duplicateMenus)) stop(paste(gettextRcmdr("Duplicate menu names:"),
        menuNames[duplicateMenus]))
    .Menus <- menus <- list()
    menuItems <- 0
    oldMenu <- ncol(Menus) == 6
    setOption("suppress.menus", FALSE)
    if (RExcelSupported()) # contributed by Erich Neuwirth
        putRExcel(".rexcel.menu.dataframe", Menus)
    modelClasses <- scan(file.path(etc, "model-classes.txt"), what="", quiet=TRUE, comment.char="#") # default recognized models
    
    # process plug-ins
    for (plugin in Plugins){
        description <- readLines(file.path(path.package(package=plugin)[1], "DESCRIPTION"))
        addModels <- description[grep("Models:", description)]
        addModels <- gsub(" ", "", sub("^Models:", "", addModels))
        addModels <- unlist(strsplit(addModels, ","))
        addRcmdrModels <- description[grep("RcmdrModels:", description)]
        addRcmdrModels <- gsub(" ", "", sub("^RcmdrModels:", "", addRcmdrModels))
        addRcmdrModels <- unlist(strsplit(addRcmdrModels, ","))
        if (length(addModels) > 0) modelClasses <- c(modelClasses, addModels)
        if (length(addRcmdrModels) > 0) modelClasses <- c(modelClasses, addRcmdrModels)
    }
    putRcmdr("modelClasses", modelClasses)
    
    # data-set edit
    onEdit <- function(){
        if (activeDataSet() == FALSE) {
            tkfocus(CommanderWindow())
            return()
        }
        dsnameValue <- ActiveDataSet()
        size <- eval(parse(text=paste("prod(dim(", dsnameValue, "))", sep=""))) #  prod(dim(save.dataset))
        if (size < 1 || size > getRcmdr("editDataset.threshold")){
            save.dataset <- get(dsnameValue, envir=.GlobalEnv)
            command <- paste("fix(", dsnameValue, ")", sep="")
            result <- justDoIt(command)
            if (class(result)[1] !=  "try-error"){ 			
                if (nrow(get(dsnameValue)) == 0){
                    errorCondition(window=NULL, message=gettextRcmdr("empty data set."))
                    justDoIt(paste(dsnameValue, "<- save.dataset"))
                    return()
                }
                else{
                    logger(command, rmd=FALSE)
                    activeDataSet(dsnameValue)
                }
            }
            else{
                errorCondition(window=NULL, message=gettextRcmdr("data set edit error."))
                return()
            }
        }
        else {
            command <- paste("editDataset(", dsnameValue, ")", sep="")
            result <- justDoIt(command)
            if (class(result)[1] !=  "try-error"){
                logger(command, rmd=FALSE)
            }
            else{
                errorCondition(window=NULL, message=gettextRcmdr("data set edit error."))
                return()
            }
        }
        tkwm.deiconify(CommanderWindow())
        tkfocus(CommanderWindow())
    }
    
    # data-set view
    onView <- function(){
#        if (packageAvailable("relimp")) Library("relimp", rmd=FALSE)
        if (activeDataSet() == FALSE) {
            tkfocus(CommanderWindow())
            return()
        }
        suppress <- if(getRcmdr("suppress.X11.warnings")) ", suppress.X11.warnings=FALSE" else ""
        view.height <- max(as.numeric(output.height) + as.numeric(log.height), 10)
        ncols <- ncol(get(ActiveDataSet()))
        command <- if (ncols <= getRcmdr("showData.threshold")){
            paste("showData(", ActiveDataSet(), ", placement='-20+200', font=getRcmdr('logFont'), maxwidth=",
                log.width, ", maxheight=", view.height, suppress, ")", sep="")
        }
        else paste("View(", ActiveDataSet(), ")", sep="")
        window <- justDoIt(command)
        if (!is.null(window)){
          open.showData.windows <- getRcmdr("open.showData.windows")
          open.window <- open.showData.windows[[ActiveDataSet()]]
          if (!is.null(open.window)) tkdestroy(open.window)
          open.showData.windows[[ActiveDataSet()]] <- window
          putRcmdr("open.showData.windows", open.showData.windows)
        }
    }
    
    # submit command in script tab or compile .Rmd file in markdown tab or compile .Rnw file in knitr tab
    onSubmit <- function(){
        .log <- LogWindow()
        .rmd <- RmdWindow()
        .rnw <- RnwWindow()
        if (as.character(tkselect(notebook)) == logFrame$ID) {
            selection <- strsplit(tclvalue(tktag.ranges(.log, "sel")), " ")[[1]]
            if (is.na(selection[1])) {
                tktag.add(.log, "currentLine", "insert linestart", "insert lineend")
                selection <- strsplit(tclvalue(tktag.ranges(.log,"currentLine")), " ")[[1]]
                tktag.delete(.log, "currentLine")
                if (is.na(selection[1])) {
                    Message(message=gettextRcmdr("Nothing is selected."),
                        type="error")
                    tkfocus(CommanderWindow())
                    return()
                }
            }
            lines <- tclvalue(tkget(.log, selection[1], selection[2]))
            lines <- strsplit(lines, "\n")[[1]]
            .console.output <- getRcmdr("console.output")
            .output <- OutputWindow()
            iline <- 1
            nlines <- length(lines)
            while (iline <= nlines){
                while (nchar(lines[iline])==0) iline <- iline + 1
                if (iline > nlines) break
                current.line <- lines[iline]
                if (.console.output) cat(paste("\n", getRcmdr("prefixes")[1], current.line,"\n", sep=""))
                else{
                    tkinsert(.output, "end", paste("\n> ", current.line,"\n", sep="")) 
                    tktag.add(.output, "currentLine", "end - 2 lines linestart", "end - 2 lines lineend")
                    tktag.configure(.output, "currentLine", foreground=getRcmdr("command.text.color"))
                }
                jline <- iline + 1
                while (jline <= nlines){
                    if (class(try(parse(text=current.line),silent=TRUE))!="try-error") break
                    if (.console.output)cat(paste(getRcmdr("prefixes")[2], lines[jline],"\n", sep=""))
                    else{
                        tkinsert(.output, "end", paste("+ ", lines[jline],"\n", sep=""))
                        tktag.add(.output, "currentLine", "end - 2 lines linestart", "end - 2 lines lineend")
                        tktag.configure(.output, "currentLine", foreground=getRcmdr("command.text.color"))
                    }
                    current.line <- paste(current.line, lines[jline],sep="\n")
                    jline <- jline + 1
                    iline <- iline + 1
                }
                if (!(is.null(current.line) || is.na(current.line))) doItAndPrint(current.line, log=FALSE, rmd=TRUE)
                iline <- iline + 1
                tkyview.moveto(.output, 1)
                tkfocus(.log)
            }
            if (length(as.character(tksearch(.log, "-regexp", "-forward",  "--", "\\n\\n$", "1.0"))
) == 0){
                tkinsert(.log, "end", "\n")
            }
            cursor.line.posn <- 1 + floor(as.numeric(tkindex(.log, "insert")))
            tkmark.set(.log, "insert", paste(cursor.line.posn, ".0", sep=""))
            tktag.remove(.log, "sel", "1.0", "end")
        }
        else if (as.character(tkselect(notebook)) == RmdFrame$ID) {
            compileRmd()
        }
        else{ 
            compileRnw()
        }
    }
    
    # right-click context menus
    contextMenuLog <- function(){
        .log <- LogWindow()
        contextMenu <- tkmenu(tkmenu(.log), tearoff=FALSE)
        tkadd(contextMenu, "command", label=gettextRcmdr("Submit"), command=onSubmit)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
        tkadd(contextMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
        tkadd(contextMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
        tkadd(contextMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Find..."), command=onFind)
        tkadd(contextMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
        tkadd(contextMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)
        tkpopup(contextMenu, tkwinfo("pointerx", .log), tkwinfo("pointery", .log))
    }
    contextMenuRmd <- function(){
        .rmd <- RmdWindow()
        contextMenu <- tkmenu(tkmenu(.rmd), tearoff=FALSE)
        tkadd(contextMenu, "command", label=gettextRcmdr("Generate report"), command=onSubmit)
        tkadd(contextMenu, "command", label=gettextRcmdr("Edit R Markdown document"), command=editMarkdown)
        tkadd(contextMenu, "command", label=gettextRcmdr("Remove last Markdown command block"), command=removeLastRmdBlock)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
        tkadd(contextMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
        tkadd(contextMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
        tkadd(contextMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
        tkadd(contextMenu, "separator")
#        tkadd(contextMenu, "command", label=gettextRcmdr("Find..."), command=onFind)  # doesn't work FIXME
        tkadd(contextMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
        tkadd(contextMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)
        tkpopup(contextMenu, tkwinfo("pointerx", .rmd), tkwinfo("pointery", .rmd))
    }
    contextMenuRnw <- function(){
        .rnw <- RnwWindow()
        contextMenu <- tkmenu(tkmenu(.rnw), tearoff=FALSE)
        tkadd(contextMenu, "command", label=gettextRcmdr("Generate PDF report"), command=onSubmit)
        tkadd(contextMenu, "command", label=gettextRcmdr("Edit knitr document"), command=editKnitr)
        tkadd(contextMenu, "command", label=gettextRcmdr("Remove last knitr command block"), command=removeLastRnwBlock)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
        tkadd(contextMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
        tkadd(contextMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
        tkadd(contextMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
        tkadd(contextMenu, "separator")
#        tkadd(contextMenu, "command", label=gettextRcmdr("Find..."), command=onFind) # doesn't work FIXME
        tkadd(contextMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
        tkadd(contextMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)
        tkpopup(contextMenu, tkwinfo("pointerx", .rnw), tkwinfo("pointery", .rnw))
    }
    contextMenuOutput <- function(){
        .output <- OutputWindow()
        contextMenu <- tkmenu(tkmenu(.output), tearoff=FALSE)
        tkadd(contextMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
        tkadd(contextMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
        tkadd(contextMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
        tkadd(contextMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Find..."), command=onFind)
        tkadd(contextMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
        tkadd(contextMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)
        tkpopup(contextMenu, tkwinfo("pointerx", .output), tkwinfo("pointery", .output))
    }
    contextMenuMessages <- function(){
        .messages <- MessagesWindow()
        contextMenu <- tkmenu(tkmenu(.messages), tearoff=FALSE)
        tkadd(contextMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
        tkadd(contextMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
        tkadd(contextMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
        tkadd(contextMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Find..."), command=onFind)
        tkadd(contextMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
        tkadd(contextMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)
        tkpopup(contextMenu, tkwinfo("pointerx", .messages), tkwinfo("pointery", .messages))
    }
    
    # main Commander window
    if (getRcmdr("crisp.dialogs")) tclServiceMode(on=FALSE)
    putRcmdr("commanderWindow", tktoplevel(class="Rcommander"))
    .commander <- CommanderWindow()
    tcl("wm", "iconphoto", .commander, "-default", "::image::RlogoIcon")
    tkwm.geometry(.commander, placement)
    tkwm.title(.commander, gettextRcmdr("R Commander"))
    tkwm.protocol(.commander, "WM_DELETE_WINDOW", 
        if (getRcmdr("quit.R.on.close")) closeCommanderAndR else CloseCommander)
    topMenu <- tkmenu(.commander)
    tkconfigure(.commander, menu=topMenu)
    position <- numeric(0)
    
    # install menus
    if (!getRcmdr("suppress.menus")){
        for (m in 1:nrow(Menus)){
            install <- if (oldMenu) "" else Menus[m, 7]
            if ((install != "") && (!eval(parse(text=install)))) next
            if (Menus[m, 1] == "menu") {
                position[Menus[m, 2]] <- 0
                assign(Menus[m, 2], tkmenu(get(Menus[m, 3]), tearoff=FALSE))
                menus[[Menus[m, 2]]] <- list(ID=get(Menus[m, 2])$ID, position=0)
            }
            else if (Menus[m, 1] == "item") {
                position[Menus[m, 2]] <- position[Menus[m, 2]] + 1
                if (Menus[m, 3] == "command"){
                    if (Menus[m, 6] == "")
                        tkadd(get(Menus[m, 2]), "command", label=gettextMenus(Menus[m, 4]),
                            command=get(Menus[m, 5]))
                    else {
                        tkadd(get(Menus[m, 2]), "command", label=gettextMenus(Menus[m, 4]),
                            command=get(Menus[m, 5]), state="disabled")
                        menuItems <- menuItems + 1
                        menus[[Menus[m, 2]]]$position <- position[Menus[m, 2]]
                        .Menus[[menuItems]] <- list(ID=menus[[Menus[m, 2]]]$ID, position=position[Menus[m, 2]],
                            activation=eval(parse(text=paste("function()", Menus[m, 6]))))
                    }
                }
                else if (Menus[m, 3] == "cascade")
                    tkadd(get(Menus[m, 2]), "cascade", label=gettextMenus(Menus[m, 4]),
                        menu=get(Menus[m, 5]))
                else if (Menus[m, 3] == "separator")
                    tkadd(get(Menus[m, 2]), "separator")
                else stop(paste(gettextRcmdr("menu definition error:"), Menus[m, ], collapse=" "),
                    domain=NA)
            }
            else stop(paste(gettextRcmdr("menu definition error:"), Menus[m, ], collapse=" "),
                domain=NA)
        }
    }
    putRcmdr("Menus", .Menus)
    putRcmdr("autoRestart", FALSE)
    activateMenus()
    
    # toolbar
    controlsFrame <- tkframe(CommanderWindow())
    editButton <- buttonRcmdr(controlsFrame, text=gettextRcmdr("Edit data set"), command=onEdit, 
        image="::image::editIcon", compound="left")
    viewButton <- buttonRcmdr(controlsFrame, text=gettextRcmdr("View data set"), command=onView,
        image="::image::viewIcon", compound="left")
    putRcmdr("dataSetName", tclVar(gettextRcmdr("<No active dataset>")))
    putRcmdr("dataSetLabel", tkbutton(controlsFrame, textvariable=getRcmdr("dataSetName"), foreground="red",
        relief="groove", command=selectActiveDataSet, image="::image::dataIcon", compound="left"))
    
    # script and markdown tabs
    notebook <- ttknotebook(CommanderWindow())
    logFrame <- ttkframe(CommanderWindow())    
    putRcmdr("logWindow", tktext(logFrame, bg="white", foreground=getRcmdr("log.text.color"),
        font=getRcmdr("logFont"), height=log.height, width=log.width, wrap="none", undo=TRUE))
    .log <- LogWindow()
    logXscroll <- ttkscrollbar(logFrame, orient="horizontal",
        command=function(...) tkxview(.log, ...))
    logYscroll <- ttkscrollbar(logFrame,
        command=function(...) tkyview(.log, ...))
    tkconfigure(.log, xscrollcommand=function(...) tkset(logXscroll, ...))
    tkconfigure(.log, yscrollcommand=function(...) tkset(logYscroll, ...))
    RmdFrame <- ttkframe(CommanderWindow())
    putRcmdr("RmdWindow", tktext(RmdFrame, bg="#FAFAFA", foreground=getRcmdr("log.text.color"),
        font=getRcmdr("logFont"), height=log.height, width=log.width, wrap="none", undo=TRUE))
    .rmd <- RmdWindow()
    rmd.template <- setOption("rmd.template", 
        system.file("etc", if (getRcmdr("capabilities")$pandoc) "Rcmdr-RMarkdown-Template.Rmd"
            else "Rcmdr-Markdown-Template.Rmd", package="Rcmdr"))
    template <- paste(readLines(rmd.template), collapse="\n")
    if (getRcmdr("use.rgl")) template <- paste0(template, 
      "\n\n```{r echo=FALSE}\n# include this code chunk as-is to enable 3D graphs\nlibrary(rgl)\nknitr::knit_hooks$set(webgl = hook_webgl)\n```\n\n")
    tkinsert(.rmd, "end", template)
    putRcmdr("markdown.output", FALSE)
    RmdXscroll <- ttkscrollbar(RmdFrame, orient="horizontal",
        command=function(...) tkxview(.rmd, ...))
    RmdYscroll <- ttkscrollbar(RmdFrame,
        command=function(...) tkyview(.rmd, ...))
    tkconfigure(.rmd, xscrollcommand=function(...) tkset(RmdXscroll, ...))
    tkconfigure(.rmd, yscrollcommand=function(...) tkset(RmdYscroll, ...))    
    
    RnwFrame <- ttkframe(CommanderWindow())
    putRcmdr("RnwWindow", tktext(RnwFrame, bg="#FAFAFA", foreground=getRcmdr("log.text.color"),
        font=getRcmdr("logFont"), height=log.height, width=log.width, wrap="none", undo=TRUE))
    .rnw <- RnwWindow()
    rnw.template <- setOption("rnw.template", 
        system.file("etc", "Rcmdr-knitr-Template.Rnw", package="Rcmdr"))
    template <- paste(readLines(rnw.template), collapse="\n")
    tkinsert(.rnw, "end", template)
    putRcmdr("knitr.output", FALSE)
    RnwXscroll <- ttkscrollbar(RnwFrame, orient="horizontal",
        command=function(...) tkxview(.rnw, ...))
    RnwYscroll <- ttkscrollbar(RnwFrame,
        command=function(...) tkyview(.rnw, ...))
    tkconfigure(.rnw, xscrollcommand=function(...) tkset(RnwXscroll, ...))
    tkconfigure(.rnw, yscrollcommand=function(...) tkset(RnwYscroll, ...))    
    
    outputFrame <- tkframe(.commander) 
    submitButtonLabel <- tclVar(gettextRcmdr("Submit"))
    submitButton <- if (getRcmdr("console.output"))
        buttonRcmdr(CommanderWindow(), textvariable=submitButtonLabel, borderwidth="2", command=onSubmit,
            image="::image::submitIcon", compound="left")
    else buttonRcmdr(outputFrame, textvariable=submitButtonLabel, borderwidth="2", command=onSubmit, 
        image="::image::submitIcon", compound="left")
    
    tkbind(CommanderWindow(), "<Button-1>", function() {
        if (as.character(tkselect(notebook)) == logFrame$ID) tclvalue(submitButtonLabel) <- gettextRcmdr("Submit")
        if (as.character(tkselect(notebook)) == RmdFrame$ID) tclvalue(submitButtonLabel) <- gettextRcmdr("Generate report")
        if (as.character(tkselect(notebook)) == RnwFrame$ID) tclvalue(submitButtonLabel) <- gettextRcmdr("Generate PDF report")
    })
    putRcmdr("outputWindow", tktext(outputFrame, bg="white", foreground=getRcmdr("output.text.color"),
        font=getRcmdr("logFont"), height=output.height, width=log.width, wrap="none", undo=TRUE))
    .output <- OutputWindow()
    outputXscroll <- ttkscrollbar(outputFrame, orient="horizontal",
        command=function(...) tkxview(.output, ...))
    outputYscroll <- ttkscrollbar(outputFrame,
        command=function(...) tkyview(.output, ...))
    tkconfigure(.output, xscrollcommand=function(...) tkset(outputXscroll, ...))
    tkconfigure(.output, yscrollcommand=function(...) tkset(outputYscroll, ...))
    # messages window
    messagesFrame <- tkframe(.commander)
    putRcmdr("messagesWindow", tktext(messagesFrame, bg="lightgray",
        font=getRcmdr("logFont"), height=messages.height, width=log.width, wrap="none", undo=TRUE))
    .messages <- MessagesWindow()
    messagesXscroll <- ttkscrollbar(messagesFrame, orient="horizontal",
        command=function(...) tkxview(.messages, ...))
    messagesYscroll <- ttkscrollbar(messagesFrame,
        command=function(...) tkyview(.messages, ...))
    tkconfigure(.messages, xscrollcommand=function(...) tkset(messagesXscroll, ...))
    tkconfigure(.messages, yscrollcommand=function(...) tkset(messagesYscroll, ...))
    
    # configure toolbar, etc., install various windows and widgets
    putRcmdr("modelName", tclVar(gettextRcmdr("<No active model>")))
    putRcmdr("modelLabel", tkbutton(controlsFrame, textvariable=getRcmdr("modelName"), foreground="red",
        relief="groove", command=selectActiveModel, image="::image::modelIcon", compound="left"))
    show.edit.button <- options("Rcmdr")[[1]]$show.edit.button
    show.edit.button <- if (is.null(show.edit.button)) TRUE else show.edit.button
    if (!getRcmdr("suppress.menus")){
        tkgrid(labelRcmdr(controlsFrame, image="::image::RlogoIcon", compound="left"),
            labelRcmdr(controlsFrame, text=gettextRcmdr("   Data set:")), getRcmdr("dataSetLabel"),
            if(show.edit.button) editButton, viewButton,
            labelRcmdr(controlsFrame, text=gettextRcmdr("Model:")), getRcmdr("modelLabel"), sticky="w", pady=c(3, 3))
        tkgrid(controlsFrame, sticky="w")
        tkgrid.configure(getRcmdr("dataSetLabel"), padx=c(2, 5))
        tkgrid.configure(getRcmdr("modelLabel"), padx=c(2, 10))
        tkgrid.configure(editButton, padx=c(10, 1))
        if (show.edit.button) tkgrid.configure(viewButton, padx=c(1, 15))
        else tkgrid.configure(viewButton, padx=c(10, 15))
    }
    .log.commands <-  getRcmdr("log.commands")
    .console.output <- getRcmdr("console.output")
    if (.log.commands) {
        tkgrid(.log, logYscroll, sticky="news", columnspan=2)
        tkgrid(logXscroll)
        tkgrid(logFrame, sticky="news", padx=10, pady=0, columnspan=2)
        tkgrid(.rmd, RmdYscroll, sticky="news", columnspan=2)
        tkgrid(RmdXscroll)
        tkgrid(.rnw, RnwYscroll, sticky="news", columnspan=2)
        tkgrid(RnwXscroll)
        if (getRcmdr("use.markdown")) tkgrid(RmdFrame, sticky="news", padx=10, pady=0, columnspan=2)
        if (getRcmdr("use.knitr")) tkgrid(RnwFrame, sticky="news", padx=10, pady=0, columnspan=2)
    }
    tkadd(notebook, logFrame, text=gettextRcmdr("R Script"), padding=6)
    if (getRcmdr("use.markdown")) tkadd(notebook, RmdFrame, text=gettextRcmdr("R Markdown"), padding=6)
    if (getRcmdr("use.knitr")) tkadd(notebook, RnwFrame, text=gettextRcmdr("knitr Document"), padding=6)
    tkgrid(notebook, sticky="news")
    if (.log.commands && .console.output) tkgrid(submitButton, sticky="w", pady=c(0, 6))
    tkgrid(labelRcmdr(outputFrame, text=gettextRcmdr("Output"), font="RcmdrOutputMessagesFont", foreground=title.color),
        if (.log.commands && !.console.output) submitButton, sticky="sw", pady=c(6, 6))
    tkgrid(.output, outputYscroll, sticky="news", columnspan=2)
    tkgrid(outputXscroll, columnspan=1 + (.log.commands && !.console.output))
    if (!.console.output) tkgrid(outputFrame, sticky="news", padx=10, pady=0, columnspan=2)
    tkgrid(labelRcmdr(messagesFrame, text=gettextRcmdr("Messages"), font="RcmdrOutputMessagesFont", foreground=title.color), 
           sticky="w", pady=c(6, 6))
    tkgrid(.messages, messagesYscroll, sticky="news", columnspan=2)
    tkgrid(messagesXscroll)
    if (!.console.output) tkgrid(messagesFrame, sticky="news", padx=10, pady=0, columnspan=2) ##rmh & J. Fox
    tkgrid.configure(logYscroll, sticky="ns")
    tkgrid.configure(logXscroll, sticky="ew")
    tkgrid.configure(RmdYscroll, sticky="ns")
    tkgrid.configure(RmdXscroll, sticky="ew")
    tkgrid.configure(RnwYscroll, sticky="ns")
    tkgrid.configure(RnwXscroll, sticky="ew")
    tkgrid.configure(outputYscroll, sticky="ns")
    tkgrid.configure(outputXscroll, sticky="ew")
    tkgrid.configure(messagesYscroll, sticky="ns")
    tkgrid.configure(messagesXscroll, sticky="ew")
    .commander <- CommanderWindow()
    tkgrid.rowconfigure(.commander, 0, weight=0)
    tkgrid.rowconfigure(.commander, 1, weight=1)
    tkgrid.rowconfigure(.commander, 2, weight=1)
    tkgrid.columnconfigure(.commander, 0, weight=1)
    tkgrid.columnconfigure(.commander, 1, weight=0)
    if (.log.commands){
        tkgrid.rowconfigure(logFrame, 0, weight=1)
        tkgrid.rowconfigure(logFrame, 1, weight=0)
        tkgrid.columnconfigure(logFrame, 0, weight=1)
        tkgrid.columnconfigure(logFrame, 1, weight=0)
        if (getRcmdr("use.markdown")){
            tkgrid.rowconfigure(RmdFrame, 0, weight=1)
            tkgrid.rowconfigure(RmdFrame, 1, weight=0)
            tkgrid.columnconfigure(RmdFrame, 0, weight=1)
            tkgrid.columnconfigure(RmdFrame, 1, weight=0)
        }
        if (getRcmdr("use.knitr")){
            tkgrid.rowconfigure(RnwFrame, 0, weight=1)
            tkgrid.rowconfigure(RnwFrame, 1, weight=0)
            tkgrid.columnconfigure(RnwFrame, 0, weight=1)
            tkgrid.columnconfigure(RnwFrame, 1, weight=0)
        }
    }
    if (!.console.output){
        tkgrid.rowconfigure(outputFrame, 0, weight=0)
        tkgrid.rowconfigure(outputFrame, 1, weight=1)
        tkgrid.rowconfigure(outputFrame, 2, weight=0)
        tkgrid.columnconfigure(outputFrame, 0, weight=1)
        tkgrid.columnconfigure(outputFrame, 1, weight=0)
    }
    tkgrid.rowconfigure(messagesFrame, 0, weight=0)
    tkgrid.rowconfigure(messagesFrame, 1, weight=0)
    tkgrid.rowconfigure(messagesFrame, 2, weight=0)
    tkgrid.columnconfigure(messagesFrame, 0, weight=1)
    tkgrid.columnconfigure(messagesFrame, 1, weight=0)
    .Tcl("update idletasks")
    tkbind(.commander, "<Control-x>", onCut)
    tkbind(.commander, "<Control-X>", onCut)
    tkbind(.commander, "<Control-c>", onCopy)
    tkbind(.commander, "<Control-C>", onCopy)
    tkbind(.commander, "<Control-r>", onSubmit)
    tkbind(.commander, "<Control-R>", onSubmit)
    tkbind(.commander, "<Control-Tab>", onSubmit)
    tkbind(.commander, "<Control-f>", onFind)
    tkbind(.commander, "<Control-F>", onFind)
    tkbind(.commander, "<F3>", onFind)
    tkbind(.commander, "<Control-s>", saveLog)
    tkbind(.commander, "<Control-S>", saveLog)
    tkbind(.commander, "<Control-a>", onSelectAll)
    tkbind(.commander, "<Control-A>", onSelectAll)
    tkbind(.commander, "<Control-w>", onRedo)
    tkbind(.commander, "<Control-W>", onRedo)
    tkbind(.commander, "<Alt-BackSpace>", onUndo)
    tkbind(.log, "<ButtonPress-3>", contextMenuLog)
    tkbind(.rmd, "<ButtonPress-3>", contextMenuRmd)
    tkbind(.rnw, "<ButtonPress-3>", contextMenuRnw)
    tkbind(.output, "<ButtonPress-3>", contextMenuOutput)
    tkbind(.messages, "<ButtonPress-3>", contextMenuMessages)
    tkbind(.log, "<Control-ButtonPress-1>", contextMenuLog)
    tkbind(.rmd, "<Control-ButtonPress-1>", contextMenuRmd)
    tkbind(.rnw, "<Control-ButtonPress-1>", contextMenuRnw)
    tkbind(.output, "<Control-ButtonPress-1>", contextMenuOutput)
    tkbind(.messages, "<Control-ButtonPress-1>", contextMenuMessages)
    tkbind(.rmd, "<Control-e>", editMarkdown)
    tkbind(.rmd, "<Control-E>", editMarkdown)
    tkbind(.rnw, "<Control-e>", editKnitr)
    tkbind(.rnw, "<Control-E>", editKnitr)
    if (MacOSXP()){
        tkbind(.commander, "<Meta-x>", onCut)
        tkbind(.commander, "<Meta-X>", onCut)
        tkbind(.commander, "<Meta-c>", onCopy)
        tkbind(.commander, "<Meta-C>", onCopy)
        tkbind(.commander, "<Meta-v>", onPaste)
        tkbind(.commander, "<Meta-V>", onPaste)
        tkbind(.commander, "<Meta-r>", onSubmit)
        tkbind(.commander, "<Meta-R>", onSubmit)
        tkbind(.commander, "<Meta-Tab>", onSubmit)
        tkbind(.commander, "<Meta-f>", onFind)
        tkbind(.commander, "<Meta-F>", onFind)
        tkbind(.commander, "<Meta-s>", saveLog)
        tkbind(.commander, "<Meta-S>", saveLog)
        tkbind(.commander, "<Meta-a>", onSelectAll)
        tkbind(.commander, "<Meta-A>", onSelectAll)
        tkbind(.commander, "<Meta-w>", onRedo)
        tkbind(.commander, "<Meta-W>", onRedo)
        tkbind(.log, "<Meta-ButtonPress-1>", contextMenuLog)
        tkbind(.rmd, "<Meta-ButtonPress-1>", contextMenuRmd)
        tkbind(.rnw, "<Meta-ButtonPress-1>", contextMenuRnw)
        tkbind(.output, "<Meta-ButtonPress-1>", contextMenuOutput)
        tkbind(.messages, "<Meta-ButtonPress-1>", contextMenuMessages)
        tkbind(.rmd, "<Meta-e>", editMarkdown)
        tkbind(.rmd, "<Meta-E>", editMarkdown)
        tkbind(.rnw, "<Meta-e>", editKnitr)
        tkbind(.rnw, "<Meta-E>", editKnitr)
    }
    tkwm.deiconify(.commander)
    tkfocus(.commander)
    if (getRcmdr("crisp.dialogs")) tclServiceMode(on=TRUE)
    tkwait.commander <- options("Rcmdr")[[1]]$tkwait.commander  # to address problem in Debian Linux
    if ((!is.null(tkwait.commander)) && tkwait.commander) {
        putRcmdr(".commander.done", tclVar("0"))
        tkwait.variable(getRcmdr(".commander.done"))
    }
    Message(paste(gettextRcmdr("R Commander Version "), " ", getRcmdr("RcmdrVersion"), ": ", date(), sep=""))
    if (.Platform$GUI == "Rgui"  && ismdi()) Message(gettextRcmdr(
        "The Windows version of the R Commander works best under\nRGui with the single-document interface (SDI); see ?Commander."),
        type="warning")
    if (RappP()  && mavericksP() && appnap() == "on") Message(gettextRcmdr(
      "The Mac OS X version of the R Commander works best under R.app\nwith app nap turned off. See ?Commander and the Tools menu."),
      type="warning")
}

# put commands in script, markdown, and knitr tabs
logger <- function(command, rmd=TRUE){
    pushCommand(command)
    .log <- LogWindow()
    .rmd <- RmdWindow()
    .rnw <- RnwWindow()
    .output <- OutputWindow()
    Rmd <- rmd && is.null(attr(command, "suppressRmd")) && (getRcmdr("use.markdown") || getRcmdr("use.knitr"))
    command <- splitCmd(command)
    if (getRcmdr("log.commands")) {
        last2 <- tclvalue(tkget(.log, "end -2 chars", "end"))
        if (last2 != "\n\n") tkinsert(.log, "end", "\n")
        tkinsert(.log, "end", paste(command,"\n", sep=""))
        tkyview.moveto(.log, 1)
        if (Rmd){
            if (getRcmdr("use.markdown")){
                if (getRcmdr("startNewCommandBlock")){
                    beginRmdBlock()
                    tkinsert(.rmd, "end", paste(command, "\n", sep=""))
                    tkyview.moveto(.rmd, 1)
                    putRcmdr("markdown.output", TRUE)
                    endRmdBlock()
                }
                else{
                    tkinsert(.rmd, "end", paste(command, "\n", sep=""))
                    tkyview.moveto(.rmd, 1)
                    putRcmdr("markdown.output", TRUE)
                    putRcmdr("rmd.generated", TRUE)
                }
            }
            if (getRcmdr("use.knitr")){
                if (getRcmdr("startNewKnitrCommandBlock")){
                    beginRnwBlock()
                    tkinsert(.rnw, "end", paste(command, "\n", sep=""))
                    tkyview.moveto(.rnw, 1)
                    putRcmdr("knitr.output", TRUE)
                    endRnwBlock()
                }
                else{
                    tkinsert(.rnw, "end", paste(command, "\n", sep=""))
                    tkyview.moveto(.rnw, 1)
                    putRcmdr("knitr.output", TRUE)
                    putRcmdr("rnw.generated", TRUE)
                }
            }
        }
        
    }
    lines <- strsplit(command, "\n")[[1]]
    tkinsert(.output, "end", "\n")
    if (getRcmdr("console.output")) {
        for (line in seq(along.with=lines)) {
            prompt <- ifelse (line==1, paste("\n", getRcmdr("prefixes")[1], sep=""), paste("\n", getRcmdr("prefixes")[2], sep=""))
            cat(paste(prompt, lines[line]))
        }
        cat("\n")
    }
    else {
        for (line in  seq(along.with=lines)) {
            prompt <- ifelse(line==1, "> ", "+ ")
            tkinsert(.output, "end", paste(prompt, lines[line], "\n", sep=""))
            tktag.add(.output, "currentLine", "end - 2 lines linestart", "end - 2 lines lineend")
            tktag.configure(.output, "currentLine", foreground=getRcmdr("command.text.color"))
            tkyview.moveto(.output, 1)
        }
    }
    command
}

justDoIt <- function(command) {
    command <- enc2native(command)
    Message()
    if (!getRcmdr("suppress.X11.warnings")){
        messages.connection <- file(open="w+")
        sink(messages.connection, type="message")
        on.exit({
            sink(type="message")
            close(messages.connection)
        })
    }
    else messages.connection <- getRcmdr("messages.connection")
    capture.output(result <- try(eval(parse(text=command), envir=.GlobalEnv), silent=TRUE))
    if (class(result)[1] ==  "try-error"){
        Message(message=paste(strsplit(result, ":")[[1]][2]), type="error")
        tkfocus(CommanderWindow())
        return(result)
    }
    checkWarnings(readLines(messages.connection))
    if (getRcmdr("RStudio")) Sys.sleep(0)
    result
}

# execute commands, save commands and output
doItAndPrint <- function(command, log=TRUE, rmd=log) {
    command <- enc2native(command)
    Message()
    .console.output <- getRcmdr("console.output")
    .output <- OutputWindow()
    if (!.console.output) {
        width <- (as.numeric(tkwinfo("width", .output)) - 2*as.numeric(tkcget(.output, borderwidth=NULL)) - 2)/
            as.numeric(tkfont.measure(tkcget(.output, font=NULL), "0"))
        eval(parse(text=paste("options(width=", floor(width), ")", sep="")))
    }
    if (!getRcmdr("suppress.X11.warnings")){
        messages.connection <- file(open="w+")
        sink(messages.connection, type="message")
        on.exit({
            sink(type="message")
            close(messages.connection)
        })
    }
    else messages.connection <- getRcmdr("messages.connection")
    output.connection <- file(open="w+")
    sink(output.connection, type="output")
    on.exit({
        if (!.console.output) sink(type="output") # if .console.output, output connection already closed
        close(output.connection)
    }, add=TRUE)
    if (log) logger(command, rmd=rmd) 
    else {
        pushCommand(command)
        if (rmd) {
            if (getRcmdr("use.markdown")) enterMarkdown(command)
            if (getRcmdr("use.knitr")) enterKnitr(command)
        }
    }
    result <- try(parse(text=paste(command)), silent=TRUE)
    if (class(result)[1] == "try-error"){
        if (rmd) {
            if (getRcmdr("use.markdown")) {
                removeLastRmdBlock()
                putRcmdr("startNewCommandBlock", TRUE)
            }
            if (getRcmdr("use.knitr")) {
                removeLastRnwBlock()
                putRcmdr("startNewKnitrCommandBlock", TRUE)
            }
        }
        Message(message=paste(strsplit(result, ":")[[1]][2]), type="error")
        if (.console.output) sink(type="output")
        tkfocus(CommanderWindow())
        return(result)
    } else {
        exprs <- result
        result <- NULL
    }
    for (i in seq_along(exprs)) {
        ei <- exprs[i]
        tcl("update")
        result <-  try(withVisible(eval(ei, envir=.GlobalEnv)), silent=TRUE)
        if (class(result)[1] ==  "try-error"){
            if (rmd) {
                if (getRcmdr("use.markdown")) {
                    removeLastRmdBlock()
                    putRcmdr("startNewCommandBlock", TRUE)
                }
                if (getRcmdr("use.knitr")) {
                    removeLastRnwBlock()
                    putRcmdr("startNewKnitrCommandBlock", TRUE)
                }
            }
            Message(message=paste(strsplit(result, ":")[[1]][2]), type="error")
            if (.console.output) sink(type="output")
            tkfocus(CommanderWindow())
            return(result)
        }
        result <- if (result$visible == FALSE) NULL else result$value
        if (!is.null(result)) pushOutput(result)
        if (isS4object(result)) show(result) else print(result)
        .Output <- readLines(output.connection)
        if (length(.Output) > 0 && .Output[length(.Output)] == "NULL")
            .Output <- .Output[-length(.Output)] # suppress "NULL" line at end of output
        if (length(.Output) != 0) {  # is there output to print?
            if (.console.output) {
                out <- .Output
                sink(type="output")
                for (line in out) cat(paste(line, "\n", sep=""))
            }
            else{
                for (line in .Output) tkinsert(.output, "end", paste(line, "\n", sep=""))
                tkyview.moveto(.output, 1)
            }
        }
        else if (.console.output) sink(type="output")
        if (RExcelSupported()) # added by Erich Neuwirth
            putRExcel(".rexcel.last.output",.Output)
        # errors already intercepted, display any warnings
        checkWarnings(readLines(messages.connection))
    }
    if (getRcmdr("RStudio")) Sys.sleep(0)
    result
}

checkWarnings <- function(messages){
    if (getRcmdr("suppress.X11.warnings")){
        X11.warning <- grep("X11 protocol error|Warning in structure", messages)
        if (length(X11.warning) > 0){
            messages <- messages[-X11.warning]
        }
        if (length(messages) == 0) Message()
        else if (length(messages) > 10) {
            messages <- c(paste(length(messages), "warnings."),
                gettextRcmdr("First and last 5 warnings:"),
                head(messages,5), ". . .", tail(messages, 5))
            Message(message=paste(messages, collapse="\n"), type="warning")
        }
        else {
            if (length(grep("warning", messages, ignore.case=TRUE)) > 0)
                Message(message=paste(messages, collapse="\n"), type="warning")
            else Message(message=paste(messages, collapse="\n"), type="note")
        }
    }
    else{
        if (length(messages) == 0) Message()
        else if (length(messages) > 10){
            messages <- c(paste(length(messages), "warnings."),
                gettextRcmdr("First and last 5 warnings:"),
                head(messages, 5), ". . .", tail(messages, 5))
            Message(message=paste(messages, collapse="\n"), type="warning")
        }
        else {
            if (length(grep("warning", messages, ignore.case=TRUE)) > 0)
                Message(message=paste(messages, collapse="\n"), type="warning")
            else Message(message=paste(messages, collapse="\n"), type="note")
        }
    }
    tkfocus(CommanderWindow())
}

pause <- function(seconds = 1){
    if (seconds <= 0) stop("seconds must be positive")
    start <- proc.time()[3]
    while (as.numeric(elapsed <- (proc.time()[3] - start)) < seconds) {}
    elapsed
}

Message <- function(message, type=c("note", "error", "warning")){
    tcl("update") 
    .message <- MessagesWindow()
    type <- match.arg(type)
    if (type != "note") tkbell()
    if (getRcmdr("retain.messages")) {
        if (missing(message) && !is.null(getRcmdr("last.message"))) {
            putRcmdr("last.message", NULL)
            tkyview.moveto(.message, 1.0)
        }
    }
    else if (type == "note"){
        lastMessage <- tclvalue(tkget(MessagesWindow(),  "end - 2 lines", "end"))
        if (length(c(grep(gettextRcmdr("ERROR:"), lastMessage), grep(gettextRcmdr("WARNING:"), lastMessage))) == 0)
            tkdelete(.message, "1.0", "end")
    }
    else tkdelete(.message, "1.0", "end")
    col <- if (type == "error") getRcmdr("error.text.color")
    else if (type == "warning") getRcmdr("warning.text.color")
    else getRcmdr("output.text.color")
    prefix <- switch(type, error=gettextRcmdr("ERROR"), warning=gettextRcmdr("WARNING"), note=gettextRcmdr("NOTE"))
    if (missing(message)){
        return()
    }
    putRcmdr("last.message", type)
    message <- paste(prefix, ": ", message, sep="")
    if (getRcmdr("retain.messages") && getRcmdr("number.messages")) {
        messageNumber <- getRcmdr("messageNumber") + 1
        putRcmdr("messageNumber", messageNumber)
        message <- paste("[", messageNumber, "] ", message, sep="")
    }
    if (RExcelSupported()) # added by Erich Neuwirth
        putRExcel(".rexcel.last.message",message)
    lines <- strsplit(message, "\n")[[1]]
    console.output <- getRcmdr("console.output")
    if (!console.output){
        width <- (as.numeric(tkwinfo("width", .message)) - 2*as.numeric(tkcget(.message, borderwidth=NULL)) - 2)/
            as.numeric(tkfont.measure(tkcget(.message, font=NULL), "0"))
        eval(parse(text=paste("options(width=", floor(width), ")", sep="")))
    }
    lines <- strwrap(lines)
    if (console.output) {
        if (sink.number() != 0) sink()
        for (jline in seq(along.with=lines)) {
            Header <- if (jline==1) getRcmdr("prefixes")[3] else getRcmdr("prefixes")[4]
            cat(paste(Header, lines[jline], "\n", sep=""))
        } 
    }
    else
        for (line in lines){
            tagName <- messageTag()
            tkinsert(.message, "end", paste(line, "\n", sep=""))
            tktag.add(.message, tagName, "end - 2 lines linestart", "end - 2 lines lineend")
            tktag.configure(.message, tagName, foreground=col)
            tkyview.moveto(.message, 1.0)
        }
}

messageTag <- function(reset=FALSE){
    if (reset){
        putRcmdr("tagNumber", 0)
        return()
    }
    tagNumber <- getRcmdr("tagNumber") + 1
    putRcmdr("tagNumber", tagNumber)
    paste("message", tagNumber, sep="")
}

pushOutput <- function(element) {
    stack <- getRcmdr("outputStack")
    stack <- c(list(element), stack[-getRcmdr("length.output.stack")])
    putRcmdr("outputStack", stack)
}

popOutput <- function(keep=FALSE){
    stack <- getRcmdr("outputStack")
    lastOutput <- stack[[1]]
    if (!keep) putRcmdr("outputStack", c(stack[-1], NA))
    lastOutput
}

pushCommand <- function(element) {
    stack <- getRcmdr("commandStack")
    stack <- c(list(element), stack[-getRcmdr("length.command.stack")])
    putRcmdr("commandStack", stack)
}

popCommand <- function(keep=FALSE){
    stack <- getRcmdr("commandStack")
    lastCommand <- stack[[1]]
    if (!keep) putRcmdr("commandStack", c(stack[-1], NA))
    lastCommand
}
