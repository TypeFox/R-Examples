if (getRversion() >= '2.15.1') globalVariables(c(
    "tl", "last.table", "top", "dtm", "corpus", "corpusList", "corpusCa", "buttonsFrame",
    "clustDtm", "corpusClust", "clusters", "corpusVars", "diss",
    "freqTerms", "sourceVariable", "lowercaseVariable",
    "punctuationVariable", "numbersVariable", "stopwordsVariable",
    "stemmingVariable", "customStemmingVariable", "sourceFrame", "processingFrame",
    "foreign", "channel", "corpusDataset", "messages",
    "messageBox", "whatVariable", "whatFrame",
    "wordsDtm", "docLabelsVariable", "termLabelsVariable",
    "varLabelsVariable", "docPointsVariable", "termPointsVariable",
    "varPointsVariable", "labelsFrame", "pointsFrame",
    "preventDoubleClick", "termFreqs", "termSeries",
    "coocs", "specTerms", "varFreqs", "docSeries", "docButton",
    "globalButton", "unitVariable", "totaltVariable", "uniqueVariable",
    "hapaxVariable", "totalwVariable", "longVariable", "vlongVariable",
    "longavgVariable", "voc", "unitFrame", "digitsVariable",
    "exclRetweetsVariable", "removeNamesVariable", "removeHashtagsVariable",
    "optionsFrame", "twitCred", "odbcDataSources", "sqlTables"
    ))

.titleLabel <- function(...) labelRcmdr(..., font="RcmdrTitleFont",
                                        fg=getRcmdr("title.color"))

# Private environment inspired from Rcmdr's
.env <- new.env(parent=emptyenv())

.putEnv <- function(x, value) assign(x, value, envir=.env)

.getEnv <- function(x, mode="any", fail=TRUE){
    if ((!fail) && (!exists(x, mode=mode, envir=.env, inherits=FALSE))) return(NULL)
    get(x, envir=.env, mode=mode, inherits=FALSE)
}

.getCorpusWindow <- function() {
    if(!is.null(.getEnv("corpusTxt", fail=FALSE))) {
        window <- .getEnv("corpusWindow")
        txt <- .getEnv("corpusTxt")
        listbox <- .getEnv("corpusList")
        tkdelete(txt, "0.0", "end")
        tkdelete(listbox, 0, "end")
    }
    else {
        window <- tktoplevel(class="Rcommander")
        tkwm.geometry(window, "-0+20")
        scr1 <- tkscrollbar(window, repeatinterval=5,
                           command=function(...) tkyview(txt,...))
        txt <- tktext(window, bg="white", font="courier 11", wrap="word",
                      width=getOption("width", 80),
                      yscrollcommand=function(...) tkset(scr1, ...))

        tktag.configure(txt, "body", font="times")
        tktag.configure(txt, "heading", font="sans 13 bold")
        tktag.configure(txt, "articlehead", font="sans 12 bold")
        tktag.configure(txt, "details", font="sans 10 italic")
        tktag.configure(txt, "small", font="sans 5")
        tktag.configure(txt, "fixed", font="courier 11")

        menu <- tkmenu(txt, tearoff=FALSE)
        tkadd(menu, "command", label=.gettext("Select All"),
              command=function() tktag.add(txt, "sel", "1.0", "end"))
        # "break" is needed to prevent default bindings from being fired
        tkbind(txt, "<Control-Key-a>", expression(tktag.add(txt, "sel", "1.0", "end"), break))
        tkadd(menu, "command", label=.gettext("Copy"),
              command=function() tkevent.generate(txt, "<<Copy>>"))
        tkbind(txt, "<Control-Key-c>", function() tkevent.generate(txt, "<<Copy>>"))
        tkbind(txt, "<Button-3>", function(x, y)
               tkpopup(menu,
                       as.integer(x) + as.integer(tkwinfo("rootx", txt)),
                       as.integer(y) + as.integer(tkwinfo("rooty", txt))))

        tkpack(txt, side="left", fill="both", expand=TRUE)
        tkpack(scr1, side="left", fill="y")

        scr2 <- tkscrollbar(window, repeatinterval=5,
                            command=function(...) tkyview(listbox, ...))
        listbox <- tklistbox(window, selectmode="single",
                             yscrollcommand=function(...) tkset(scr2,...))
        tkpack(listbox, side="left", fill="y")
        tkpack(scr2, side="left", fill="y")

        tkbind(listbox, "<<ListboxSelect>>", function() {
            tkyview(txt, paste("mark", tkcurselection(listbox), sep=""))
        })

        .putEnv("corpusWindow", window)
        .putEnv("corpusTxt", txt)
        .putEnv("corpusList", listbox)

	tkwm.protocol(window, "WM_DELETE_WINDOW", function() {
            tkdestroy(.getEnv("corpusWindow"))
            .putEnv("corpusWindow", NULL)
            .putEnv("corpusTxt", NULL)
            .putEnv("corpusList", NULL)
        })
    }

    list(window=window, txt=txt, listbox=listbox)
}

.checkAndInstall <- function(package, message) {
    if(!package %in% rownames(installed.packages())) {
            # Create a function because dialog does not close until function returns
            msgbox <- function() tkmessageBox(title=.gettext("Package required"), message=message,
                                              icon="question", type="yesno")

            if (tclvalue(msgbox()) != "yes")
                return(FALSE)

            setBusyCursor()
            on.exit(setIdleCursor())

            if(package %in% available.packages()[,1]) {
                install.packages(package)
            }
            else {
                tkmessageBox(title=.gettext("Package not available"),
                             message=sprintf(.gettext("Package %s is not available. Please check your Internet connection, restart R and try again."),
                                             package),
                             icon="error", type="ok")
                return(FALSE)
            }
    }

    if(!require(package, character.only=TRUE)) {
        tkmessageBox(title=.gettext("Could not load package"),
                     message=sprintf(.gettext("Package %s could not be loaded. See errors in the \"Messages\" area."), package),
                     icon="error", type="ok")

        return(FALSE)
    }

    return(TRUE)
}

.Message <- function(message, type=c("info", "error", "warning"), parent=NULL) {
    type <- match.arg(type)

    caption <- switch(type,
                      info=.gettext("Information"),
                      error=.gettext("Error"),
                      warning=.gettext("Warning"))

    if(is.null(parent))
        tk_messageBox("ok", message, caption, icon=type)
    else
        tk_messageBox("ok", message, caption, icon=type, parent=parent)
}

# Adapted from Rcmdr's OKCancelHelp()
# Copyright John Fox. License: GPL >= 2.
.customCloseHelp <- defmacro(window=top, helpSubject=NULL, model=FALSE,
                             reset=NULL, apply=NULL, custom.button="OK",
    expr={
        memory <- getRcmdr("retain.selections")
        button.strings <- c(custom.button, "Close", 
                            if (!is.null(helpSubject)) "Help", 
                            if (!is.null(reset) && memory) "Reset", 
                            if (!is.null(apply)) "Apply")
        width <- max(nchar(c(gettextRcmdr("Help", "Reset", "Apply"), .gettext("Close"), custom.button)))
        if (WindowsP()) width <- width + 2
        buttonsFrame <- tkframe(window)
        leftButtonsBox <- tkframe(buttonsFrame)
        rightButtonsBox <- tkframe(buttonsFrame)
        
        customButton <- buttonRcmdr(rightButtonsBox, text=custom.button, foreground="darkgreen", width=width, command=onCustom, default="active",
            image="::image::okIcon", compound="left")
        
        onClose <- function() {
            if (exists(".exit")){
                result <- .exit()
                if (result == "abort") return()
            }
            putRcmdr("restoreTab", FALSE)
            if (model) putRcmdr("modelNumber", getRcmdr("modelNumber") - 1)
            if (GrabFocus()) tkgrab.release(window)
            tkdestroy(window)
            tkfocus(CommanderWindow())
        }
        
        closeButton <- buttonRcmdr(rightButtonsBox, text=.gettext("Close"), foreground="red", width=width, command=onClose, # borderwidth=3,
            image="::image::cancelIcon", compound="left")
        
        if (!is.null(helpSubject)){
            onHelp <- function() {
                if (GrabFocus() && (!WindowsP())) tkgrab.release(window)
                if (as.numeric(R.Version()$major) >= 2) print(help(helpSubject))
                else help(helpSubject)
            }
            helpButton <- buttonRcmdr(leftButtonsBox, text=gettextRcmdr("Help"), width=width, command=onHelp, # borderwidth=3,
                image="::image::helpIcon", compound="left")
        }
        
        if (!is.null(reset) && memory){
            onReset <- function(){
                ID <- window$ID
                putRcmdr("cancelDialogReopen", TRUE)
                putRcmdr("open.dialog.here", as.character(.Tcl(paste("winfo geometry", ID))))
                if (model) putRcmdr("modelNumber", getRcmdr("modelNumber") - 1)
                putDialog(reset, NULL)
                putDialog(reset, NULL, resettable=FALSE)
                closeDialog()
                eval(parse(text=paste(reset, "()")))
                putRcmdr("open.dialog.here", NULL)
                putRcmdr("restoreTab", FALSE)
            }
            resetButton <- buttonRcmdr(leftButtonsBox, text=gettextRcmdr("Reset"), width=width, command=onReset,
                image="::image::resetIcon", compound="left")
        }
        
        if (!is.null(apply)){
            onApply <- function(){
                putRcmdr("restoreTab", TRUE)
                putRcmdr("cancelDialogReopen", FALSE)
                ID <- window$ID
                putRcmdr("open.dialog.here", as.character(.Tcl(paste("winfo geometry", ID))))
                if (getRcmdr("use.markdown")) {
                    putRcmdr("startNewCommandBlock", FALSE)
                    beginRmdBlock()
                }
                if (getRcmdr("use.knitr")) {
                    putRcmdr("startNewKnitrCommandBlock", FALSE)
                    beginRnwBlock()
                }
                setBusyCursor()
                on.exit(setIdleCursor())
                onCustom()
                if (getRcmdr("use.markdown")){
                    removeNullRmdBlocks()
                    putRcmdr("startNewCommandBlock", TRUE)
                    if (getRcmdr("rmd.generated")) {
                        endRmdBlock()
                        putRcmdr("rmd.generated", FALSE)
                    }
                    removeNullRmdBlocks()
                }
                if (getRcmdr("use.knitr")){
                    removeNullRnwBlocks()
                    putRcmdr("startNewKnitrCommandBlock", TRUE)
                    if (getRcmdr("rnw.generated")) {
                        endRnwBlock()
                        putRcmdr("rnw.generated", FALSE)
                    }
                    removeNullRnwBlocks()
                }
                if (getRcmdr("cancelDialogReopen")){
                    putRcmdr("cancelDialogReopen", FALSE)
                }
                else{
                    eval(parse(text=paste(apply, "()")))
                    putRcmdr("open.dialog.here", NULL)
                }
            }
            applyButton <- buttonRcmdr(rightButtonsBox, text=gettextRcmdr("Apply"), foreground="yellow", width=width, command=onApply,
                image="::image::applyIcon", compound="left")
        }
        
        if(!WindowsP()) {
            if (!is.null(apply)){
                tkgrid(applyButton, closeButton, customButton, sticky="w")
                tkgrid.configure(customButton, padx=c(6, 0))
            }
            else{
                tkgrid(closeButton, customButton, sticky="w")
            }
            tkgrid.configure(closeButton, padx=c(6, 6))
        }
        else {
            if (!is.null(apply)){
                tkgrid(customButton, closeButton, applyButton, sticky="w")
                tkgrid.configure(applyButton, padx=c(6, 0))
            }
            else{
                tkgrid(customButton, closeButton, sticky="w")
            }
            tkgrid.configure(customButton, padx=c(6, 6))
        }
        if (!is.null(reset) && memory) {
            if (! is.null(helpSubject)){
                tkgrid (helpButton, resetButton, pady=6)
            }
            else tkgrid (resetButton, pady=6)
            if (!WindowsP()) tkgrid.configure(resetButton, padx=c(0, 6))
        }
        else if (! is.null(helpSubject)){
            tkgrid(helpButton, pady=6)
        }
        tkgrid(leftButtonsBox, rightButtonsBox, pady=6, sticky="ew")
        if (!is.null(helpSubject)) tkgrid.configure(helpButton, padx=c(0, 18))
        else if (!is.null(reset) && memory) tkgrid.configure(resetButton, padx=c(0, 18))
        tkgrid.columnconfigure(buttonsFrame, 0, weight=1)
        tkgrid.columnconfigure(buttonsFrame, 1, weight=1)
        tkgrid.configure(leftButtonsBox, sticky="w")
        tkgrid.configure(rightButtonsBox, sticky="e")
    })

.validate.unum <- function(P, ..., fun=NULL) {
    # Empty string must be allowed so that the user can remove
    # all chars before typing a new value
    if(P == "") {
        tcl("expr", "true")
    }
    else if(suppressWarnings(!is.na(as.numeric(P))) && as.numeric(P) > 0) {
        if(!is.null(fun)) fun(as.numeric(P))
        tcl("expr", "true")
    }
    else {
        tcl("expr", "false")
    }
}

.validate.uint <- function(P, ..., fun=NULL) {
    # Empty string must be allowed so that the user can remove
    # all chars before typing a new value
    if(P == "") {
        tcl("expr", "true")
    }
    else if(suppressWarnings(!is.na(as.integer(P))) && as.integer(P) > 0) {
        if(!is.null(fun)) fun(as.integer(P))
        tcl("expr", "true")
    }
    else {
        tcl("expr", "false")
    }
}
