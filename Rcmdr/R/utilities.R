# last modified 2016-03-31 by J. Fox

# utility functions

# listing objects etc.

listDataSets <- function(envir=.GlobalEnv, ...) {
    Vars <- ls(envir = envir, all.names = TRUE) # + PhG
    if (length(Vars) == 0) return(Vars) # + PhG
    
    names(which(sapply(Vars, function(.x) is.data.frame(get(.x, envir=envir)))))
}

listLinearModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
        function(.x) "lm" == (class(get(.x, envir=envir))[1]))]
}

listAOVModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
        function(.x) "aov" == (class(get(.x, envir=envir))[1]))]
}

listGeneralizedLinearModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
        function(.x) "glm" == (class(get(.x, envir=envir))[1]))]
}

listMultinomialLogitModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
        function(.x) "multinom" == (class(get(.x, envir=envir))[1]))]
}

listProportionalOddsModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
        function(.x) "polr" == (class(get(.x, envir=envir))[1]))]
}

listAllModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
        function(.x) (class(get(.x, envir=envir))[1])) %in% getRcmdr("modelClasses")]
}

activeDataSet <- function(dsname, flushModel=TRUE, flushDialogMemory=TRUE){
    .activeDataSet <- ActiveDataSet()
    if (missing(dsname)) {
        if (is.null(.activeDataSet)){
            Message(message=gettextRcmdr("There is no active data set."), type="error")
            return(FALSE)
        }
        else return(.activeDataSet)
    }
    if (!is.data.frame(ds <- get(dsname, envir=.GlobalEnv))){
        if (!exists.method("as.data.frame", ds, default=FALSE)){
            Message(message=paste(dsname, gettextRcmdr(" is not a data frame and cannot be attached."),
                sep=""), type="error")
            tkfocus(CommanderWindow())
            return()
        }
        command <- paste(dsname, " <- as.data.frame(", dsname, ")", sep="")
        justDoIt(command)
        logger(command)
        Message(message=paste(dsname, gettextRcmdr(" has been coerced to a data frame"), sep=""),
            type="warning")
    }
    varnames <- names(get(dsname, envir=.GlobalEnv))
    newnames <- make.names(varnames)
    badnames <- varnames != newnames
    if (any(badnames)){
        command <- paste("names(", dsname, ") <- make.names(names(",
            dsname, "))", sep="")
        doItAndPrint(command)
    }
    if (!is.null(.activeDataSet) && getRcmdr("attach.data.set")
        && (length(grep(.activeDataSet, search())) !=0)) {
        detach(pos = match(.activeDataSet, search()))
        logger(paste("detach(", .activeDataSet, ")", sep=""))
    }
    if (flushModel) {
        putRcmdr(".activeModel", NULL)
        RcmdrTclSet("modelName", gettextRcmdr("<No active model>"))
        tkconfigure(getRcmdr("modelLabel"), foreground="red")
    }
    if (flushDialogMemory) putRcmdr("dialog.values", list())
    ActiveDataSet(dsname)
    Message(sprintf(gettextRcmdr("The dataset %s has %d rows and %d columns."), dsname,
        nrow(get(dsname, envir=.GlobalEnv)), ncol(get(dsname, envir=.GlobalEnv))), type="note")
    if (any(badnames)) Message(message=paste(dsname, gettextRcmdr(" contains non-standard variable names:\n"),
        paste(varnames[badnames], collapse=", "),
        gettextRcmdr("\nThese have been changed to:\n"), paste(newnames[badnames], collapse=", "),
        sep=""), type="warning")
    RcmdrTclSet("dataSetName", paste(" ", dsname, " "))
    tkconfigure(getRcmdr("dataSetLabel"), foreground="blue")
    activateMenus()
    dsname
}


activeModel <- function(model){
    if (missing(model)) {
        .activeModel <- ActiveModel()
        if (is.null(.activeModel)){
            Message(message=gettextRcmdr("There is no active model."), type="error")
            return(FALSE)
        }
        else return(.activeModel)
    }
    ActiveModel(model)
    RcmdrTclSet("modelName", paste(" ", model, " "))
    tkconfigure(getRcmdr("modelLabel"), foreground="blue")
    activateMenus()
    model
}

listVariables <- function(dataSet=ActiveDataSet()) {
    if(missing(dataSet)) {
        Variables()
    }
    else {
        vars <- names(get(dataSet, envir=.GlobalEnv))
        if (getRcmdr("sort.names")) sortVarNames(vars) else vars
    }
}

listFactors <- function(dataSet=ActiveDataSet()) {
  if(missing(dataSet)) {
    Factors()
  }
  else {
    variables <- listVariables(dataSet)
    variables[sapply(variables, function(.x){
      .v <- eval(parse(text=.x), envir=get(dataSet, envir=.GlobalEnv))
      is.factor(.v) || is.logical(.v) || is.character(.v)
    })]
  }
}

listTwoLevelFactors <- function(dataSet=ActiveDataSet()){
  if(missing(dataSet)) {
    TwoLevelFactors()
  }
  else {
    factors <- listFactors(dataSet)
    if(length(factors) == 0) return(NULL)
    factors[sapply(factors, function(.x){
      .v <- eval(parse(text=.x), envir=get(dataSet, envir=.GlobalEnv))
      2 == length(levels(.v)) || length(unique(.v)) == 2
    })]
  }
}

listNumeric <- function(dataSet=ActiveDataSet()) {
    if(missing(dataSet)) {
        Numeric()
    }
    else {
        variables <- listVariables(dataSet)
        variables[sapply(variables,function(.x)
            is.numeric(eval(parse(text=.x), envir=get(dataSet, envir=.GlobalEnv))))]
    }
}

trim.blanks <- function(text){
    gsub("^\ *", "", gsub("\ *$", "", text))
}

is.valid.name <- function(x){
    length(x) == 1 && is.character(x) && x == make.names(x)
}


# statistical

Coef <- function(object, ...) UseMethod("Coef")

Coef.default <- function(object, ...) coef(object, ...)

Coef.multinom <- function (object, ...) {
    # the following adapted from nnet:
    cf <- function (object, ...) 
    {
        r <- length(object$vcoefnames)
        if (length(object$lev) == 2L) {
            coef <- object$wts[1L + (1L:r)]
            names(coef) <- object$vcoefnames
        }
        else {
            coef <- matrix(object$wts, nrow = object$n[3L], byrow = TRUE)[, 
                                                                          1L + (1L:r), drop = FALSE]
            if (length(object$lev)) 
                dimnames(coef) <- list(object$lev, object$vcoefnames)
            if (length(object$lab)) 
                dimnames(coef) <- list(object$lab, object$vcoefnames)
            coef <- coef[-1L, , drop = FALSE]
        }
        coef
    }
    b <- cf(object, ...)
    cn <- colnames(b)
    rn <- rownames(b)
    b <- as.vector(t(b))
    names(b) <- as.vector(outer(cn, rn, function(c, r) paste(r, c, sep = ":")))
    b
}


Confint <- function(object, parm, level=0.95, ...) UseMethod("Confint")

Confint.default <- function(object, parm, level = 0.95, ...) {
    ci <- confint(object, parm, level, ...)
    ci <- cbind(Coef(object)[parm], ci)
    colnames(ci)[1] <- "Estimate"
    ci
}

Confint.glm <- function (object, parm, level=0.95, type=c("LR", "Wald"), ...){
    # adapted from stats:::confint.lm
    type <- match.arg(type)
    cf <- coef(object)
    pnames <- names(cf)
    if (type == "LR") 
        ci <- confint(object, parm, level, ...)
    else {
        if (missing(parm))
            parm <- seq(along = pnames)
        else if (is.character(parm))
            parm <- match(parm, pnames, nomatch = 0)
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        pct <- paste(round(100 * a, 1), "%")
        ci <- array(NA, dim = c(length(parm), 2), dimnames = list(pnames[parm],
                                                                  pct))
        ses <- sqrt(diag(vcov(object)))[parm]
        fac <- qnorm(a)
        ci[] <- cf[parm] + ses %o% fac
    }
    ci <- cbind(cf[parm], ci)
    colnames(ci)[1] <- "Estimate"
    fam <- family(object)
    if (((fam$family == "binomial" || fam$family == "quasibinomial")  && fam$link == "logit")
      || ((fam$family == "poisson" || fam$family == "quasipoisson")  && fam$link == "log"))
      {
        expci <- exp(ci)
        colnames(expci)[1] <- "exp(Estimate)"
        ci <- cbind(ci, expci)
    }
    ci
}

Confint.polr <- function (object, parm, level=0.95, ...){
    # adapted from stats:::confint.lm
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm))
        parm <- seq(along = pnames)
    else if (is.character(parm))
        parm <- match(parm, pnames, nomatch = 0)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(round(100 * a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(pnames[parm],
                                                              pct))
    ses <- sqrt(diag(vcov(object)))[parm]
    fac <- qnorm(a)
    ci[] <- cf[parm] + ses %o% fac
    ci <- cbind(cf[parm], ci)
    colnames(ci)[1] <- "Estimate"
    ci
}

# confint.multinom <- function (object, parm, level=0.95, ...){
#     # adapted from stats:::confint.lm
#     cf <- coef(object)
#     if (is.vector(cf)) cf <- matrix(cf, nrow=1,
#         dimnames=list(object$lev[2], names(cf)))
#     pnames <- colnames(cf)
#     if (missing(parm))
#         parm <- seq(along = pnames)
#     else if (is.character(parm))
#         parm <- match(parm, pnames, nomatch = 0)
#     a <- (1 - level)/2
#     a <- c(a, 1 - a)
#     ses <- matrix(sqrt(diag(vcov(object))),
#         ncol=ncol(cf), byrow=TRUE)[,parm, drop=FALSE]
#     cf <- cf[,parm, drop=FALSE]
#     fac <- qnorm(a)
#     ci <- abind::abind(cf + fac[1]*ses, cf + fac[2]*ses, along=3)
#     dimnames(ci)[[3]] <- paste(round(100 * a, 1), "%")
#     aperm(ci, c(2,3,1))[,,1]
# }

Confint.multinom <- function(object, parm, level = 0.95, ...) {
    # adapted from stats:::confint.lm
    cf <- Coef(object)
    if (is.vector(cf)) cf <- matrix(cf, nrow=1,
                                    dimnames=list(object$lev[2], names(cf)))
    pnames <- colnames(cf)
    if (missing(parm))
        parm <- seq(along = pnames)
    else if (is.character(parm))
        parm <- match(parm, pnames, nomatch = 0)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    ses <- matrix(sqrt(diag(vcov(object))),
                  ncol=ncol(cf), byrow=TRUE)[,parm, drop=FALSE]
    cf <- cf[,parm, drop=FALSE]
    fac <- qnorm(a)
    ci <- abind::abind(cf + fac[1]*ses, cf + fac[2]*ses, along=3)
    dimnames(ci)[[3]] <- paste(round(100 * a, 1), "%")
    ci <- aperm(ci, c(2,3,1))[,,1]
    ci <- cbind(cf[parm], ci)
    colnames(ci)[1] <- "Estimate"
    ci
}

# Pager

# this is slightly modified from tkpager to use the Rcmdr monospaced font
#   and a white background

RcmdrPager <- function (file, header, title, delete.file)
{
    title <- paste(title, header)
    for (i in seq(along = file)) {
        zfile <- file[[i]]
        tt <- tktoplevel()
        if (WindowsP()) tkwm.iconbitmap(tt, system.file("etc", "R-logo.ico", package="Rcmdr"))
        tkwm.title(tt, if (length(title))
            title[(i - 1)%%length(title) + 1]
            else "")
        txt <- tktext(tt, bg = "white", font = getRcmdr("logFont"))
        scr <- ttkscrollbar(tt, command = function(...) tkyview(txt,
            ...))
        tkconfigure(txt, yscrollcommand = function(...) tkset(scr,
            ...))
        tkpack(txt, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
        chn <- tcl("open", zfile)
        tkinsert(txt, "end", gsub("_\b", "", tclvalue(tcl("read",
            chn))))
        tcl("close", chn)
        tkconfigure(txt, state = "disabled")
        tkmark.set(txt, "insert", "0.0")
        tkfocus(txt)
        if (delete.file)
            tcl("file", "delete", zfile)
    }
}

# help functions

helpCommander <- function() {
    PDF <- file.access(paste(file.path(path.package(package="Rcmdr")[1], "doc"), 
        "/", gettextRcmdr("Commander"), ".pdf", sep=""), mode=4)
    if (PDF == 0){
        browseURL(paste(file.path(path.package(package="Rcmdr")[1], "doc"),
            "/", gettextRcmdr("Commander"), ".pdf", sep=""))
    } 
    else if (as.numeric(R.Version()$major) >= 2) print(help(gettextRcmdr("Commander")))
    else help(gettextRcmdr("Commander"))
}

helpAboutCommander <- function() {
    if (as.numeric(R.Version()$major) >= 2) print(help("Rcmdr"))
    else help("Rcmdr")
}

browseManual <- function() {
    browseURL(paste(file.path(path.package(package="Rcmdr")[1], "doc"),
        "/", gettextRcmdr("Getting-Started-with-the-Rcmdr"), ".pdf", sep=""))
}

browseEnglishManual <- function() {
    browseURL(paste(file.path(path.package(package="Rcmdr")[1], "doc"),
        "/Getting-Started-with-the-Rcmdr.pdf", sep=""))
}

manualTranslationP <- function(){
    gettextRcmdr("Getting-Started-with-the-Rcmdr") != "Getting-Started-with-the-Rcmdr"
}

browseRcmdrWebsite <- function() browseURL("http://socserv.socsci.mcmaster.ca/jfox/Misc/Rcmdr/")

browseRWebsite <- function() browseURL("http://www.r-project.org/")

browseRMarkdown <- function() browseURL("http://rmarkdown.rstudio.com/")



# functions for building dialog boxes

# the following function is slightly modified, with permission, from Thomas Lumley, 
#   "Programmer's Niche: Macros in R," R-News, Sept. 2001, Vol. 1, No. 3, pp.11-13.
defmacro <- function(..., expr){
    expr <- substitute(expr)
    len <- length(expr)
    expr[3:(len+1)] <- expr[2:len]
    ## delete "macro" variables starting in ..
    expr[[2]] <- quote(on.exit(remove(list=objects(pattern="^\\.\\.", all.names=TRUE))))
    a <- substitute(list(...))[-1]
    ## process the argument list
    nn <- names(a)
    if (is.null(nn)) nn <- rep("", length(a))
    for (i in seq(length.out=length(a))){
        if (nn[i] == "") {
            nn[i] <- paste(a[[i]])
            msg <- paste(a[[i]], gettext("not supplied", domain="R-Rcmdr"))
            a[[i]] <- substitute(stop(foo), list(foo = msg))
        }
    }
    names(a) <- nn
    a <- as.list(a)
    ff <- eval(substitute(
        function(){
            tmp <- substitute(body)
            eval(tmp, parent.frame())
        },
        list(body = expr)))
    ## add the argument list
    formals(ff) <- a
    ## create a fake source attribute
    mm <- match.call()
    mm$expr <- NULL
    mm[[1]] <- as.name("macro")
    expr[[2]] <- NULL # get "local" variable removal out of source
    attr(ff, "source") <- c(deparse(mm), deparse(expr))
    ## return the macro
    ff
}

OKCancelHelp <- defmacro(window=top, helpSubject=NULL,  model=FALSE, reset=NULL, apply=NULL, helpPackage=NULL,
    expr={
        memory <- getRcmdr("retain.selections")
        button.strings <- c("OK", "Cancel", 
                            if (!is.null(helpSubject)) "Help", 
                            if (!is.null(reset) && memory) "Reset", 
                            if (!is.null(apply)) "Apply")
        width <- max(nchar(gettextRcmdr(button.strings)))
        if (WindowsP()) width <- width + 2
        buttonsFrame <- tkframe(window)
        leftButtonsBox <- tkframe(buttonsFrame)
        rightButtonsBox <- tkframe(buttonsFrame)
        
        OnOK <- function(){
            putRcmdr("restoreTab", FALSE)
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
            onOK()
            if (model) putDialog ("effectPlots", NULL)
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
            putRcmdr("rgl.command", FALSE)
        }
        
        OKbutton <- buttonRcmdr(rightButtonsBox, text=gettextRcmdr("OK"), foreground="darkgreen", width=width, command=OnOK, default="active",
            image="::image::okIcon", compound="left")
        
        onCancel <- function() {
            if (exists(".exit")){
                result <- .exit()
                if (result == "abort") return()
            }
            putRcmdr("restoreTab", FALSE)
            if (model) putRcmdr("modelNumber", getRcmdr("modelNumber") - 1)
            if (GrabFocus()) tkgrab.release(window)
            tkdestroy(window)
            putRcmdr("rgl.command", FALSE)
            tkfocus(CommanderWindow())
        }
        
        cancelButton <- buttonRcmdr(rightButtonsBox, text=gettextRcmdr("Cancel"), foreground="red", width=width, command=onCancel, # borderwidth=3,
            image="::image::cancelIcon", compound="left")
        
        if (!is.null(helpSubject)){
            onHelp <- function() {
                if (GrabFocus() && (!WindowsP())) tkgrab.release(window)
                if (as.numeric(R.Version()$major) >= 2) print(help(helpSubject, package=helpPackage))
                else help(helpSubject, package=helpPackage)
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
                onOK()
                putRcmdr("rgl.command", FALSE)
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
                tkgrid(applyButton, cancelButton, OKbutton, sticky="w")
                tkgrid.configure(OKbutton, padx=c(6, 0))
            }
            else{
                tkgrid(cancelButton, OKbutton, sticky="w")
            }
            tkgrid.configure(cancelButton, padx=c(6, 6))
        }
        else {
            if (!is.null(apply)){
                tkgrid(OKbutton, cancelButton, applyButton, sticky="w")
                tkgrid.configure(applyButton, padx=c(6, 0))
            }
            else{
                tkgrid(OKbutton, cancelButton, sticky="w")
            }
            tkgrid.configure(OKbutton, padx=c(6, 6))
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

subOKCancelHelp <- defmacro(window=subdialog, helpSubject=NULL,
    expr={
        
        button.strings <- c("OK", "Cancel", 
                            if (!is.null(helpSubject)) "Help")
        width <- max(nchar(gettextRcmdr(button.strings)))
        if (WindowsP()) width <- width + 2
        subButtonsFrame <- tkframe(window)
        subLeftButtonsBox <- tkframe(subButtonsFrame)
        subRightButtonsBox <- tkframe(subButtonsFrame)
        subOKbutton <- buttonRcmdr(subRightButtonsBox, text=gettextRcmdr("OK"), foreground="darkgreen", width=width, command=onOKsub, default="active",
            image="::image::okIcon", compound="left")
        onCancelSub <- function() {
          if (exists(".subexit")){
            .subexit()
          }
          if (GrabFocus()) tkgrab.release(window)
          tkdestroy(window)
          tkfocus(CommanderWindow())
        }
        subCancelButton <- buttonRcmdr(subRightButtonsBox, text=gettextRcmdr("Cancel"), foreground="red", width=width, command=onCancelSub,
            image="::image::cancelIcon", compound="left") # borderwidth=3, 
        if (!is.null(helpSubject)){
            onHelpSub <- function(){
                if (GrabFocus() && (!WindowsP())) tkgrab.release(window)
                if (as.numeric(R.Version()$major) >= 2) print(help(helpSubject))
                else help(helpSubject)
            }
            subHelpButton <- buttonRcmdr(subLeftButtonsBox, text=gettextRcmdr("Help"), width=width, command=onHelpSub, 
                image="::image::helpIcon", compound="left")
        }
        if(!WindowsP()) {
            tkgrid(subCancelButton, subOKbutton, sticky="w")
            tkgrid.configure(subOKbutton, padx=c(6, 0))
        }
        else {
            tkgrid(subOKbutton, subCancelButton, sticky="w")
            tkgrid.configure(subOKbutton, padx=c(0, 6))
        }
        if (! is.null(helpSubject)){
            tkgrid(subHelpButton, pady=6, padx=c(0, 18))
        }
        tkgrid(subLeftButtonsBox, subRightButtonsBox, pady=6, sticky="ew")
        tkgrid.columnconfigure(subButtonsFrame, 0, weight=1)
        tkgrid.columnconfigure(subButtonsFrame, 1, weight=1)
        tkgrid.configure(subLeftButtonsBox, sticky="w")
        tkgrid.configure(subRightButtonsBox, sticky="e")
    })

checkActiveDataSet <- function(){
    if (activeDataSet() == FALSE) {
        tkfocus(CommanderWindow())
        FALSE
    }
    else TRUE
}

checkActiveModel <- function(){
    if (activeModel() == FALSE) {
        tkfocus(CommanderWindow())
        FALSE
    }
    else TRUE
}

checkFactors <- function(n=1){
    if (length(Factors()) < n){
        if (n > 1)
            Message(message=sprintf(gettextRcmdr("There fewer than %d factors in the active data set."), n),
                type="error")
        else Message(message=gettextRcmdr("There are no factors in the active data set."),
            type="error")
        tkfocus(CommanderWindow())
        FALSE
    }
    else TRUE
}

checkTwoLevelFactors <- function(n=1){
    if (length(TwoLevelFactors()) < n){
        if (n > 1)
            Message(message=sprintf(gettextRcmdr("There fewer than %d two-level factors in the active data set."), n),
                type="error")
        else Message(message=gettextRcmdr("There are no two-level factors in the active data set."),
            type="error")
        tkfocus(CommanderWindow())
        FALSE
    }
    else TRUE
}

checkNumeric <- function(n=1){
    if (length(Numeric()) < n){
        if (n > 1)
            Message(message=sprintf(gettextRcmdr("There fewer than %d numeric variables in the active data set."), n),
                type="error")
        else Message(message=gettextRcmdr("There are no numeric variables in the active data set."),
            type="error")
        tkfocus(CommanderWindow())
        FALSE
    }
    else TRUE
}

checkVariables <- function(n=1){
    if (length(Variables()) < n){
        if (n > 1)
            Message(message=sprintf(gettextRcmdr("There fewer than %d variables in the active data set."), n),
                type="error")
        else Message(message=gettextRcmdr("There are no variables in the active data set."),
            type="error")
        tkfocus(CommanderWindow())
        FALSE
    }
    else TRUE
}

commanderPosition <- function (){
    ID <- CommanderWindow()$ID
    as.numeric(c(tclvalue(.Tcl(paste("winfo rootx", ID))),
        tclvalue(.Tcl(paste("winfo rooty", ID)))))
}

initializeDialog <- defmacro(window=top, title="", offset=10, preventCrisp, 
    use.tabs=FALSE, notebook=notebook, tabs=c("dataTab", "optionsTab"),
    suppress.window.resize.buttons=TRUE,
    expr={
        if (getRcmdr("crisp.dialogs")) tclServiceMode(on=FALSE)
        window <- tktoplevel(borderwidth=10)
        if (use.tabs){
            notebook <- ttknotebook(window)
            for (tab in tabs) assign(tab, tkframe(window))
        }
        tkwm.title(window, title)
        location <- getRcmdr("open.dialog.here")
        position <- if (!is.null(location)) location
        else {
            pos <- offset + commanderPosition() 
            if (any(pos < 0)) "-50+50"
            else paste("+", paste(pos, collapse="+"), sep="")
        }
        tkwm.geometry(window, position)
        if (suppress.window.resize.buttons) tkwm.transient(window, CommanderWindow())
    }
)

closeDialog <- defmacro(window=top, release=TRUE,
    expr={
        if (release && GrabFocus()) tkgrab.release(window)
        tkdestroy(window)
    }
)

dialogSuffix <- defmacro(window=top, onOK=onOK, onCancel=onCancel, rows, columns, focus=top,
    bindReturn=TRUE, preventGrabFocus=FALSE, preventDoubleClick=FALSE,
    preventCrisp, 
    use.tabs=FALSE, notebook=notebook, tabs=c("dataTab", "optionsTab"), tab.names=c("Data", "Options"),
    grid.buttons=FALSE, resizable=FALSE, force.wait=FALSE,
    expr={
        if (use.tabs){
            for (i in 1:length(tabs)){
                tkadd(notebook, get(tabs[i]), text=gettextRcmdr(tab.names[i]), padding=6, sticky="nsew")
            }
            tkgrid(notebook, sticky="nsew")
        }
        if (grid.buttons) tkgrid(buttonsFrame, sticky = "ew")
        if (use.tabs && exists("dialog.values") && !is.null(dialog.values$initial.tab) && getRcmdr("restoreTab")) 
            tkselect(notebook, dialog.values$initial.tab)
        .Tcl("update idletasks")
        tkwm.resizable(window, as.numeric(resizable), as.numeric(resizable))
        if (bindReturn) tkbind(window, "<Return>", onOK)
        tkbind(window, "<Escape>", onCancel)
        if (getRcmdr("double.click") && (!preventDoubleClick)) tkbind(window, "<Double-ButtonPress-1>", onOK)
        tkwm.deiconify(window)
        # focus grabs appear to cause problems for some dialogs
        if (GrabFocus() && (!preventGrabFocus)) tkgrab.set(window)
        tkfocus(focus)
        if (getRcmdr("tkwait.dialog") || force.wait) tkwait.window(window)
        if (getRcmdr("crisp.dialogs")) tclServiceMode(on=TRUE)
    }
)

variableListBox <- function(parentWindow, variableList=Variables(), bg="white",
    selectmode="single", export="FALSE", initialSelection=NULL, listHeight=getRcmdr("variable.list.height"), title){
    if (selectmode == "multiple") selectmode <- getRcmdr("multiple.select.mode")
    if (length(variableList) == 1 && is.null(initialSelection)) initialSelection <- 0
    frame <- tkframe(parentWindow)
    minmax <- getRcmdr("variable.list.width")
    listbox <- tklistbox(frame, height=min(listHeight, length(variableList)),
        selectmode=selectmode, background=bg, exportselection=export, 
        width=min(max(minmax[1], 2 + nchar(variableList)), minmax[2]))
    scrollbar <- ttkscrollbar(frame, command=function(...) tkyview(listbox, ...))
    tkconfigure(listbox, yscrollcommand=function(...) tkset(scrollbar, ...))
    for (var in variableList) tkinsert(listbox, "end", var)
    if (is.numeric(initialSelection)) for (sel in initialSelection) tkselection.set(listbox, sel)
    firstChar <- tolower(substr(variableList, 1, 1))
    len <- length(variableList)
    onLetter <- function(letter){
        letter <- tolower(letter)
        current <- 1 + round(as.numeric(unlist(strsplit(tclvalue(tkyview(listbox) ), " "))[1])*len)
        mat <- match(letter, firstChar[-(1:current)])
        if (is.na(mat)) return()
        tkyview.scroll(listbox, mat, "units")
    }
    onA <- function() onLetter("a")
    onB <- function() onLetter("b")
    onC <- function() onLetter("c")
    onD <- function() onLetter("d")
    onE <- function() onLetter("e")
    onF <- function() onLetter("f")
    onG <- function() onLetter("g")
    onH <- function() onLetter("h")
    onI <- function() onLetter("i")
    onJ <- function() onLetter("j")
    onK <- function() onLetter("k")
    onL <- function() onLetter("l")
    onM <- function() onLetter("m")
    onN <- function() onLetter("n")
    onO <- function() onLetter("o")
    onP <- function() onLetter("p")
    onQ <- function() onLetter("q")
    onR <- function() onLetter("r")
    onS <- function() onLetter("s")
    onT <- function() onLetter("t")
    onU <- function() onLetter("u")
    onV <- function() onLetter("v")
    onW <- function() onLetter("w")
    onX <- function() onLetter("x")
    onY <- function() onLetter("y")
    onZ <- function() onLetter("z")
    for (letter in c(letters, LETTERS)){
        tkbind(listbox, paste("<", letter, ">", sep=""),
            get(paste("on", toupper(letter), sep="")))
    }
    onClick <- function() tkfocus(listbox)
    toggleSelection <- function(){
        active <- tclvalue(tkindex(listbox, "active"))
        selected <- tclvalue(tkcurselection(listbox))
        if (selected == active) tkselection.clear(listbox, "active") else tkselection.set(listbox, "active")
    }
    tkbind(listbox, "<ButtonPress-1>", onClick)
    if (selectmode == "single") tkbind(listbox, "<Control-ButtonPress-1>", toggleSelection)
    tkgrid(labelRcmdr(frame, text=title, fg=getRcmdr("title.color"), font="RcmdrTitleFont"), columnspan=2, sticky="w")
    tkgrid(listbox, scrollbar, sticky="nw")
    tkgrid.configure(scrollbar, sticky="wns")
    tkgrid.configure(listbox, sticky="ewns")
    result <- list(frame=frame, listbox=listbox, scrollbar=scrollbar,
        selectmode=selectmode, varlist=variableList)
    class(result) <- "listbox"
    result
}

getSelection <- function(object) UseMethod("getSelection")

getSelection.listbox <- function(object){
    object$varlist[as.numeric(tkcurselection(object$listbox)) + 1]
}

getFrame <- function(object) UseMethod("getFrame")

getFrame.listbox <- function(object){
    object$frame
}

variableComboBox <- function(parentWindow, variableList=Variables(),
                             export="FALSE", state="readonly",
                             initialSelection=gettextRcmdr("<no variable selected>"),
                             title=""){
  variableList <- c(gettextRcmdr("<no variable selected>"), variableList)
  frame <- tkframe(parentWindow)
  combovar <- tclVar()
  tclvalue(combovar) <- initialSelection
  combobox <- ttkcombobox(frame, values=variableList, textvariable=combovar, 
                          state=state, export=export)
  firstChar <- tolower(substr(variableList, 1, 1))
  onLetter <- function(letter){
    letter <- tolower(letter)
    current <- as.numeric(tcl(combobox, "current"))
    current <- if (current == -1) 1 else current + 1
    mat <- match(letter, firstChar[-(1:current)])
    if (is.na(mat)) return()
    tcl(combobox, "current", current + mat - 1)
  }
  onA <- function() onLetter("a")
  onB <- function() onLetter("b")
  onC <- function() onLetter("c")
  onD <- function() onLetter("d")
  onE <- function() onLetter("e")
  onF <- function() onLetter("f")
  onG <- function() onLetter("g")
  onH <- function() onLetter("h")
  onI <- function() onLetter("i")
  onJ <- function() onLetter("j")
  onK <- function() onLetter("k")
  onL <- function() onLetter("l")
  onM <- function() onLetter("m")
  onN <- function() onLetter("n")
  onO <- function() onLetter("o")
  onP <- function() onLetter("p")
  onQ <- function() onLetter("q")
  onR <- function() onLetter("r")
  onS <- function() onLetter("s")
  onT <- function() onLetter("t")
  onU <- function() onLetter("u")
  onV <- function() onLetter("v")
  onW <- function() onLetter("w")
  onX <- function() onLetter("x")
  onY <- function() onLetter("y")
  onZ <- function() onLetter("z")
  for (letter in c(letters, LETTERS)){
    tkbind(combobox, paste("<", letter, ">", sep=""),
           get(paste("on", toupper(letter), sep="")))
  }
  tkgrid(labelRcmdr(frame, text=title, fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w") # , columnspan=2
  tkgrid(combobox, sticky="nw")
  result <- list(frame=frame, combobox=combobox, varlist=variableList, combovar=combovar)
  class(result) <- "combobox"
  result
}

getSelection.combobox <- function(object){
  tclvalue(object$combovar)
}

getFrame.combobox <- function(object){
  object$frame
}

# This function modified based on code by Liviu Andronic (13 Dec 09) and on code by Milan Bouchet-Valat (29 Jun 12):
radioButtons <- defmacro(window=top, name, buttons, values=NULL, initialValue=..values[1], labels, 
    title="", title.color=getRcmdr("title.color"), right.buttons=FALSE, command=function(){},
    expr={
        ..values <- if (is.null(values)) buttons else values
        ..frame <- paste(name, "Frame", sep="")
        assign(..frame, tkframe(window))
        ..variable <- paste(name, "Variable", sep="")
        assign(..variable, tclVar(initialValue))
        if(title != ""){
            tkgrid(labelRcmdr(eval(parse(text=..frame)), text=title, foreground=title.color, font="RcmdrTitleFont"), columnspan=2, sticky="w")
        }
        for (i in 1:length(buttons)) {
            ..button <- paste(buttons[i], "Button", sep="")
            if (right.buttons) {
                assign(..button, ttkradiobutton(eval(parse(text=..frame)), variable=eval(parse(text=..variable)), 
                    value=..values[i], command=command))
                tkgrid(labelRcmdr(eval(parse(text=..frame)), text=labels[i], justify="left"), eval(parse(text=..button)), sticky="w")
            }
            else{
                assign(..button, ttkradiobutton(eval(parse(text=..frame)), variable=eval(parse(text=..variable)), 
                    value=..values[i], text=labels[i], command=command))
                tkgrid(eval(parse(text=..button)), sticky="w")
            }
        }
    }
)


checkBoxes <- defmacro(window=top, frame, boxes, initialValues=NULL, labels, title=NULL, ttk=FALSE,
    expr={
        ..initialValues <- if (is.null(initialValues)) rep("1", length(boxes)) else initialValues
        assign(frame, if (ttk) ttklabelframe(window, labelwidget=tklabel(window, text=title, 
                                          font="RcmdrTitleFont", foreground=getRcmdr("title.color"))) else tkframe(window))
        if (!is.null(title) && !ttk) tkgrid(labelRcmdr(eval(parse(text=frame)), text=title, fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
        ..variables <- paste(boxes, "Variable", sep="")
        for (i in 1:length(boxes)) {
            assign(..variables[i], tclVar(..initialValues[i]))
            ..checkBox <- paste(boxes[i], "CheckBox", sep="")
            assign(..checkBox,
                ttkcheckbutton(eval(parse(text=frame)), variable=eval(parse(text=..variables[i])), text=labels[i]))
            tkgrid(eval(parse(text=..checkBox)), sticky="w")
        }
    }
)

checkReplace <- function(name, type=gettextRcmdr("Variable")){
    RcmdrTkmessageBox(message=sprintf(gettextRcmdr("%s %s already exists.\nOverwrite %s?"),
        type, name, tolower(type)), icon="warning", type="yesno", default="no")
}

errorCondition <- defmacro(window=top, recall=NULL, message, model=FALSE,
    expr={
        putRcmdr("cancelDialogReopen", TRUE)
        if (model) putRcmdr("modelNumber", getRcmdr("modelNumber") - 1)
        if (!is.null(window)){
            if (GrabFocus()) tkgrab.release(window)
            tkdestroy(window)
        }
        Message(message=message, type="error")
        if (!is.null(recall)) recall()
        else tkfocus(CommanderWindow())
    })

subsetBox <- defmacro(window=top, subset.expression=NULL, model=FALSE,
    expr={
        subsetVariable <- if (!is.null(subset.expression)) tclVar(gettextRcmdr(subset.expression))
        else if (model){
            if (currentModel && currentFields$subset != "")
                tclVar(currentFields$subset) else tclVar(gettextRcmdr("<all valid cases>"))
        }
        else tclVar(gettextRcmdr("<all valid cases>"))
        subsetFrame <- tkframe(window)
        subsetEntry <- ttkentry(subsetFrame, width="20", textvariable=subsetVariable)
        subsetScroll <- ttkscrollbar(subsetFrame, orient="horizontal",
            command=function(...) tkxview(subsetEntry, ...))
        tkconfigure(subsetEntry, xscrollcommand=function(...) tkset(subsetScroll, ...))
        tkgrid(labelRcmdr(subsetFrame, text=gettextRcmdr("Subset expression"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
        tkgrid(subsetEntry, sticky="ew")
        tkgrid(subsetScroll, sticky="ew")
        tkgrid.columnconfigure(subsetFrame, 0, weight=1)
    })


groupsBox <- defmacro(recall=NULL, label=gettextRcmdr("Plot by:"), initialLabel=gettextRcmdr("Plot by groups"),
                      errorText=gettextRcmdr("There are no factors in the active data set."),
                      variables=Factors(), plotLinesByGroup=FALSE, positionLegend=FALSE, plotLinesByGroupsText=gettextRcmdr("Plot lines by group"),
                      initialGroup=NULL, initialLinesByGroup=1, window=top,
                      expr={
                          env <- environment()
                          .groups <- if (is.null(initialGroup)) FALSE else initialGroup
                          .linesByGroup <- initialLinesByGroup == 1
                          .groupsLabel <- tclVar(if (!is.null(initialGroup)) initialLabel else paste(initialLabel, "...", sep=""))
                          .factors <- variables
                          onGroups <- function(){
                              if (length(.factors) == 0){
                                  errorCondition(recall=recall, message=errorText)
                                  return()
                              }
                              initializeDialog(subdialog, title=gettextRcmdr("Groups"))
                              groupsBox <- variableListBox(subdialog, .factors, title=gettextRcmdr("Groups variable (pick one)"),
                                                           initialSelection=varPosn(initialGroup, "factor"))
                              if (plotLinesByGroup){
                                  linesByGroupFrame <- tkframe(subdialog)
                                  linesByGroup <- tclVar(if(initialLinesByGroup == 1) "1" else "0")
                                  linesCheckBox <- ttkcheckbutton(linesByGroupFrame, variable=linesByGroup)
                                  tkgrid(labelRcmdr(linesByGroupFrame, text=plotLinesByGroupsText), linesCheckBox, sticky="w")
                              }
                              onOKsub <- function() {
                                  groups <- getSelection(groupsBox)
                                  if (length(groups) == 0){
                                      assign(".groups", FALSE, envir=env)
                                      tclvalue(.groupsLabel) <- paste(gettextRcmdr("Plot by groups"), "...", sep="")
                                      tkconfigure(groupsButton, foreground="black")
                                      if (GrabFocus()) tkgrab.release(subdialog)
                                      tkdestroy(subdialog)
                                      tkwm.deiconify(top)
                                      if (GrabFocus()) tkgrab.set(top)
                                      tkfocus(top)
                                      tkwait.window(top)
                                      return()
                                  }
                                  assign(".groups", groups, envir=env)
                                  tclvalue(.groupsLabel) <- paste(label, groups)
                                  tkconfigure(groupsButton, foreground=getRcmdr("title.color"))
                                  tkconfigure(groupsButton, font="RcmdrTitleFont")
                                  if (plotLinesByGroup) {
                                      lines <- as.character("1" == tclvalue(linesByGroup))
                                      assign(".linesByGroup", lines, envir=env)
                                  }
                                  if (GrabFocus()) tkgrab.release(subdialog)
                                  tkdestroy(subdialog)
                                  tkwm.deiconify(top)
                                  if (GrabFocus()) tkgrab.set(top)
                                  tkfocus(top)
                                  tkwait.window(top)
                              }
                              subOKCancelHelp()
                              tkgrid(getFrame(groupsBox), sticky="nw")
                              if (plotLinesByGroup) tkgrid(linesByGroupFrame, sticky="w")
                              tkgrid(subButtonsFrame, sticky="ew")
                              if (positionLegend) tkgrid(labelRcmdr(subdialog, text=gettextRcmdr("Position legend with mouse click"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"))
                              dialogSuffix(subdialog, onOK=onOKsub, focus=subdialog, force.wait=TRUE)
                          }
                          groupsFrame <- tkframe(window)
                          groupsButton <- tkbutton(groupsFrame, textvariable=.groupsLabel, command=onGroups)
                          if (!is.null(initialGroup)) tkconfigure(groupsButton, foreground=getRcmdr("title.color"), font="RcmdrTitleFont")
                          tkgrid(groupsButton, sticky="we")
                          tkgrid.columnconfigure(groupsFrame, 0, weight=1)
                      })


groupsLabel <- defmacro(frame=top, groupsBox=groupsBox, columnspan=1, initialText=NULL,
    expr={
        initial.label <- if (exists("dialog.values")) dialog.values$initial.label else NULL
        if  (is.null(initial.label)) {
            group <- getSelection(groupsBox)
            initial.label <- if (length(group) == 0) NULL 
            else {
                levels <- eval(parse(text = paste("levels(", ActiveDataSet(), 
                    "$", group, ")", sep = "")))
                paste(levels[1], "-", levels[2])
            }
        }
        groupsFrame <- tkframe(frame)
        .groupsLabel <- if (!is.null(initialText)) initialText 
        else if (is.null(initial.label)) gettextRcmdr("<No groups selected>") 
        else initial.label
        groupsLabel <- labelRcmdr(groupsFrame, text=.groupsLabel)
        tkgrid(labelRcmdr(groupsFrame, text=gettextRcmdr("Difference: "), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), groupsLabel, sticky="w")
        tkgrid(groupsFrame, sticky="w", columnspan=columnspan)
        onSelect <- function(){
            group <- getSelection(groupsBox)
            if (length(group) == 0) {
                .groupsLabel <<- gettextRcmdr("<No groups selected>") 
            }
            else {
                levels <- eval(parse(text=paste("levels(", ActiveDataSet(), "$", group, ")", sep="")))
                .groupsLabel <<- paste(levels[1], "-", levels[2])
            }
            tkconfigure(groupsLabel, text=.groupsLabel)
        }
        tkbind(groupsBox$listbox, "<ButtonRelease-1>", onSelect)
    })

modelFormula <- defmacro(frame=top, hasLhs=TRUE, rhsExtras=NULL, formulaLabel=gettextRcmdr("Model Formula"),
                         expr={
  .rhsExtras <- if (is.null(rhsExtras)) hasLhs else rhsExtras
  checkAddOperator <- function(rhs){
    rhs.chars <- rev(strsplit(rhs, "")[[1]])
    if (length(rhs.chars) < 1) return(FALSE)
    check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
      rhs.chars[1] else rhs.chars[2]
    !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
  }
  .variables <- Variables()
  word <- paste("\\[", gettextRcmdr("factor"), "\\]", sep="")
  variables <- paste(.variables,
                     ifelse(is.element(.variables, Factors()), paste("[", gettextRcmdr("factor"), "]", sep=""), ""))
  xBox <- variableListBox(frame, variables, selectmode="multiple", title=gettextRcmdr("Variables (double-click to formula)"))
  onDoubleClick <- if (!hasLhs){
    function(){
      var <- getSelection(xBox)
      tkselection.clear(xBox$listbox, "0", "end")            		
      if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
      tkfocus(rhsEntry)
      rhs <- tclvalue(rhsVariable)
      rhs.chars <- rev(strsplit(rhs, "")[[1]])
      check.char <- if (length(rhs.chars) > 0){
        if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
          rhs.chars[1] else rhs.chars[2]
      }
      else ""
      tclvalue(rhsVariable) <- if (rhs == "" ||
                                   is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
        paste(rhs, var, sep="")
      else paste(rhs, "+", var)
      tkicursor(rhsEntry, "end")
      tkxview.moveto(rhsEntry, "1")
    }
  }
  else{
    function(){
      var <- getSelection(xBox)
      which <- tkcurselection(xBox$listbox)
      tkselection.clear(xBox$listbox, "0", "end")
      if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
      lhs <- tclvalue(lhsVariable)
      if (lhs == "" || tclvalue(tkselection.present(lhsEntry)) == "1"){
        tclvalue(lhsVariable) <- var
        tkselection.clear(lhsEntry)
        tkfocus(rhsEntry)
      }
      else {
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        check.char <- if (length(rhs.chars) > 0){
          if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
            rhs.chars[1] else rhs.chars[2]
        }
        else ""
        tclvalue(rhsVariable) <- if (rhs == "" ||
                                     is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
          paste(rhs, var, sep="")
        else paste(rhs, "+", var)
      }
      tkicursor(rhsEntry, "end")
      tkxview.moveto(rhsEntry, "1")
    }
  }
  tkbind(xBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
  onPlus <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")										
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      if (length(var) > 1) var <- paste(var, collapse=" + ")
    }
    tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onTimes <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")						
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      var <- trim.blanks(var)
      if (length(var) > 1) var <- paste(var, collapse="*")
      tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    }
    else tclvalue(rhsVariable) <- paste(rhs, if (!check) "*", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onColon <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")						
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      var <- trim.blanks(var)
      if (length(var) > 1) var <- paste(var, collapse=":")
      tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    }
    else tclvalue(rhsVariable) <- paste(rhs, if (!check) ":", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onSlash <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onIn <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "%in% ")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onMinus <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "- ")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onPower <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onLeftParen <- function(){
    tkfocus(rhsEntry)
    rhs <- tclvalue(rhsVariable)
    tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onRightParen <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  outerOperatorsFrame <- tkframe(frame)
  operatorsFrame <- tkframe(outerOperatorsFrame)
  splinePolyFrame <- tkframe(outerOperatorsFrame)
  plusButton <- buttonRcmdr(operatorsFrame, text="+", width="3", command=onPlus)
  timesButton <- buttonRcmdr(operatorsFrame, text="*", width="3", command=onTimes)
  colonButton <- buttonRcmdr(operatorsFrame, text=":", width="3", command=onColon)
  slashButton <- buttonRcmdr(operatorsFrame, text="/", width="3", command=onSlash)
  inButton <- buttonRcmdr(operatorsFrame, text="%in%", width="5", command=onIn)
  minusButton <- buttonRcmdr(operatorsFrame, text="-", width="3", command=onMinus)
  powerButton <- buttonRcmdr(operatorsFrame, text="^", width="3", command=onPower)
  leftParenButton <- buttonRcmdr(operatorsFrame, text="(", width="3", command=onLeftParen)
  rightParenButton <- buttonRcmdr(operatorsFrame, text=")", width="3", command=onRightParen)
  onBSpline <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("spline requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + bs(", var, ", df=", tclvalue(dfSplineVar), ")", sep="") 
                                   else paste(" bs(", var, ", df=", tclvalue(dfSplineVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onNatSline <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("spline requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + ns(", var, ", df=", tclvalue(dfSplineVar), ")", sep="") 
                                   else paste(" ns(", var, ", df=", tclvalue(dfSplineVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onPoly <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("polynomial requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + poly(", var, ", degree=", tclvalue(degPolyVar), ")", sep="") 
                                   else paste(" poly(", var, ", degree=", tclvalue(degPolyVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onRawPoly <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("polynomial requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + poly(", var, ", degree=", tclvalue(degPolyVar), ", raw=TRUE)", sep="") 
                                   else paste(" poly(", var, ", degree=", tclvalue(degPolyVar), ", raw=TRUE)", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  bsplineButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("B-spline\n"), width="10", command=onBSpline)
  nsplineButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("natural\nspline"), width="10", command=onNatSline)
  polyButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("orthogonal\npolynomial"), width="10", command=onPoly)
  RawPolyButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("raw\npolynomial"), width="10", command=onRawPoly)
  dfSplineVar <- tclVar("5")
  degPolyVar <- tclVar("2")
  dfDegFrame <- tkframe(outerOperatorsFrame)
  dfSplineSpin <- tkspinbox(dfDegFrame, textvariable=dfSplineVar, state="readonly", from=2, to=10, width=2)
  degPolySpin <- tkspinbox(dfDegFrame, textvariable=degPolyVar, state="readonly", from=2, to=5, width=2)
  tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, minusButton,
         powerButton, leftParenButton, rightParenButton, sticky="w")
  tkgrid(labelRcmdr(dfDegFrame, text=gettextRcmdr("df for splines: ")), dfSplineSpin,  sticky="se")
  tkgrid(labelRcmdr(dfDegFrame, text=gettextRcmdr("deg. for polynomials: ")), degPolySpin, sticky="se")
  formulaFrame <- tkframe(frame)
  formulaFrameMain <- tkframe(formulaFrame)
  onFormulaHelp <- function () print(help("formula"))
  formulaHelpButton <- buttonRcmdr(formulaFrame, text=gettextRcmdr("Model formula\nhelp"), command=onFormulaHelp,
                                   image="::image::helpIcon", compound="left")
  if (hasLhs){
    tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Model Formula"), 
                      fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(labelRcmdr(outerOperatorsFrame, text="Operators (click to formula):  "), operatorsFrame, sticky="nw")
    if (.rhsExtras){
      tkgrid(bsplineButton, nsplineButton, polyButton, RawPolyButton, sticky="nw")
      tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Splines/Polynomials:\n(select variable and click)")), 
             splinePolyFrame, dfDegFrame, sticky="nw")
    }
    lhsVariable <- if (currentModel) tclVar(currentFields$lhs) else tclVar("")
    rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
    rhsEntry <- ttkentry(formulaFrameMain, width="75", textvariable=rhsVariable)
    rhsXscroll <- ttkscrollbar(formulaFrameMain,
                               orient="horizontal", command=function(...) tkxview(rhsEntry, ...))
    tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
    lhsEntry <- ttkentry(formulaFrameMain, width="10", textvariable=lhsVariable)
    lhsScroll <- ttkscrollbar(formulaFrameMain,
                              orient="horizontal", command=function(...) tkxview(lhsEntry, ...))
    tkconfigure(lhsEntry, xscrollcommand=function(...) tkset(lhsScroll, ...))
    tkgrid(lhsEntry, labelRcmdr(formulaFrameMain, text=" ~    "), rhsEntry, sticky="w")
    tkgrid(lhsScroll, labelRcmdr(formulaFrameMain, text=""), rhsXscroll, sticky="w")
    tkgrid.configure(lhsScroll, sticky="ew")
    tkgrid(formulaFrameMain, labelRcmdr(formulaFrame, text="  "), formulaHelpButton, sticky="nw")
  }
  else{
    if (.rhsExtras){
      tkgrid(labelRcmdr(outerOperatorsFrame, text=formulaLabel, 
                        fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
      tkgrid(labelRcmdr(outerOperatorsFrame, text="Operators (click to formula):  "), operatorsFrame, sticky="nw")
      tkgrid(bsplineButton, nsplineButton, polyButton, RawPolyButton, sticky="nw")
      tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Splines/Polynomials:\n(select variable and click)")), 
             splinePolyFrame, dfDegFrame, sticky="nw")
    }
    rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
    rhsEntry <- ttkentry(formulaFrameMain, width="75", textvariable=rhsVariable)
    rhsXscroll <- ttkscrollbar(formulaFrameMain,
                               orient="horizontal", command=function(...) tkxview(rhsEntry, ...))
    tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
    tkgrid(labelRcmdr(formulaFrameMain, text="   ~ "), rhsEntry, sticky="w")
    tkgrid(labelRcmdr(formulaFrameMain, text=""), rhsXscroll, sticky="w")
    tkgrid(formulaFrameMain, labelRcmdr(formulaFrame, text="  "), formulaHelpButton, sticky="nw")
  }
  tkgrid.configure(rhsXscroll, sticky="ew")
})

exists.method <- function(generic, object, default=TRUE, strict=FALSE){
    classes <- class(object)
    if (default) classes <- c(classes, "default")
    if (strict) classes <- classes[1]
    any(paste(generic, ".", classes, sep="") %in%
            as.character(methods(generic)))
}

checkMethod <- defmacro(generic, object, message=NULL, default=FALSE, strict=FALSE, reportError=TRUE,
    expr={
        msg <- if (is.null(message)) sprintf(gettextRcmdr("No appropriate %s method exists\nfor a model of this class."), generic)
        else message
        method <- exists.method(generic, get(object), default=default, strict=strict)
        if ((!method) && reportError) Message(message=msg, type="error")
        method
    }
)

checkClass <- defmacro(object, class, message=NULL,
    expr={
        msg <- if (is.null(message)) sprintf(gettextRcmdr('The model is not of class "%s".'), class)
        else message
        properClass <- class(get(object))[1] == class
        if (!properClass) Message(message=msg, type="error")
        properClass
    }
)

# the following function is from John Chambers (plus new test for R 2.4.0)

isS4object <- function(object) {
    if (getRversion() < "2.4.0"){
        if (length(attr(object, "class"))!= 1)
            return(FALSE)
        !isVirtualClass(getClass(class(object), TRUE))
    }
    else isS4(object)
}

.RcmdrEnv <- new.env(parent=emptyenv())

# putRcmdr <- function(x, value) assign(x, value, envir=.RcmdrEnv)
# 
# getRcmdr <- function(x, mode="any") get(x, envir=.RcmdrEnv, mode=mode, inherits=FALSE)

RcmdrEnv <- function() .RcmdrEnv

putRcmdr <- function(x, value) assign(x, value, envir=RcmdrEnv())

# getRcmdr <- function(x, mode="any") get(x, envir=RcmdrEnv(), mode=mode, inherits=FALSE)

getRcmdr <- function(x, mode="any", fail=TRUE){
    if ((!fail) && (!exists(x, mode=mode, envir=RcmdrEnv(), inherits=FALSE))) return(NULL)
    get(x, envir=RcmdrEnv(), mode=mode, inherits=FALSE)
}


RcmdrTclSet <- function(name, value){
    name <- ls(unclass(getRcmdr(name))$env)
    tcl("set", name, value)
}

# functions to store or retrieve Rcmdr state information

Variables <- function(names){
    if (missing(names)) getRcmdr("variables")
    else putRcmdr("variables", names)
}

Numeric <- function(names){
    if (missing(names)) getRcmdr("numeric")
    else putRcmdr("numeric", names)
}

Factors <- function(names){
    if (missing(names)) getRcmdr("factors")
    else putRcmdr("factors", names)
}

TwoLevelFactors <- function(names){
    if (missing(names)) getRcmdr("twoLevelFactors")
    else putRcmdr("twoLevelFactors", names)
}

# The following two functions were modified by Erich Neuwrith
#  and subsequently by John Fox (23 July 07)
#  and Milan Bouchet-Valat (27 March 14)

ActiveDataSet <- function(name){
    if (missing(name)) {
        temp <- getRcmdr(".activeDataSet")
        if (is.null(temp))
            return(NULL)
        else
            if (!exists(temp) || !is.data.frame(get(temp,envir=.GlobalEnv))) {
                Message(sprintf(gettextRcmdr("the dataset %s is no longer available"),
                    temp), type="error")
                putRcmdr(".activeDataSet", NULL)
                Variables(NULL)
                Numeric(NULL)
                Factors(NULL)
                TwoLevelFactors(NULL)
                RcmdrTclSet("dataSetName", gettextRcmdr("<No active dataset>"))
                putRcmdr(".activeModel", NULL)
                RcmdrTclSet("modelName", gettextRcmdr("<No active model>"))
                tkconfigure(getRcmdr("dataSetLabel"), foreground="red") 
                tkconfigure(getRcmdr("modelLabel"), foreground="red") 
                activateMenus()
                if (getRcmdr("suppress.menus") && RExcelSupported()) return(NULL)
            }
        return(temp)
    }
    else {
        putRcmdr(".activeDataSet", name)

      if(!is.null(name)) {
        Variables(listVariables(name))
        Numeric(listNumeric(name))
        Factors(listFactors(name))
        TwoLevelFactors(listTwoLevelFactors(name))
        open.showData.windows <- getRcmdr("open.showData.windows")
        if (!is.null(open.showData.windows) && name %in% names(open.showData.windows)){
          ID <- open.showData.windows[[name]]$ID
          posn <- as.numeric(c(tclvalue(.Tcl(paste("winfo x", ID))),
                       tclvalue(.Tcl(paste("winfo y", ID)))))
          posn <- paste("+", paste(posn, collapse = "+"), sep = "")
          tkdestroy(open.showData.windows[[name]])
          suppress <- if(getRcmdr("suppress.X11.warnings")) ", suppress.X11.warnings=FALSE" else ""
          view.height <- max(as.numeric(getRcmdr("output.height")) + as.numeric(getRcmdr("log.height")), 10)
          command <- paste("showData(", name, ", placement='", posn, "', font=getRcmdr('logFont'), maxwidth=",
                           getRcmdr("log.width"), ", maxheight=", view.height, suppress, ")", sep="")
          window <- justDoIt(command)
          open.showData.windows[[ActiveDataSet()]] <- window
          putRcmdr("open.showData.windows", open.showData.windows)
        }
        
      }
        else {
            Variables(NULL)
            Numeric(NULL)
            Factors(NULL)
            TwoLevelFactors(NULL)
            RcmdrTclSet("dataSetName", gettextRcmdr("<No active dataset>"))
            putRcmdr(".activeModel", NULL)
            RcmdrTclSet("modelName", gettextRcmdr("<No active model>"))
            tkconfigure(getRcmdr("dataSetLabel"), foreground="red") 
            tkconfigure(getRcmdr("modelLabel"), foreground="red") 
            activateMenus()
            if (getRcmdr("suppress.menus") && RExcelSupported()) return(NULL)
        }
    }
}

ActiveModel <- function(name){
    if (missing(name)) {
        temp <- getRcmdr(".activeModel")
        if (is.null(temp))
            return(NULL)
        else
            if (!exists(temp) || !is.model(get(temp,envir=.GlobalEnv))) {
                Message(sprintf(gettextRcmdr("the model %s is no longer available"),
                    temp), type="error")
                putRcmdr(".activeModel", NULL)
                RcmdrTclSet("modelName", gettextRcmdr("<No active model>"))
                tkconfigure(getRcmdr("modelLabel"), foreground="red")
                activateMenus()
                return(NULL)
            }
        else return(temp)
    }
    else putRcmdr(".activeModel", name)
}

GrabFocus <- function(value){
    if (missing(value)) getRcmdr("grab.focus")
    else putRcmdr("grab.focus", value)
}

UpdateModelNumber <- function(increment=1){
    modelNumber <- getRcmdr("modelNumber")
    modelNumber <- modelNumber + increment
    if (modelNumber < 1) modelNumber <- 1 # sanity check
    putRcmdr("modelNumber", modelNumber)
}

CommanderWindow <- function() getRcmdr("commanderWindow")

LogWindow <- function() getRcmdr("logWindow")

RmdWindow <- function() getRcmdr("RmdWindow")

RnwWindow <- function() getRcmdr("RnwWindow")

OutputWindow <- function() getRcmdr("outputWindow")

MessagesWindow <- function() getRcmdr("messagesWindow")

# some predicates for the menu system

activeDataSetP <- function() !is.null(ActiveDataSet())

dataSetsP <- function(n=1){
    datasets <- listDataSets()
    (!is.null(datasets)) && length(datasets) >= n
}

numericP <- function(n=1) activeDataSetP() && length(listNumeric()) >= n

factorsP <- function(n=1) activeDataSetP() && length(listFactors()) >= n

twoLevelFactorsP <- function(n=1) activeDataSetP() && length(listTwoLevelFactors()) >= n

modelsP <- function(n=1) activeDataSetP() && length(listAllModels()) >= n

activeModelP <- function() !is.null(ActiveModel())

lmP <- function() activeModelP() && any(class(get(ActiveModel()))[1] == c('lm', 'aov'))

glmP <- function() activeModelP() && class(get(ActiveModel()))[1] == 'glm'

aicP <- function() activeModelP() && exists.method("extractAIC", get(ActiveModel()))

polrP <- function() activeModelP() && class(get(ActiveModel()))[1] == 'polr'

multinomP <- function() activeModelP() && class(get(ActiveModel()))[1] == 'multinom'

hclustSolutionsP <- function() length(listHclustSolutions()) > 0

MacOSXP <- function(release) {
    sys <- Sys.info()
    OSX <- !is.null(sys) && length(grep("[Dd]arwin", sys["sysname"]) > 0)
    if (missing(release)) OSX
    else (OSX && release <= sys["release"])
}

packageAvailable <- function(name) 0 != length(find.package(name, quiet=TRUE))

rglLoaded <- function() 0 != length(grep("^rgl", loadedNamespaces()))

activateMenus <- function(){
    if (getRcmdr("suppress.menus")) return()
    for (item in getRcmdr("Menus")){
        if (item$activation()) .Tcl(paste(item$ID, " entryconfigure ", item$position - 1," -state normal", sep=""))
        else .Tcl(paste(item$ID, " entryconfigure ", item$position - 1," -state disabled", sep=""))
    }
}

# for internationalization

gettextRcmdr <- function(...) gettext(..., domain="R-Rcmdr")

gettextMenus <- function(...){
    text <- gettextRcmdr(...)
    plugins <- getOption("Rcmdr")$plugins
    if (is.null(plugins)) return(text)
    plugins <- paste("R-", plugins, sep="")
    for (plugin in plugins){
        text <- gettext(text, domain=plugin)
    }
    text
}

English <- function() {
    env <- Sys.getenv()
    names(env) <- toupper(names(env))
    LANG <- env["LANGUAGE"]
    LC_CTYPE <- Sys.getlocale("LC_CTYPE")
    if (!is.na(LANG)) length(grep("^en", LANG, ignore.case=TRUE)) > 0
    else LC_CTYPE == "C" || length(grep("^en", LC_CTYPE, ignore.case=TRUE)) > 0
}

# to replace tkmessageBox on non-English Windows systems,
#  to allow for translation of button text

RcmdrTkmessageBox <- function(message, icon=c("info", "question", "warning",
    "error"), type=c("okcancel", "yesno", "ok"), default, title="") {
    if ( (English()) || (!WindowsP()) ){
        if (missing(default)){
            default <- switch(type,
                okcancel="ok",
                yesno="yes",
                ok="ok")}
        return(tkmessageBox(message=message, icon=icon, type=type,
            default=default, title=title))
    }
    icon <- match.arg(icon)
    type <- match.arg(type)
    initializeDialog(messageBox, title=title)
    messageFrame <- tkframe(messageBox, borderwidth=5)
    buttonFrame <- tkframe(messageBox,  borderwidth=5)
    if (icon != "question") tkbell()
    result <- tclVar()
    iconColor <- switch(icon, info=getRcmdr("title.color"), question=getRcmdr("title.color"), warning="black",
        error="red")
    onOK <- function() {
        if (GrabFocus()) tkgrab.release(messageBox)
        tkdestroy(messageBox)
        tkfocus(CommanderWindow())
        tclvalue(result) <- "ok"
    }
    OKbutton <- buttonRcmdr(buttonFrame, text=gettextRcmdr("OK"),
        foreground="darkgreen", width="12", command=onOK, borderwidth=3,
        default=if (missing(default)) "active"
        else if (default == "ok") "active" else "normal")
    onCancel <- function() {
        if (GrabFocus()) tkgrab.release(messageBox)
        tkdestroy(messageBox)
        tkfocus(CommanderWindow())
        tclvalue(result) <- "cancel"
    }
    cancelButton <- buttonRcmdr(buttonFrame, text=gettextRcmdr("Cancel"),
        foreground="red", width="12", command=onCancel, borderwidth=3,
        default=if (missing(default)) "normal"
        else if (default == "cancel") "active" else "normal")
    onYes <- function() {
        if (GrabFocus()) tkgrab.release(messageBox)
        tkdestroy(messageBox)
        tkfocus(CommanderWindow())
        tclvalue(result) <- "yes"
    }
    yesButton <- buttonRcmdr(buttonFrame, text=gettextRcmdr("Yes"),
        foreground="darkgreen", width="12", command=onYes, borderwidth=3,
        default=if (missing(default)) "active"
        else if (default == "yes") "active" else "normal")
    onNo <- function() {
        if (GrabFocus()) tkgrab.release(messageBox)
        tkdestroy(messageBox)
        tkfocus(CommanderWindow())
        tclvalue(result) <- "no"
    }
    noButton <- buttonRcmdr(buttonFrame, text=gettextRcmdr("No"),
        foreground="red", width="12", command=onNo, borderwidth=3,
        default=if (missing(default)) "normal"
        else if (default == "no") "active" else "normal")
    ## FIXME -- left in old style
    tkgrid(tklabel(messageFrame, bitmap=icon, fg=iconColor),
        tklabel(messageFrame, text="    "),
        tklabel(messageFrame, text=message))
    tkgrid(messageFrame)
    switch(type,
        okcancel = {
            tkgrid(OKbutton, labelRcmdr(buttonFrame, text="    "), cancelButton)
            if (missing(default) || default == "ok") tkbind(messageBox, "<Return>",
                onOK)
            else if (default == "cancel") tkbind(messageBox, "<Return>", onCancel)
        },
        yesno =  {
            tkgrid(yesButton, labelRcmdr(buttonFrame, text="    "), noButton)
            if (missing(default) || default == "yes") tkbind(messageBox, "<Return>",
                onYes)
            else if (default == "no") tkbind(messageBox, "<Return>", onNo)
        },
        ok = {
            tkgrid(OKbutton)
            if (missing(default) || default == "ok") tkbind(messageBox, "<Return>",
                onOK)
        }
    )
    tkgrid(buttonFrame)
    dialogSuffix(messageBox, focus=messageBox, bindReturn=FALSE, force.wait=TRUE)
    result
}

# The following function was contributed by Matthieu Lesnoff (added 20 July 06)

trim.col.na <- function(dat){
    # Remove variables with only missing values (occurs sometimes with modified Excel file)
    colsup <- NULL
    for (i in 1:ncol(dat))
    {
        if (length(dat[is.na(dat[,i])==T,i]) ==length(dat[,i]))
            colsup <- c(colsup,i)
    }
    if (length(colsup) > 0)
        dat <- dat[,-colsup]
    dat
}

# check whether packages are available

packagesAvailable <- function(packages){
    sapply(sapply(packages, find.package, quiet=TRUE),
        function(x) length(x) != 0)
}

# insert a row (or rows) in a matrix or data frame

insertRows <- function(object1, object2, where=NULL, ...){
    if (ncol(object1) != ncol(object2))
        stop(gettextRcmdr("objects have different numbers of columns"))
    if (!(TRUE == all.equal(colnames(object1), colnames(object2))))
        stop(gettextRcmdr("objects have different column names"))
    n <- nrow(object1)
    if (is.null(where) || where >= n) rbind(object1, object2)
    else if (where < 1) rbind(object2, object1)
    else rbind(object1[1:floor(where),], object2,
        object1[(floor(where) + 1):n,])
}

# functions for handling Rcmdr plug-in packages

# the following function based on a suggestion by Brian Ripley

listPlugins <- function(loaded=FALSE){
    plugins <- unlist(lapply(.libPaths(),
        function(x) Sys.glob(file.path(x, "*/etc/menus.txt"))))
    plugins <- sub(".*/([^/]*)/etc/menus.txt", "\\1", plugins)
    if (loaded) plugins else sort(setdiff(plugins, .packages()))
}


loadPlugins <- function(){
    plugins <- listPlugins()
    initializeDialog(title=gettextRcmdr("Load Plug-ins"))
    packagesBox <- variableListBox(top, plugins, title=gettextRcmdr("Plug-ins (pick one or more)"),
        selectmode="multiple", listHeight=10)
    onOK <- function(){
        plugins <- getSelection(packagesBox)
        closeDialog(top)
        if (length(plugins) == 0){
            errorCondition(recall=loadPlugins, message=gettextRcmdr("You must select at least one plug-in."))
            return()
        }
        opts <- options("Rcmdr")
        opts$Rcmdr$plugins <- c(plugins, opts$Rcmdr$plugins)
        options(opts)
        for (plugin in plugins) {
            command <- paste('library("', plugin, '", character.only=TRUE)', sep="")
            justDoIt(command)
        }
        Message(paste(gettextRcmdr("Plug-ins loaded:"), paste(plugins, collapse=", ")), type="note")
        response <- tkmessageBox(message=paste(gettextRcmdr(
            "The plug-in(s) will not be available until the Commander is restarted.\nRestart now?")),
            icon="question", type="yesno")
        if (tclvalue(response) == "yes") {
            putRcmdr("autoRestart", TRUE)
            closeCommander(ask=FALSE)
            Commander()
        }
    }
    OKCancelHelp(helpSubject="Plugins")
    tkgrid(getFrame(packagesBox), sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}

# the following two functions contributed by Erich Neuwirth (added 22 July 07)

whitespaceonly <- function(str) sub('[[:space:]]+$', '', str) == ''

is.model <- function(object) {
    any(class(object) %in% getRcmdr("modelClasses"))
}

# the following lines, adding support for ttk widgets, adapted from code by Brian Ripley
if (!(as.character(tcl("info", "tclversion")) >= "8.5" && getRversion() >= "2.7.0")){
    buttonRcmdr <- tkbutton
    labelRcmdr <- tklabel
    ttkentry <- function(parent, ...) tkentry(parent, ...)
    ttkframe <- tkframe
    ttkradiobutton <- tkradiobutton
    ttkscrollbar <- function(...) tkscrollbar(..., repeatinterval=5)
} else {
    buttonRcmdr <- function(..., borderwidth, fg, foreground, relief) ttkbutton(...)
    labelRcmdr <- function(..., fg)
        if(missing(fg)) ttklabel(...) else ttklabel(..., foreground=fg)
}

# Label looking like that of a TtkLabelFrame
titleLabel <- function(...) labelRcmdr(..., font="RcmdrTitleFont", fg=getRcmdr("title.color"))

# the following function alters the default behaviour of tclvalue() by trimming leading and trailing blanks

tclvalue <- function(x) trim.blanks(tcltk::tclvalue(x))

# the following function splits a character string at blanks and commas according to width

splitCmd <- function(cmd, width=getOption("width") - 4, at="[ ,]"){
  if (length(grep("\n", cmd)) >0 ){
    cmds <- strsplit(cmd, "\n")[[1]]
    allcmds <- character(length(cmds))
    for (i in 1:length(cmds))
      allcmds[i] <- splitCmd(cmds[i], width=width, at=at)
    return(paste(allcmds, collapse="\n"))
  }
  if (nchar(cmd) <= width) return(cmd)
  where <- gregexpr(at, cmd)[[1]]
  if (where[1] < 0) return(cmd)
  singleQuotes <- gregexpr("'", cmd)[[1]]
  doubleQuotes <- gregexpr('"', cmd)[[1]]
  comment <- regexpr("#", cmd)
  if (singleQuotes[1] > 0 && (singleQuotes[1] < doubleQuotes[1] || doubleQuotes[1] < 0 ) && (singleQuotes[1] < comment[1] || comment[1] < 0 )){
    nquotes <- length(singleQuotes)
    if (nquotes < 2) stop("unbalanced quotes")
    for(i in seq(nquotes/2))
      where[(where > singleQuotes[2 * i - 1]) & (where < singleQuotes[2 * i])] <- NA
    where <- na.omit(where)
  }  
  else if (doubleQuotes[1] > 0 && (doubleQuotes[1] < singleQuotes[1] || singleQuotes[1] < 0) && (doubleQuotes[1] < comment[1] || comment[1] < 0 )){
    nquotes <- length(doubleQuotes)
    if (nquotes < 2) stop("unbalanced quotes")
    for(i in seq(nquotes/2))
      where[(where > doubleQuotes[2 * i - 1]) & (where < doubleQuotes[2 * i])] <- NA
    where <- na.omit(where)
  }
  else if (comment > 0){
    where[where > comment] <- NA
    where <- na.omit(where)
  }
  if (length(where) == 0) return(cmd)
  where2 <- where[where <= width]
  where2 <- if (length(where2) == 0) where[1]
  else where2[length(where2)]
  paste(substr(cmd, 1, where2), "\n  ", 
        Recall(substr(cmd, where2 + 1, nchar(cmd)), width, at), sep="")
} 

# the following function sorts names containing numerals "more naturally" than does sort()

sortVarNames <- function(x){
    sort.helper <- function(x){
        prefix <- strsplit(x, "[0-9]+")
        prefix <- sapply(prefix, "[", 1)
        prefix[is.na(prefix)] <- ""
        suffix <- strsplit(x, "[^0-9]+")
        suffix <- as.numeric(sapply(suffix, "[", 2))
        suffix[is.na(suffix)] <- -Inf
        remainder <- sub("[^0-9]+", "", x)
        remainder <- sub("[0-9]+", "", remainder)
        if (all (remainder == "")) list(prefix, suffix)
        else c(list(prefix, suffix), Recall(remainder))
    }
    ord <- do.call("order", sort.helper(x))
    x[ord]
}

# to load packages

Library <- function(package, pos=length(search()), rmd=TRUE){
    dependencies <- tools::package_dependencies(package, db=getRcmdr("installed.packages"), which="Depends")
    loaded <- search()
    loaded <- loaded[grep("^package:", loaded)]
    loaded <- sub("^package:", "", loaded)
    if (!getRcmdr("suppress.X11.warnings")){
        messages.connection <- file(open="w+")
        sink(messages.connection, type="message")
        on.exit({
            sink(type="message")
            close(messages.connection)
        })
    }
    if (!(package %in% loaded)){
        for (pkg in dependencies[[package]]){
            Library(pkg, pos=pos, rmd=rmd)
        }
        command <- paste("library(", package, ", pos=", pos, ")", sep="")
        logger(command, rmd=rmd)
        result <- try(eval(parse(text=command), envir=.GlobalEnv), silent=TRUE)
        if (class(result)[1] ==  "try-error"){
            Message(message=paste(strsplit(result, ":")[[1]][2]), type="error")
            tkfocus(CommanderWindow())
            return("error")
        }
        return(package)
    }
    else return(invisible(NULL))
}

# Library <- function(package, pos=4, rmd=TRUE){
#     loaded <- search()
#     loaded <- loaded[grep("^package:", loaded)]
#     loaded <- sub("^package:", "", loaded)
#     if (!getRcmdr("suppress.X11.warnings")){
#         messages.connection <- file(open="w+")
#         sink(messages.connection, type="message")
#         on.exit({
#             sink(type="message")
#             close(messages.connection)
#         })
#     }
#     if (!(package %in% loaded)){
#         command <- paste("library(", package, ", pos=", pos, ")", sep="")
#         logger(command, rmd=rmd)
#         result <- try(eval(parse(text=command), envir=.GlobalEnv), silent=TRUE)
#         if (class(result)[1] ==  "try-error"){
#             Message(message=paste(strsplit(result, ":")[[1]][2]), type="error")
#             tkfocus(CommanderWindow())
#             return("error")
#         }
#         return(package)
#     }
#     else return(invisible(NULL))
# }

# start help system

startHelp <- function(){
    Sys.sleep(2)
    help.start()
}

# dialog memory support

putDialog <- function (dialog, values=NULL, resettable=TRUE){
    if (resettable){
        dialog.values <- getRcmdr("dialog.values")
        dialog.values[[dialog]] <- values
        putRcmdr("dialog.values", dialog.values)
    }
    else{
        dialog.values <- getRcmdr("dialog.values.noreset")
        dialog.values[[dialog]] <- values
        putRcmdr("dialog.values.noreset", dialog.values)
    }
}

getDialog <- function(dialog, defaults=NULL){
    values <- getRcmdr("dialog.values.noreset")[[dialog]]
    if (getRcmdr("retain.selections") && !is.null(values)) return(values)
    values <- getRcmdr("dialog.values")[[dialog]]
    if (!getRcmdr("retain.selections") || is.null(values)) return(defaults)
    else return (values)
}

varPosn <- function(variables, 
    type=c("all", "factor", "numeric", "nonfactor", "twoLevelFactor"), vars=NULL){
    if (is.null(variables)) return(NULL)
    type <- match.arg(type)
    if (is.null(vars)) vars <- switch(type,
        all = Variables(),
        factor = Factors(),
        numeric = Numeric(),
        nonfactor = setdiff(Variables(), Factors()),
        twoLevelFactor = TwoLevelFactors()
    )
    if (any(!variables %in% vars)) NULL
    else apply(outer(variables, vars, "=="), 1, which) - 1
}

flushDialogMemory <- function(what){
    if (missing(what)) putRcmdr("dialog.values", list())
    else{
        dialog.values <- getRcmdr("dialog.values")
        dialog.values.noreset <- getRcmdr("dialog.values.noreset")
        for (dialog in what){
            dialog.values[dialog] <- NULL
            dialog.values.noreset[dialog] <- NULL
        }
        putRcmdr("dialog.values", dialog.values)
        putRcmdr("dialog.values.noreset", dialog.values.noreset)
    }
}

# for assignments to the global environment

gassign <- function(x, value){
    if (!(is.valid.name(x))) stop("argument x not a valid R name")
    G <- .GlobalEnv
    assign(x, value, envir=G)
}


tkfocus <- function(...) tcl("focus", ...)

tkspinbox <- function(parent, ...) tkwidget(parent, "spinbox", ...)

# the following two functions adapted from Milan Bouchet-Valat

WindowsP <- function() {
    .Platform$OS.type == "windows"
}

X11P <- function(){
    .Platform$GUI == "X11"
}

# the following functions to support R Markdown

trimTrailingNewLines <- function(string){
  repeat{
    where <- regexpr("\n\n[ ]*$", string)
    if (where == -1) break
    string <- paste0(substr(string, 1, where - 1), substr(string, where + 2, nchar(string)))
  }
  paste0(string, "\n")
}

suppressMarkdown <- function(command){
    attr(command, "suppressRmd") <- TRUE
    command
}

beginRmdBlock <- function(){
    .rmd <- RmdWindow()
    last2 <- tclvalue(tkget(.rmd, "end -2 chars", "end"))
    if (last2 != "\n\n") tkinsert(.rmd, "end", "\n")
    tkinsert(.rmd, "end", "\n")
    if (getRcmdr("rgl.command") && getRcmdr("use.rgl")) tkinsert(.rmd, "end", "```{r, webgl=TRUE}\n")
      else tkinsert(.rmd, "end", "```{r}\n")
}

endRmdBlock <- function(){
    .rmd <- RmdWindow()
    rmd <- tclvalue(tkget(.rmd, "1.0", "end"))
    rmd <- paste(substring(rmd, 1, nchar(rmd) - 1), "```\n", sep="")
    rmd <- trimHangingEndRmdBlock(rmd)
    rmd <- trimTrailingNewLines(rmd)
    tkdelete(.rmd, "1.0", "end")
    tkinsert(.rmd, "end", rmd)
    tkyview.moveto(.rmd, 1)
}

removeNullRmdBlocks <- function(){
    .rmd <- RmdWindow()
    rmd <- tclvalue(tkget(.rmd, "1.0", "end"))
    rmd <- gsub("\n+$", "\n", rmd)
    rmd <- gsub("```\\{r\\}\n$", "", rmd)
    rmd <- gsub("```\\{r\\}\n```\n$", "", rmd)
    tkdelete(.rmd, "1.0", "end")
    tkinsert(.rmd, "end", rmd)
    tkyview.moveto(.rmd, 1)
}

removeStrayRmdBlocks <- function(){
    .rmd <- RmdWindow()
    rmd <- tclvalue(tkget(.rmd, "1.0", "end"))
    rmd <- strsplit(rmd, "\\n")[[1]]
    starts <- grep("^```\\{r.*\\}$", rmd)
    ends  <- grep("^```$", rmd)
    n.ends <- length(ends)
    j <- 1
    if (length(starts) > 1){
        for (i in 1:(length(starts) - 1)){
            if (j > n.ends || ends[j] > starts[i + 1]) {
                rmd[starts[i]] <- ""
            }
            else {
                j <- j + 1
                next
            }
        }
    }
    else return()
    rmd <- paste(rmd, collapse="\n")
    tkdelete(.rmd, "1.0", "end")
    tkinsert(.rmd, "end", rmd)
    tkyview.moveto(.rmd, 1)
}

enterMarkdown <- function(command){
    if (!getRcmdr("use.markdown")) return()
    .rmd <- RmdWindow()
    command <- splitCmd(command)
    beginRmdBlock()
    tkinsert(.rmd, "end", paste(command, "\n", sep=""))
    tkyview.moveto(.rmd, 1)
    putRcmdr("markdown.output", TRUE)
    endRmdBlock()
    command
}

trimHangingEndRmdBlock <- function(string){
    loc.ends <- gregexpr("```\n", string)[[1]]
    n.ends <- length(loc.ends)
    if (n.ends > 1){
        substr <- substring(string, loc.ends[n.ends - 1], loc.ends[n.ends])
        if (!grepl("```\\{r\\}|```\\{r, webgl=TRUE\\}", substr)){
            string <- cutstring(string, loc.ends[n.ends], loc.ends[n.ends] + 3)
        }
    }
    string
}

removeLastRmdBlock <- function(){
    .rmd <- RmdWindow()    
    rmd <- tclvalue(tkget(.rmd, "1.0", "end"))
    start <- gregexpr("```\\{r\\}\n|```\\{r, webgl=TRUE\\}\n", rmd)
    if (start[[1]][1] > 0){
        start <- start[[1]]
        start <- start[length(start)]
        tail <- substring(rmd, start, nchar(rmd))
        end <- gregexpr("```\n", tail)
        end <- if (end[[1]][1] > 0) end[[1]][1] + 3 else nchar(tail)
        rmd <- cutstring(rmd, start, start + end)
        rmd <- trimTrailingNewLines(rmd)
        tkdelete(.rmd, "1.0", "end")
        tkinsert(.rmd, "end", rmd)
        tkyview.moveto(.rmd, 1)
    }
}

removeRglRmdBlocks <- function(string){
  repeat{
    match <- regexpr("```\\{r, webgl=TRUE\\}\n", string)
    if (match == -1) return(trimTrailingNewLines(string))
    substring <- cutstring(string, end=match)
    match.end <- regexpr("```\n", substring)
    string <- cutstring(string, match, match + match.end + 3)
  }
}

cutstring <- function(x, start=1, end=nchar(x)){
    one <- if (start > 1) substr(x, 1, start - 1) else ""
    two <- if (end < nchar(x)) substr(x, end + 1, nchar(x)) else ""
    paste0(one, two)
}

MarkdownP <- function(){
    getRcmdr("log.commands") && getRcmdr("use.markdown")
}

compileRmd <- function() {
    ChooseOutputFormat <- function(){
        initializeDialog(title=gettextRcmdr("Select Output Format"))
        format <- getRcmdr("rmd.output.format")
        putRcmdr("abort.compile.rmd", TRUE)
        hasLatex <- getRcmdr("capabilities")$pdflatex
        radioButtons(name="formatButtons", 
            buttons=c("html", if (hasLatex) "pdf", "docx"), 
            initialValue=format,
            labels=c(gettextRcmdr(".html (web page)"), 
                if (hasLatex) gettextRcmdr(".pdf (PDF file)"), gettextRcmdr(".docx (Word file)")))
        onOK <- function(){
            putRcmdr("abort.compile.rmd", FALSE)
            format <- tclvalue(formatButtonsVariable)
            putRcmdr("rmd.output.format", format)
            closeDialog()
        }
        OKCancelHelp()
        tkgrid(formatButtonsFrame, sticky="w")
        dialogSuffix(force.wait=TRUE, grid.buttons=TRUE)
    }
    .RmdFile <- getRcmdr("RmdFileName")
    rmdDir <- dirname(.RmdFile)
    saveDir <- setwd(rmdDir)
    on.exit(setwd(saveDir))
    fig.files <- list.files("./figure")
    fig.files <- fig.files[grep("^unnamed-chunk-[0-9]*\\..*$", fig.files)]
    if (length(fig.files) != 0) {
        response <- tkmessageBox(message = gettextRcmdr("Delete previously created R Markdown\ngraphics files (recommended)?"),
            icon = "question", type = "okcancel", default = "ok")
        if (tclvalue(response) == "ok") unlink(paste("./figure/", fig.files, sep=""))
    }
    removeStrayRmdBlocks()
    lines <- tclvalue(tkget(RmdWindow(), "1.0", "end"))
    lines <- sub("date: \"AUTOMATIC\"", paste("date: \"", as.character(Sys.time()), "\"", sep=""), lines)
    .filename <- sub("\\.Rmd$", "", trim.blanks(.RmdFile))
    writeLines(lines, .RmdFile)
    if (getRcmdr("capabilities")$pandoc){
        ChooseOutputFormat()
        if (getRcmdr("abort.compile.rmd")){
            putRcmdr("abort.compile.rmd", NULL)
            return()
        }
        else putRcmdr("abort.compile.rmd", NULL)
        format <- getRcmdr("rmd.output.format")
        switch(format,
            html = {
                rmarkdown::render(.RmdFile, rmarkdown::html_document())
                .html.file <- paste(.filename, ".html", sep="")
                .html.file.location <- paste("file:///", normalizePath(.html.file), sep="")
                Message(paste(gettextRcmdr("HTML file written to:"), normalizePath(.html.file)), type="note")
                browseURL(.html.file.location)
            },
            pdf = {
                lines <- removeRglRmdBlocks(lines)
                writeLines(lines, .RmdFile)
                rmarkdown::render(.RmdFile, rmarkdown::pdf_document())
                .pdf.file <- paste(.filename, ".pdf", sep="")
                .pdf.file.location <- paste("file:///", normalizePath(.pdf.file), sep="")
                Message(paste(gettextRcmdr("PDF file written to:"), normalizePath(.pdf.file)), type="note")
                browseURL(.pdf.file.location)
            },
            docx = {
              lines <- removeRglRmdBlocks(lines)
              writeLines(lines, .RmdFile)
                rmarkdown::render(.RmdFile, rmarkdown::word_document())
                .docx.file <- paste(.filename, ".docx", sep="")
                Message(paste(gettextRcmdr("Word file written to:"), normalizePath(.docx.file)), type="note")
            }
        )
    }
    else{
        knitr::knit(.RmdFile, paste(.filename, ".md", sep=""), quiet=TRUE)
        .html.file <- paste(.filename, ".html", sep="")
        markdown::markdownToHTML(paste(.filename, ".md", sep=""), .html.file)
        .html.file.location <- paste("file:///", normalizePath(.html.file), sep="")
        Message(paste(gettextRcmdr("HTML file written to:"), normalizePath(.html.file)), type="note")
        browseURL(.html.file.location)
    }
}

# the following functions to support knitr

beginRnwBlock <- function(){
    .rnw <- RnwWindow()
    last2 <- tclvalue(tkget(.rnw, "end -2 chars", "end"))
    if (last2 != "\n\n") tkinsert(.rnw, "end", "\n")
    tkinsert(.rnw, "end", "\n")
    tkinsert(.rnw, "end", "\\newpage\n")
    tkinsert(.rnw, "end", "<<>>=\n")
}

endRnwBlock <- function(){
    .rnw <- RnwWindow()
    rnw <- tclvalue(tkget(.rnw, "1.0", "end"))
    rnw <- paste(substring(rnw, 1, nchar(rnw) - 1), "@\n", sep="")
    rnw <- trimHangingEndRnwBlock(rnw)
    rmw <- trimTrailingNewLines(rnw)
    tkdelete(.rnw, "1.0", "end")
    tkinsert(.rnw, "end", rnw)
    tkyview.moveto(.rnw, 1)    
}

removeNullRnwBlocks <- function(){
    .rnw <- RnwWindow()
    rnw <- tclvalue(tkget(.rnw, "1.0", "end"))
    rnw <- gsub("\n+$", "\n", rnw)
    rnw <- gsub("<<>>=\n$", "", rnw)
    rnw <- gsub("<<>>=\n@\n$", "", rnw)
    rnw <- gsub("\\\\newpage\n*$", "", rnw)
    rnw <- gsub("\\\\newpage\n*$", "", rnw)
    tkdelete(.rnw, "1.0", "end")
    tkinsert(.rnw, "end", rnw)
    tkyview.moveto(.rnw, 1)
}

removeStrayRnwBlocks <- function(){
    .rnw <- RnwWindow()
    rnw <- tclvalue(tkget(.rnw, "1.0", "end"))
    rnw <- strsplit(rnw, "\\n")[[1]]
    starts <- grep("^<<.*>>=$", rnw)
    ends  <- grep("^@$", rnw)
    n.ends <- length(ends)
    j <- 1
    if (length(starts) > 1){
        for (i in (length(starts) - 1)){
            if (j > n.ends || ends[j] > starts[i + 1]) {
                rnw[starts[i]] <- ""
            }
            else {
                j <- j + 1
                next
            }
        }
    }
    else return()
    rnw <- paste(rnw, collapse="\n")
    tkdelete(.rnw, "1.0", "end")
    tkinsert(.rnw, "end", rnw)
    tkyview.moveto(.rnw, 1)
}

enterKnitr <- function(command){
    .rnw <- RnwWindow()
    if (!getRcmdr("use.knitr")) return()
    command <- splitCmd(command)
    beginRnwBlock()
    tkinsert(.rnw, "end", paste(command, "\n", sep=""))
    tkyview.moveto(.rnw, 1)
    putRcmdr("knitr.output", TRUE)
    endRnwBlock()
    command
}

trimHangingEndRnwBlock <- function(string){
    loc.ats <- gregexpr("@\n", string)[[1]]
    n.ats <- length(loc.ats)
    if (n.ats > 1){
        substr <- substring(string, loc.ats[n.ats - 1], loc.ats[n.ats])
        if (!grepl("<<>>=", substr)){
            string <- cutstring(string, loc.ats[n.ats], loc.ats[n.ats] + 1)
        }
    }
    string
}

removeLastRnwBlock <- function(){
    .rnw <- RnwWindow()
    rnw <- tclvalue(tkget(.rnw, "1.0", "end"))
    start <- gregexpr("\\\\newpage\n<<>>=\n", rnw)
    if (start[[1]][1] > 0){
        start <- start[[1]]
        start <- start[length(start)]
        tail <- substring(rnw, start, nchar(rnw))
        end <- gregexpr("@\n", tail)
        end <- if (end[[1]][1] > 0) end[[1]][1] + 1 else nchar(tail)
        rnw <- cutstring(rnw, start, start + end)
        rnw <- trimTrailingNewLines(rnw)
        tkdelete(.rnw, "1.0", "end")
        tkinsert(.rnw, "end", rnw)
        tkyview.moveto(.rnw, 1)
    }
}

compileRnw <- function(){
    .RnwFile <- getRcmdr("RnwFileName")
    rnwDir <- dirname(.RnwFile)
    saveDir <- setwd(rnwDir)
    on.exit(setwd(saveDir))
    fig.files <- list.files("./figure")
    fig.files <- fig.files[grep("^unnamed-chunk-[0-9]*\\..*$", fig.files)]
    if (length(fig.files) != 0) {
        response <- tkmessageBox(message = gettextRcmdr("Delete previously created knitr\ngraphics files (recommended)?"),
            icon = "question", type = "okcancel", default = "ok")
        if (tclvalue(response) == "ok") unlink(paste("./figure/", fig.files, sep=""))
    }
    removeStrayRnwBlocks()
    lines <- tclvalue(tkget(RnwWindow(), "1.0", "end"))
    lines <- paste(lines, "\n\\end{document}\n")
    .filename <- sub("\\.Rnw$", "", trim.blanks(.RnwFile))
    writeLines(lines, .RnwFile)
    knitr::knit2pdf(.RnwFile)
    .pdf.file <- paste(.filename, ".pdf", sep="")
    .pdf.file.location <- paste("file:///", normalizePath(.pdf.file), sep="")
    browseURL(.pdf.file.location)
}


knitrP <- function(){
    getRcmdr("log.commands") && getRcmdr("use.knitr")
}

# editor for R Markdown and knitr documents

RcmdrEditor <- function(buffer, title="R Commander Editor", ok,
                        help=NULL, file.menu=NULL, edit.menu=NULL, context.menu=NULL, toolbar.buttons=NULL){
  contextMenu <- function(){
    contextMenu <- tkmenu(tkmenu(editor), tearoff=FALSE)
    if (!is.null(context.menu)){
      for (item in context.menu){
        tkadd(contextMenu, "command", label=gettextRcmdr(item$label), command=item$command)
      }
      tkadd(contextMenu, "separator")
    }
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
    tkpopup(contextMenu, tkwinfo("pointerx", editor), tkwinfo("pointery", editor))
  }
  onCopy <- function(){
    selection <- strsplit(tclvalue(tktag.ranges(editor, "sel")), " ")[[1]]
    if (is.na(selection[1])) return()
    text <- tclvalue(tkget(editor, selection[1], selection[2]))
    tkclipboard.clear()
    tkclipboard.append(text)
  }
  onDelete <- function(){
    selection <- strsplit(tclvalue(tktag.ranges(editor, "sel")), " ")[[1]]
    if (is.na(selection[1])) return()
    tkdelete(editor, selection[1], selection[2])
  }
  onCut <- function(){
    onCopy()
    onDelete()
  }
  onPaste <- function(){
    onDelete()
    text <- tclvalue(.Tcl("selection get -selection CLIPBOARD"))
    if (length(text) == 0) return()
    tkinsert(editor, "insert", text)
  }
  onFind <- function(){
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
      where.txt <- if (case) tksearch(editor, type, direction, "--", text, "insert", stop)
      else tksearch(editor, type, direction, "-nocase", "--", text, "insert", stop)
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
      tkfocus(editor)
      tkmark.set(editor, "insert", where.txt)
      tksee(editor, where.txt)
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
    tktag.add(editor, "sel", "1.0", "end")
    tkfocus(editor)
  }
  onClear <- function(){
    onSelectAll()
    onDelete()
  }
  onUndo <- function(){
    tcl(editor, "edit", "undo")
  }
  onRedo <- function(){
    tcl(editor, "edit", "redo")
  }
  initializeDialog(title = gettextRcmdr(title), suppress.window.resize.buttons=FALSE)
  toolbarFrame <- tkframe(top) 
  cutButton <- buttonRcmdr(toolbarFrame, image="::image::cutIcon", command=onCut)
  copyButton <- buttonRcmdr(toolbarFrame, image="::image::copyIcon", command=onCopy)
  pasteButton <- buttonRcmdr(toolbarFrame, image="::image::pasteIcon", command=onPaste)
  deleteButton <- buttonRcmdr(toolbarFrame, image="::image::deleteIcon", command=onDelete)
  undoButton <- buttonRcmdr(toolbarFrame, image="::image::undoIcon", command=onUndo)
  redoButton <- buttonRcmdr(toolbarFrame, image="::image::redoIcon", command=onRedo)
  findButton <- buttonRcmdr(toolbarFrame, image="::image::findIcon", command=onFind)
  if (!is.null(toolbar.buttons)){
    for (i in 1:length(toolbar.buttons)){
      tool <- toolbar.buttons[[i]]
      assign(paste("var", i, sep=""), tclVar(gettextRcmdr(tool$label)))
      assign(paste("button", i, sep=""), buttonRcmdr(toolbarFrame, textvariable=eval(parse(text=paste("var", i, sep=""))), 
                                                     borderwidth="2", command=tool$command, image=tool$image, compound="left"))
    }
  }
  editorFrame <- tkframe(top)
  screenheight <- as.numeric(.Tcl(paste("winfo screenheight", top$ID)))
  char.size <- as.numeric(.Tcl(paste("font metrics", getRcmdr('logFont'))))[6]
  width <- as.numeric(tkcget(LogWindow(), "-width")) + 5
  height <- max(floor(screenheight/(2.5*char.size)), 25)   
  editor <- tktext(editorFrame, bg = "white", font = getRcmdr("logFont"), 
                   height = height, width = width, wrap = "none", undo=TRUE)
  putRcmdr("editor.text", editor)
  editorXscroll <- ttkscrollbar(editorFrame, orient = "horizontal", 
                                command = function(...) tkxview(editor, ...))
  editorYscroll <- ttkscrollbar(editorFrame, command = function(...) tkyview(editor, 
                                                                             ...))
  tkconfigure(editor, xscrollcommand = function(...) tkset(editorXscroll, 
                                                           ...))
  tkconfigure(editor, yscrollcommand = function(...) tkset(editorYscroll, 
                                                           ...))
  tkinsert(editor, "1.0", buffer)
#  putRcmdr("buffer", NULL)
  onOK <- function(){
#    putRcmdr("buffer", tclvalue(tkget(editor, "1.0", "end")))
    ok()
    closeDialog()
  }
  .exit <- function(){
    answer <- RcmdrTkmessageBox("Discard edits?", icon="question", type="yesno")
    if (as.character(answer) == "no") "abort" else ""
  }
  OKCancelHelp(helpSubject = "ScriptEditor")
  editorMenu <- tkmenu(top)
  tkconfigure(top, menu = editorMenu)
  fileMenu <- tkmenu(editorMenu, tearoff=FALSE)
  if (!is.null(file.menu)){
    for (item in file.menu){
      tkadd(fileMenu, "command", label=gettextRcmdr(item$label), command=item$command)
    }
    tkadd(fileMenu, "separator")
  }
  tkadd(fileMenu, "command", label=gettextRcmdr("Exit editor"), command=onOK)
  tkadd(fileMenu, "command", label=gettextRcmdr("Cancel"), command=onCancel)
  tkadd(editorMenu, "cascade", label=gettextRcmdr("File"), menu=fileMenu)
  editMenu <- tkmenu(editorMenu, tearoff=FALSE)
  if (!is.null(edit.menu)){
    for (item in edit.menu){
      tkadd(editMenu, "command", label=gettextRcmdr(item$label), command=item$command)
    }
    tkadd(editMenu, "separator")
  }
  tkadd(editMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
  tkadd(editMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
  tkadd(editMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
  tkadd(editMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
  tkadd(editMenu, "separator")
  tkadd(editMenu, "command", label=gettextRcmdr("Find..."), command=onFind)
  tkadd(editMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
  tkadd(editMenu, "separator")
  tkadd(editMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
  tkadd(editMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
  tkadd(editMenu, "separator")
  tkadd(editMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)
  tkadd(editorMenu, "cascade", label=gettextRcmdr("Edit"), menu=editMenu)
  helpMenu <- tkmenu(editorMenu, tearoff=FALSE)
  onEditorHelp <- function() print(help("ScriptEditor", package="Rcmdr"))
  tkadd(helpMenu, "command", label=gettextRcmdr("Editor help"), command=onEditorHelp)
  if (!is.null(help)){
    tkadd(helpMenu, "command", label=gettextRcmdr(help$label), command=help$command)
  }
  tkadd(editorMenu, "cascade", label=gettextRcmdr("Help"), menu=helpMenu)
  tkgrid(editor, editorYscroll, sticky = "nsew")
  tkgrid(editorXscroll)
  if (!is.null(toolbar.buttons)){
    for (i in 1:length(toolbar.buttons)){
      tkgrid(eval(parse(text=paste("button", i, sep=""))), sticky="w", row=0, column=i - 1,
             padx=c(3, 3), pady=c(0, 8))
    }
  }
  else i <- 0
  tkgrid(cutButton, sticky="w", row=0, column=i, padx=c(3, 3), pady=c(0, 8))
  tkgrid(copyButton, sticky="w", row=0, column=i + 1, padx=c(3, 3), pady=c(0, 8))
  tkgrid(pasteButton, sticky="w", row=0, column=i + 2, padx=c(3, 3), pady=c(0, 8))
  tkgrid(deleteButton, sticky="w", row=0, column=i + 3, padx=c(3, 3), pady=c(0, 8))
  tkgrid(undoButton, sticky="w", row=0, column=i + 4, padx=c(3, 3), pady=c(0, 8))
  tkgrid(redoButton, sticky="w", row=0, column=i + 5, padx=c(3, 3), pady=c(0, 8))
  tkgrid(findButton, sticky="w", row=0, column=i + 6, padx=c(3, 3), pady=c(0, 8))
  tkgrid(toolbarFrame, sticky="w")
  tk2tip(cutButton, gettextRcmdr("Cut"))
  tk2tip(copyButton, gettextRcmdr("Copy"))
  tk2tip(pasteButton, gettextRcmdr("Paste"))
  tk2tip(deleteButton, gettextRcmdr("Delete"))
  tk2tip(undoButton, gettextRcmdr("Undo"))
  tk2tip(redoButton, gettextRcmdr("Redo"))
  tk2tip(findButton, gettextRcmdr("Find"))
  tkgrid(editorFrame, sticky = "nsew")
  tkgrid.configure(editorXscroll, sticky = "ew")
  tkgrid.configure(editorYscroll, sticky = "ns")
  tkgrid.configure(editor, sticky = "nsew")
  tkgrid.configure(editorFrame, sticky = "nsew")
  tkgrid(buttonsFrame, sticky = "ew")
  tkbind(top, "<ButtonPress-3>", contextMenu)
  tkbind(top, "<Control-x>", onCut)
  tkbind(top, "<Control-X>", onCut)
  tkbind(top, "<Control-c>", onCopy)
  tkbind(top, "<Control-C>", onCopy)
  tkbind(top, "<Control-f>", onFind)
  tkbind(top, "<Control-F>", onFind)
  tkbind(top, "<F3>", onFind)
  tkbind(top, "<Control-a>", onSelectAll)
  tkbind(top, "<Control-A>", onSelectAll)
  tkbind(top, "<Control-w>", onRedo)
  tkbind(top, "<Control-W>", onRedo)
  tkbind(top, "<Alt-BackSpace>", onUndo)
  if (MacOSXP()){
    tkbind(top, "<Meta-x>", onCut)
    tkbind(top, "<Meta-X>", onCut)
    tkbind(top, "<Meta-c>", onCopy)
    tkbind(top, "<Meta-C>", onCopy)
    tkbind(top, "<Meta-v>", onPaste)
    tkbind(top, "<Meta-V>", onPaste)
    tkbind(top, "<Meta-f>", onFind)
    tkbind(top, "<Meta-F>", onFind)
    tkbind(top, "<Meta-a>", onSelectAll)
    tkbind(top, "<Meta-A>", onSelectAll)
    tkbind(top, "<Meta-w>", onRedo)
    tkbind(top, "<Meta-W>", onRedo)
  }
  tkwm.protocol(top, "WM_DELETE_WINDOW", onCancel)
  dialogSuffix(bindReturn = FALSE, resizable=TRUE, focus=editor)
  tkgrid.rowconfigure(top, 0, weight = 0)
  tkgrid.rowconfigure(top, 1, weight = 1)
  tkgrid.rowconfigure(top, 2, weight = 0)
  tkgrid.columnconfigure(top, 0, weight=1)
  tkgrid.rowconfigure(editorFrame, 1, weight=0)
  tkgrid.rowconfigure(editorFrame, 0, weight=1)
  tkgrid.columnconfigure(editorFrame, 0, weight=1)
  tkgrid.columnconfigure(editorFrame, 1, weight=0)
}

# the rgb2col function translates #RRGGBB colors to names if a named color exists or otherwise a "close" color (not exported)
#  uses code from r-help adapted from Kevin Wright

rgb2col <- local({
    all.names <- colors(distinct=TRUE)
    all.lab <- t(convertColor(t(col2rgb(all.names)), from = "sRGB", 
        to = "Lab", scale.in = 255))
    findNear <- function(x.lab) {
        sq.dist <- colSums((all.lab - x.lab)^2)
        rbind(all.names[which.min(sq.dist)], min(sq.dist))
    }
    function(cols.hex, near = 15) { # near = 2.3 is nominally the JND
        cols.lab <- t(convertColor(t(col2rgb(cols.hex)), from = "sRGB", 
            to = "Lab", scale.in = 255))
        cols.near <- apply(cols.lab, 2, findNear)
        ifelse(as.numeric(cols.near[2, ]) < near^2, cols.near[1, ], toupper(cols.hex))
    }
})

# the following function is for plug-ins that test for SciViews (which is no longer supported)

is.SciViews <- function() FALSE

# the following two functions from Milan Bouchet-Valat

setBusyCursor <- function() {
    .commander <- CommanderWindow()
    .menu <- tkcget(.commander, menu=NULL)
    .log <- LogWindow()
    .output <- OutputWindow()
    .messages <- MessagesWindow()
    
    tkconfigure(.commander, cursor="watch")
    tkconfigure(.menu, cursor="watch")
    tkconfigure(.log, cursor="watch")
    tkconfigure(.output, cursor="watch")
    tkconfigure(.messages, cursor="watch")
}

setIdleCursor <- function() {
    .commander <- CommanderWindow()
    .menu <- tkcget(.commander, menu=NULL)
    .log <- LogWindow()
    .output <- OutputWindow()
    .messages <- MessagesWindow()
    
    tkconfigure(.commander, cursor="")
    tkconfigure(.menu, cursor="")
    tkconfigure(.log, cursor="xterm")
    tkconfigure(.output, cursor="xterm")
    tkconfigure(.messages, cursor="xterm")
}


# hasJava <- function(){
#   getRcmdr("capabilities")$java
# }

# setupHelp <- function(){
#   if (MacOSXP() && .Platform$GUI == "AQUA"){
#     current <- system("defaults read org.R-project.R", intern=TRUE)
#     use.external.help <- grep("use.external.help", current)
#     if (length(use.external.help) < 1 || 
#           length(grep("YES", current[use.external.help])) < 1){
#       system("defaults write org.R-project.R use.external.help YES")
#       putRcmdr("restore.use.external.help", TRUE)
#     }
#   }
# }

# Rcmdr data editor

editDataset <- function(data, dsname){
  putRcmdr("dataset.modified", FALSE)
  if (missing(data)){
    if (missing(dsname)) dsname <- "Dataset"
    data <- data.frame(V1="NA")
  }
  else {
    if (!inherits(data, "data.frame")) stop ("data argument must be a data frame")
    if (missing(dsname)) dsname <- deparse(substitute(data))
  }
  if (getRcmdr("crisp.dialogs")) tclServiceMode(on=FALSE)
  top <- tktoplevel(borderwidth = 10)
  tkwm.title(top, paste(gettextRcmdr("Data Editor"), ": ", dsname, sep=""))
  location <- getRcmdr("open.dialog.here")
  pos <- 10 + commanderPosition()
  position <- if (any(pos < 0)) "-50+50" 
  else paste("+", paste(pos, collapse = "+"), sep = "")
  tkwm.geometry(top, position)
  #  tkwm.geometry(top, '-20+200')
  tcl.array <- tclArray()
  nr <- nrow(data)
  nc <- ncol(data)
  for (j in 1:nc){
    data[, j] <- as.character(data[, j])
  }
  colnames <- colnames(data)
  rownames <- rownames(data)
  putRcmdr("data.dim", list(nr=nr, nc=nc, NR=nr, NC=nc))
  # NR, NC not decremented on row/column deletion
  #   to avoid possibly duplicate 
  #   auto-generated row/column names
  for (i in 1:nr) {
    tcl.array[[i + 1, 0]] <- i
    tcl.array[[i + 1, 1]] <- rownames[i]
  }
  for (j in 1:nc){
    tcl.array[[0, j + 1]] <- j
    tcl.array[[1, j + 1]] <- colnames[j]
  }
  tcl.array[[1, 1]] <- "rowname"
  for (i in 1:nr){
    for (j in 1:nc){
      tcl.array[[i + 1, j + 1]] <- data[i, j]
    }
  }
  tableFrame <- tkframe(top)
  data.table <- tk2table(tableFrame, rows=nr + 2, cols=nc + 2, 
                         titlerows=1, titlecols=1,
                         width=nc + 2, height=nr + 2, sparsearray=0,
                         cache=1, flashmode=1, autoclear=1, wrap=1, 
                         colstretchmode="all", rowstretchmode="all",
                         font=getRcmdr('logFont'), anchor="e", padx=6, 
                         resizeborders="both", drawmode="slow",
                         xscrollcommand=function(...) tkset(xscroll,...),
                         yscrollcommand=function(...) tkset(yscroll,...))
  tcl(data.table, "width", 0, max(max(nchar(as.character(nr))), 3))
  tcl(data.table, "width", 1, max(max(nchar(c(rownames, "rowname"))), 3))
  for (j in 1:nc){
    tcl(data.table, "width", j + 1, 
        max(max(nchar(c(colnames[j], data[, j]))), 8))
  }
  xscroll <- ttkscrollbar(tableFrame, orient="horizontal", 
                          command=function(...) tkxview(data.table,...))
  yscroll <- ttkscrollbar(tableFrame,
                          command=function(...) tkyview(data.table,...))
  deleteCell <- function() {
    result <- try(tkdelete(data.table, "active", "0", "end"), silent=TRUE)
    if (inherits(result, "try-error")) return()
    tkinsert(data.table, "active", "0", "NA")
  }
  copyCell <- function(){
    text <- try(tclvalue(tkget(data.table, "active")), silent=TRUE)
    if (inherits(text, "try-error")) return()
    tkclipboard.clear()
    tkclipboard.append(text)
  }
  pasteCell <- function(){
    text <- tclvalue(.Tcl("selection get -selection CLIPBOARD"))
    if (length(text) == 0) return()
    result <- try(tkdelete(data.table, "active", "0", "end"), silent=TRUE)
    if (inherits(result, "try-error")) return()
    tkinsert(data.table, "active", "0", text)
  }
  cutCell <- function(){
    copyCell()
    deleteCell()
  }
  addRow <- function(){
    dims <- getRcmdr("data.dim")
    nr <- dims$nr + 2
    nc <- dims$nc + 1
    NR <- dims$NR + 1
    NC <- dims$NC
    tkinsert(data.table, "row", nr + 1, 1)
    putRcmdr("data.dim", list(nr=nr - 1, nc=nc - 1, NR=NR, NC=NC))
    tcl.array[[nr, 0]] <- NR
    tcl.array[[nr, 1]] <- NR
    for (j in 1:nc) tcl.array[[nr, j + 1]] <- "NA"
    tkconfigure(data.table, width=nc + 1, height=nr + 1)
    tkactivate(data.table, paste0(nr, ",", 2))
    tcl(data.table, "yview", nr)
    tcl(data.table, "xview", 0)
  }
  addCol <- function(){
    dims <- getRcmdr("data.dim")
    nr <- dims$nr + 1
    nc <- dims$nc + 2
    NR <- dims$NR
    NC <- dims$NC + 1
    tkinsert(data.table, "cols", nc + 1, 1)
    putRcmdr("data.dim", list(nr=nr - 1, nc=nc - 1, NR=NR, NC=NC))
    tcl.array[[0, nc]] <- NC
    tcl.array[[1, nc]] <- paste("V", NC, sep="")
    for (i in 1:nr) tcl.array[[i + 1, nc]] <- "NA"
    tkconfigure(data.table, width=nc + 1, height=nr + 1)
    tkactivate(data.table, paste0(2, ",", nc))
    tcl(data.table, "xview", nc)
    tcl(data.table, "yview", 0)
  }
  deleteRow <- function(){
    result <- try(tkdelete(data.table , "rows",
                           tclvalue(tkindex(data.table, "active" ,"row")), 1),
                  silent=TRUE)
    if (inherits(result, "try-error")) return()
    dims <- getRcmdr("data.dim")
    nr <- dims$nr - 1
    nc <- dims$nc
    NR <- dims$NR
    NC <- dims$NC
    putRcmdr("data.dim", list(nr=nr, nc=nc, NR=NR, NC=NC))
  }
  deleteCol <- function(){
    result <- try(tkdelete(data.table , "cols",
                           tclvalue(tkindex(data.table, "active" ,"col")), 1),
                  silent=TRUE)
    if (inherits(result, "try-error")) return()
    dims <- getRcmdr("data.dim")
    nr <- dims$nr
    nc <- dims$nc - 1
    NR <- dims$NR
    NC <- dims$NC
    putRcmdr("data.dim", list(nr=nr, nc=nc, NR=NR, NC=NC))
  }
  onContextMenu <- function(){
    contextMenu <- tkmenu(tkmenu(data.table), tearoff=FALSE)
    tkadd(contextMenu, "command", label=gettextRcmdr("Delete current row"), 
          command=deleteRow)
    tkadd(contextMenu, "command", label=gettextRcmdr("Delete current column"),
          command=deleteCol)
    tkadd(contextMenu, "command", label=gettextRcmdr("Delete cell"), 
          command=deleteCell)
    tkadd(contextMenu, "command", label=gettextRcmdr("Cut cell"), 
          command=cutCell)
    tkadd(contextMenu, "command", label=gettextRcmdr("Copy cell"), 
          command=copyCell)
    tkadd(contextMenu, "command", label=gettextRcmdr("Paste cell"), 
          command=pasteCell)
    tkpopup(contextMenu, tkwinfo("pointerx", data.table), 
            tkwinfo("pointery", data.table))
  }
  onOK <- function(){
    closeDialog()
    dims <- getRcmdr("data.dim")
    nr <- dims$nr + 1
    nc <- dims$nc + 1
    data <- matrix("", nc + 1, nr)
    for (i in 1:nr){
      for (j in 1:nc){
        data[j, i] <- tclvalue(tcl.array[[i, j]])
        if (trim.blanks(data[j, i]) == "") data[j, i] <- "NA"
      }
      data[nc + 1, i] <- "\n"
    }
    data <- paste(data[-1], collapse=" ")
    Data <- read.table(textConnection(data), header=TRUE)
    gassign(dsname, Data)
    activeDataSet(dsname)
    putRcmdr("dataset.modified", TRUE)
  }
  onReturn <- function(){
    location <- try(as.numeric(unlist(strsplit(tclvalue(tkindex(data.table, "active")), ","))), 
                    silent=TRUE)
    if (inherits(location, "try-error")) return()
    text <- tclvalue(tcl.array[[location[1], location[2]]])
    on.exit(tcl.array[[location[1], location[2]]] <- sub("\n", "", text))
    addRow()
  }
  .exit <- function(){
    answer <- RcmdrTkmessageBox("Discard edits?", icon="question", type="yesno", default="no")
    if (as.character(answer) == "no") "abort" else ""
  }
  OKCancelHelp(helpSubject="editDataset")
  editorMenu <- tkmenu(top)
  tkconfigure(top, menu = editorMenu)
  fileMenu <- tkmenu(editorMenu, tearoff=FALSE)
  tkadd(fileMenu, "command", label=gettextRcmdr("Exit and save"), command=onOK)
  tkadd(fileMenu, "command", label=gettextRcmdr("Cancel"), command=onCancel)
  tkadd(editorMenu, "cascade", label=gettextRcmdr("File"), menu=fileMenu)   
  editMenu <- tkmenu(editorMenu, tearoff=FALSE)
  tkadd(editMenu, "command", label=gettextRcmdr("Delete current row"), 
        command=deleteRow)
  tkadd(editMenu, "command", label=gettextRcmdr("Delete current column"), 
        command=deleteCol)
  tkadd(editMenu, "command", label=gettextRcmdr("Add row"), command=addRow)
  tkadd(editMenu, "command", label=gettextRcmdr("Add column"), command=addCol)
  tkadd(editMenu, "command", label=gettextRcmdr("Cut cell"), command=cutCell)
  tkadd(editMenu, "command", label=gettextRcmdr("Copy cell"), command=copyCell)
  tkadd(editMenu, "command", label=gettextRcmdr("Paste cell"), 
        command=pasteCell)
  tkadd(editorMenu, "cascade", label=gettextRcmdr("Edit"), menu=editMenu)   
  helpMenu <- tkmenu(editorMenu, tearoff=FALSE)
  onEditorHelp <- function() print(help("editDataset"))
  tkadd(helpMenu, "command", label=gettextRcmdr("Editor help"), 
        command=onEditorHelp)
  tkadd(editorMenu, "cascade", label=gettextRcmdr("Help"), menu=helpMenu)    
  tkbind(data.table, "<Control-x>", cutCell) # FIXME!
  tkbind(data.table, "<Control-X>", cutCell) #  doesn't work -- source of error unclear
  tkbind(data.table, "<Control-c>", copyCell)
  tkbind(data.table, "<Control-C>", copyCell)
  tkbind(data.table, "<Control-v>", pasteCell)
  tkbind(data.table, "<Control-V>", pasteCell) 
  tkbind(data.table, "<ButtonPress-3>", onContextMenu)
  tkbind(data.table, "<Control-ButtonPress-1>", onContextMenu)
  tkbind(data.table, "<Double-Button-1>", deleteCell)
  if (MacOSXP()){
    tkbind(data.table, "<Meta-x>", cutCell) # FIXME!
    tkbind(data.table, "<Meta-X>", cutCell) #  doesn't work -- source of error unclear
    tkbind(data.table, "<Meta-c>", copyCell)
    tkbind(data.table, "<Meta-C>", copyCell)
    tkbind(data.table, "<Meta-v>", pasteCell)
    tkbind(data.table, "<Meta-V>", pasteCell) 
    tkbind(data.table, "<Meta-ButtonPress-1>", onContextMenu)
  }
  buttonsAddFrame <- tkframe(top)
  addRowButton <- ttkbutton(buttonsAddFrame, command=addRow, 
                            text=gettextRcmdr("Add row"))
  addColButton <- ttkbutton(buttonsAddFrame, command=addCol, 
                            text=gettextRcmdr("Add column"))
  tkgrid(addRowButton, addColButton, sticky="w")
  tkgrid(buttonsAddFrame, sticky="w")
  tkgrid(data.table, yscroll, sticky="news")
  tkgrid.configure(yscroll, sticky="ns")
  tkgrid(xscroll, sticky="ew")
  tkconfigure(data.table, variable=tcl.array, background="lightgray", 
              selectmode="extended")
  tktag.configure(data.table, "active", fg="black", bg="white")
  tktag.configure(data.table, "flash", fg="white", bg="gray")
  tcl(data.table, "tag", "col", "rownos", 0)
  tktag.configure(data.table, "rownos", anchor="e")  
  warn <- options(warn=-1)
  on.exit(warn)
  row.numbers <- !any(is.na(as.numeric(rownames)))
  tcl(data.table, "tag", "col", "rownames", 1)
  tktag.configure(data.table, "rownames", 
                  anchor=if (row.numbers) "e" else "w", bg="darkgray")  
  tcl(data.table, "tag", "row", "colnames", 1)
  tktag.configure(data.table, "colnames", bg="darkgray")  
  tkgrid(tableFrame, sticky="news")
  tkgrid(buttonsFrame, sticky="w")
  tkwm.protocol(top, "WM_DELETE_WINDOW", onCancel)
  dialogSuffix(resizable=TRUE)
  tkgrid.rowconfigure(top, 0, weight = 0)
  tkgrid.rowconfigure(top, 1, weight = 1)
  tkgrid.rowconfigure(top, 2, weight = 0)
  tkgrid.columnconfigure(top, 0, weight = 1)
  tkgrid.rowconfigure(tableFrame, 0, weight = 1)
  tkgrid.rowconfigure(tableFrame, 1, weight = 0)
  tkgrid.columnconfigure(tableFrame, 0, weight = 1)
  tkgrid.columnconfigure(tableFrame, 1, weight = 0)
  tkconfigure(data.table, selectmode = "extended", rowseparator = "\"\n\"", colseparator = "\"\t\"")
  tkconfigure(data.table, multiline = FALSE)
  tkbind(top, "<Key-Return>", onReturn)
  tkbind(top, "<Key-Tab>", addCol)
  tkwait.window(top)
}

# editDataset <- function(data, dsname){
#     putRcmdr("dataset.modified", FALSE)
#     if (missing(data)){
#         if (missing(dsname)) dsname <- "Dataset"
#         data <- data.frame(V1="NA")
#     }
#     else {
#         if (!inherits(data, "data.frame")) stop ("data argument must be a data frame")
#         if (missing(dsname)) dsname <- deparse(substitute(data))
#     }
#     if (getRcmdr("crisp.dialogs")) tclServiceMode(on=FALSE)
#     top <- tktoplevel(borderwidth = 10)
#     tkwm.title(top, paste(gettextRcmdr("Data Editor"), ": ", dsname, sep=""))
# #     location <- getRcmdr("open.dialog.here")  # FIXME!
# #     pos <- 10 + commanderPosition()           # Don't do this because window doesn't stay on top
# #     position <- if (any(pos < 0)) "-50+50" 
# #     else paste("+", paste(pos, collapse = "+"), sep = "")
# #     tkwm.geometry(top, position)
#     tkwm.geometry(top, '-20+200')
#     tcl.array <- tclArray()
#     nr <- nrow(data)
#     nc <- ncol(data)
#     for (j in 1:nc){
#         data[, j] <- as.character(data[, j])
#     }
#     colnames <- colnames(data)
#     rownames <- rownames(data)
#     putRcmdr("data.dim", list(nr=nr, nc=nc, NR=nr, NC=nc))
#     # NR, NC not decremented on row/column deletion
#     #   to avoid possibly duplicate 
#     #   auto-generated row/column names
#     for (i in 1:nr) {
#         tcl.array[[i + 1, 0]] <- i
#         tcl.array[[i + 1, 1]] <- rownames[i]
#     }
#     for (j in 1:nc){
#         tcl.array[[0, j + 1]] <- j
#         tcl.array[[1, j + 1]] <- colnames[j]
#     }
#     tcl.array[[1, 1]] <- "rowname"
#     for (i in 1:nr){
#         for (j in 1:nc){
#             tcl.array[[i + 1, j + 1]] <- data[i, j]
#         }
#     }
#     tableFrame <- tkframe(top)
#     data.table <- tk2table(tableFrame, rows=nr + 2, cols=nc + 2, 
#         titlerows=1, titlecols=1,
#         width=nc + 2, height=nr + 2, sparsearray=0,
#         cache=1, flashmode=1, autoclear=1, wrap=1, 
#         colstretchmode="all", rowstretchmode="all",
#         font=getRcmdr('logFont'), anchor="e", padx=6, 
#         resizeborders="both",
#         xscrollcommand=function(...) tkset(xscroll,...),
#         yscrollcommand=function(...) tkset(yscroll,...))
#     tcl(data.table, "width", 0, max(max(nchar(as.character(nr))), 3))
#     tcl(data.table, "width", 1, max(max(nchar(c(rownames, "rowname"))), 3))
#     for (j in 1:nc){
#         tcl(data.table, "width", j + 1, 
#             max(max(nchar(c(colnames[j], data[, j]))), 8))
#     }
#     xscroll <- ttkscrollbar(tableFrame, orient="horizontal", 
#         command=function(...) tkxview(data.table,...))
#     yscroll <- ttkscrollbar(tableFrame,
#         command=function(...) tkyview(data.table,...))
#     deleteCell <- function() {
#         result <- try(tkdelete(data.table, "active", "0", "end"), silent=TRUE)
#         if (inherits(text, "try-error")) return()
#         tkinsert(data.table, "active", "0", "NA")
#     }
#     copyCell <- function(){
#         text <- try(tclvalue(tkget(data.table, "active")), silent=TRUE)
#         if (inherits(text, "try-error")) return()
#         tkclipboard.clear()
#         tkclipboard.append(text)
#     }
#     pasteCell <- function(){
#         text <- tclvalue(.Tcl("selection get -selection CLIPBOARD"))
#         if (length(text) == 0) return()
#         result <- try(tkdelete(data.table, "active", "0", "end"), silent=TRUE)
#         if (inherits(text, "try-error")) return()
#         tkinsert(data.table, "active", "0", text)
#     }
#     cutCell <- function(){
#         copyCell()
#         deleteCell()
#     }
#     addRow <- function(){
#         dims <- getRcmdr("data.dim")
#         nr <- dims$nr + 2
#         nc <- dims$nc + 1
#         NR <- dims$NR + 1
#         NC <- dims$NC
#         tkinsert(data.table, "row", nr + 1, 1)
#         putRcmdr("data.dim", list(nr=nr - 1, nc=nc - 1, NR=NR, NC=NC))
#         tcl.array[[nr, 0]] <- NR
#         tcl.array[[nr, 1]] <- NR
#         for (j in 1:nc) tcl.array[[nr, j + 1]] <- "NA"
#         tkconfigure(data.table, width=nc + 1, height=nr + 1)
#     }
#     addCol <- function(){
#         dims <- getRcmdr("data.dim")
#         nr <- dims$nr + 1
#         nc <- dims$nc + 2
#         NR <- dims$NR
#         NC <- dims$NC + 1
#         tkinsert(data.table, "cols", nc + 1, 1)
#         putRcmdr("data.dim", list(nr=nr - 1, nc=nc - 1, NR=NR, NC=NC))
#         tcl.array[[0, nc]] <- NC
#         tcl.array[[1, nc]] <- paste("V", NC, sep="")
#         for (i in 1:nr) tcl.array[[i + 1, nc]] <- "NA"
#         tkconfigure(data.table, width=nc + 1, height=nr + 1)
#     }
#     deleteRow <- function(){
#         result <- try(tkdelete(data.table , "rows",
#             tclvalue(tkindex(data.table, "active" ,"row")), 1),
#             silent=TRUE)
#         if (inherits(result, "try-error")) return()
#         dims <- getRcmdr("data.dim")
#         nr <- dims$nr - 1
#         nc <- dims$nc
#         NR <- dims$NR
#         NC <- dims$NC
#         putRcmdr("data.dim", list(nr=nr, nc=nc, NR=NR, NC=NC))
#     }
#     deleteCol <- function(){
#         result <- try(tkdelete(data.table , "cols",
#             tclvalue(tkindex(data.table, "active" ,"col")), 1),
#             silent=TRUE)
#         if (inherits(result, "try-error")) return()
#         dims <- getRcmdr("data.dim")
#         nr <- dims$nr
#         nc <- dims$nc - 1
#         NR <- dims$NR
#         NC <- dims$NC
#         putRcmdr("data.dim", list(nr=nr, nc=nc, NR=NR, NC=NC))
#     }
#     onContextMenu <- function(){
#         contextMenu <- tkmenu(tkmenu(data.table), tearoff=FALSE)
#         tkadd(contextMenu, "command", label=gettextRcmdr("Delete current row"), 
#             command=deleteRow)
#         tkadd(contextMenu, "command", label=gettextRcmdr("Delete current column"),
#             command=deleteCol)
#         tkadd(contextMenu, "command", label=gettextRcmdr("Delete cell"), 
#             command=deleteCell)
#         tkadd(contextMenu, "command", label=gettextRcmdr("Cut cell"), 
#             command=cutCell)
#         tkadd(contextMenu, "command", label=gettextRcmdr("Copy cell"), 
#             command=copyCell)
#         tkadd(contextMenu, "command", label=gettextRcmdr("Paste cell"), 
#             command=pasteCell)
#         tkpopup(contextMenu, tkwinfo("pointerx", data.table), 
#             tkwinfo("pointery", data.table))
#     }
#     onOK <- function(){
#         closeDialog()
#         dims <- getRcmdr("data.dim")
#         nr <- dims$nr + 1
#         nc <- dims$nc + 1
#         data <- matrix("", nc + 1, nr)
#         for (i in 1:nr){
#             for (j in 1:nc){
#                 data[j, i] <- tclvalue(tcl.array[[i, j]])
#             }
#             data[nc + 1, i] <- "\n"
#         }
#         data <- paste(data[-1], collapse=" ")
#         Data <- read.table(textConnection(data), header=TRUE)
#         gassign(dsname, Data)
#         activeDataSet(dsname)
#         putRcmdr("dataset.modified", TRUE)
#     }
#     .exit <- function(){
#         answer <- RcmdrTkmessageBox("Discard edits?", icon="question", type="yesno")
#         if (as.character(answer) == "no") "abort" else ""
#     }
#     OKCancelHelp(helpSubject="editDataset")
#     editorMenu <- tkmenu(top)
#     tkconfigure(top, menu = editorMenu)
#     fileMenu <- tkmenu(editorMenu, tearoff=FALSE)
#     tkadd(fileMenu, "command", label=gettextRcmdr("Exit and save"), command=onOK)
#     tkadd(fileMenu, "command", label=gettextRcmdr("Cancel"), command=onCancel)
#     tkadd(editorMenu, "cascade", label=gettextRcmdr("File"), menu=fileMenu)   
#     editMenu <- tkmenu(editorMenu, tearoff=FALSE)
#     tkadd(editMenu, "command", label=gettextRcmdr("Delete current row"), 
#         command=deleteRow)
#     tkadd(editMenu, "command", label=gettextRcmdr("Delete current column"), 
#         command=deleteCol)
#     tkadd(editMenu, "command", label=gettextRcmdr("Add row"), command=addRow)
#     tkadd(editMenu, "command", label=gettextRcmdr("Add column"), command=addCol)
#     tkadd(editMenu, "command", label=gettextRcmdr("Cut cell"), command=cutCell)
#     tkadd(editMenu, "command", label=gettextRcmdr("Copy cell"), command=copyCell)
#     tkadd(editMenu, "command", label=gettextRcmdr("Paste cell"), 
#         command=pasteCell)
#     tkadd(editorMenu, "cascade", label=gettextRcmdr("Edit"), menu=editMenu)   
#     helpMenu <- tkmenu(editorMenu, tearoff=FALSE)
#     onEditorHelp <- function() print(help("editDataset"))
#     tkadd(helpMenu, "command", label=gettextRcmdr("Editor help"), 
#         command=onEditorHelp)
#     tkadd(editorMenu, "cascade", label=gettextRcmdr("Help"), menu=helpMenu)    
#     #    tkbind(data.table, "<Control-x>", cutCell) # FIXME!
#     #    tkbind(data.table, "<Control-X>", cutCell) #  doesn't work -- source of error unclear
#     tkbind(data.table, "<Control-c>", copyCell)
#     tkbind(data.table, "<Control-C>", copyCell)
#     tkbind(data.table, "<Control-v>", pasteCell)
#     tkbind(data.table, "<Control-V>", pasteCell) 
#     tkbind(data.table, "<ButtonPress-3>", onContextMenu)
#     tkbind(data.table, "<Control-ButtonPress-1>", onContextMenu)
#     tkbind(data.table, "<Double-Button-1>", deleteCell)
#     if (MacOSXP()){
#         tkbind(data.table, "<Meta-c>", copyCell)
#         tkbind(data.table, "<Meta-C>", copyCell)
#         tkbind(data.table, "<Meta-v>", pasteCell)
#         tkbind(data.table, "<Meta-V>", pasteCell) 
#         tkbind(data.table, "<Meta-ButtonPress-1>", onContextMenu)
#     }
#     buttonsAddFrame <- tkframe(top)
#     addRowButton <- ttkbutton(buttonsAddFrame, command=addRow, 
#         text=gettextRcmdr("Add row"))
#     addColButton <- ttkbutton(buttonsAddFrame, command=addCol, 
#         text=gettextRcmdr("Add column"))
#     tkgrid(addRowButton, addColButton, sticky="w")
#     tkgrid(buttonsAddFrame, sticky="w")
#     tkgrid(data.table, yscroll, sticky="news")
#     tkgrid.configure(yscroll, sticky="ns")
#     tkgrid(xscroll, sticky="ew")
#     tkconfigure(data.table, variable=tcl.array, background="lightgray", 
#         selectmode="extended")
#     tktag.configure(data.table, "active", fg="black", bg="white")
#     tktag.configure(data.table, "flash", fg="white", bg="gray")
#     tcl(data.table, "tag", "col", "rownos", 0)
#     tktag.configure(data.table, "rownos", anchor="e")  
#     warn <- options(warn=-1)
#     on.exit(warn)
#     row.numbers <- !any(is.na(as.numeric(rownames)))
#     tcl(data.table, "tag", "col", "rownames", 1)
#     tktag.configure(data.table, "rownames", 
#         anchor=if (row.numbers) "e" else "w", bg="darkgray")  
#     tcl(data.table, "tag", "row", "colnames", 1)
#     tktag.configure(data.table, "colnames", bg="darkgray")  
#     tkgrid(tableFrame, sticky="news")
#     tkgrid(buttonsFrame, sticky="w")
#     tkwm.protocol(top, "WM_DELETE_WINDOW", onCancel)
#     dialogSuffix(resizable=TRUE)
#     tkgrid.rowconfigure(top, 0, weight = 0)
#     tkgrid.rowconfigure(top, 1, weight = 1)
#     tkgrid.rowconfigure(top, 2, weight = 0)
#     tkgrid.columnconfigure(top, 0, weight = 1)
#     tkgrid.rowconfigure(tableFrame, 0, weight = 1)
#     tkgrid.rowconfigure(tableFrame, 1, weight = 0)
#     tkgrid.columnconfigure(tableFrame, 0, weight = 1)
#     tkgrid.columnconfigure(tableFrame, 1, weight = 0)
#     tkwait.window(top)
# }

# some Mac OS X related functions

RappP <- function() .Platform$GUI == "AQUA"

mavericksP <- function(){
  info <- Sys.info()
  info["sysname"] == "Darwin" && info["release"] >= "13.0.0"
}

appnap <- function(state=c("on", "off", "delete")){
  if (!mavericksP()) stop("requires OS X 10.9 or greater")
  save <- options(warn = -1)
  on.exit(options(save))
  if (missing(state)){
    res <- system("defaults read org.R-project.R NSAppSleepDisabled", 
                  intern=TRUE, ignore.stderr=TRUE)
    return(c("on", "off")[1 + (length(res) > 0 && res == "1")])
  }
  state <- match.arg(state)
  switch(state,
         delete = system("defaults delete org.R-project.R NSAppSleepDisabled", ignore.stderr=TRUE),
         off = system("defaults write org.R-project.R NSAppSleepDisabled -bool YES"),
         on = system("defaults write org.R-project.R NSAppSleepDisabled -bool NO")
  )
  return(state)
}

# replacement for standard tkmenu() to play better with ttk themes
#  courtesy of Philippe Grosjean

tkmenu <- function (parent, activebackground, activeforeground, ...) {
  if (!is.ttk()) 
    stop("Tcl/Tk >= 8.5 is required")
  w <- tkwidget(parent, "menu", ...)
  if (missing(activebackground)) activebackground <- tk2style("tk2button", "selectbackground")
  if (activebackground == "") activebackground = "darkblue" # Default value
  if (missing(activeforeground)) activeforeground <- tk2style("tk2button", "selectforeground")
  if (activeforeground == "") activeforeground = "white" # Default value
  tkconfigure(w, activebackground = activebackground, activeforeground = activeforeground)
  class(w) <- c("tk2menu", "tk2widget", class(w))
  return(w)
}

hasProgram <- function(program, version, prefix="--", line=1, compare=`>=`){
    # Args:
    #   program: quoted name of program
    #   version: quoted version number (numerals . -)
    #   prefix:  for version switch
    #   line:    output line containing version number
    #   compare: comparison operator for version
    # Example: hasProgram("pandoc", version="1.12")
    path <- Sys.which(program)
    present <-  path != ""
    if (missing(version) || !present) return(as.vector(present))
    result <- try(system2(path, args=paste(prefix, "version", sep=""), stderr=TRUE, stdout=TRUE),
        silent=TRUE)
    if (inherits(result, "try-error")){
        warning("could not execute '", path, " --version'")
        return(FALSE)
    }
    result <- result[line]
    match <- regexpr("((?:(\\d+)[\\.-])?)*(\\*|\\d+)", result)
    if (match == -1) {
        warning("could not locate version number in\n", result)
        return(FALSE)
    }
    result <- substr(result, match, match + attr(match, "match.length") - 1)
    compare(result, version)
}

RcmdrCapabilities <- function(check=list(c("pdflatex"), c("pandoc", version="1.12.3"))){
    result <- vector(length(check), mode="list")
    names(result) <- sapply(check, function(x) x[1])
    for (i in 1:length(check)){
        result[[i]] <- do.call(hasProgram, as.list(check[[i]]))
    }
    result
}

browsePDF <- function(file) {
    if (WindowsP()) shell.exec(file)
    else if (MacOSXP()) system(paste("open -a Preview", shQuote(file)))
    else system(paste(shQuote(getOption("pdfviewer")), shQuote(file)), wait=FALSE)
}

# function to insure that "levels" of character variables are returned

levels.character <- function(x) sort(unique(x))

