wdEqn <-
function (eqtext, bookmark = NULL, iknow=FALSE, waitsec=2,wdapp = .R2wd, paragraph = TRUE){
    if (!iknow) stop("Beware: this function is dangerous (it works using sendkeys). If you wish to use it anyway you must set the argument iknow=TRUE")
    wsshell<-try(COMCreate("WScript.Shell"),silent=TRUE)
    if(class(wsshell)=="try-error") {
        wsshell<-comGetObject("WScript.Shell")
        if (is.null(wsshell)) wsshell<-comCreateObject("WScript.Shell")
    }
    eqtext<-gsub("\\("," \\{\\(\\}",eqtext)
    eqtext<-gsub("\\)"," \\{\\)\\}",eqtext)
    eqtext<-gsub("\\+"," \\{\\+\\}",eqtext)
    eqtext<-gsub("\\*"," \\{\\*\\}",eqtext)
    eqtext<-gsub("\\^"," \\{\\^\\}",eqtext)
    eqtext<-gsub("="," \\{=\\}",eqtext)
    eqtext<-gsub("~"," \\{~\\}",eqtext)
    wdsel <- wdapp[["Selection"]]
    objrange <- wdsel[["Range"]]
    objrange[["Text"]] <- " "
    wdsel[["OMaths"]]$Add(objrange)
    wdsel[["OMaths"]]$Item(1)$BuildUp()
    wdsel$EndKey()
    sendtext<-ifelse(paragraph,paste(eqtext,"{RIGHT}~"),paste(eqtext,"{RIGHT} "))
    wsshell$AppActivate("Microsoft Word")
    wsshell$SendKeys(sendtext)
    Sys.sleep(waitsec)
    wsshell<-NULL
    return()
}
