loadworkspace <- function () 
{
    file <- tclvalue(tkgetOpenFile(filetypes = gettextRcmdr("{\"R Data Files\" {\".RData\" \".rda\" \".Rda\" \".RDA\"}} {\"All Files\" {\"*\"}}")))
    if (file == "") 
        return()
    command <- paste("load(\"", file, "\")", sep = "")
    dsname <- justDoIt(command)
    logger(command)
    if (length(dsname)==1){ 
        if (class(dsname)[1] != "try-error") 
        activeDataSet(dsname)
    }
    else{
       dsname <- dsname[sapply(dsname, function(obj) "data.frame" %in% class(get(obj)))]
       dsname <- intersect(dsname, listDesigns())
       if (length(dsname)>=1) 
          activeDataSet(dsname[1])
    }
    tkfocus(CommanderWindow())
}