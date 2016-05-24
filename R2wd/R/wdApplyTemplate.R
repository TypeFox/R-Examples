wdApplyTemplate <-
function(filename,wdapp=.R2wd){
    wddoc<-wdapp[['ActiveDocument']]
    wddoc[["UpdateStylesOnOpen"]]<-TRUE
    wddoc[["AttachedTemplate"]]<-filename
    return()
}

