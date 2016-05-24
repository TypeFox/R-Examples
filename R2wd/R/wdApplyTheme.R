wdApplyTheme <-
function(theme="None",wdapp=.R2wd){
    wddoc<-wdapp[['ActiveDocument']]
    if (theme=="None") wddoc$RemoveTheme() else wddoc$ApplyTheme(theme)
    return()
}

