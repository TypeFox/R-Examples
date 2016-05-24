wdUndo <-
function (times=1, wdapp = .R2wd)
{
    wddoc <- wdapp[["ActiveDocument"]]
    wddoc$Undo(times)
}

