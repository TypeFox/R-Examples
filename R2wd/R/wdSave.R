wdSave <-
function (Name = NULL, wdapp = .R2wd)
{
    wddoc <- wdapp[["ActiveDocument"]]
    if (is.null(Name))
        wddoc$Save()
    else wddoc$SaveAs(Name)
}

