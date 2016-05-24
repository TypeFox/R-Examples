wdEnumerate <-
function (Gallery=2,Template=1, wdapp = .R2wd)
{
    wdsel <- wdapp[['Selection']]
    ttt<-(wdapp[["ListGalleries"]]$Item(Gallery))[["ListTemplates"]]$Item(Template)
    wdsel[["Range"]][["ListFormat"]]$ApplyListTemplateWithLevel(ttt)
}

