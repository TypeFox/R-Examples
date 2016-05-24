wdSetFont <-
function (fontname=NULL,fontsize=NULL, bold=NULL,italic=NULL,wdapp = .R2wd)
{
    wdsel <- wdapp[['Selection']]
    if (!is.null(fontname)) wdsel[['Font']][['Name']] <- fontname
    if (!is.null(fontsize)) wdsel[['Font']][['Size']] <- fontsize
    if (!is.null(bold)) wdsel[['Font']][['Bold']]<-bold
    if (!is.null(italic)) wdsel[['Font']][['Italic']]<-italic
    invisible()
}

