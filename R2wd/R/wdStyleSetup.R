wdStyleSetup <-
    function (style="Normal",fontsize=11,align=3,wdapp = .R2wd)
{
    handle<-.R2wd[["ActiveDocument"]][["Styles"]]$Item(style)
    handle[["Font"]][["Size"]]<-fontsize
    handle[["ParagraphFormat"]][["Alignment"]]<-align
    return()
}
