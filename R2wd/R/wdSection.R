wdSection <-
function (title, label = gsub("[.,-:?!@#+* ]", "_", paste("sec", title, sep = "_")),
    newpage = FALSE, wdapp = .R2wd)
{
    wdSectionBreak(continuous = !newpage,wdapp=wdapp)
    wdInsertBookmark(label,wdapp=wdapp)
    wdHeading(1,title,wdapp=wdapp)
}

