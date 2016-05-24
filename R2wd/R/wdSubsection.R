wdSubsection <-
function (title, label = gsub("[.,-:?!@#+* ]", "_", paste("subsec", title,
    sep = "_")), newpage = FALSE, wdapp = .R2wd)
{
    wdSectionBreak(continuous = !newpage,wdapp=wdapp)
    wdInsertBookmark(label,wdapp=wdapp)
    wdHeading(2,title,wdapp=wdapp)
}

