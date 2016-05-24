wdSectionBreak <-
function (continuous = TRUE, bookmark = NULL,wdapp = .R2wd)
{
    wdsel <- wdapp[["Selection"]]
    if (continuous)
        wdsel$InsertBreak(3)
    else wdsel$InsertBreak(2)
    if (!is.null(bookmark)) wdInsertBookmark(bookmark)
}

