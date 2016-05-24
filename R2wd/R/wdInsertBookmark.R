wdInsertBookmark <-
function (text, wdapp = .R2wd)
{
    wdapp[['ActiveDocument']][['bookmarks']]$Add(text)
}

