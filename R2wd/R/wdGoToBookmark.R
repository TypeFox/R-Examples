wdGoToBookmark <-
function (bookmark, wdapp = .R2wd)
{
    bookmarks<-wdapp[['ActiveDocument']][['Bookmarks']]
    ## locate the bookmark
    ## this is very crude, in principle it should be possible to navigate
    ## directly to the bookmark but I haven't found the right
    ## way of doing it.
    for (i in 1:bookmarks[['Count']]) if (bookmarks$Item(i)[['Name']]==bookmark) break
    bookmarks$Item(i)$Select()
    invisible()
}

