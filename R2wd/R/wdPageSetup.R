wdPageSetup<-
function (orientation="portrait",margins=rep(1,4),scope="section", wdapp = .R2wd) {
    handle<-switch(scope,
           "section"=wdapp[['Selection']],
           "all"=wdapp[['Selection']],
           stop("scope has to be either section or all"))[['Sections']][['PageSetup']]
    if (orientation=="portrait") handle[['Orientation']]<-0
    if (orientation=="landscape") handle[['Orientation']]<-1
    handle[['Bottommargin']]<-margins[1]*72
    handle[['Leftmargin']]<-margins[2]*72
    handle[['Topmargin']]<-margins[3]*72
    handle[['Rightmargin']]<-margins[4]*72
    return()
}
