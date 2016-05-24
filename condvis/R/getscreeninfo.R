getscreeninfo <-
# This is old and likely to be deprecated soon
function (plotobject)
{
    UseMethod("getscreeninfo")
}

getscreeninfo.ceplot <-
function (plotobject)
{
    xcplots.index <- seq_along(plotobject$xcplots)
    screen.number <- vapply(xcplots.index, 
        function(i) plotobject$xcplots[[i]]$screen, integer(1L))
    coords <- t(vapply(screen.number, 
        function(i){ screen(i, new = F); par("fig") }, numeric(4L)))
        colnames(coords) <- c("xleft", "xright", "ybottom", "ytop")
    data.frame(xcplots.index, screen.number, coords)
}

getscreeninfo.conditionselector <-
function (plotobject)
{
    xcplots.index <- seq_along(plotobject$xcplots)
    screen.number <- plotobject$screens
    coords <- t(vapply(screen.number, 
        function(i){ screen(i, new = F); par("fig") }, numeric(4L)))
        colnames(coords) <- c("xleft", "xright", "ybottom", "ytop")
    data.frame(xcplots.index, screen.number, coords)
}