

qdataview <- function(x, ...)
{
    UseMethod("qdataview")
}

qdataview.matrix <- qdataview.table <- qdataview.array <-
    function(x, ...)
{
    cdim <- dim(x)
    cdimnames <- dimnames(x)
    if (length(cdim) > 2)
        stop("Arrays of more than two dimensions not supported yet")
    if (length(cdim) == 1)
    {
        cdim <- dim(x) <- c(cdim, 1L)
        ## FIXME: need to change cdimnames?
    }
    ans <- Qt$QTableWidget(cdim[1], cdim[2])
    sx <- as.character(x)
    rowx <- row(x) - 1L
    colx <- col(x) - 1L
    for (i in seq_len(cdim[1] * cdim[2]))
    {
        item <- Qt$QTableWidgetItem(sx[i])
        item$setFlags(33L) ## selectable(1) | can interact(32)
        ans$setItem(rowx[i], colx[i], item)
    }
    if (!is.null(cdimnames))
    {
        if (!is.null(cdimnames[[1]]))
            ans$setVerticalHeaderLabels(cdimnames[[1]])
        if (length(cdimnames) > 1 && !is.null(cdimnames[[2]]))
            ans$setHorizontalHeaderLabels(cdimnames[[2]])
    }
    ans$resize(600, 400)
    ans
}

qdataview.data.frame <- function(x, ...)
{
    ## qdataview.array(do.call(cbind, lapply(x, as.character)))
    m <- qdataFrameModel(x, ...)
    v <- Qt$QTableView()
    v$setModel(m)
    v$resize(600, 400)
    v
}



