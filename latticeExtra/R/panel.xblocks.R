##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


panel.xblocks <- function(x, ...)
    UseMethod("panel.xblocks")

panel.xblocks.default <-
    function (x, y, ..., col = NULL, border = NA, 
              height = unit(1, "npc"),
              block.y = unit(0, "npc"), vjust = 0,
              name = "xblocks", gaps = FALSE,
              last.step = median(diff(tail(x))))
{
    if (is.function(y))
        y <- y(x)
    x <- as.numeric(x)
    if (length(x) == 0) return()
    if (is.unsorted(x, na.rm = TRUE))
        stop("'x' should be ordered (increasing)")
    if (is.na(last.step))
        last.step <- 0
    if (gaps) {
        .Deprecated(msg = "The 'gaps' argument is deprecated; use panel.xblocks(time(z), is.na(z))")
        y <- is.na(y)
    }
    ## Three cases:
    ## (1) If y is character, assume it gives the block colours
    ## -- unless 'col' is given, which over-rides it.
    ## (2) If y is logical, show blocks of TRUE values.
    ## (3) If y is numeric, show blocks of non-NA values.
    if (is.logical(y)) {
        y <- y
    } else if (is.numeric(y)) {
        y <- !is.na(y)
    } else {
        ## this will convert factor, Date, etc to character:
        y <- as.character(y)
    }
    ## Note: rle treats each NA as unique (does not combine runs of NAs)
    ## so we need to replace NAs with a temporary value.
    NAval <-
        if (is.character(y)) "" else FALSE
    y[is.na(y)] <- NAval
    ## find blocks (runs of constant values)
    yrle <- rle(y)
    ## substitute NA values back in
    blockCol <- yrle$values
    blockCol[blockCol == NAval] <- NA
    ## for logical series, col default comes from current theme
    if (is.logical(y) && is.null(col))
        col <- trellis.par.get("plot.line")$col
    ## set block colours if 'col' given
    if (length(col) > 0) {
        if (is.character(col))
            col[col == ""] <- NA
        ok <- !is.na(blockCol)
        blockCol[ok] <- rep(col, length = sum(ok)) ## rep to avoid warnings
    }
    ## work out block geometry
    idxBounds <- cumsum(c(1, yrle$lengths))
    idxStart <- head(idxBounds, -1)
    idxEnd <- tail(idxBounds, -1)
    idxEnd[length(idxEnd)] <- length(y)
    blockStart <- x[idxStart]
    blockEnd <- x[idxEnd]
    blockEnd[length(blockEnd)] <- tail(blockEnd, 1) + last.step
    blockWidth <- blockEnd - blockStart
    ## draw it
    grid::grid.rect(x = blockStart, width = blockWidth,
                    y = block.y, height = height,
                    hjust = 0, vjust = vjust,
                    default.units = "native", name = name,
                    gp = gpar(fill = blockCol, col = border, ...))
}

panel.xblocks.ts <-
    function(x, y = x, ...)
{
    if (!is.function(y)) y <- as.vector(y)
    panel.xblocks(as.vector(time(x)), y, ...)
}

panel.xblocks.zoo <-
    function(x, y = x, ...)
{
    if (!is.function(y)) y <- zoo::coredata(y)
    panel.xblocks(zoo::index(x), y, ...)
}

