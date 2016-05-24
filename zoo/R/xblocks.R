
xblocks <- function(x, ...)
    UseMethod("xblocks")

xblocks.default <-
    function (x, y, ..., col = NULL, border = NA, 
              ybottom = par("usr")[3], ytop = ybottom + height,
              height = diff(par("usr")[3:4]),
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
    ## for logical series, col default comes from palette()
    if (is.logical(y) && is.null(col))
        col <- palette()[1]
    ## set block colours if 'col' given
    if (length(col) > 0) {
        if (is.character(col))
            col[col == ""] <- NA
        ok <- !is.na(blockCol)
        blockCol[ok] <- rep(col, length.out = sum(ok)) ## rep to avoid warnings
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
    ## adjust for log scales
    if (par("ylog")) {
        ybottom <- 10^ybottom
        ytop <- 10^ytop
    }
    ## draw it
    rect(xleft = blockStart, xright = blockEnd,
         ybottom = ybottom, ytop = ytop,
         col = blockCol, border = border, ...)
}

xblocks.zoo <-
xblocks.ts <-
    function(x, y = x, ...)
{
    if (!is.function(y))
        y <- coredata(y)
    xblocks(index(x), y, ...)
}
