## necessary only because filter() is not generic

tisFilter <- function(x, ...)  naWindow(tis(filter(x, ...), start = start(x)))

