plotbg <-
function (map = "kola.background", which.map = c(1, 2, 3, 4),
    map.col = c(5, 1, 3, 4), map.lwd = c(2, 1, 2, 1), add.plot = FALSE, ...)
{
    all = get(map)
    xrange = c(min(all[[1]][, 1], na.rm = TRUE), max(all[[1]][,
        1], na.rm = TRUE))
    yrange = c(min(all[[1]][, 2], na.rm = TRUE), max(all[[1]][,
        2], na.rm = TRUE))
    if (!add.plot) {
        plot(1, 1, xlim = xrange, ylim = yrange, ...)
    }
    for (i in 1:length(which.map)) {
        lines(all[[which.map[i]]], col = map.col[i], lwd = map.lwd[i])
    }
}

