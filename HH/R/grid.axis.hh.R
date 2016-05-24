"grid.yaxis.hh" <-
function (at = NULL, label = TRUE, main = TRUE, gp = gpar(),
    draw = TRUE, vp = NULL, labels)
{
    if (is.null(at))
        if (is.null(vp)) {
            major <- NULL
            ticks <- NULL
            labels <- NULL
        }
        else at <- grid.pretty(vp$yscale)
    if (!is.null(at)) {
        ## major <- grid:::make.yaxis.major(at, main)
        ## ticks <- grid:::make.yaxis.ticks(at, main)
        major <- grid.make.yaxis.major(at, main)
        ticks <- grid.make.yaxis.ticks(at, main)
        if (label) {
          if(missing(labels))
            labels <- make.yaxis.hh.labels(at, main)
          else
            labels <- make.yaxis.hh.labels(at, main, labels)
        }
        else labels <- NULL
    }
    grob(list(at = at, major = major, ticks = ticks, labels = labels,
        label = label, gp = gp, main = main, vp = vp), c("yaxis",
        "axis"), draw)
}

"make.yaxis.hh.labels" <-
function (at, main, labels=at)
{
    if (main) {
        hjust <- "right"
        label.x <- unit(-1, "lines")
    }
    else {
        hjust <- "left"
        label.x <- unit(1, "npc") + unit(1, "lines")
    }
    just <- c(hjust, "centre")
    grid.text(as.character(labels), label.x, unit(at, "native"),
        just = just, rot = 0, check.overlap = TRUE, draw = FALSE)
}

"grid.xaxis.hh" <-
function (at = NULL, label = TRUE, main = TRUE, gp = gpar(),
    draw = TRUE, vp = NULL, labels)
{
    if (is.null(at))
        if (is.null(vp)) {
            major <- NULL
            ticks <- NULL
            labels <- NULL
        }
        else at <- grid.pretty(vp$xscale)
    if (!is.null(at)) {
        ## major <- grid:::make.xaxis.major(at, main)
        ## ticks <- grid:::make.xaxis.ticks(at, main)
        major <- grid.make.xaxis.major(at, main)
        ticks <- grid.make.xaxis.ticks(at, main)
        if (label) {
          if(missing(labels))
            labels <- make.xaxis.hh.labels(at, main)
          else
            labels <- make.xaxis.hh.labels(at, main, labels)
        }
        else labels <- NULL
    }
    grob(list(at = at, major = major, ticks = ticks, labels = labels,
        label = label, gp = gp, main = main, vp = vp), c("xaxis",
        "axis"), draw)
}

"make.xaxis.hh.labels" <-
function (at, main, labels=at)
{
    if (main)
        label.y <- unit(-1.5, "lines")
    else label.y <- unit(1, "npc") + unit(1.5, "lines")
    grid.text(as.character(labels), unit(at, "native"), label.y,
        just = "centre", rot = 0, check.overlap = TRUE, draw = FALSE)
}
