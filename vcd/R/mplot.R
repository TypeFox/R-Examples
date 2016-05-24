mplot <- function(...,
                  .list = list(),
                  layout = NULL, cex = NULL,
                  main = NULL, gp_main = gpar(fontsize = 20),
                  sub = NULL, gp_sub = gpar(fontsize = 15),
                  keep_aspect_ratio = TRUE) {
    l <- c(list(...), .list)
    ll <- length(l)
    m <- !is.null(main)
    s <- !is.null(sub)

    ## calculate layout
    if (is.null(layout))
        layout <- c(trunc(sqrt(ll)), ceiling(ll / trunc(sqrt(ll))))

    ## push base layout
    grid.newpage()
    hts = unit(1 - 0.1 * m - 0.1 * s, "null")
    if (m)
        hts <- c(unit(0.1, "null"), hts)
    if (s)
        hts <- c(hts, unit(0.1, "null"))
    pushViewport(viewport(layout =
                              grid.layout(nrow = 1 + m + s,
                                          ncol = 1, heights = hts)
                          )
                 )

    ## push main, if any
    if (!is.null(main)) {
        pushViewport(viewport(layout.pos.row = 1,
                              layout.pos.col = NULL))

        grid.text(main, gp = gp_main)
        popViewport(1)
    }

    ## push strucplots
    if (is.null(cex))
        cex <- sqrt(1/layout[1])
    pushViewport(viewport(layout.pos.row = 1 + m,
                          layout.pos.col = NULL))
    pushViewport(viewport(layout = grid.layout(nrow = layout[1],
                          ncol = layout[2]),
                          gp = gpar(cex = cex)
                          )
                 )
    count <- 1
    for (i in seq_len(layout[1]))
        for (j in seq_len(layout[2]))
            if(count <= ll) {
                pushViewport(viewport(layout.pos.row = i,
                                      layout.pos.col = j))
                pushViewport(viewport(width = 1,
                                      height = 1,
                                      default.units = if (keep_aspect_ratio) "snpc" else "npc"))
                if (inherits(l[[count]], "grob"))
                    grid.draw(l[[count]])
                else
                    if (!is.null(attr(l[[count]], "grob")))
                        grid.draw(attr(l[[count]], "grob"))
                popViewport(2)
                count <- count + 1
            }
    popViewport(2)

    ## push sub, if any
    if (!is.null(sub)) {
        pushViewport(viewport(layout.pos.row = 1 + m + s,
                              layout.pos.col = NULL))
        grid.text(sub, gp = gp_sub)
        popViewport()
    }

    popViewport(1)
}

