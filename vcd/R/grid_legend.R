grid_legend <- function (x, y, pch = NA, col = par('col'), labels, frame = TRUE, hgap = unit(0.8, "lines"), vgap = unit(0.8, "lines"), default_units = "lines", gp = gpar(), draw = TRUE, title = NULL, just = 'center', lwd = NA, lty = NA, size = 1, gp_title = NULL, gp_labels = NULL, gp_frame = gpar(fill = "transparent"), inset = c(0, 0))

{
    inset <- rep(inset, length.out = 2)

    if((length(x) > 1) && missing(y)) {
        y <- x[2]
        x <- x[1]
    }

    if(is.character(x))
        switch(x,
               left = {x = unit(0 + inset[1],'npc'); y = unit(0.5 + inset[2],'npc'); just = c("left","center")},
               topleft = {x = unit(0 + inset[1],'npc'); y = unit(1 - inset[2],'npc'); just = c(0,1)},
               top = {x = unit(0.5 + inset[1],'npc'); y = unit(1 - inset[2],'npc'); just = c("center", "top")},
               topright = {x = unit(1 - inset[1],'npc'); y = unit(1 - inset[2],'npc'); just = c(1,1)},
               center = {x = unit(0.5 + inset[1],'npc'); y = unit(0.5 + inset[2],'npc'); just = c("center","center")},
               bottom = {x = unit(0.5 - inset[1],'npc'); y = unit(0 + inset[2],'npc'); just = c("center","bottom")},
               bottomright = {x = unit(1 - inset[1],'npc'); y = unit(0 + inset[2],'npc'); just = c(1,0)},
               right = {x = unit(1 - inset[1],'npc'); y = unit(0.5 + inset[2],'npc'); just = c("right","center")},
               bottomleft = {x = unit(0 + inset[1],'npc'); y = unit(0 + inset[2],'npc'); just = c(0,0)})

    labels <- as.character(labels)
    nlabs <- length(labels)

    if(length(pch) == 1)
        pch <- rep(pch, nlabs)
    if(length(lwd) == 1)
        lwd <- rep(lwd, nlabs)
    if(length(lty) == 1)
        lty <- rep(lty, nlabs)
    if(length(col) == 1)
        col <- rep(col, nlabs)
    if(length(gp_labels) == 1)
        gp_labels <- rep(list(gp_labels), nlabs)


    if (is.logical(title) && !title)
        title <- NULL
    if(is.null(title))
        tit <- 0
    else
        tit <- 1

    if (!is.unit(hgap))
        hgap <- unit(hgap, default_units)
    if (length(hgap) != 1)
        stop("hgap must be single unit")
    if (!is.unit(vgap))
        vgap <- unit(vgap, default_units)
    if (length(vgap) != 1)
        stop("vgap must be single unit")

    if(tit)
        legend.layout <- grid.layout(nlabs + tit, 3,
                                     widths = unit.c(unit(2, "lines"),
                                     max(unit(rep(1, nlabs), "strwidth", as.list(c(labels))),
                                         unit(1, "strwidth", title) - unit(2, "lines")), hgap),
                                     heights = unit.pmax(unit(1, "lines"),
                                     vgap + unit(rep(1, nlabs + tit ),
                                                 "strheight", as.list(c(labels,title)))))
    else
        legend.layout <- grid.layout(nlabs, 3,
                                     widths = unit.c(unit(2, "lines"),
                                     max(unit(rep(1, nlabs), "strwidth", as.list(labels))), hgap),
                                     heights = unit.pmax(unit(1, "lines"),
                                     vgap + unit(rep(1, nlabs), "strheight", as.list(labels))))

    fg <- frameGrob(layout = legend.layout, gp = gp)

    if (frame)
        fg <- placeGrob(fg, rectGrob(gp = gp_frame))

    if (tit)
        fg <- placeGrob(fg, textGrob(title, x = .2, y = 0.5, just = c("left", "center"), gp = gp_title), col = 1, row = 1)

    for (i in 1:nlabs) {
        if(!is.na(pch[i]))
            fg <- placeGrob(fg, pointsGrob(0.5, 0.5, pch = pch[i],
                                           size = unit(size, "char"),
                                           gp = gpar(col = col[i])), col = 1, row = i + tit)
        else if(!is.na(lwd[i]) || !is.na(lty[i]))
            fg <- placeGrob(fg, linesGrob( unit(c(0.2, .8), "npc"),  unit(c(.5), "npc"),
                                          gp = gpar(col = col[i], lwd = lwd[i], lty=lty[i])), col = 1, row = i + tit)

        fg <- placeGrob(fg, textGrob(labels[i], x = .1, y = 0.5, just = c("left", "center"), gp = gp_labels[[i]]), col = 2, row = i + tit)
    }

    pushViewport(viewport(x, y, height = grobHeight(fg), width = grobWidth(fg), just = just ))

    if (draw)
        grid.draw(fg)
    popViewport(1)
    invisible(fg)
}

