# yes.no.R

draw.yes.no <- function(type, draw.shadow,
    xflip, left, branch, xlim, ylim, node.xy, lwd,
    yesno.yshift, split.boxes, split.cex, split.box.col, split.border.col,
    split.shadow.col, split.shadow.offset,
    nn.cex, nn.font, nn.family, nn.col, nn.box.col,
    nn.border.col, nn.lty, nn.round, bg)
{
    if(is.null(nn.cex)) # auto cex?
        nn.cex <- .8 * min(split.cex)
    # 1.5 for white space aound the text
    height1 <- 1.5 * my.strheight("yes", nn.cex, nn.font, nn.family)
    width1  <- 1.5 * my.strwidth("yes", nn.cex, nn.font, nn.family)
    xleft  <- split.boxes$x1[1] # horiz posn of left text
    xright <- split.boxes$x2[1] # horiz posn of right text

    # get vertical position of text
    y <- if(type == TYPE.all.under)
            .5 * (split.boxes$y1[1] + split.boxes$y2[1]) # center of split box
        else if(branch >= .5)   # just below center of top split box
            .5 * (split.boxes$y1[1] + split.boxes$y2[1]) - .4 * height1
        else                    # near the bottom of the top box
            max(split.boxes$y1[1], node.xy$y[1])

    y <- y + yesno.yshift * height1

    # draw the boxes

    round <-  nn.round[1] * height1
    border.col <- if(is.box.invisible(split.box.col[1], split.border.col[1], bg))
                      nn.border.col
                  else
                      c(bg, bg)

    if(is.invisible(border.col, bg) && is.invisible(split.shadow.col, bg)) { # no border?
        split.shadow.col <- 0    # no shadow if no border
        height1 <- height1 / 1.3 # just small space around box so branch lines not blocked
        width1  <- width1 / 1.3
    }
    if(draw.shadow) {
        if(!is.invisible(split.shadow.col, bg)) {
            draw.shadow(xleft - 1.2 * width1, y - .6 * height1,
                        xleft - .2 * width1,  y + .6 * height1,
                        xlim, ylim, round,
                        split.shadow.col[1], split.shadow.offset[1])

            draw.shadow(xright + .2 * width1,  y - .6 * height1,
                        xright + 1.2 * width1, y + .6 * height1,
                        xlim, ylim, round,
                        split.shadow.col[1], split.shadow.offset[1])
        }
    } else {
        rounded.rect(xleft - 1.2 * width1, y - .6 * height1,
                     xleft - .2 * width1,  y + .6 * height1,
                     xlim, ylim, round, nn.box.col[1], border.col[1],
                     nn.lty[1], lwd[1])

        rounded.rect(xright + .2 * width1,  y - .6 * height1,
                     xright + 1.2 * width1, y + .6 * height1,
                     xlim, ylim, round, nn.box.col[1], border.col[1],
                     nn.lty[1], lwd[1])
    }
    # draw the text

    xcenter <- (xleft - 1.2 * width1 + xleft - .2 * width1) / 2
    if(xflip == left) {
        left.msg <- "no"
        right.msg <- "yes"
    } else {
        left.msg <- "yes"
        right.msg <- "no"
    }
    # The .15 height adjustment is because the driver centers the "yes" accounting for
    # the descender on "y", so must adjust to get "yes" and "no" at approx. same height.
    height <- my.strheight("no", nn.cex, nn.font, nn.family)
    text(xcenter, y - .15 * height, left.msg,
         cex=nn.cex[1], font=nn.font[1], family=nn.family[1], col=nn.col[1])
    xcenter <- (xright + .2 * width1 + xright + 1.2 * width1) / 2
    text(xcenter, y, right.msg,
         cex=nn.cex[1], font=nn.font[1], family=nn.family, col=nn.col[1])
}
