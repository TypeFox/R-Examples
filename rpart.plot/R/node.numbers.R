# node.numbers.R

draw.node.numbers <- function(nn, ni, draw.shadow, type, branch,
    Margin, xflip, yflip, cex,
    main, sub, col.main, cex.main, col.sub, cex.sub,
    xlim, ylim, node.xy, is.leaf, nodes,
    node.labs, font,  family, box.col, border.col, shadow.col,
    under.cex, under.font, under.ygap, ygap,
    split.labs, split.cex, split.font, split.family, split.box.col,
    split.border.col, split.shadow.col,
    nn.cex, nn.font, nn.family, nn.col, nn.box.col,
    nn.border.col, nn.lty, nn.lwd, nn.round,
    split.adj, split.space, split.yspace, split.yshift,
    yshift, adj, space, yspace, shadow.offset,
    nn.adj, nn.yshift, nn.space, nn.yspace, bg)
{
    merge1 <- function(vec, split.vec)
    {
        vec       <- recycle(vec, nodes)
        split.vec <- recycle(split.vec, nodes)
        split.vec[is.leaf] <- vec[is.leaf]
        split.vec
    }
    #--- draw.node.numbers starts here ---
    # The node numbers go on top of the node box or split box, whichever
    # is higher.  So we have to get the positions of those boxes first.
    if(is.fancy(type) || type == TYPE.all.under) {
        all.labs   <- node.labs
        all.yshift <- yshift
        all.cex    <- cex
        all.font   <- font
        all.family <- family
        all.adj    <- adj
        all.space  <- space
        all.yspace  <- yspace
    } else {
        all.labs   <- merge1(node.labs, split.labs)
        all.yshift <- merge1(yshift, split.yshift)
        all.cex    <- merge1(cex, split.cex)
        all.font   <- merge1(font, split.font)
        all.family <- merge1(family, split.family)
        all.adj    <- merge1(adj, split.adj)
        all.space  <- merge1(space, split.space)
        all.yspace <- merge1(yspace, split.yspace)
    }
    if(!is.invisible(shadow.col, bg)) {
        # following prevents shadows on node numbers when not on their split boxes
        want.interior.node.shadows <-
            identical(shadow.col, split.shadow.col) || type == TYPE.fancy.all || type == TYPE.all.under
        shadow.col <- recycle(shadow.col, nodes)
        if(!want.interior.node.shadows)
            shadow.col[!is.leaf] <- if(is.character(col)) par("bg") else 0
    }
    box <- get.boxes("default", all.labs, node.xy$x, node.xy$y,
                     xlim, ylim, nodes, branch,
                     Margin, xflip, yflip, main, sub,
                     col.main, cex.main, col.sub, cex.sub,
                     all.cex, all.font, all.family, all.adj,
                     all.yshift, box.col, border.col,
                     all.space, all.yspace,
                     ygap, bg,
                     do.init.plot=FALSE)

    if(is.null(nn.cex)) # auto cex?
        nn.cex <- .7 * min(all.cex)
    strheight1 <- my.strheight("M", cex=nn.cex, font=nn.font, family=nn.family)
    x <- (box$x2 + box$x1) / 2
    if(ni) # ni's look better with a bit more whitespace
        nn.yspace <- 1.5 * nn.yspace
    nn.yshift <- nn.yshift - nn.yspace / 4 # TODO revisit but works for default case
    y <- box$y2 + (nn.yshift + nn.yspace / 2) * strheight1 + strheight1 / 2
    if(ni && nn)
        nodes <- sprintf("[%g] %g", 1:length(nodes), nodes)
    else if(ni)
        nodes <- sprintf("[%g]", 1:length(nodes))
    else
        nodes <- sprintf("%g", nodes)
    nn.adj    <- recycle(nn.adj, nodes)
    nn.yshift <- recycle(nn.yshift, nodes)
    nn.space  <- recycle(nn.space, nodes)
    nn.yspace <- recycle(nn.yspace, nodes)
    boxes <- draw.boxes("default", draw.shadow,
               nodes, list(x=x, y=y),
               xlim, ylim, nodes, branch,
               Margin, xflip, yflip, main, sub,
               col.main, cex.main, col.sub, cex.sub,
               nn.cex, nn.font, nn.family, nn.adj, nn.yshift,
               nn.box.col, nn.border.col,
               nn.lty, nn.lwd,
               nn.space, nn.yspace, nn.round * strheight1,
               under.cex, under.font, under.ygap, ygap,
               shadow.col, shadow.offset, bg)
    xy <- get.box.centers(boxes)
    text(xy$x, xy$y,
         nodes, cex=nn.cex, font=nn.font, family=nn.family, col=nn.col)
}
