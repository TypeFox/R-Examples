# branches.R

# get.branches: Each column in the returned branch.xy is the line
#               from a node to its parent.
#
# This function is similar to rpart.branch, but the returned arrays have
# nrow(frame) columns and are indexed on rows in frame --- so indexing of
# prp's branch.col etc. args is same as col, box.col, etc. args.
#
# An example branch.xy:
#
#    $x
#         [,1]  [,2]   [,3]  [,4]  [,5]  [,6]   # index in frame
# [1,]   0.432 0.127 0.0524 0.202 0.152 0.252   # bottom of branch line
# [2,]      NA 0.127 0.0524 0.202 0.152 0.252   # shoulder of branch line
# [3,]      NA 0.432 0.1271 0.127 0.202 0.202   # top of branch line
#
#    $y
#        [,1]  [,2]   [,3]  [,4]  [,5]   [,6]
# [1,]  0.974 0.786 0.0211 0.746 0.095 0.0211
# [2,]     NA 0.974 0.7858 0.786 0.746 0.7456
# [3,]     NA 0.974 0.7858 0.786 0.746 0.7456

get.branches <- function(x, y, nodes, branch)
{
    parent <- match(nodes %/% 2, nodes)
    sibling <- match(ifelse(is.left(nodes), nodes+1, nodes-1), nodes)
    dist.to.sibling <- x[sibling] - x
    list(x=rbind(x, x + (1-branch) * dist.to.sibling / 2, x[parent]),
         y=rbind(y, y[parent],                            y[parent]))
}
get.branch.widths <- function(obj, branch.type)
{
    if(is.function(branch.type)) {
        branch.type <- check.func.args(branch.type, "branch.type", function(x) NA)
        branch.type(x=obj)
    } else {
        frame <- obj$frame
        width <- switch(branch.type,
            frame$dev,                      # 1
            sqrt(frame$dev),                # 2
            frame$dev / frame$wt,           # 3
            sqrt(frame$dev / frame$wt),     # 4
            frame$wt,                       # 5
            frame$complexity,               # 6
            abs(frame$yval),                # 7
            frame$yval - min(frame$yval),   # 8
            rep(1, nrow(frame)))            # 9

        if(is.null(width)) # needed because switch has no default mechanism
            stop0("branch.type=", branch.type, " is illegal")

        width
    }
}
draw.branches <- function(obj, branch.type, branch.col,
        branch.lty, branch.lwd, branch.fill, branch.tweak,
        node.labs, split.labs, node.xy, strheight,
        type, branch, xflip, yflip, Margin, space, yspace,
        cex, font, family, adj, box.col, border.col,
        under.cex, under.font, under.ygap,
        split.cex, split.font, split.family, split.adj, split.yshift,
        split.box.col, split.border.col, split.space, split.yspace,
        main, sub, col.main, cex.main, col.sub, cex.sub,
        xlim, ylim, yshift, ygap, bg,
        min.branch.width)
{
    check.returned.width <- function(width)
    {
        print.and.stop <- function(...)
        {
            cat("width:\n")
            print(width)
            cat("\n")
            stop0("the call to the branch.type function returned a bad result: ", ...)
        }
        if(length(width) == 0)
            print.and.stop("length(width) == 0")
        if(!is.numeric(width))
            print.and.stop("!is.numeric(width)")
        if(any(width < 0))
          print.and.stop("widths less than zero")
        if(length(width) != nrow(obj$frame))
            print.and.stop("\nthe length ", length(width),
                " of the returned width is not equal to the number of rows in frame ",
                nrow(obj$frame))
    }
    get.branch.widths.wrapper <- function()
    {
        width <- get.branch.widths(obj, branch.type)
        check.returned.width(width)

        # normalize branch widths so widest width is one fifth the plot width
        sibling <- match(ifelse(is.left(nodes), nodes+1, nodes-1), nodes)
        max.width <- max(width + width[sibling], na.rm=TRUE)
        xrange <- xlim[2] - xlim[1]
        width <- branch.tweak * width * .2 * xrange / max.width

        # clamp at min.branch.width
        min.width <- min.branch.width * xrange
        width[width < min.width] <- min.width

        width
    }
    # Adjust branch vertical position for the actual positions of the boxes, and
    # for branch.type.  Several somewhat arbitrary aesthetic judgments here.
    get.branch.y <- function(y)
    {
        stopifnot(branch.type == 0)
        if(type == TYPE.default) {
            y[1,] <- .667 * split.boxes$y1 + .333 * split.boxes$y2
            y[1, is.leaf] <- node.boxes$y2[is.leaf]
        } else if(type == TYPE.all) {
            y[1,] <- .667 * split.boxes$y1 + .333 * split.boxes$y2
            y[1,is.leaf] <- node.boxes$y2[is.leaf]
            y[2,] <- y[3,] <- ((node.boxes$y1 + node.boxes$y2) / 2)[parent]
        } else if(type == TYPE.all.under) {
            y[1,] <- node.boxes$y2
            y[2,] <- y[3,] <- ((split.boxes$y1 + split.boxes$y2) / 2)[parent]
        } else if (type == TYPE.fancy.noall) {
            y[1, is.leaf] <- node.boxes$y2[is.leaf]
        } else if (type == TYPE.fancy.all) {
            y[1,] <- node.boxes$y2
            y[2,] <- y[3,] <- ((node.boxes$y1 + node.boxes$y2) / 2)[parent]
        } else
            stop("illegal type ", type) # programming error

        # ensure the bottom of branch is below top (want no upward sloping branches)
        y[1,] <- pmin(y[1,], y[2,], na.rm=TRUE)

        y
    }
    # return the branch shapes and a modified node.xy
    # shape rows:  1=bot_left  2=top_left  3=top_right  4=bot_right
    get.wide.branches <- function(branch.xy)
    {
        stopifnot(wide.branches)
        x <- branch.xy$x
        y <- branch.xy$y

        # set top of branch line to bottom of parent box
        # and bottom of branch line to top of box
        y[1,] <- pmax(split.boxes$y2, node.boxes$y2, na.rm=TRUE)
        y[2,] <- y[3,] <- pmin(split.boxes$y1, node.boxes$y1, na.rm=TRUE)[parent]

        shape.x <- shape.y <- matrix(nrow=4, ncol=ncol(x))
        width <- get.branch.widths.wrapper()
        is.left <- is.left(nodes)
        sibling <- match(ifelse(is.left, nodes+1, nodes-1), nodes)
        stopifnot(branch == 0 || branch == 1)
        if(branch == 1) { # square shoulders?
            shape.x[1,] <- shape.x[2,] <- x[1,] - width / 2     # tl bl
            shape.x[3,] <- shape.x[4,] <- x[1,] + width / 2     # tr br
        } else {          # v shaped branch lines
            shape.x[1,] <- x[1,] - width / 2    # bl
            shape.x[4,] <- x[1,] + width / 2    # br
            for(i in 1:nrow(frame)) {
                if(is.left[i]) { # left branch
                    shape.x[2, i] <- x[2,i] - (width[i] + width[sibling[i]]) / 2 # tl
                    shape.x[3, i] <- shape.x[2, i] + width[i]                    # tr
                } else {         # right branch
                    shape.x[3, i] <- x[2,i] + (width[i] + width[sibling[i]]) / 2 # tr
                    shape.x[2, i] <- shape.x[3, i] - width[i]                    # tl
                }
            }
        }
        shape.y[1,] <- shape.y[4,] <- y[1,]
        shape.y[2,] <- shape.y[3,] <- y[2,]
        x[1, is.left]  <- shape.x[2, is.left]
        x[1, !is.left] <- shape.x[3, !is.left]
        x[2,] <- x[3,]
        x[3,] <- NA
        y[1,] <- shape.y[2,]
        y[2,] <- y[3,]
        y[3,] <- NA
        if(branch == 1) { # square shoulders?
            # Put top of branch a fraction below the bottom of the box,
            # else slight overlap of (possibly invisible) box blanks out
            # parts of the horizontal branch, looks odd.
            fudge <- 0.003 * (ylim[2] - ylim[1]) # TODO .003 determined empirically
            shape.y[2,] <- shape.y[2,] - fudge
            shape.y[3,] <- shape.y[3,] - fudge
            y[1,] <- y[1,] - fudge
            y[2,] <- y[2,] - fudge
        }
        list(branch.xy=list(x=x, y=y), shape.xy=list(x=shape.x, y=shape.y))
    }
    draw.branch.lines <- function(branch.xy, branch.col, branch.lty, branch.lwd)
    {
        # lines() doesn't recycle the matrix branch.col the way we want, so need loop
        # lend=1 else line ends visible when branch.type != 0
        for(i in 1:length(nodes))
            lines(branch.xy$x[,i], branch.xy$y[,i],
                  col=branch.col[i], lty=branch.lty[i], lwd=branch.lwd[i], lend=1)
    }
    #--- draw.branches starts here ---
    frame <- obj$frame
    nodes <- as.numeric(row.names(frame))
    parent <- match(nodes %/% 2, nodes) # index of node's parent in frame
    is.leaf <- is.leaf(obj$frame)
    stopifnot(is.function(branch.type) ||
                (is.numeric(branch.type) && length(branch.type) == 1))
    wide.branches <- !is.na.or.zero(branch.type)
    if(wide.branches &&
            type != TYPE.default && type != TYPE.all && type != TYPE.all.under)
        stop0("branch.type=", branch.type, " is not yet supported with type=", type)

    branch.xy <- get.branches(node.xy$x, node.xy$y, nodes, branch)
    if(length(nodes) == 1)  # single node tree?
        return(branch.xy)   # NOTE: return

    branch.col  <- recycle(branch.col, nodes)
    branch.lty  <- recycle(branch.lty, nodes)
    branch.lwd  <- recycle(branch.lwd, nodes)
    branch.fill <- recycle(branch.fill, nodes)

    # we need the boxes for adjusting the branch positions
    node.boxes <- get.boxes("default",
        node.labs, node.xy$x, node.xy$y, xlim, ylim, nodes, branch,
        Margin, xflip, yflip, main, sub, col.main, cex.main, col.sub, cex.sub,
        cex, font, family, adj, yshift, box.col, border.col, space, yspace,
        ygap, bg, do.init.plot=FALSE, box.around.all.text=wide.branches)

    split.boxes <- get.boxes("default",
        split.labs, node.xy$x, node.xy$y, xlim, ylim, nodes, branch,
        Margin, xflip, yflip, main, sub, col.main, cex.main, col.sub, cex.sub,
        split.cex * cex, split.font, split.family, split.adj, split.yshift,
        split.box.col, split.border.col, split.space, split.yspace,
        ygap, bg, do.init.plot=FALSE, box.around.all.text=wide.branches)

    if(!wide.branches) {
        branch.xy$y <- get.branch.y(branch.xy$y)
        draw.branch.lines(branch.xy, branch.col, branch.lty, branch.lwd)
    } else {
        temp <- get.wide.branches(branch.xy)
        draw.branch.lines(temp$branch.xy, branch.col, branch.lty, branch.lwd)

        # Omit shape border lines, they artificially expand the size of the polygon
        # slightly, messing up the proportions especially for small widths.
        # But can only do that if border col is same as fill col,
        # else user won't see the border col she specified.
        if(branch.col[2] == branch.fill[2])
            branch.col <- recycle(NA, nodes)

        # draw the branch shape
        for(i in 1:length(nodes))
            polygon(temp$shape.xy$x[,i], temp$shape.xy$y[,i],
                    col=branch.fill[i], border=branch.col[i],
                    lty=branch.lty[i], lwd=branch.lwd[i])
    }
    if(type == TYPE.fancy.noall) { # draw small vertical line at top split?
        lines(c(node.xy$x[1], node.xy$x[1]),
              c(node.xy$y[1], node.xy$y[1] + strheight),
              col=branch.col[1], lty=branch.lty[1], lwd=branch.lwd[1])
    }
    branch.xy
}
