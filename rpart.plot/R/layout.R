# layout.R

get.layout <- function(obj, type, nn, fallen.leaves, branch,
    uniform, Margin, cex, auto.cex, compress, ycompress,
    trace, main, sub,
    node.labs, node.font, node.family, box.col, border.col,
    under.font, under.cex,
    split.labs, right.split.labs, split.cex, split.font, split.family,
    split.box.col, split.border.col,
    nspace, minbranch, node.adj, node.yshift, node.space, node.yspace,
    split.adj, split.yshift, split.space, split.yspace,
    gap, ygap, under.ygap, xcompact, ycompact, xcompact.ratio,  min.inter.height,
    max.auto.cex, min.auto.cex, ycompress.cex, accept.cex,
    shift.amounts, Fallen.yspace, bg)
{
    printf2 <- function(...) if(trace >= 2) cat0(sprintf(..., sep=""))
    printf3 <- function(...) if(trace >= 4) cat0("        ", sprintf(..., sep=""))

    merge1 <- function(vec, split.vec)
    {
        vec       <- recycle(vec, nodes)
        split.vec <- recycle(split.vec, nodes)
        split.vec[is.leaf] <- vec[is.leaf]
        split.vec
    }
    get.strheight <- function(x, y, s, Margin, cex, font, family)
    {
        init.plot(x, y, Margin, FALSE, FALSE, main, sub, 0, 1, 0, 1)
        strwidth <- my.strheight(s, cex, font, family)
        par(new=TRUE)
        strwidth
    }
    # Return the box(es) at each node.  Returns one row per box.
    # If !is.fancy, we merge each split box and corresponding node box
    # into one large box, and return nbr.nodes boxes.
    # If is.fancy, we keep each node box, its left split box, and its
    # right split box separate, and return 3 * nbr.nodes boxes

    get.combined.boxes <- function(x, y, xmax, ymax, scale, type, split.yshift)
    {
        interleave3 <- function(a, b, c) # a b c must all have the same length
        {
            x <- double(3 * length(a))
            x[seq(from=1, to=length(x), by=3)] <- a
            x[seq(from=2, to=length(x), by=3)] <- b
            x[seq(from=3, to=length(x), by=3)] <- c
            x
        }
        #--- get.combined.boxes starts here ---
        if(is.fancy(type)) {
            node.boxes <- get.boxes("default", node.labs, x, y,
                xmax, ymax, nodes, branch,
                Margin, FALSE, FALSE, main, sub, 0, 1, 0, 1,
                scale * node.cex, node.font, node.family, node.adj,
                node.yshift, box.col, border.col,
                node.space + gap/2, node.yspace + ygap/2,
                ygap, bg)

            left.boxes <- get.boxes("left", split.labs, x, y,
                xmax, ymax, nodes, branch,
                Margin, FALSE, FALSE, main, sub, 0, 1, 0, 1,
                scale * node.cex * split.cex, split.font, split.family, split.adj,
                split.yshift, split.box.col, split.border.col,
                split.space + gap/2, split.yspace + ygap/2,
                ygap, bg,
                do.init.plot=FALSE) # did init.plot in above call, so save time

            right.boxes <- get.boxes("right", right.split.labs, x, y,
                xmax, ymax, nodes, branch,
                Margin, FALSE, FALSE, main, sub, 0, 1, 0, 1,
                scale * node.cex * split.cex, split.font, split.family, split.adj,
                split.yshift, split.box.col, split.border.col,
                split.space + gap/2, split.yspace + ygap/2,
                ygap, bg,
                do.init.plot=FALSE)

            # interleave:
            # for each node, first the node box, then the left box, then the right box
            combined.boxes <- NULL
            combined.boxes$x1 <- interleave3(node.boxes$x1, left.boxes$x1, right.boxes$x1)
            combined.boxes$y1 <- interleave3(node.boxes$y1, left.boxes$y1, right.boxes$y1)
            combined.boxes$x2 <- interleave3(node.boxes$x2, left.boxes$x2, right.boxes$x2)
            combined.boxes$y2 <- interleave3(node.boxes$y2, left.boxes$y2, right.boxes$y2)
        } else {
            node.boxes <- get.boxes("default", node.labs, x, y,
                xmax, ymax, nodes, branch,
                Margin, FALSE, FALSE, main, sub, 0, 1, 0, 1,
                scale * node.cex, node.font, node.family, node.adj,
                node.yshift, box.col, border.col,
                node.space + gap/2, node.yspace + ygap/2,
                ygap, bg,
                do.init.plot=TRUE)
            split.boxes <- get.boxes(
                # extra space under split if type==TYPE.all.under so can see branch lines
                if(type == TYPE.all.under) "undersplit" else "default",
                split.labs, x, y,
                xmax, ymax, nodes, branch,
                Margin, FALSE, FALSE, main, sub, 0, 1, 0, 1,
                scale * node.cex * split.cex, split.font, split.family, split.adj,
                split.yshift, split.box.col, split.border.col,
                split.space + gap/2, split.yspace + ygap/2,
                ygap, bg,
                do.init.plot=FALSE) # did init.plot in above call, so save time

            combined.boxes <- node.boxes
            combined.boxes$x1 <- pmin(node.boxes$x1, split.boxes$x1, na.rm=TRUE)
            combined.boxes$y1 <- pmin(node.boxes$y1, split.boxes$y1, na.rm=TRUE)
            combined.boxes$x2 <- pmax(node.boxes$x2, split.boxes$x2, na.rm=TRUE)
            combined.boxes$y2 <- pmax(node.boxes$y2, split.boxes$y2, na.rm=TRUE)
        }
        # Shift and scale so the leftmost box edge is at 0, the rightmost at 1.
        # TODO calculation of min and max here is not right, must look at
        # how much edge labs jut out and use that to calculate scale.
        # (to easily reproduce, use a stump with very long box labs).
        xmin <- min(combined.boxes$x1, na.rm=TRUE)
        xmax <- max(combined.boxes$x2, na.rm=TRUE)
        new.x <- (x - xmin) / (xmax - xmin)
        delta.x <- (new.x - x)
        if(is.fancy(type)) {
            # For split boxes use the child's offset, if there is a child.
            # TODO should really take into account branch argument
            left.child <- match(nodes * 2, nodes)
            delta.left.child <- delta.x[left.child]
            no.child <- is.na(delta.left.child)
            delta.left.child[no.child] <- delta.x[no.child]

            right.child <- match(nodes * 2 + 1, nodes)
            delta.right.child <- delta.x[right.child]
            no.child <- is.na(delta.right.child)
            delta.right.child[no.child] <- delta.x[no.child]

            delta.x <- interleave3(delta.x, delta.left.child, delta.right.child)
        }
        combined.boxes$x1 <- combined.boxes$x1 + delta.x
        combined.boxes$x2 <- combined.boxes$x2 + delta.x

        # Do likewise vertically.
        ymin <- min(combined.boxes$y1, na.rm=TRUE)
        nn.height <- 0
        if(nn) # estimate height of the node number box # TODO incorporate nn.yshift, etc.?
            nn.height <- 1.2 * my.strheight("M", scale * node.cex, node.font, node.family)
        ymax <- max(combined.boxes$y2, na.rm=TRUE) + nn.height
        new.y <- (y - ymin) / (ymax - ymin)
        delta.y <- (new.y - y)
        if(is.fancy(type))
            delta.y <- interleave3(delta.y, delta.y, delta.y)
        combined.boxes$y1 <- combined.boxes$y1 + delta.y
        combined.boxes$y2 <- combined.boxes$y2 + delta.y
        list(boxes=combined.boxes, x=new.x, y=new.y)
    }
    # Get the amount needed to scale each node by to get a gap of exactly
    # "gap" between the label on the node and the label in its nearest neighbor,
    # where neighbor is the right neighbor or the neighbor below, whichever is closest.
    # Will be NA for nodes that have no neighbor.

    get.scales <- function(x, y, xmax, ymax, scale, split.yshift)
    {
        # Return the first node whose center is to the right of my center
        # and whose label vertically overlaps mine.
        # Return 0 if none, i.e. clear space to the right.
        get.right.neighbor <- function(i)
        {
            y1 <- by1[i]
            y2 <- by2[i]
            # overlaps is TRUE for all nodes that vertically overlap me
            # (usually these will be nodes on my level)
            overlaps <- (y1 >= by1 & y1 <= by2) | # -o
                        (y2 >= by1 & y2 <= by2) | # _o
                        (y1 <= by1 & y2 >= by2) | # o-
                        (y1 >= by1 & y2 <= by2)   # o_
            # exclude me and nodes whose center is to the left of my center
            overlaps[(xcenters <= xcenters[i]) | (is.na(overlaps))] <- FALSE
            # get indices of overlapping nodes
            overlaps <- ((1:length(overlaps))[overlaps])
            # of the overlapping nodes, return the one with the leftmost left edge
            if(length(overlaps) == 0) NA else overlaps[which.min(boxes$x1[overlaps])]
        }
        # Return the first node whose center is below my center
        # and whose label horizontally overlaps mine.
        # Return 0 if none, i.e. clear space below me
        get.lower.neighbor <- function(i)
        {
            x1 <- bx1[i]
            x2 <- bx2[i]
            # overlaps is TRUE for all nodes that horizontally overlap me
            overlaps <- (x1 >= bx1 & x1 <= bx2) |
                        (x2 >= bx1 & x2 <= bx2) |
                        (x1 <= bx1 & x2 >= bx2) |
                        (x1 >= bx1 & x2 <= bx2)
            # exclude me and nodes whose center is above my center
            overlaps[(ycenters >= ycenters[i]) | (is.na(overlaps))] <- FALSE
            # get indices of overlapping nodes
            overlaps <- ((1:length(overlaps))[overlaps])
            # of the overlapping nodes, return the one with the highest upper edge
            if(length(overlaps) == 0) NA else overlaps[which.max(boxes$y2[overlaps])]
        }
        #--- get.scales starts here ---
        boxes <- get.combined.boxes(x, y, xmax, ymax, scale, type, split.yshift)$boxes
        bx1 <- boxes$x1; bx2 <- boxes$x2; by1 <- boxes$y1; by2 <- boxes$y2
        xcenters <- abs(bx2 + bx1) / 2    # the centers of the boxes
        ycenters <- abs(by2 + by1) / 2
        widths2  <- abs(bx2 - bx1) / 2    # the box widths divided by 2
        heights2 <- abs(by2 - by1) / 2
        neighbors <- yneighbors <- double(length(bx1))
        for(i in 1:length(bx1)) # TODO could vectorize?
            neighbors[i] <- get.right.neighbor(i)
        xscales <- (xcenters[neighbors] - xcenters) / (widths2[neighbors] + widths2)
        for(i in 1:length(bx1))
            yneighbors[i] <- get.lower.neighbor(i)
        yscales <- (ycenters - ycenters[yneighbors]) / (heights2[yneighbors] + heights2)
        # required scale is the max of scale in horiz or vert direction
        # except that if a scale is 1 or more then we must use the minimum
        both.scales <- pmax(xscales, yscales, na.rm=TRUE)
        which <- xscales >= 1 | yscales >= 1
        which[is.na(which)] <- FALSE
        both.scales[which] <- pmin(xscales, yscales, na.rm=TRUE)[which]

        # TODO fix this, currently yneighbors not included in neighbors, so trace msg can be wrong
        #       if(trace) {
        #           which <- yscales > both.scales
        #           which[is.na(which)] <- FALSE
        #           neighbors[which] <- yneighbors[which]
        #       }

        if(any(both.scales < min.auto.cex, na.rm=TRUE))
            both.scales[both.scales < min.auto.cex] <- min.auto.cex

        if(any(is.infinite(both.scales), na.rm=TRUE)) {
            # TODO should never get here, legacy of implementation before min.auto.cex
            both.scales[is.infinite(both.scales)] <- 1
            warning0("setting infinite scales to 1")
        }
        list(scales=both.scales, neighbors=neighbors)
    }
    # Get the amount we need to scale xmax by to get a gap of
    # exactly gap between the labs of the worst two nodes.
    # (So scale could be less than or greater than 1.
    # A scale of 1 or greater means that all boxes fit
    # with at least gap space between them.)
    # This gives a result which ignores that changing scale can change
    # which nodes are neighbors (as nodes shift vertically), so
    # will need to be adjusted in get.actual.scale.

    get.scale <- function(x, y, xmax, ymax, scale, split.yshift)
    {
        temp <- get.scales(x, y, xmax, ymax, scale, split.yshift)
        imin <- which.min(temp$scales) # index of worst scale
        if(length(imin) == 0) # TODO look into this
            0
        else
            temp$scales[imin]
    }
    # get.scale gives a result which ignores that changing scale can change
    # which nodes are neighbors.
    # It also does not take into account that font sizes are discrete,
    # so the cex you get may not be the cex you asked for.
    # This function takes care of all of that, using a binary search.

    get.actual.scale <- function(x, y, split.yshift)
    {
        scale <- .8 # initial guess
        lower <- 0
        upper <- 5  # never need to scale up by more than this (used in do.xcompact)
        best.scale <- -Inf # needed because relative.scale not always monotonic with scale
        while(1) {
            relative.scale <- get.scale(x, y, xmax, ymax, scale, split.yshift)
            if(relative.scale >= 1) {    # boxes fit?
                lower <- scale
                best.scale <- max(scale, best.scale)
                if(best.scale > max.auto.cex) # TODO correct? ok for do.xcompact?
                    break           # good enough, don't need more detail
            } else {                # boxes don't fit
                upper <- scale
                if(upper <= min.auto.cex) {
                    best.scale <- min.auto.cex
                    break
                }
            }
            printf3(
"get.actual.scale: scale %4.2f relative.scale %4.2f best.scale %4.2f upper-lower %5.3f\n",
                scale, relative.scale, best.scale, upper - lower)
            if(upper - lower < .02)
                break
            scale <- (lower + upper) / 2
        }
        best.scale
    }
    # shift nodes vertically looking for a better cex
    shifter <- function(start.scale)
    {
        issue.shifter.msg <- function() # called only if trace >= 2
        {
            printf(
"Shifter: cex improvement %.3g best.shift.amount %g best.split.yshift.amount %g%s\n",
                best.scale.after.shifting / start.scale, best.shift.amount,
                best.split.yshift[2] - split.yshift[2],
                if(best.scale.after.shifting / scale >= accept.cex)
                    " (will be used)"
                else
                    " (will not be used)")

            if(best.scale.after.shifting / scale >= accept.cex) {
                improvement <- best.scale.after.shifting / scale
                if(ycompress && scale <= ycompress.cex)
                    if(auto.cex)
                        printf("ycompress increased cex by %.2f\n", improvement)
                    else
                        printf("ycompress increased available space by %.2f\n", improvement)
                else {
                    stopifnot(is.fancy(type))
                    printf("Shifting split labels increased available space by %.2f\n",
                           improvement)
                }
            }
        }
        # init is.shift (bool vec of nodes to be shifted)
        init.is.shift <- function()
        {
            is.shift <- NULL
            if(is.fancy(type)) {
                if(fallen.leaves)
                    is.shift <- add.fallen.leaves(rep(FALSE, 3 * nnodes))
            } else {
                is.shift <- nodes %% 2 == 1 # alternate all nodes
                if(fallen.leaves)
                    is.shift <- add.fallen.leaves(is.shift)
            }
            is.shift
        }
        # return is.shift but with leaves that are odd (in the
        # display) set FALSE and even leaves set TRUE
        add.fallen.leaves <- function(is.shift)
        {
            even <- FALSE
            for(i in 1:length(nodes))
                if(is.leaf[i]) {
                    is.shift[i] <- even
                    even <- !even
                }
            is.shift
        }
        get.shifted.y <- function(shift.amount, ref.shift, nnodes)
        {
            is.shift <- is.shift[1:nnodes] # needed when is.fancy
            y[is.shift] <- y[is.shift] + shift.amount * ref.shift
            y
        }
        search.for.best.shift <- function()
        {
            for(i in 1:length(shift.amounts)) {
                shift.amount <- shift.amounts[i]
                shifted.y <- get.shifted.y(shift.amount, ref.shift, nnodes)
                # check that a shift doesn't move nodes above the nodes for the level above
                if(any(shifted.y > shifted.y[iparents], na.rm=TRUE)) {
                    printf2("    Node shifter: skipping invalid      shift.amount %-4.1f\n", shift.amount)
                    next
                }
                scale.after.shifting <- get.actual.scale(x.org, shifted.y, split.yshift)
                printf2("    Node shifter: cex improvement %-5.3g shift.amount %-4.1f ",
                        scale.after.shifting / start.scale, shift.amount)
                # Note use of >= versus > below, so will use 1st shift unless actually worse.
                # We do want to move enough to allow some room for expansion, but
                # don't want to move labels too far up towards the fancy split labels.
                if((i <= 1 && scale.after.shifting >= best.scale.after.shifting) ||
                            scale.after.shifting > best.scale.after.shifting) {
                    best.scale.after.shifting <- scale.after.shifting
                    best.shift.amount <- shift.amount
                    if(trace >= 2)
                        printf("<new best")
                }
                printf2("\n")
            }
            list(best.scale.after.shifting=best.scale.after.shifting,
                 best.shift.amount=best.shift.amount)
        }
        #--- shifter starts here ---
        best.scale.after.shifting <- start.scale
        best.shift.amount <- ref.shift <- 0
        best.split.yshift <- split.yshift
        nnodes <- length(nodes)
        if(ycompress && scale <= ycompress.cex) { # try shifting nodes vertically?
            is.shift <- init.is.shift() # init is.shift (bool vec of nodes to be shifted)

            # Use current box heights as an estimate of amount to shift (in get.shifted.y).
            # This will not be correct as scales change, but is just an estimate.
            # ref.shift is the min amount we must move for any box to clear its neighbor.
            # The [1:nnodes] is necessary when is.fancy(type) and length(y) == 3 * nnodes
            # TODO Conservative approach for now, use max height of all boxes.

            boxes <- get.combined.boxes(x.org, y.org, xmax, ymax, start.scale,
                                        type, split.yshift)$boxes
            ref.shift <- max(abs(boxes$y2 - boxes$y1)[1:nnodes], na.rm=TRUE)
            if(auto.cex && (start.scale <= min.auto.cex)) {
                # Text is too small to display properly.
                # Force shift to be accepted, which usually improves legibility.
                # Note that the scale returned by get.actual.scale is clamped at
                # min.auto.cex and thus can't be used to compare different configurations.

                printf2("Shifter: forcing shift, because start.scale %.3g <= min.auto.cex %.3g\n",
                        start.scale, min.auto.cex)
                if(accept.cex > .98)
                    accept.cex <- .98 # force shift to be accepted, unless actually worse
            }
            if(!is.null(is.shift)) {
                temp <- search.for.best.shift()
                best.scale.after.shifting <- temp$best.scale.after.shifting
                best.shift.amount <- temp$best.shift.amount
                y <- get.shifted.y(best.shift.amount, ref.shift, nnodes)
            }
        }
        if(trace >= 2)
            issue.shifter.msg()
        list(y=y,
             scale=best.scale.after.shifting,
             split.yshift=best.split.yshift,
             accept.cex=accept.cex)
    }
    do.xcompact <- function()
    {
        new.scale <- get.actual.scale(x, y, split.yshift)
        xmax <- round(new.scale * xcompact.ratio - .05, 1) # round down to one digit after point
        xmax <- max(xmax, 1) # never expand horizontally
        if(xmax != 1)
            printf2("Compacted horizontally, new xlim is c(0, %.3g)\n", xmax)
        xmax
    }
    do.ycompact <- function(scale)
    {
        if(min.inter.height < .25)
            min.inter.height <- .25
        if(nn) # estimate size of node-number boxes TODO improve and move into get.combined.boxes
            min.inter.height <- min.inter.height + 1.5
        max.cex <- max(scale * merge1(node.cex, split.cex * node.cex))
        font1 <- node.font[1]
        family1 <- node.family[1]
        if(fallen.leaves)
            tree.depth[is.leaf] <- max(tree.depth)
        # look first for the ymax that causes vertical gap to be closer than min.inter.height
        ymax <- new.ymax <- 1
        while(TRUE) {
            strheight1 <- get.strheight(xmax, new.ymax, "M",
                                        Margin, max.cex, font1, family1)
            # TODO fix type1 handling, although works ok most of the time
            type1 <- type
            if(type == TYPE.fancy.all)
                type1 <- TYPE.all
            else if(type == TYPE.fancy.noall)
                type1 <- TYPE.default
            boxes <- get.combined.boxes(x.org, y.org, xmax, ymax, scale, type1, split.yshift)$boxes
            gap <- Inf # min diff between this ceiling and previous floor
            # TODO this doesn't work with fallen leaves
            for(depth in 1:max.tree.depth) # first depth is 0, start one beyond that
                gap <- min(gap, min(boxes$y1[tree.depth == depth-1], na.rm=TRUE) -
                                max(boxes$y2[tree.depth == depth], na.rm=TRUE))
            if(gap / strheight1 < min.inter.height || # gap too small?
                    new.ymax >= 5 || # compacted by more than 5?
                    get.scale(x.org, y.org, xmax, new.ymax, scale, split.yshift) < 1) # touching?
                break # gone too far
            ymax <- new.ymax
            new.ymax <- new.ymax + .1 # TODO could use an intelligent bump here for speed
        }
        ymax <- max(1, round(ymax - .05, 1)) # round down to one digit after point
        if(ymax > 1)
            printf2("Compacted vertically, new ylim is c(0, %.3g)\n", ymax)
        ymax
    }
    tree.depth <- function (nodes) # lifted from rpart::tree.depth.R
    {
        depth <- floor(log(nodes, base = 2) + 1e-7)
        depth - min(depth) # doesn't seem to need as.vector.
    }
    #--- get.layout starts here ---
    frame <- obj$frame
    nodes <- as.numeric(row.names(frame))
    is.leaf <- is.leaf(frame)
    iparents <- match(nodes %/% 2, nodes) # row index of parent node in frame
    node.cex <- 1
    auto.gap <- is.null(gap)
    if(auto.gap)
        gap <- ygap <- .5
    node.xy <- get.node.coords(obj, uniform, branch, compress,
                               nspace, minbranch, fallen.leaves, Fallen.yspace)
    if(length(node.xy$x) == 1) {
        # tree is just a root
        return(list(node.xy=list(x=.5, y=.5), cex=cex, # NOTE: return
                    xlim=c(0, 1),
                    ylim=c(0, 1),
                    split.yshift=0,
                    gap=gap, ygap=ygap,
                    node.space=node.space,
                    node.yspace=node.yspace))
    }
    x.org <- x <- node.xy$x
    y.org <- y <- node.xy$y
    tree.depth <- tree.depth(nodes)
    max.tree.depth <- max(tree.depth)
    scale <- ymax <- xmax <- 1
    scale <- get.actual.scale(x.org, y.org, split.yshift)

    printf2("Initial scale %.3g\n", scale)
    shifted <- shifter(scale)
    if(shifted$scale / scale >= shifted$accept.cex) {
        y <- shifted$y
        scale <- shifted$scale
        split.yshift <- shifted$split.yshift
    }
    scale.before.clip <- scale
    if(auto.cex && scale > max.auto.cex) {
        printf2("Clipped scale %.3g to max.auto.cex %.3g\n", scale, max.auto.cex)
        scale <- max.auto.cex
    }
    # get node xy taking into account xmin, xmax etc. adjustment for labs on edges
    xy <- get.combined.boxes(x.org, y, xmax, ymax, scale, type, split.yshift)
    x <- xy$x
    y <- xy$y
    if(scale <= min.auto.cex) {
        warning0(sprintf(
            "labs do not fit even at cex %.3g, there may be some overplotting",
            min.auto.cex))
        scale <- min.auto.cex
    }
    if(auto.cex) { # TODO change compact code to deal with manual cex too
        if(xcompact && scale.before.clip >= 1)
            xmax <- do.xcompact()
        if(ycompact)
            ymax <- do.ycompact(scale)
    }
    if(auto.cex)
        cex  <- scale * cex
    list(node.xy=list(x=x, y=y), cex=cex,
         xlim=c(-xmax / 2 + .5, xmax / 2 + .5), # center graph horizontally
         ylim=c(-ymax / 2 + .5, ymax / 2 + .5), # center graph vertically
         split.yshift=split.yshift,
         gap=gap, ygap=ygap,
         node.space=node.space,
         node.yspace=node.yspace)
}
