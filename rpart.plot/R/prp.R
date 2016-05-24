# prp.R: expanded version of plot.rpart and related functions
#
# This is derived from the rpart code by written by
#   Terry M Therneau and Beth Atkinson:
#   http://mayoresearch.mayo.edu/mayo/research/biostat/splusfunctions.cfm
# and the R port and modifications of that code by Brian Ripley:
#   www.stats.ox.ac.uk/~ripley/
#
# Stephen Milborrow, Gardens, Dec 2010
#------------------------------------------------------------------------------
#
# TODO
#
# Much could be cleaned up.
# Revisit the ycompress algorithm.  Rework, or at least tidy up.
# Incorporate node number boxes into graph layout algorithm.
# Add code to undo unnecessary shifts when ycompress=TRUE.
# Add newdata arg or similar, allowing stats be displayed for test data and explain
#    how to use the predicted values at a node, as per Hastie et al. ESL Figure 9.5.
# Test device.supports.alpha on linux boxes. How best to determine device.supports.alpha?
# Extra stuff as in party::plot.mob and Zhan and Singer's book?
#    For example, if fallen.leaves==2 then call earth::plotd or similar for each leaf.
# Support for ellipses like text.rpart fancy?
# Allow lty etc. to optionally index on node numbers instead of indices in frame?
# Top of plot is slightly off the screen with is.fancy, all=F (need split.yshift=-2 too?)
# shift.amounts was empirically determined, could be improved?
# What is required of the obj to be displayable by prp? inherits(x, "rpart") is too restrictive?
# Allow plotting of selected tree(s) in randomForest, and objects produced by tree().
# Allow use of par srt?
# Allow other values for uniform, like branch.type.
# Option to print deviance and normalized deviance.
# Option to use expression in split.labs as requested by Jeffrey Evans.
#
# Code assumes that combined split and label boxes with type=1 are rectangular.
# This can causes smaller text than necessary esp with unif=F.
# To reproduce (note no touching boxes even with tweak):
#   a <- rpart(survived~.,data=ptitanic, method="anova") # anova for short yval strings
#   pdf(); par(mfrow=c(1,2))
#   prp(a, type=1, unif=F, boxes.include=T, split.border=1)
#   prp(a, type=1, unif=F, boxes.include=T, split.border=1, tweak=1.3, main="tweak=1.3")
#   dev.off()
#
#------------------------------------------------------------------------------

TYPE.default     <- 0 # allowable values of prp's type argument
TYPE.all         <- 1
TYPE.all.under   <- 2
TYPE.fancy.noall <- 3
TYPE.fancy.all   <- 4

rpart.plot <- function(x=stop("no 'x' arg"),
    type=0, extra=0, under=FALSE, clip.right.labs=TRUE,
    fallen.leaves=FALSE, branch=if(fallen.leaves) 1 else .2,
    uniform=TRUE,
    digits=2, varlen=-8, faclen=3,
    cex=NULL, tweak=1,
    compress=TRUE, ycompress=uniform,
    snip=FALSE,
    ...)
{
    prp(x,
        type=type, extra=extra, under=under, clip.right.labs=clip.right.labs,
        fallen.leaves=fallen.leaves, branch=branch,
        uniform=uniform,
        digits=digits, varlen=varlen, faclen=faclen,
        cex=cex, tweak=tweak,
        compress=compress, ycompress=ycompress,
        snip=snip,
        ...)
}

prp <- function(x=stop("no 'x' arg"),
    type=0, extra=0, under=FALSE, clip.right.labs=TRUE,
    nn=FALSE, ni=FALSE, yesno=TRUE,
    fallen.leaves=FALSE, branch=if(fallen.leaves) 1 else .2,
    uniform=TRUE, left=TRUE, xflip=FALSE, yflip=FALSE, Margin=0, space=1, gap=NULL,
    digits=2, varlen=-8, faclen=3,
    cex=NULL, tweak=1,
    compress=TRUE, ycompress=uniform,
    trace=FALSE, snip=FALSE, snip.fun=NULL,

    box.col=0, border.col=col,
    round=NULL, leaf.round=NULL,
    shadow.col=0, prefix="", suffix="", xsep=NULL,

    under.font=font, under.col=1, under.cex=.8,

    split.cex=1, split.font=2, split.family=family, split.col=1,
    split.box.col=0, split.border.col=0,
    split.lty=1, split.lwd=NULL, split.round=0,
    split.shadow.col=0,
    split.prefix="", right.split.prefix=NULL,
    split.suffix="", right.split.suffix=NULL,
    facsep=",", eq=" = ", lt=" < ", ge=" >= ",

    branch.col=if(identical(branch.type, 0)) 1 else "gray",
    branch.lty=1, branch.lwd=NULL,
    branch.type=0, branch.tweak=1,
    min.branch.width=.002, branch.fill=branch.col,

    nn.cex=NULL, nn.font=3, nn.family="", nn.col=1,
    nn.box.col=0, nn.border.col=nn.col,
    nn.lty=1, nn.lwd=NULL, nn.round=.3,

    node.fun=internal.node.labs,
    split.fun=internal.split.labs,
    FUN="text",

    nspace=branch, minbranch=.3, do.par=TRUE,
    add.labs=TRUE, clip.left.labs=FALSE, fam.main="",
    yshift=0, yspace=space, shadow.offset=.4,

    split.adj=NULL, split.yshift=0, split.space=space,
    split.yspace=yspace, split.shadow.offset=shadow.offset,

    nn.adj=.5, nn.yshift=0, nn.space=.8, nn.yspace=.5,

    ygap=gap/2, under.ygap=.5, yesno.yshift=0,
    xcompact=TRUE, ycompact=uniform, xcompact.ratio=.8, min.inter.height=4,

    max.auto.cex=1, min.auto.cex=.15, ycompress.cex=.7, accept.cex=1.1,
    shift.amounts=c(1.5, 2),
    Fallen.yspace=.1,
    boxes.include.gap=FALSE,
    ...)
{
    check.dots <- function(dots) # check dots arguments, if any
    {
        legal.dots.args <- # they are legal if we have code to process them later
            c("adj", "cex.main", "cex.sub", "col", "col.main", "col.sub", "family",
              "font", "lty", "lwd", "main", "mar", "sub", "xlim", "xpd", "ylim")

        if(length(dots) > 0) {
            names <- names(dots)
            pmatch <- pmatch(names, legal.dots.args, duplicates.ok=TRUE)
            if(any(is.na(pmatch))) {
                # report the first illegal arg
                ibad <- (1:length(dots))[is.na(pmatch)]
                stop0("prp: illegal argument \"", names[ibad][1], "\"")
            }
            duplicated <- duplicated(pmatch)
            if(any(duplicated))
                stop0("prp: duplicated argument \"", names[duplicated][1], "\"")
        }
    }
    merge1 <- function(vec, split.vec)
    {
        split.vec <- recycle(split.vec, nodes)
        split.vec[is.leaf] <- recycle(vec, nodes)[is.leaf]
        split.vec
    }
    draw.labs <- function(draw.shadows1, draw.split.shadows1)
    {
        # put the labels on the screen, text after \n\n (if any) goes under the box
        draw.labs1 <- function(labs, boxes, yspace, cex, font, family, col,
                               draw.shadows1, make.space.for.shadows, shadow.col, round)
        {
            draw.under.text <- function() { # draw the text under the box, and its shadow
                height1 <- my.strheight("M", cex, font, family)
                cex <- under.cex * cex
                under.height <- my.strheight(sep.labs$under.box, cex, under.font, family)
                x <- xy$x
                y <- boxes$y1 - under.ygap * height1 - .5 * under.height
                width  <- .5 * my.strwidth(sep.labs$under.box,  cex, under.font, family)
                height <- .5 * my.strheight(sep.labs$under.box, cex, under.font, family)
                # the magic numbers 1.4 and 1.2 seem about right visually
                if(make.space.for.shadows)
                    height <- 1.4 * height
                if(draw.shadows1)
                    draw.shadow(x - 1.2 * width, y - height,
                         x + 1.2 * width, y + height,
                         xlim, ylim, 0, shadow.col, shadow.offset)
                else {
                    # draw a white box to hide anything underneath, then draw the text
                    rect(x - 1.2 * width, y - height,
                         x + 1.2 * width, y + height, col=bg, border=0)
                    FUN(x, y, sep.labs$under.box,
                        cex=cex, font=under.font, family=family, col=under.col)
                }
            }
            #--- draw.labs1 starts here ---
            FUN <- check.func.args(FUN, "FUN argument to the prp", graphics::text)
            sep.labs <- separate.labs(labs)
            xy <- get.box.centers(boxes)
            # draw the text after \n\n if any under the box
            if(!all(nchar(sep.labs$under.box) == 0))
                draw.under.text()
            # draw the text before \n\n in the box
            if(!draw.shadows1)
                FUN(xy$x, xy$y, sep.labs$in.box,
                    cex=cex, font=font, family=family, col=col)
        }
        #--- draw.labs starts here ---
        if(boxes.include.gap) {
            # For debugging: make get.boxes expand the boxes to include
            # what would normally be the gap between the boxes.
            # With optimum cex, at least one pair of boxes will just touch.
            printf("boxes.include.gap is TRUE\n")
            split.space <- split.space + gap/2
            split.yspace <- split.yspace + ygap/2
            space <- space + gap/2
            yspace <- yspace + ygap/2
            gap <- ygap <- 0
        }
        # With type==TYPE.all.under and no visible split box, branch lines
        # look better if just a small space around labs.
        small.underspace <- type == TYPE.all.under &&
            is.box.invisible(split.box.col, split.border.col, bg)

        split.boxes <-
            draw.boxes(if(is.fancy(type)) "left" else "default", draw.split.shadows1,
                   split.labs, node.xy, xlim, ylim,
                   nodes, branch,
                   Margin, xflip, yflip, main, sub,
                   col.main, cex.main, col.sub, cex.sub,
                   split.cex * cex, split.font, split.family, split.adj, split.yshift,
                   split.box.col, split.border.col,
                   split.lty, split.lwd, split.space, split.yspace,
                   split.round * split.strheight,
                   under.cex, under.font, under.ygap, ygap,
                   split.shadow.col, split.shadow.offset, bg,
                   small.underspace, split.strwidth, split.strheight)

        if(!draw.split.shadows1)
            draw.labs1(split.labs, split.boxes, split.yspace,
                split.cex * cex, split.font, split.family, split.col,
                draw.split.shadows1, draw.split.shadows, split.shadow.col, split.round)

        if(is.fancy(type)) { # right hand boxes
            right.split.boxes <-
                draw.boxes("right", draw.split.shadows1,
                   right.split.labs, node.xy, xlim, ylim,
                   nodes, branch,
                   Margin, xflip, yflip, main, sub,
                   col.main, cex.main, col.sub, cex.sub,
                   split.cex * cex, split.font, split.family, split.adj, split.yshift,
                   split.box.col, split.border.col,
                   split.lty, split.lwd,
                   split.space, split.yspace,
                   split.round * split.strheight,
                   under.cex, under.font, under.ygap, ygap,
                   split.shadow.col, split.shadow.offset, bg)

            if(!draw.split.shadows1)
                draw.labs1(right.split.labs, right.split.boxes, split.yspace,
                    split.cex * cex, split.font, split.family, split.col,
                    draw.split.shadows1, draw.split.shadows, split.shadow.col, split.round)
        }
        node.boxes <- draw.boxes("default", draw.shadows1,
                   node.labs, node.xy, xlim, ylim,
                   nodes, branch,
                   Margin, xflip, yflip, main, sub,
                   col.main, cex.main, col.sub, cex.sub,
                   cex, font, family, adj, yshift,
                   box.col, border.col,
                   lty, lwd,
                   space, yspace,
                   round * strheight,
                   under.cex, under.font, under.ygap, ygap,
                   shadow.col, shadow.offset, bg)

        draw.labs1(node.labs, node.boxes, yspace,
            cex, font, family, col,
            draw.shadows1, draw.shadows, shadow.col, round)

        if(yesno && !is.fancy(type) && !snip) # draw "yes" and "no" at root?
            draw.yes.no(type, draw.shadows1,
                    xflip, left, branch, xlim, ylim, node.xy, lwd,
                    yesno.yshift,
                    split.boxes, split.cex * cex, split.box.col, split.border.col,
                    split.shadow.col, split.shadow.offset,
                    nn.cex, nn.font, nn.family, nn.col, nn.box.col,
                    nn.border.col, nn.lty, nn.round, bg)

        if(nn || ni)
            draw.node.numbers(nn, ni, draw.shadows1, type, branch,
                    Margin, xflip, yflip, cex,
                    main, sub, col.main, cex.main, col.sub, cex.sub,
                    xlim, ylim, node.xy, is.leaf, nodes,
                    node.labs, font,  family, box.col, border.col, shadow.col,
                    under.cex, under.font, under.ygap, ygap,
                    split.labs, split.cex * cex, split.font, split.family, split.box.col,
                    split.border.col, split.shadow.col,
                    nn.cex, nn.font, nn.family, nn.col, nn.box.col,
                    nn.border.col, nn.lty, nn.lwd, nn.round,
                    split.adj, split.space, split.yspace, split.yshift,
                    yshift, adj, space, yspace, shadow.offset,
                    nn.adj, nn.yshift, nn.space, nn.yspace, bg)

        list(node.boxes=node.boxes, split.boxes=split.boxes)
    }
    do.bg <- function(col) # set elems of col that are 0 etc. to bg
    {
        if(is.null(col))
            col <- bg
        else
            col[which(col == 0) | is.na(col)] <- bg
        col
    }
    #--- prp starts here ---
    if(!inherits(x, "rpart"))
        stop0("the object passed to prp is not an rpart object")

    # Process dots args.  The call to eval.parent is necessary to evaluate the
    # call to say "c" when user does something like "xlim=c(0,2)". Note
    # also that we allow abbreviations by using say "dots$fo" instead of "dots$font".
    # TODO Is there a better way? This approach is fragile, we have to be
    #       extremely careful that abbreviation doesn't alias with other args.

    dots <- match.call(expand.dots=FALSE)$...
    check.dots(dots)
    adj      <- eval.parent(dots$adj);   if(is.null(adj))      adj      <- par("adj")
    cex.main <- eval.parent(dots$cex.m)
    cex.sub  <- eval.parent(dots$cex.s)
    col      <- eval.parent(dots$col);   if(is.null(col))      col      <- par("col")
    col.main <- eval.parent(dots$col.m); if(is.null(col.main)) col.main <- par("col.main")
    col.sub  <- eval.parent(dots$col.s); if(is.null(col.sub))  col.sub  <- par("col.sub")
    family   <- eval.parent(dots$fam);   if(is.null(family))   family   <- par("family")
    font     <- eval.parent(dots$fo);    if(is.null(font))     font     <- par("font")
    lty      <- eval.parent(dots$lt);    if(is.null(lty))      lty      <- par("lty")
    lwd      <- eval.parent(dots$lw);    if(is.null(lwd))      lwd      <- par("lwd")
    main     <- eval.parent(dots$mai)
    mar      <- eval.parent(dots$mar)
    sub      <- eval.parent(dots$sub);
    xlim     <- eval.parent(dots$xl)
    xpd      <- eval.parent(dots$xp)
    ylim     <- eval.parent(dots$yl)

    if(is.null(under.col))  under.col  <- col
    if(is.null(border.col)) border.col <- col
    if(is.null(branch.lwd)) branch.lwd <- lwd
    if(is.null(split.lwd))  split.lwd  <- lwd
    if(is.null(nn.lwd))     nn.lwd     <- lwd
    if(is.null(split.adj))  split.adj  <- adj

    # Set bg to the background color or "white" if transparent.
    # The idea is that we want a color that is opaque but matches background.
    bg <- par("bg") # TODO this incorrectly returns transparent with mfrow
    if(bg[1] == "transparent")
        bg <- "white"
    box.col          <- do.bg(box.col)
    border.col       <- do.bg(border.col)
    shadow.col       <- do.bg(shadow.col)
    under.col        <- do.bg(under.col)
    split.col        <- do.bg(split.col)
    split.box.col    <- do.bg(split.box.col)
    split.shadow.col <- do.bg(split.shadow.col)
    nn.col           <- do.bg(nn.col)
    nn.box.col       <- do.bg(nn.box.col)
    nn.border.col    <- do.bg(nn.border.col)

    # The idea with the following  argument checking is to catch user
    # errors here where possible before they cause an obscure message
    # later on.  But it is impossible to be exhaustive.

    stopifnot(is.numeric(type) && length(type) == 1 && floor(type) == type)
    if(type < TYPE.default || type > TYPE.fancy.all)
        stop0("type must be ", TYPE.default, "...",
              TYPE.fancy.all, ", you have type=", type)
    stopifnot.boolean(under)
    stopifnot.boolean(clip.left.labs[1])
    stopifnot.boolean(clip.right.labs[1])
    stopifnot.boolean(nn)
    stopifnot.boolean(ni)
    stopifnot.boolean(yesno)
    stopifnot.boolean(fallen.leaves)
    stopifnot.boolean(uniform)
    stopifnot.boolean(left)
    stopifnot.boolean(xflip)
    stopifnot.boolean(yflip)
    stopifnot.boolean(do.par)
    stopifnot.boolean(snip)
    stopifnot.boolean(compress)
    stopifnot.boolean(ycompress)
    stopifnot.boolean(xcompact)
    stopifnot.boolean(ycompact)
    stopifnot.boolean(add.labs)
    stopifnot.boolean(boxes.include.gap)
    stopifnot(all(split.round >= 0))
    stopifnot(all(nn.round >= 0))
    stopifnot(tweak > 0 && tweak <= 10) # upper limit is arb
    stopifnot(max.auto.cex >= 1)
    stopifnot(min(shift.amounts) >= 0 && max(shift.amounts) <= 10) # 10 is arb
    stopifnot(xcompact.ratio > 0 && xcompact.ratio <= 2) # 2 is arb, max useful is probably 1
    stopifnot(min.auto.cex >= 0 && min.auto.cex <= 1)
    stopifnot(branch >= 0 && branch <= 1)
    if(!is.null(snip.fun))
        check.func.args(snip.fun, "snip.fun", function(tree) NULL)
    if(length(family) != 1 || length(split.family) != 1 || length(nn.family) != 1)
        stop0("prp: family argument must be length 1 (family cannot be vectorized)")
    stopifnot(is.numeric(digits) && length(digits) == 1 &&
              floor(digits) == digits && digits >= 0)
    if(digits == 0)
        digits <- getOption("digits")
    if(!is.na.or.zero(branch.type)) {
        branch <- if(branch > .5) 1 else 0
        ycompact <- FALSE # want branches to be as vertical as possible
    }
    auto.cex <- FALSE
    if(is.null(cex)) {
        auto.cex <- TRUE    # automatically calculate cex
        cex <- 1
    }
    if(is.null(split.cex))
        split.cex <- 1
    if(fallen.leaves)
        compress <- FALSE
    obj <- x
    if(!is.null(obj$frame$splits))
        stop0("Old-style rpart object?  (frame$splits is NULL)")
    frame <- obj$frame
    is.leaf <- is.leaf(frame)
    nodes <- as.numeric(row.names(frame))
    if(do.par) {
        # Make the side margins small.
        # Retain the top edge for the main title but only if necessary.
        # Likewise the bottom edge for the subtitle.
        # Note that family may change in my.strheight and init.plot, so we on.exit it here.
        init.plot(1, 1, Margin, xflip, yflip, main, sub,
                  col.main, cex.main, col.sub, cex.sub)
        par <- par("mar", "xpd", "family")
        on.exit(par(par))
        if(is.null(mar)) { # user did not explictly set mar when invoking prp?
            mar <- par$mar
            if(is.null(sub))  mar[1] <- 1
            if(is.null(main)) mar[3] <- 1
            mar[2] <- mar[4] <- 1
        }
        if(is.null(xpd)) # user did not explicitly set xpd when invoking prp?
            xpd <- NA
        par(mar=mar, xpd=xpd)
        par(new=TRUE) # done par for now, start next plot on the same page
    }
    node.labs <- internal.node.labs(obj, node.fun, deparse(substitute(node.fun)),
                                    type, extra, under, xsep, digits, varlen,
                                    prefix, suffix)

    split.labs <- split.labs.wrapper(x, split.fun,
                deparse(substitute(split.fun)),
                split.prefix, split.suffix,
                right.split.prefix, right.split.suffix,
                type, clip.left.labs, clip.right.labs, xflip, digits,
                varlen, faclen, facsep, eq, lt, ge)

    if(is.fancy(type)) {
        right.split.labs <- split.labs[match(2 * nodes+1, nodes)]
        split.labs <- split.labs[match(2 * nodes, nodes)]
        if(!left) # following msg assumes hard coded TYPE.fancy, TYPE.fancy.all
            stop0("left=FALSE is not yet supported with type=3 or 4")
    } else {
        if(left != xflip)   # default, set right labs to NA
            split.labs <- split.labs[match(2 * nodes, nodes)]
        else                # set left labs to NA
            split.labs <- split.labs[match(2 * nodes+1, nodes)]
    }
    draw.shadows <- !is.invisible(shadow.col, bg)
    draw.split.shadows <- !is.invisible(shadow.col, bg)

    # Recycle stuff that doesn't get recyled automatically.  It's more efficient
    # to recycle it here once rather than over and over in get.boxes etc.
    adj                 <- recycle(adj, nodes)
    space               <- recycle(space, nodes)
    yspace              <- recycle(yspace, nodes)
    shadow.offset       <- recycle(shadow.offset, nodes)
    under.cex           <- recycle(under.cex, nodes)
    under.ygap          <- recycle(under.ygap, nodes)
    split.adj           <- recycle(adj, nodes)
    split.space         <- recycle(split.space, nodes)
    split.yspace        <- recycle(split.yspace, nodes)
    split.shadow.offset <- recycle(split.shadow.offset, nodes)
    nn.adj              <- recycle(nn.adj, nodes)
    nn.space            <- recycle(nn.space, nodes)
    nn.yspace           <- recycle(nn.yspace, nodes)

    temp <- get.yshift(type, nodes, is.leaf,
                       cex, node.labs, yshift, yspace, under.cex,
                       split.labs, split.cex, split.yshift, split.yspace, ygap)
    yshift       <- temp$yshift
    split.yshift <- temp$split.yshift

    layout <- get.layout(obj, type, nn, fallen.leaves, branch,
        uniform, Margin, cex, auto.cex, compress, ycompress,
        trace, main, sub,
        node.labs, font, family, box.col, border.col,
        under.font, under.cex,
        split.labs, right.split.labs, split.cex, split.font, split.family,
        split.box.col, split.border.col,
        nspace, minbranch, adj, yshift, space, yspace,
        split.adj, split.yshift, split.space, split.yspace,
        gap, ygap, under.ygap, xcompact, ycompact, xcompact.ratio,  min.inter.height,
        max.auto.cex, min.auto.cex, ycompress.cex, accept.cex,
        shift.amounts, Fallen.yspace, bg)

    cex <- layout$cex
    gap <- layout$gap
    ygap <- layout$ygap
    # we use pmax here so there is always a little space even if it causes overlapping
    space <- pmax(.25, layout$node.space)
    yspace <- pmax(.25, layout$node.yspace)
    if(is.null(xlim))
        xlim <- layout$xlim
    stopifnot(is.numeric(xlim) && length(xlim) == 2)
    if(is.null(ylim))
        ylim <- layout$ylim
    stopifnot(is.numeric(ylim) && length(ylim) == 2)
    split.yshift <- layout$split.yshift
    if(trace) {
        tweak.msg <- if(tweak == 1) "" else sprintf(" (before applying tweak %g)", tweak)
        printf("cex %.3g%s   xlim c(%.3g, %.3g)   ylim c(%.3g, %.3g)\n",
               cex[1], tweak.msg, xlim[1], xlim[2], ylim[1], ylim[2])
    }
    if(!auto.cex && tweak != 1)
        warning0("cex and tweak both specified, applying both")
    cex <- tweak * cex
    all.cex <- merge1(cex, split.cex * cex)

    split.lwd  <- recycle(cex * split.lwd,  nodes)
    branch.lwd <- recycle(cex * branch.lwd, nodes)
    nn.lwd     <- recycle(cex * nn.lwd,     nodes)
    # do this last because split.lwd etc. above use lwd as the default
    lwd <- recycle(cex * lwd, nodes)

    node.xy <- layout$node.xy
    init.plot(xlim, ylim, Margin, xflip, yflip, main, sub,
              col.main, cex.main, col.sub, cex.sub,
              fam.main=fam.main, cex=cex[1], trace=trace, hide.title=FALSE)
    split.strwidth  <- my.strwidth("M", split.cex * cex, split.font, split.family)
    strheight <- my.strheight("M", cex, font, family)
    split.strheight <- my.strheight("M", split.cex * cex, split.font, split.family)
    node.boxes <- split.boxes <- NA
    if(add.labs) {
        if(is.null(round))
            round <- max(1, 2 * min(space, yspace))
        stopifnot(all(round >= 0))
        round <- recycle(round, nodes)
        if(is.null(leaf.round))
            leaf.round <- round
        stopifnot(all(leaf.round >= 0))
        leaf.round <- recycle(leaf.round, nodes)
        round[is.leaf] <- leaf.round[is.leaf]
        # draw shadows first, if any, so boxes and lines are over shadows
        if(draw.shadows || draw.split.shadows)
            draw.labs(draw.shadows, draw.split.shadows)
    }
    # draw branch lines now so they are over the shadows, under the boxes
    branch.xy <- draw.branches(obj, branch.type, branch.col,
                    branch.lty, branch.lwd,  branch.fill, branch.tweak,
                    node.labs, split.labs, node.xy, strheight,
                    type, branch, xflip, yflip, Margin, space, yspace,
                    cex, font, family, adj, box.col, border.col,
                    under.cex, under.font, under.ygap,
                    split.cex, split.font, split.family, split.adj, split.yshift,
                    split.box.col, split.border.col, split.space, split.yspace,
                    main, sub, col.main, cex.main, col.sub, cex.sub,
                    xlim, ylim, yshift, ygap, bg,
                    min.branch.width)
    if(add.labs) {
        temp <- draw.labs(FALSE, FALSE)
        node.boxes <- temp$node.boxes
        split.boxes <- temp$split.boxes
    }
    snipped.nodes <- NULL
    if(snip) {
        temp <- do.snip(obj, nodes, split.labs, node.xy, branch.xy,
                        branch.lty, branch.lwd, xlim, ylim, digits,
                        snip.fun)
        obj <- temp$obj
        snipped.nodes <- temp$snipped.nodes
    }
    invisible(list(obj=obj, snipped.nodes=snipped.nodes,
              xlim=xlim, ylim=ylim,
              x=node.xy$x, y=node.xy$y,
              branch.x=branch.xy$x, branch.y=branch.xy$y,
              labs=node.labs, cex=cex, boxes=node.boxes,
              split.labs="", split.cex=split.cex, split.box=split.boxes))
}
init.plot <- function(x, y,
                      Margin, xflip, yflip, main, sub,
                      col.main, cex.main, col.sub, cex.sub,
                      fam.main="", cex=1, trace=0, hide.title=TRUE)
{
    if(length(x) == 1)
        x <- c(0, x)
    if(length(y) == 1)
        y <- c(0, y)
    xlim <- range(x) + diff(range(x)) * c(-Margin, Margin)
    if(xflip)
        xlim <- rev(xlim)
    ylim <- range(y) + diff(range(y)) * c(-Margin, Margin)
    if(yflip)
        ylim <- rev(ylim)
    old.family <- NA
    if(hide.title) {
        # need to plot main and sub to get screen layout that accounts
        # for their allotted space, but want them to be invisible
        col.main <- col.sub <- 0
    } else {
        if(is.null(cex.main))
            cex.main <- max(.8, min(1.5 * cex, 1.2))
        if(is.null(cex.sub))
            cex.sub <- max(.8, min(1.5 * cex, 1.2))
        if(!is.null(sub)) # hack to make subtitle visible with our do.par()
            sub <- paste(sub, "\n")
        if(!identical(fam.main, ""))
            old.family <- par(family=fam.main) # on.exit for this already set up in prp()
    }
    plot(0, 0, xlim=xlim, ylim=ylim, type="n", axes=FALSE, xlab="", ylab="",
         main=main, sub=sub,
         col.main=col.main, cex.main=cex.main, col.sub=col.sub, cex.sub=cex.sub)
    if(!is.na(old.family))
        par(family=old.family)
    if(trace >= 2) { # draw the grid and region boxes
        col <- "palegreen"
        # set xpd so grid lines stay in our region
        old.xpd <- par(xpd=FALSE)
        on.exit(par(xpd=old.xpd))
        grid(col=col, lwd=cex)
        axis(1, col=col, col.axis=col)
        axis(2, col=col, col.axis=col)
        rect(xlim[1], ylim[1], xlim[2], ylim[2], col=NA, border=col, lty=1, lwd=cex)
        text((xlim[1] + xlim[2]) / 2, ylim[2], "xlim ylim", col=col)
        usr <- par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], col=NA, border=col, lty=1, lwd=cex)
        text((usr[1] + usr[2]) / 2, usr[4], "usr", col=col, xpd=NA)
    }
}
get.yshift <- function(type, nodes, is.leaf,
                       cex, node.labs, yshift, yspace, under.cex,
                       split.labs, split.cex, split.yshift, split.yspace, ygap)
{
    # Return number of lines (separated by \n) in each lab in labs
    # A \n\n counts as one \n (it should really equal split.yshift)
    # TODO there must be a simpler way to do this
    get.nlines <- function(labs)
    {
        labs <- gsub("\n\n", "\n", labs) # replace \n\n with \n
        sapply(strsplit(labs, "\n"), function(x) length(x))
    }
    #--- get.yshift starts here ---
    yshift       <- recycle(yshift, nodes)
    split.yshift <- recycle(split.yshift, nodes)

    # always want _some_ ygap else cex has to be very small to prevent overplotting
    if(length(ygap) == 0)
        ygap <- .5

    # We use the number of lines to estimate the vert space taken by the label
    # We discount the first line of the inbox text.

    node.sep.labs <- separate.labs(node.labs)
    node.nlines  <- get.nlines(node.sep.labs$in.box) - 1 +
                    under.cex * get.nlines(node.sep.labs$under.box)

    split.sep.labs <- separate.labs(split.labs)
    split.nlines  <- get.nlines(split.sep.labs$in.box) - 1 +
                    under.cex * get.nlines(split.sep.labs$under.box)

    # Note that get.boxes uses split.yshift with split.cex, but need yshift
    # w.r.t. node.cex --- hence the node.to.split conversion ratio.

    ratio <- cex / split.cex

    # The following calculations must match the y1 and y2 calculations in get.boxes
    # We have to calculate the position of the _center_ of the boxes.
    # TODO investigate .9 etc. below (empirically determined) --- space for newline?

    node.shift  <- -yshift + .6 + .9 * node.nlines + .5 * yspace
    split.shift <- .5 * yspace + .6 + .9 * split.nlines

    if(is.fancy(type)) {
        # Want the node box on the node, and the top of left split box below
        # the bottom of node box, with a little vertical space.
        # The right split box position is calculated in get.boxes, once we
        # have a known cex and calculate the absolute box sizes.
        # The node box, left split and right split get treated as three
        # separate boxes in get.layout.

        split.yshift <- split.yshift -
                        ratio * node.shift - max(ygap, .25) - split.shift

    } else if(type == TYPE.all) {
        # Want the split box on the node, and the top of the node
        # box just above the bottom of the split box (a slight overlap).
        # These get combined and treated as one large box in get.layout.

        new.yshift <- -node.shift - split.shift / ratio

        # following needed when user explicitly sets yshift
        split.yshift <- split.yshift + yshift
        # we only want to shift node labels with a split label sitting on them
        yshift[!is.leaf] <- new.yshift[!is.leaf]
        # except that we do move leaf nodes down a bit to roughly match their split brothers
        yshift[is.leaf] <- yshift[is.leaf] - node.nlines[is.leaf]

    } else if(type == TYPE.all.under) {
        # Want the node box on the node, and the top of the split
        # box just below the bottom of the node box.
        # These get combined and treated as one large box in get.layout.
        split.yshift <- split.yshift -
                        ratio * node.shift - split.shift
    }
    list(yshift=yshift, split.yshift=split.yshift)
}
get.node.coords <- function(obj, uniform, branch, compress,
                            nspace, minbranch, fallen.leaves, Fallen.yspace)
{
    if(NROW(obj$frame) <= 1)    # tree is just a root?
        compress <- FALSE       # prevent rpartco from issuing an error
    # initialize nspace for rpartco()
    if(!compress)
        nspace <- -1    # magic value for rpartco meaning no compression
    if(is.null(nspace))
        nspace <- branch
    xy <- my.rpartco(obj, uniform, nspace, minbranch)
    x <- xy$x
    y <- xy$y

    # scale x to 0 ... 1
    x <- x - min(x)
    if(length(x) > 1)   # needed when tree is just a root
        x <- x  / max(x)

    # scale y to 0 ... 1 so top node is always at 1, bottom leaves at 0
    y <- y - min(y)
    if(length(y) > 1)
        y <- y  / max(y)

    if(fallen.leaves) {
        is.leaf <- is.leaf(obj$frame)
        y[is.leaf] <- 0
        # Make more space above leaves by shifting all other nodes up to make to
        # more likely that we will be able to shift fallen leaves for more space.
        # It usually looks a little better too.
        y[!is.leaf] <- (y[!is.leaf] + Fallen.yspace) / (1 + Fallen.yspace)
    }
    list(x=x, y=y)
}
# Get the box coords, a row for each box.
# Use do.init.plot=FALSE if want to use char sizes etc. of the existing plot.

get.boxes <- function(boxtype,  # one of "default", "left", "right", "undersplit"
    labs, x, y, xlim, ylim, nodes, branch,
    Margin, xflip, yflip, main, sub, col.main, cex.main, col.sub,  cex.sub,
    cex, font, family, adj, yshift, box.col, border.col, space, yspace,
    ygap, bg, do.init.plot=TRUE,
    box.around.all.text=TRUE)   # else box only around "in box" text i.e. text before \n\n
{                               # TRUE when figuring out box spacing, FALSE when drawing boxes
    if(do.init.plot)
        init.plot(xlim, ylim, Margin, xflip, yflip, main, sub,
                  col.main, cex.main, col.sub, cex.sub)

    # to minimize blanking out parts of the branch lines, we want only a
    # small white space around letters when fancy and non-visible boxes
    if((boxtype == "left" || boxtype == "right") && is.box.invisible(box.col, border.col, bg))
        space <- recycle(min(.2, space), labs)

    # sanity check that variables are already expanded correctly for recycling
    stopifnot(length(adj)    == length(labs) &&
              length(yshift) == length(labs) &&
              length(space)  == length(labs) &&
              length(yspace) == length(labs))

    # TODO simplistic approach for now, assumes \n approx equal to under.ygap
    stripped.labs <- gsub("\n\n", "\n", labs) # replace \n\n with \n
    sep.labs <- separate.labs(labs)
    in.box.labs <- sep.labs$in.box

    height1         <- recycle(my.strheight("M", cex, font, family), labs)
    stripped.height <- my.strheight(stripped.labs, cex, font, family)
    width1          <- recycle(my.strwidth("M", cex, font, family), labs)
    width <- if(box.around.all.text)
                my.strwidth(labs,        cex, font, family)
             else
                my.strwidth(in.box.labs, cex, font, family)
    if(do.init.plot)
        par(new=TRUE) # done strwidth and strheight, ensure we stay on same page
    xy.to.calc.xshift <- list(x=x, y=y) # grab these before they change

    # The text function (which draws the labels elsewhere), centers the labels
    # vertically on the point, a fact that is used to calculate box positions below.
    # The following calculations must match the calculations in get.yshift.

    x[is.na(labs)] <- y[is.na(labs)] <- NA
    x1 <- x - adj     * width - space/2 * width1                      # left edge
    x2 <- x + (1-adj) * width + space/2 * width1                      # right edge
    y2 <- y + yshift * height1 + stripped.height/2 + yspace * height1 # top of box
    if(box.around.all.text)
        y1 <- y2 - stripped.height - yspace * height1
    else {
        in.box.height <- my.strheight(in.box.labs, cex, font, family)
        y1 <- y2 - in.box.height - yspace * height1
    }
    xshift <- 0 # stays at zero unless left or right split
    if(length(x1) > 1) { # tree is not just a root?
        if(boxtype == "undersplit") {
            # splits are under the node boxes,
            # force a little extra space under the split box so branch line is visible
            stopifnot(box.around.all.text) # should use this only in get.layout
            y2 <- y2 - height1
        } else if(boxtype == "left") {
            child <- match(2 * nodes, nodes)        # left child
        } else if(boxtype == "right") {
            child <- match(2 * nodes + 1, nodes)    # right child
            # lower the right splits relative to the left splits
            box.heights <- y2 - y1
            adjust <- box.heights + (max(ygap, .4) * height1)
            y1 <- y1 - adjust
            y2 <- y2 - adjust
            yshift <- yshift - min(box.heights / height1, na.rm=TRUE)
        }
        if(boxtype == "left" || boxtype == "right") {
            # adjust x coords so labels are centered on the branch lines
            branch.xy <- get.branches(xy.to.calc.xshift$x, xy.to.calc.xshift$y, nodes, branch)
            x <- branch.xy$x[, child]
            y <- branch.xy$y[, child]
            xshift <- (x[2, ] - x[3, ]) + # 1.3 below to exaggerate the separation, looks better
                      1.3 * yshift * height1 * (x[2, ] - x[1, ]) / (y[2, ] - y[1, ])
        }
    }
    list(x1=x1 + xshift, y1=y1, x2=x2 + xshift, y2=y2)
}
draw.boxes <- function(fancy.style, draw.shadow, labs, xy,
                       xlim, ylim, nodes, branch,
                       Margin, xflip, yflip, main, sub,
                       col.main, cex.main, col.sub, cex.sub,
                       cex, font, family, adj, yshift,
                       box.col, border.col,
                       lty, lwd, space, yspace, r,
                       under.cex, under.font, under.ygap, ygap,
                       shadow.col, shadow.offset, bg,
                       small.underspace=FALSE, split.strwidth=0, split.strheight=0)
{
    box <- get.boxes(fancy.style, labs, xy$x, xy$y, xlim, ylim, nodes, branch,
                     Margin, xflip, yflip, main, sub,
                     col.main, cex.main, col.sub, cex.sub,
                     cex, font, family, adj,
                     yshift, box.col, border.col, space, yspace,
                     ygap, bg,
                     do.init.plot=FALSE,
                     box.around.all.text=FALSE)

    new.box <- box
    if(small.underspace) {
        # Splits are under the node boxes.  Reduce sides and bottom of box slightly so
        # just a little white space below the split box so branch line is more visible.
        add.space  <- pmin(space, .6)
        add.yspace <- pmin(yspace, .8)
        new.box$x1 <- new.box$x1 + (space - add.space)   * split.strwidth
        new.box$x2 <- new.box$x2 - (space - add.space)   * split.strwidth
        new.box$y1 <- new.box$y1 + (yspace - add.yspace) * split.strheight
    }
    if(!draw.shadow)
        rounded.rect(new.box$x1, new.box$y1, new.box$x2, new.box$y2,
                     xlim, ylim, r, box.col, border.col, lty, lwd)
    else if(!is.invisible(shadow.col, bg))
        draw.shadow(new.box$x1, new.box$y1, new.box$x2, new.box$y2,
                    xlim, ylim, r, shadow.col, shadow.offset)
    box
}
is.fancy <- function(type)
{
    type == TYPE.fancy.noall || type == TYPE.fancy.all
}
# text before \n\n goes in the box
# text after \n\n if any goes under the box
separate.labs <- function(labs)
{
    labs <- strsplit(labs, "\n\n")
    list(in.box    = sapply(labs, function(x) x[1]),
         under.box = sapply(labs, function(x) paste(x[-1], collapse="\n")))
}
get.box.centers <- function(box)
{
    list(x=(box$x1 + box$x2)/2, y=(box$y1 + box$y2)/2)
}
