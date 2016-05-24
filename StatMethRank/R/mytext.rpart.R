mytext.rpart <- function (x, splits = TRUE, label, FUN = text, all = FALSE, pretty = NULL, 
    digits = getOption("digits") - 3L, use.n = FALSE, fancy = FALSE, 
    fwidth = 0.2, fheight = 0.2, bg = par("bg"), minlength = 1L, parms) 
{
    if (nrow(x$frame) <= 1L) 
        stop("fit is not a tree, just a root")
    frame <- x$frame
    if (!missing(label)) 
        warning("argument 'label' is no longer used")
    col <- names(frame)
    ylevels <- attr(x, "ylevels")
    if (!is.null(ylevels <- attr(x, "ylevels"))) 
        col <- c(col, ylevels)
    cxy <- par("cxy")
    xy <- myrpartco(x, parms)
    node <- as.numeric(row.names(frame))
    is.left <- (node%%2L == 0L)
    node.left <- node[is.left]
    parent <- match(node.left/2L, node)
    if (splits) { # default = TRUE
        left.child <- match(2L * node, node)
        right.child <- match(node * 2L + 1L, node)
        rows <- if (!missing(pretty) && missing(minlength)) 
            mylabels.rpart(x, pretty = pretty)
        else mylabels.rpart(x, minlength = minlength)
        if (fancy) {
            xytmp <- myrpart.branch(x = xy$x, y = xy$y, node = node, branch = parms$branch)
            leftptx <- (xytmp$x[2L, ] + xytmp$x[1L, ])/2
            leftpty <- (xytmp$y[2L, ] + xytmp$y[1L, ])/2
            rightptx <- (xytmp$x[3L, ] + xytmp$x[4L, ])/2
            rightpty <- (xytmp$y[3L, ] + xytmp$y[4L, ])/2
            FUN(leftptx, leftpty + 0.52 * cxy[2L], rows[left.child[!is.na(left.child)]]) 
            FUN(rightptx, rightpty - 0.52 * cxy[2L], rows[right.child[!is.na(right.child)]])
        }
        else FUN(xy$x, xy$y + 0.5 * cxy[2L], rows[left.child])
    }
    leaves <- if (all)
        rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    stat <- mytext(yval = frame$yval[leaves], digits = digits, n = frame$n[leaves], use.n = use.n)
    if (fancy) {
        if (col2rgb(bg, alpha = TRUE)[4L, 1L] < 255) 
            bg <- "white"
        oval <- function(middlex, middley, a, b) {
            theta <- seq(0, 2 * pi, pi/30)
            newx <- middlex + a * cos(theta)
            newy <- middley + b * sin(theta)
            polygon(newx, newy, border = TRUE, col = bg)
        }
        rectangle <- function(middlex, middley, a, b) {
            newx <- middlex + c(a, a, -a, -a)
            newy <- middley + c(b, -b, -b, b)
            polygon(newx, newy, border = TRUE, col = bg)
        }
        maxlen <- max(mystring.bounding.box(stat)$columns) + 1L 
        maxht <- max(mystring.bounding.box(stat)$rows) + 1L
        a.length <- if (fwidth < 1) 
            fwidth * cxy[1L] * maxlen
        else fwidth * cxy[1L]
        b.length <- if (fheight < 1) 
            fheight * cxy[2L] * maxht
        else fheight * cxy[2L]
        for (i in parent) oval(xy$x[i], xy$y[i], sqrt(2) * a.length/2, 
            sqrt(2) * b.length/2)
        child <- match(node[frame$var == "<leaf>"], node)
        for (i in child) rectangle(xy$x[i], xy$y[i], a.length/2, 
            b.length/2)
    }
    if (fancy) 
        FUN(xy$x[leaves], xy$y[leaves] - 0.5 * cxy[2L], stat)
    else FUN(xy$x[leaves], xy$y[leaves] - 0.5 * cxy[2L], stat, adj = 0.5)
    invisible()
}