################################################################## 
# THIS FILE CONTAINS FUNCTIONS THAT ARE RELATED TO OBTAINING     #
# TREE-BASED MODELS WITH THE RPART PACKAGE. IT BASICALLY PROVIDES#
# SOME EXTRA FUNCTIONALITY.                                      #
# IT IS A PART OF THE PACKAGE DMwR                               #
##################################################################
# Author : Luis Torgo (ltorgo@inescporto.pt)     Date: Jan 2009  #
# License: GPL (>= 2)                                            #
##################################################################



# =====================================================
# Function that obtains a tree-based model using the
# x-SE post pruning rule of CART (Breiman et al. 1984).
# The idea is to grow an overly large tree and then post
# prune it using the internal cross validation estimates
# obtained by the initial call to rpart(), which are
# accessible through the cptable component of rpart objects.
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
rpartXse <- function(form,data,se=1,cp=0,minsplit=6,verbose=F,...) {
#   require(rpart)
   tree <- rpart(form,data,cp=cp,minsplit=minsplit,...)
   if (verbose && ncol(tree$cptable) < 5) 
     warning("No pruning will be carried out because no estimates were obtained.")
   rt.prune(tree,se,verbose)
 }

#
# Helper function to actually carry out the prunning
rt.prune <- function(tree,se=1,verbose=T,...) {
   if (ncol(tree$cptable) < 5) tree
   else {
     lin.min.err <- which.min(tree$cptable[,4])
     if (verbose && lin.min.err == nrow(tree$cptable))
       warning("Minimal Cross Validation Error is obtained 
                at the largest tree.\n  Further tree growth 
               (achievable through smaller 'cp' parameter value),\n
                could produce more accurate tree.\n")
     tol.err <- tree$cptable[lin.min.err,4] + se * tree$cptable[lin.min.err,5]
     se.lin <- which(tree$cptable[,4] <= tol.err)[1]
     prune.rpart(tree,cp=tree$cptable[se.lin,1]+1e-9)
   }
}




# =====================================================
# Function that plots an rpart tree in a nice way.
# The code is based on a slight adaptation of the original
# code of function text.rpart() by
# Terry M Therneau and Beth Atkinson. R port by Brian Ripley.
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
prettyTree <- function(t,
                     compress=F, branch=0.2, margin=0.1, uniform=T,
                     all=T, cex=0.8, font=10, use.n=T,
                     fwidth=0.5, fheight=0.45, center=0,
                     ...) { 
  if (!inherits(t, "rpart")) 
    stop("Not a valid rpart object.")
  plot(t,
       uniform = uniform, branch = branch, margin = margin,
       compress=compress,...)
  changed.text.rpart(t,
                 all = all, fancy = TRUE, digits = 3, pretty = 0,
                 cex=cex,font=font,use.n=use.n,
                 fwidth=fwidth,fheight=fheight,center=center,...)
}


# =====================================================
# This is a very slight adaptation of the original text.rpart()
# function found in package rpart by
# Terry M Therneau and Beth Atkinson. R port by Brian Ripley.
# =====================================================
# Original code by
# Terry M Therneau and Beth Atkinson. R port by Brian Ripley. (2008)
# Small modification by
# Luis Torgo, Jan 2009
# =====================================================
changed.text.rpart <- function (x,
           splits = TRUE, label, FUN = text, all = FALSE, pretty = NULL, 
           digits = getOption("digits") - 3, use.n = FALSE, fancy = FALSE, 
           fwidth = 0.65, fheight = 0.8,
           bg.color="white",center=0,...) 
{
#  require(rpart)
  if (!inherits(x, "rpart")) 
    stop("Not legitimate rpart")
  if (!is.null(x$frame$splits)) 
    x <- rpart:::rpconvert(x)
  if (nrow(x$frame) <= 1) 
    stop("fit is not a tree, just a root")
    frame <- x$frame
    if (!missing(label)) 
        warning("argument 'label' is currently unused")
    cxy <- par("cxy")
    if (!is.null(srt <- list(...)$srt) && srt == 90) 
        cxy <- rev(cxy)
    xy <- rpart:::rpartco(x)
    node <- as.numeric(row.names(x$frame))
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    if (splits) {
        left.child <- match(2 * node, node)
        right.child <- match(node * 2 + 1, node)
        rows <- labels(x, pretty = pretty)
        if (fancy) {
            xytmp <- rpart:::rpart.branch(x = xy$x, y = xy$y, node = node)
            leftptx <- (xytmp$x[2L, ] + xytmp$x[1L, ])/2
            leftpty <- (xytmp$y[2L, ] + xytmp$y[1L, ])/2
            rightptx <- (xytmp$x[3L, ] + xytmp$x[4, ])/2
            rightpty <- (xytmp$y[3L, ] + xytmp$y[4L, ])/2
            FUN(leftptx, leftpty + 0.2 * cxy[2L], rows[left.child[!is.na(left.child)]], 
                ...)
            FUN(rightptx, rightpty - 0.2 * cxy[2L], rows[right.child[!is.na(right.child)]], 
                ...)
        }
        else FUN(xy$x, xy$y + 0.5 * cxy[2L], rows[left.child], 
            ...)
    }
    leaves <- if (all) 
        rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    ylevels <- attr(x, "ylevels")
    stat <- if (is.null(frame$yval2)) 
        x$functions$text(yval = frame$yval[leaves], dev = frame$dev[leaves], 
            wt = frame$wt[leaves], ylevel = ylevels, digits = digits, 
            n = frame$n[leaves], use.n = use.n)
    else x$functions$text(yval = frame$yval2[leaves, ], dev = frame$dev[leaves], 
        wt = frame$wt[leaves], ylevel = ylevels, digits = digits, 
        n = frame$n[leaves], use.n = use.n)
    oval <- function(middlex, middley, a, b) {
        theta <- seq(0, 2 * pi, pi/30)
        newx <- middlex + a * cos(theta)
        newy <- middley + b * sin(theta)
        polygon(newx, newy, border = TRUE, col = bg.color)
    }
    rectangle <- function(middlex, middley, a, b) {
        newx <- middlex + c(a, a, -a, -a)
        newy <- middley + c(b, -b, -b, b)
        polygon(newx, newy, border = TRUE, col = bg.color)
    }
    if (fancy) {
        maxlen <- max(rpart:::string.bounding.box(stat)$columns) + 1L
        maxht <- max(rpart:::string.bounding.box(stat)$rows) + 1L
        if (fwidth < 1) 
            a.length <- fwidth * cxy[1L] * maxlen
        else a.length <- fwidth * cxy[1L]
        if (fheight < 1) 
            b.length <- fheight * cxy[2L] * maxht
        else b.length <- fheight * cxy[2L]
        for (i in parent) oval(xy$x[i], xy$y[i], a = sqrt(2) * 
            a.length/2, b = sqrt(2) * b.length/2)
        child <- match(node[frame$var == "<leaf>"], node)
        for (i in child) rectangle(xy$x[i], xy$y[i], a = a.length/2, 
            b = b.length/2)
    }
    if (fancy) 
        FUN(xy$x[leaves], xy$y[leaves] + center * cxy[2], stat, 
            ...)
    else FUN(xy$x[leaves], xy$y[leaves] - 0.5 * cxy[2], stat, 
        adj = 0.5, ...)
    invisible()
}

