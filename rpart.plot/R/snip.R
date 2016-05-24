# do.snip.R

do.snip <- function(obj, nodes, split.labels, node.xy, branch.xy,
                    branch.lty, branch.lwd, xlim, ylim, digits, snip.fun)
{
    snip.rpart1 <- function(obj, deleted.nodes)
    {
        if(length(deleted.nodes))
            my.snip.rpart(obj, deleted.nodes)
        else
            obj # no changes
    }
    do.mouse.snip <- function()
    {
        draw.quit.button <- function(col)
        {
            width <- my.strwidth("QUIT", 1, 2)
            height <- my.strheight("QUIT", 1, 2)
            x <- xlim[1] + width
            y <- node.xy$y[1] - height
            rounded.rect(x - .6 * width, y - height,
                         x + .6 * width, y + height,
                         xlim, ylim, 1, 0, col, 1, 1)
            text(x, y, "QUIT", font=2, col=col)
            # add pseudo nodes for identify(), so can recognize a
            # click on the quit button (actually anywhere near the button)
            node.xy$x <- c(node.xy$x, x-.6 * width, x, x + .6 * width)
            node.xy$y <- c(node.xy$y, y,            y, y)
            node.xy
        }
        parent <- function(node) # parent of the given node
        {
            node %/% 2
        }
        get.parents <- function(node) # path to root, including node
        {
            if(node == 1)   # root?
                node
            else            # recurse
                c(node, get.parents(parent(node)))
        }
        get.children <- function(node) # node and all its children
        {
            if(is.leaf[match(node, nodes)])
                node
            else
                c(node,
                  get.children(2 * node),     # left child
                  get.children(2 * node + 1)) # right child
        }
        show.branches <- function()
        {
            branch.col <- recycle("black", nodes)
            # want pink line to nodes with deleted parents
            dnodes <- nodes[deleted.nodes]
            branch.col[match(dnodes[parent(dnodes) %in% dnodes], nodes)] <- "pink"
            # need loop for proper recycling of lwd etc.
            for(i in 1:length(nodes))
                lines(branch.xy$x[,i], branch.xy$y[,i],
                      col=branch.col[i], lty=branch.lty, lwd=branch.lwd)
        }
        print.node.info <- function(msg, inode)
        {
            if(nchar(msg)) {
                msg1 <- sprintf("%s node %d", msg, nodes[inode])
                printf("%-18s %s\n", msg1,
                    if(is.leaf[inode]) "" else split.labels[inode])
            }
            print(obj$frame[inode, ])
            cat("\n")
            flush.console()
        }
        #--- do.mouse.snip starts here ---
        old.options <- options(width=1000, digits=digits) # so no wrap in print.node.info
        on.exit(options(width=old.options$width, digits=old.options$digits))
        cat("Click to snip ...\n")
        if(!is.null(snip.fun))
            snip.fun(obj)
        flush.console()
        is.leaf <- is.leaf(obj$frame)
        node.xy <- draw.quit.button("black")
        # don't display the shoulders (to minimize overplotting)
        branch.xy$x[3,] <- branch.xy$y[3,] <- NA
        deleted.nodes <- rep(FALSE, length(nodes))
        show.branches()
        while((inode <- identify(node.xy$x, node.xy$y, n=1, plot=FALSE)) &&
                  inode <= length(nodes)) { # not a click on QUIT?
            if(is.leaf[inode])
                print.node.info("Leaf", inode)
            else {
                if(!deleted.nodes[inode]) { # if node is not currently deleted, then delete it
                    deleted.nodes[match(get.children(nodes[inode]), nodes)] <- TRUE
                    if(is.null(snip.fun)) # reduce clutter
                        print.node.info("Delete", inode)
                    else
                        snip.fun(snip.rpart1(obj, nodes[deleted.nodes]))
                } else {
                    # Node is currently deleted, so undelete it and its children ---
                    # but not if any of its ancestors are deleted.
                    # [-1] below removes node itself from its path to the root.
                    if(any(nodes[deleted.nodes] %in% get.parents(nodes[inode])[-1])) {
                        cat("Cannot delete node", nodes[inode],
                            "because its parent is already deleted\n")
                    } else {
                        deleted.nodes[match(get.children(nodes[inode]), nodes)] <- FALSE
                        if(is.null(snip.fun))
                            print.node.info("Undelete", inode)
                        else
                            snip.fun(snip.rpart1(obj, nodes[deleted.nodes]))
                    }
                }
                show.branches()
            }
        }
        draw.quit.button("gray")
        cat("Quit\n")
        nodes[deleted.nodes]
    }
    #--- do.snip starts here ---
    if(!dev.interactive()) { # can't proceed if output device is not the screen
        warning0("ignoring snip=TRUE for ", names(dev.cur())[1], " device")
        return(list(obj=obj, snipped.nodes=NULL)) # NOTE: return
    }
    snipped.nodes <- do.mouse.snip()
    if(length(snipped.nodes) == 0)
        list(obj=obj, snipped.nodes=NULL)
    else
        list(obj=snip.rpart1(obj, snipped.nodes), snipped.nodes=snipped.nodes)
}
