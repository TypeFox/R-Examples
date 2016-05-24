subClasses <- function(Cl, directOnly = TRUE, complete = TRUE, ...)
{
    ## utility for classTree():
    if (isClassDef(Cl)) {
        cDef <- Cl
        Cl <- cDef@className
    } else { ## need getClass() can give error because sub classes can
        ## be "not defined" (?!)   -- e.g. "iMatrix"
        cDef <- if (complete) getClass(Cl) else getClassDef(Cl)
    }

    subs <- showExtends(cDef@subclasses, printTo = FALSE)
    if(directOnly) subs$what[subs$how == "directly"] else subs$what
}

numOutEdges <- function(g)
{
    ## Purpose: returns a named integer vector giving for each node in g,
    ##		the number of edges *from* the node
    ## ----------------------------------------------------------------------
    ## Arguments: g: graph
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  8 Feb 2007, 22:59
    el <- sapply(edgeL(g), `[[`, "edges")
    sapply(el, length)
}

is.leaf <- function(g) numOutEdges(g) == 0
## The graph package now defines a leaves() generic {w/ degree.dir}
##     leaves  <- function(g) nodes(g)[is.leaf(g)]


bGraph <- function(n, root = "Mom",
                   leaves = paste(l.prefix, seq(length=n), sep=""),
                   l.prefix = "D", # for 'D'aughter
                   weights = NULL,
                   mode = c("undirected", "directed"))
{
    ## Purpose: Create a "branch graph", a simple tree with root and
    ##		n branches / leaves
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: Aug 2005
    if(!missing(leaves)) {
        stopifnot(is.character(leaves))
        n <- length(leaves)
    } else stopifnot(is.numeric(n), length(n) == 1, n >= 0)

    mode <- match.arg(mode)
    ftM2graphNEL(cbind(root, leaves), W = weights, edgemode = mode)
}

## agopen() has
## layoutType = c("dot","neato","twopi", "circo", "fdp")

abbrMatrixcl <- function(clnames, level = 1) {
    ### Do "Matrixclass" name abbrevation
    doSub <- clnames != "Matrix"
    clnames[doSub] <- sub("Matrix$", "*", clnames[doSub])
    ## sparse
    iSp <- grep("sparse", clnames)
    if(level >= 2)
        clnames[iSp] <- sub("sparse\\*", ".sp.", clnames[iSp])

    ## dense
    iD <- grep("dense", clnames)
    if(level >= 2)
        clnames[iD] <- sub("dense\\*", ".D.", clnames[iD])
    list(clnames = clnames,  iSparse = iSp, iDense = iD)
}


mRagraph <- function(gr, lType,
                     fixedsize = FALSE, ## <---- this is it !
                     fill = c("lightblue", "gray90"),
                     color =c("blue3", "gray60"),
                     labcol = c("blue3","green4","purple"))
{
    ## Produce a layed out graph, an "Ragraph" -- to be plotted
    if (!validGraph(gr))
        stop("The graph to be plotted is not a valid graph structure")
    if (missing(lType))
        lType <- "dot"

    ng <- nodes(gr)
##     leaves  <- function(g) nodes(g)[is.leaf(g)]
##     nonVirtual <- leaves(gr) ## the leaves are *non*-virtual classes
    nonVirtual <- ng[is.leaf(gr)] ## the leaves are *non*-virtual classes

    r <- abbrMatrixcl(ng)

    nAtt <- makeNodeAttrs(gr, label = r$clnames, shape = "ellipse",
                          fixedsize = fixedsize,
                          fillcolor = fill[1], color = color[1],
                          fontcolor = labcol[1])

    nAtt$fontcolor[r$iSparse] <- labcol[2]
    nAtt$fontcolor[r$iDense]  <- labcol[3]

    nAtt$fillcolor[nonVirtual] <- fill[2]
    nAtt $   color[nonVirtual] <- color[2]

    ## but make one exception (for visualization):
    nAtt$fillcolor["pMatrix"] <- "thistle"
    nAtt $   color["pMatrix"] <- "turquoise"

    if(getOption("verbose")) { cat("mRagraph(): nodeAttrs: "); str(nAtt) }

    ### Returns the "layouted graph";  is +- ==  method("plot", "graph"):
    agopen(gr, name = "", layout = TRUE, layoutType = lType,
           attrs = list(), nodeAttrs = nAtt, edgeAttrs = list(),
           subGList = list(), recipEdges = "combined")
}

## plotRag() : a bit more than selectMethod("plot", "Ragraph")
##   --        but building on that
.optRagargs <- function(side = 1, adj = 0.05, cex = 0.75, line = 3)
    list(side = side, adj = adj, cex = cex, line = line)

plotRag <- function(ragr, sub, subArgs = .optRagargs(), ...)
{
    stopifnot(is(ragr, "Ragraph"))

    if(missing(sub)) {
	## nEdges <- length(unlist(edgeL(gr), use.names=FALSE))
	sub <- paste(length(ragr@AgNode), "nodes with",
		     length(ragr@AgEdge), "edges")
    }
### BUG in Rgraphviz ----> FIXME: bug report, ...
###    plot(ragr, sub = sub, ...)
### workaround {{but more flexible anyway:
    plot(ragr, ...)
    op <- par(xpd = NA) ; on.exit(par(op))
    do.call(mtext, c(list(text = sub), subArgs))
}


## Now do this recusively

classTree <- function(Cl, all = FALSE, ...)
{
    ## First a check
    if (isClassDef(Cl)) {
        cDef <- Cl
        Cl <- cDef@className
    } else cDef <- getClass(Cl)

    pkg <- cDef@package
    where <- if(pkg == ".GlobalEnv") .GlobalEnv else asNamespace(pkg)

    ## Now define a recursive function that computes the extended subtree
    ## for one class, and uses this for all sub-classes of Cl
    subtree <- function(cl, all) {
        stopifnot(isClassDef(cl))
        clN <- cl@className
        if(getOption('verbose')) cat(" ST",clN,":")
        sc <- subClasses(cl, directOnly = !all)
        if(length(sc) == 0) {
            if(getOption('verbose'))  cat(" is leaf\n")
            ## one node named 'cl':
            g <- new("graphNEL", nodes = clN, edgemode = "dir")
        }
        else {
            if(getOption('verbose'))  cat(" has leaves:\n\t")
            g <- bGraph(root = clN, leaves = sc, mode = "dir")
            for(cc in sc) {
                if(getOption('verbose'))  cat(":: ",clN,"-",cc,sep="")
                st <- subtree(getClass(cc, where = where), all = all)
                ##    -------## recursive
                if(numNodes(st) > 1)
                    g <- join(g, st)
            }
        }
        g
    }

    subtree(cDef, all = all)
}

