
# implements a wrapper to distance (similarity) computation on 
# collections of sequences. auto and cross distances can be
# computed (compare with dist in package proxy)
#
# note that 1) we can supply lists of vectors or vectors of 
#              character (strings)
#           2) operation weights are in the order of
#              insertion/deletion, equality, and replacing
#           3) the first row/column of the matrix of alphabet
#              weights are used for replacement with the empty
#              symbol (space)
#           4) include NA, etc if exclude = NULL
#           5) but the C function returns NA if NAs are encounterd
#           6) use parallel mode only if y != NULL
# 
# ceeboo 2006, 2008

sdists <- function(x,y=NULL, method="ow", weight=c(1,1,0,2),
                   exclude=c(NA,NaN,Inf,-Inf), pairwise = FALSE) {
    METHODS <- c("ow","aw","awl")
    code <- pmatch(method, METHODS)
    if (is.na(code))
       stop("invalid method")
    if (code == -1)
       stop("ambiguous method")
    if (is.character(x))
       x <- strsplit(x,"")
    if (!is.list(x))
       stop("'x' not a list")
    if (!is.null(y)) {
       if (is.character(y))
          y <- strsplit(y,"")
       if (!is.list(y))
          stop("'y' not a list")
    }
    if (code >= 2) {
       if (!is.matrix(weight))
          stop("'weight' not a matrix")
       if (dim(weight)[1] != dim(weight)[2])
          stop("'weight' not square")
       if (is.null(colnames(weight)))
          stop("'weight' no colnames")
       l <- colnames(weight)
    }
    else {
       if (length(weight) < 4)
	  stop("'weight' invalid")
       # determine common symbol set
       l <- sort(unique(c(unlist(x),unlist(y),"")),na.last=TRUE)       
    }
    x <- lapply(x,function(x) 
             factor(x,levels=l,exclude=if(is.integer(x))NA else exclude))
    if (!is.null(y)) { 
       y <- lapply(y,function(x) 
                factor(x,levels=l,exclude=if(is.integer(x))NA else exclude))
       if (pairwise && length(x) != length(y))
           stop("'pairwise', lengths of 'x' and 'y' do not conform")
    }
    if (!is.double(weight))
       storage.mode(weight) <- "double"
    obj <- .Call(R_sdists,x,y,as.integer(code),weight,pairwise)
    if (is.null(y))
       obj <- structure(obj, Size=length(x), class="dist",
                             Diag=FALSE, Upper=FALSE,
                             Labels=names(x), method=method)
    else if (!pairwise) {
       rownames(obj) <- names(x)
       colnames(obj) <- names(y)
    }
    obj
}

# as there is no unique space symbol available not 'excluding'
# NA has NA as result (see the C implementation).
#
# if graph = TRUE the vector of transcripts is transformed into 
# graph data that can be supplied to 'segments', or 'grid.segments', 
# etc. the dynmic programming table is returned as attribute
# 'table' and the traceback graph in attribute 'graph'.

sdists.trace <- 
function(x,y, method="ow", weight=c(1,1,0,2), exclude=c(NA,NaN,Inf,-Inf), graph = FALSE, partial = FALSE) {
    METHODS <- c("ow","aw","awl")
    code <- pmatch(method, METHODS)
    if (is.na(code))
       stop("invalid method")
    if (code == -1)
        stop("ambiguous method")
    if (is.character(x)) {
       if (length(x) != 1)
           stop("'x' not a scalar string")
       x <- strsplit(x,"")[[1]]
    }
    if (is.factor(x))
        x <- as.character(x)
    if (!is.vector(x))
        stop("'x' not a vector")
    if (is.character(y)) {
       if (length(y) != 1)
           stop("'y' not a scalar string")
       y <- strsplit(y,"")[[1]]
    }
    if (is.factor(y))
        y <- as.character(y)
    if (!is.vector(y))
        stop("'y' not a vector")
    if (code >= 2) {
       if (partial)
	  stop("'partial' not implemented")
       if (!is.matrix(weight))
          stop("'weight' not a matrix")
       if (is.null(colnames(weight)))
          stop("'weight' no colnames")
       l2 <- colnames(weight)
       if (is.null(rownames(weight))) {
          if (dim(weight)[1] != dim(weight)[2])
             stop("'weight' not square")
          l1 <- l2
       } else
          l1 <- rownames(weight)  
    }
    else {
       if (length(weight) < 4)
          stop("'weight' invalid")
       if (partial) {
	  if (length(weight) < 5)
	     weight <- c(weight, weight[1], 0)
	  if (length(weight) < 6)
	     weight <- c(weight, 0)
       }
       # determine symbol sets
       l1 <- l2 <- sort(unique(c(x,y,"")),na.last=TRUE)     
    }
    x <- factor(x,levels=l1,exclude=if(is.integer(x))NA else exclude)
    y <- factor(y,levels=l2,exclude=if(is.integer(y))NA else exclude)
    if (!is.double(weight))
        storage.mode(weight) <- "double"
    t <- .Call(R_sdists_transcript, x, y, as.integer(code), weight, graph, partial)
    if (is.na(t[1]))
        return(t)
    # reduce set of transcripts/paths
    if (partial) {
	z <- t
	## reduce to maximum number of trailing inserts
	k <- attr(regexpr("I+$", z), "match.length")
	z <- z[k == max(k)]
	## reduce to maximum number of matches
	k <- sapply(lapply(strsplit(z, ""), table), "[", "M")
	k <- which(k == max(k, na.rm = TRUE))
	if (length(k))
	    z <- z[k]
	## reduce to maximum number of leading inserts
	k <- attr(regexpr("^I+", z), "match.length")
	z <- z[k == max(k)]
	attributes(z) <- attributes(t)
	t <- z
    }
    if (graph) {
        dimnames(attr(t, "table")) <- 
            list(x = c("", as.character(x)), y = c("", as.character(y)))
        attr(t, "graph") <- .Call(R_sdists_graph, t)
        names(attr(t, "graph")) <- c("x0", "y0", "x1", "y1", "weight")
        names(attr(t, "pointer")) <- c("x0", "y0", "x1", "y1")
        class(t) <- "sdists.graph"
        return(t)
    }
    z <- lapply(t, function(t) .Call(R_sdists_align, x, y, t))
    names(z) <- t
    attr(z, "value") <- attr(t, "value")
    attr(z, "partial") <- attr(t, "partial")
    class(z) <- "sdists.trace"
    z
}

### experimental plot function for 
#
# idea from: http://home.uchicago.edu/~aabbott/
#
# in R 2.4.x we will fix yscale = c(ny, 0)
#
# label in grid.xaxis cannot contain "", i.e. does not
# produce output if it does.
#
# with pdf() it produces a garbage file that segfaults
# xpdf (but not acroread) :-(
#
# fixme: use another line type for prefixes or suffixes
#        of local alignments.
#
# ceeboo 2006

plot.sdists.graph <- 
function(x, circle.col = 1, graph.col = 2, circle.scale = c("mean", "max", "last", "text"), main = "", ...) {
    
    circle.scale <- match.arg(circle.scale)
    
    g <- attr(x, "graph")
    b <- attr(x, "pointer")
    t <- attr(x, "table")

    nx <- dim(t)[2]
    ny <- dim(t)[1]

    if (circle.scale == "text")
	fontsize <- 24	## FIXME
    else {    
	t <- t - min(t)
	t <- t / switch(circle.scale, mean = mean(t),
				      max  = max(t),
                                      last = t[ny,nx])
    }
    cn <- colnames(t)
    rn <- rownames(t)

    # bug?fix
    cn[cn == ""] <- " "
    rn[rn == ""] <- " "

    grid.newpage()

    grid.text(y = 0.95, label = main,
              gp = gpar(fontface = "bold"))

    vp <- viewport(xscale = c(0, nx), yscale = c(0, ny),
                   width  = nx / max(nx, ny) * 0.70, 
                   height = ny / max(nx, ny) * 0.70)

    pushViewport(vp)

    grid.grill(h  = seq(ny)-0.5, v = seq(nx)-0.5, default.units = "native")
    
    grid.xaxis(at = seq(nx)-0.5, label = cn)
    grid.yaxis(at = seq(ny)-0.5, label = rn)
    
    if (circle.scale == "text")
        mapply(grid.text,
            label = t,
            x = rep(1:nx,  each = ny) - 1/2,
            y = rep(1:ny, times = nx) - 1/2,
            MoreArgs = list(
                check.overlap = TRUE,
                default.units = "native",
                gp = gpar(col = "lightgrey", fontsize = fontsize))
        )
    else
        grid.circle(x = rep(1:nx,  each = ny) - 1/2,
                    y = rep(1:ny, times = nx) - 1/2,
                    r = t / 2,
                    default.units = "native",
                    gp = gpar(col = circle.col))

    grid.segments(x0 = b$y0 + 1/2, y0 = b$x0 + 1/2,
                  x1 = b$y1 + 1/2, y1 = b$x1 + 1/2,
                  default.units = "native",
                  gp = gpar(lty = 3))
 
    grid.segments(x0 = g$y0 + 1/2, y0 = g$x0 + 1/2,
                  x1 = g$y1 + 1/2, y1 = g$x1 + 1/2,
                  default.units = "native",
                  gp = gpar(col = graph.col, lwd = g$weight,
                            lty = (g$y1 > g$y0 & g$x1 > g$x0 & 
                                   cn[g$y1+1] == rn[g$x1+1]) + 1))

    popViewport()
}

###
