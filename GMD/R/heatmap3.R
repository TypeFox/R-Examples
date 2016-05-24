
## ************************************************************************
## FILENAME: heatmap3.R
## 
## AUTHOR: Xiaobei ZHAO <xiaobei _at_ binf.ku.dk>
##
## v0.3.3 2014-08-25 17:59:25 EDT
## v0.3.1 Sat Feb 04 16:13:09 CET 2012
## v0.3   Fri Nov 18 02:30:01 CET 2011
##
## DESCRIPTION:
## The R package \code{GMD} source code is available at
## \url{http://cran.r-project.org/web/packages/GMD/}
## under GPL license. 
##
## COMMENTS:
## Some code is adapted from `gplots:::heatmap.2'.
## ************************************************************************



##' Enhanced heatmap representation with dendrograms and partition given the \emph{elbow criterion}
##' or a desired number of clusters.\cr
##' 1) a dendrogram added to the left side and to the top, according to cluster analysis;\cr
##' 2) partitions in highlighted rectangles, according to the "elbow" rule or a desired number of clusters.\cr
##'
##' Enhanced heatmap representation with partition and summary statistics (optional).
##' This is an enhanced version of `heatmap.2' function in the Package \pkg{gplots}. The enhancement includes:
##' 1) Improved performance with optional input of precomputed \emph{\code{dist}} object
##' and \emph{\code{hclust}} object.
##' 2) Highlight of specific cells using rectangles. For instance, the cells of clusters of interests.
##' (Examples should be included in future.)
##' 3) Add-on plots in addition to the heatmap, such as cluster-wise summary plots and
##' overall clustering summary plots, to the right of or under the heatmap.
##' @title Enhanced Heatmap Representation with Dendrogram and Partition
##' @param x data matrix or data frame, or dissimilarity matrix or `dist' object
##' determined by the value of the 'diss' argument.
##' ##diss logical flag: if TRUE (default for \code{dist} or \code{dissimilarity} objects),
##' then \code{x} is assumed to be a dissimilarity matrix. If FALSE,then \code{x} is treated
##' as a matrix of observations by variables.
##' @param diss logical, whether the \code{x} is a dissimilarity matrix
##' @param Rowv one of the following: TRUE, a `dend' object, a vector or NULL/FALSE;
##' determines if and how the \emph{row} dendrogram should be reordered.
##' @param Colv one of the following: "Rowv", TRUE, a `dend' object, a vector or NULL/FALSE; 
##' determines if and how the \emph{column} dendrogram should be reordered.
##' @param dendrogram character string indicating whether to draw 'none', 'row', 'column' or
##' 'both' dendrograms.  Defaults to 'both'.
##' @param dist.row a \code{dist} object for \emph{row} observations.
##' @param dist.col a \code{dist} object for \emph{column} observations.
##' @param dist.FUN function used to compute the distance (dissimilarity) between
##' both rows and columns.  Defaults to \code{gdist}.
##' @param dist.FUN.MoreArgs a list of other arguments to be passed to \code{gdist}
##' @param hclust.row a \code{hclust} object (as produced by \code{hclust}) for \emph{row} observations.
##' @param hclust.col a \code{hclust} object (as produced by \code{hclust}) for \emph{column} observations.
##' @param hclust.FUN function used to compute the hierarchical clustering when
##' "Rowv" or "Colv" are not dendrograms. Defaults to \code{hclust}.
##' @param hclust.FUN.MoreArgs a list of other arguments to be passed to \code{hclust}.\cr
##' Defaults to \code{list(method="ward")}
##' @param scale character indicating if the values should be centered and scaled in either the row direction
##' or the column direction, or none.  The default is \code{"none"}.
##' @param na.rm logical, whether NA values will be removed when scaling.
##' @param cluster.by.row logical, whether to cluster \emph{row} observations and reorder the input accordingly.
##' @param cluster.by.col logical, whether to cluster \emph{column} observations and reorder the input accordingly.
##' @param kr numeric, number of clusters in rows; suppressed when \code{row.cluster} is specified.
##' DEFAULT: NULL.
##' @param kc numeric, number of clusters in columns; suppressed when \code{col.cluster} is specified.
##' DEFAULT: NULL.
##' @param row.clusters a numerical vector, indicating the cluster labels of \emph{row} observations.
##' @param col.clusters a numerical vector, indicating the cluster labels of \emph{column} observations.
##' @param revR logical indicating if the row order should be 'rev'ersed for plotting.
##' @param revC logical indicating if the column order should be 'rev'ersed for plotting, such that
##' e.g., for the symmetric case, the symmetry axis is as usual.
##' @param add.expr expression that will be evaluated after the call to \code{image}.
##' Can be used to add components to the plot.
##' @param breaks numeric, either a numeric vector indicating the splitting
##' points for binning \code{x} into colors, or a integer number of
##' break points to be used, in which case the break points will
##' be spaced equally between \code{range(x)}. DEFAULT: 16 when not specified.
##' @param x.center numeric, a value of \code{x} for centering colors to
##' @param color.FUN function or function name in characters, for colors in the heatmap
##' @param sepList a \code{list} of length 2 giving the row and column lines of separation.
##' @param sep.color color for lines of separation.
##' @param sep.lty line type for lines of separation.
##' @param sep.lwd line width for lines of separation.
##' @param cellnote (optional) matrix of character strings which will be placed within each color cell,
##' e.g. cell labels or p-value symbols.
##' @param cex.note relative font size of \code{cellnote}.
##' @param notecol color of \code{cellnote}.
##' @param na.color Color to use for missing value (\code{NA}). Defaults to the plot background color.
##' @param trace character string indicating whether a solid "trace" line should be drawn across
##' \code{"row"}s or down \code{"column"}s, \code{"both"} or \code{"none"}.
##' The distance of the line from the center of each color-cell is proportional to the size of
##' the measurement. Defaults to \code{"none"}.
##' @param tracecol character string giving the color for "trace" line. Defaults to "cyan";
##' @param hline Vector of values within cells where a horizontal dotted line should be drawn.
##' only plotted if 'trace' is 'row' or 'both'. Default to the median of the breaks.
##' @param vline Vector of values within cells where a vertical dotted line should be drawn;
##' only drawn if 'trace' 'column' or 'both'. \code{vline} default to the median of the breaks.
##' @param linecol the color of \code{hline} and \code{vline}. Defaults to the value of 'tracecol'.
##' @param labRow character vectors with row labels to use; defaults to \code{rownames(x)}.
##' @param labCol character vectors with column labels to use; defaults to \code{colnames(x)}.
##' @param srtRow numerical, specifying (in degrees) how row labels should be rotated.
##' See \code{help("par", package="graphics")}.
##' @param srtCol numerical, specifying (in degrees) how col labels should be rotated.
##' See \code{help("par", package="graphics")}.
##' @param sideRow 2 or 4, which side row labels display.
##' @param sideCol 1 or 3, which side row labels display.
##' @param margin.for.labRow a numerical value gives the margin to plot \code{labRow}.
##' @param margin.for.labCol a numerical value gives the margin to plot \code{labCol}.
##' @param ColIndividualColors (optional) character vector of length \code{ncol(x)} containing
##' the color names for a horizontal side bar that may be used to annotate the columns of \code{x}.
##' @param RowIndividualColors (optional) character vector of length \code{nrow(x)} containing
##' the color names for a vertical side bar that may be used to annotate the rows of \code{x}.
##' @param cexRow positive numbers, used as 'cex.axis' in for column axis labeling.
##' The default currently only uses number of columns.
##' @param cexCol positive numbers, used as 'cex.axis' in for the row axis labeling.
##' The default currently only uses number of rows.
##' @param labRow.by.group logical, whether group unique labels for rows.
##' @param labCol.by.group logical, whether group unique labels for columns.
##' @param key logical indicating whether a color-key should be shown.
##' @param key.title character, title of the color-key ["Color Key"]
##' @param key.xlab character, xlab of the color-key ["Value"]
##' @param key.ylab character, ylab of the color-key ["Count"]
##' @param keysize numeric value indicating the relative size of the key
##' @param mapsize numeric value indicating the relative size of the heatmap.
##' @param mapratio the width-to-height ratio of the heatmap.
##' @param sidesize numeric value indicating the relative size of the sidebars.
##' @param cex.key.main a numerical value giving the amount by which \code{main}-title of color-key should be
##' magnified relative to the default.
##' @param cex.key.xlab a numerical value giving the amount by which \code{xlab} of color-key should be
##' magnified relative to the default.
##' @param cex.key.ylab a numerical value giving the amount by which \code{ylab} of color-key should be
##' magnified relative to the default.
##' @param density.info character string indicating whether to superimpose a 'histogram',
##' a 'density' plot, or no plot ('none') on the color-key.
##' @param denscol character string giving the color for the density display specified by 'density.info',
##' defaults to the same value as 'tracecol'.
##' @param densadj Numeric scaling value for tuning the kernel width when a density plot is drawn on the
##' color key.  (See the 'adjust' parameter for the 'density' function for details.)  Defaults to 0.25.
##' @param main an overall title for the plot. See \code{help("title", package="graphics")}.
##' @param sub a subtitle for the plot, describing the distance and/or alignment gap (the "shift").
##' @param xlab a title for the x axis. See \code{help("title", package="graphics")}.
##' @param ylab a title for the y axis. See \code{help("title", package="graphics")}.
##' @param cex.main a numerical value giving the amount by which \code{main}-title should be
##' magnified relative to the default.
##' @param cex.sub a numerical value giving the amount by which \code{sub}-title should be
##' magnified relative to the default.
##' @param font.main An integer which specifies which font to use for \code{main}-title.
##' @param font.sub An integer which specifies which font to use for \code{sub}-title.
##' @param adj.main The value of 'adj' determines the way in which \code{main}-title strings are justified.
##' @param mgp.main the margin line (in 'mex' units) for the \code{main}-title.
##' @param mar.main a numerical vector of the form \code{c(bottom, left, top, right)} which gives the
##' number of lines of margin to be specified on the four sides of the \code{main}-title.
##' @param mar.sub a numerical vector of the form \code{c(bottom, left, top, right)} which gives the
##' number of lines of margin to be specified on the four sides of the \code{sub}-title.
##' @param if.plot logical, whether to plot. Reordered matrix is returned without graphical output if FALSE.
##' @param plot.row.partition logical, whether to plot \emph{row} partition.
##' @param plot.col.partition logical, whether to plot \emph{column} partition.
##' @param cex.partition a numerical value giving the amount by which \code{partition} should be
##' magnified relative to the default.
##' @param color.partition.box color for the \code{partition} box.
##' @param color.partition.border color for the \code{partition} border.
##' @param plot.row.individuals logical, whether to make a plot of \emph{row} observations.
##' @param plot.col.individuals logical, whether to make a plot of \emph{column} observations.
##' @param plot.row.clusters logical, whether to make a summary plot of \emph{row} clusters.
##' @param plot.col.clusters logical, whether to make a summary plot of \emph{column} clusters.
##' @param plot.row.clustering logical, whether to make a summary plot of overall \emph{row} clustering.
##' @param plot.col.clustering logical, whether to make a summary plot of overall \emph{column} clustering.
##' @param plot.row.individuals.list a list of expressions that is used to \code{plot.row.individuals}
##' @param plot.col.individuals.list a list of expressions that is used to \code{plot.col.individuals}
##' @param plot.row.clusters.list a list of expressions that is used to \code{ plot.row.clusters}
##' @param plot.col.clusters.list a list of expressions that is used to \code{plot.col.clusters}
##' @param plot.row.clustering.list a list of expressions that is used to \code{plot.row.clustering}
##' @param plot.col.clustering.list a list of expressions that is used to \code{plot.col.clustering}
##' @param row.data (optional) data used to \code{plot.row.individuals}, \code{ plot.row.clusters} or \code{plot.row.clustering}
##' @param col.data (optional) data used to \code{plot.col.individuals}, \code{ plot.col.clusters} or \code{plot.col.clustering}
##' @param if.plot.info logical, whether to plot \code{text.box}.
##' @param text.box character plotted when \code{if.plot.info} is TRUE.
##' @param cex.text a numerical value giving the amount by which \code{text.box} should be
##' magnified relative to the default.
##' @param ... arguments to be passed to method \code{heatmap.3}.\cr
##'e \code{help("image", package="graphics")}.
##' 
##' @return A reordered matrix according to \emph{row} or/and \emph{col} dendrogram(s) and
##' indices that used for reordering.
##' @examples
##' 
##' ## ------------------------------------------------------------------------
##' ## Example1: mtcars
##' ## ------------------------------------------------------------------------
##' ## load library
##' require("GMD")
##' 
##' ## load data
##' data(mtcars)
##' 
##' ## heatmap on raw data
##' x  <- as.matrix(mtcars)
##' 
##' dev.new(width=10,height=8)
##' heatmap.3(x)                               # default, with reordering and dendrogram
##' \dontrun{
##' heatmap.3(x, Rowv=FALSE, Colv=FALSE)       # no reordering and no dendrogram
##' heatmap.3(x, dendrogram="none")            # reordering without dendrogram
##' heatmap.3(x, dendrogram="row")        # row dendrogram with row (and col) reordering
##' heatmap.3(x, dendrogram="row", Colv=FALSE) # row dendrogram with only row reordering
##' heatmap.3(x, dendrogram="col")             # col dendrogram
##' heatmap.3(x, dendrogram="col", Rowv=FALSE) # col dendrogram with only col reordering
##' heatmapOut <-
##'   heatmap.3(x, scale="column")             # sacled by column
##' names(heatmapOut)                          # view the list that is returned
##' heatmap.3(x, scale="column", x.center=0)   # colors centered around 0
##' heatmap.3(x, scale="column",trace="column")  # trun "trace" on
##' }
##' 
##' ## coloring cars (row observations) by brand
##' brands <- sapply(rownames(x), function(e) strsplit(e," ")[[1]][1]) 
##' names(brands) <- c()
##' brands.index <- as.numeric(as.factor(brands))
##' RowIndividualColors <- rainbow(max(brands.index))[brands.index]
##' heatmap.3(x, scale="column", RowIndividualColors=RowIndividualColors)
##' 
##' ## coloring attributes (column features) randomly (just for a test :)
##' heatmap.3(x, scale="column", ColIndividualColors=rainbow(ncol(x)))
##' 
##' ## add a single plot for all row individuals 
##' dev.new(width=12,height=8)
##' expr1 <- list(quote(plot(row.data[rowInd,"hp"],1:nrow(row.data),
##' xlab="hp",ylab="",yaxt="n",main="Gross horsepower")),
##' quote(axis(2,1:nrow(row.data),rownames(row.data)[rowInd],las=2)))
##' heatmap.3(x, scale="column", plot.row.individuals=TRUE, row.data=x,
##'           plot.row.individuals.list=list(expr1))
##'
##' 
##' ## ------------------------------------------------------------------------
##' ## Example2: ruspini
##' ## ------------------------------------------------------------------------
##' ## load library
##' require("GMD")
##' require(cluster)
##' 
##' ## load data 
##' data(ruspini)
##' 
##' ## heatmap on a `dist' object
##' x <- gdist(ruspini)
##' main <- "Heatmap of Ruspini data"
##' dev.new(width=10,height=10)
##' heatmap.3(x, main=main, mapratio=1) # with a title and a map in square!
##' \dontrun{
##' heatmap.3(x, main=main, revC=TRUE)  # reverse column for a symmetric look
##' heatmap.3(x, main=main, kr=2, kc=2) # partition by predefined number of clusters
##' }
##' ## show partition by elbow
##' css.multi.obj <- css.hclust(x,hclust(x))
##' elbow.obj <- elbow.batch(css.multi.obj,ev.thres=0.90,inc.thres=0.05)
##' heatmap.3(x, main=main, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k)
##'
##' \dontrun{
##' ## show elbow info as subtitle
##' heatmap.3(x, main=main, sub=sub("\n"," ",attr(elbow.obj,"description")),
##' cex.sub=1.25,revC=TRUE,kr=elbow.obj$k, kc=elbow.obj$k)
##' }
heatmap.3 <-
  function(x,
           
           ## whether a dissimilarity matrix
           diss=inherits(x,"dist"),

           ## dendrogram control
           Rowv=TRUE,
           Colv=TRUE,
           dendrogram=c("both","row","column","none"),

           ## dist object
           dist.row,
           dist.col,
           dist.FUN=gdist,
           dist.FUN.MoreArgs=list(method="euclidean"),
           
           ## hclust object
           hclust.row,
           hclust.col,
           hclust.FUN=hclust,
           hclust.FUN.MoreArgs=list(method="ward"),
           
           ## data scaling
           scale=c("none","row","column"),
           na.rm=TRUE,

           ## clustering control
           cluster.by.row=FALSE,
           cluster.by.col=FALSE,
           kr=NA,      
           kc=NA,
           row.clusters=NA,
           col.clusters=NA,

           ## image plot
           revR=FALSE,
           revC=FALSE,
           add.expr,
           
           ## mapping data to colors
           breaks,
           ## centering colors to a value
           x.center,
           ## colors
           color.FUN=gplots::bluered,
           ##
           ## block sepration
           sepList=list(NULL,NULL),
           sep.color=c("gray45","gray45"),
           sep.lty=1,
           sep.lwd=2,
           

           ## cell labeling
           cellnote,
           cex.note=1.0,
           notecol="cyan",
           na.color=par("bg"),

           ## level trace
           trace=c("none","column","row","both"),
           tracecol="cyan",
           hline,
           vline,
           linecol=tracecol,

           ## Row/Column Labeling
           labRow=TRUE, ## shown by default
           labCol=TRUE, ## shown by default
           srtRow=NULL,
           srtCol=NULL,
           sideRow=4,
           sideCol=1,
           ##
           margin.for.labRow,
           margin.for.labCol,
           ColIndividualColors,
           RowIndividualColors,
           cexRow,
           cexCol,
           labRow.by.group=FALSE,
           labCol.by.group=FALSE,
           

           ## plot color key + density info
           key=TRUE,
           key.title="Color Key",
           key.xlab="Value",
           key.ylab="Count",
           
           keysize=1.5,
           mapsize=9,
           mapratio=4/3,
           sidesize=3,
           cex.key.main=0.75,
           cex.key.xlab=0.75,
           cex.key.ylab=0.75,
           density.info=c("histogram","density","none"),
           denscol=tracecol,
           densadj=0.25,

           ## plot titles/labels
           main="Heatmap",
           sub="",
           xlab="",
           ylab="",
           cex.main=2,
           cex.sub=1.5,
           font.main=2,
           font.sub=3,
           adj.main=0.5,
           mgp.main=c(1.5,0.5,0),
           mar.main=3,
           mar.sub=3,
           ## plot ##
           if.plot=TRUE,
           
           ## plot of partition (left/top of heatmap)
           plot.row.partition=FALSE,
           plot.col.partition=FALSE,
           cex.partition=1.25,
           color.partition.box="gray45",
           color.partition.border="#FFFFFF",

           ## plot of summary (right/bottom of heatmap)
           plot.row.individuals=FALSE,
           plot.col.individuals=FALSE,
           plot.row.clusters=FALSE,
           plot.col.clusters=FALSE,
           plot.row.clustering=FALSE,
           plot.col.clustering=FALSE,

           ##
           plot.row.individuals.list=FALSE,
           plot.col.individuals.list=FALSE,
           plot.row.clusters.list=FALSE,
           plot.col.clusters.list=FALSE,
           plot.row.clustering.list=FALSE,
           plot.col.clustering.list=FALSE,
           
           ## for plot of clusters - row
           row.data=FALSE,
           ## for plot of clusters - col
           col.data=FALSE,

           ##
           if.plot.info=FALSE,
           text.box,
           cex.text=1.0,
           ## extras
           ...
           )
{
  

  ## check input - take1 ##
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  x.ori <- x
  
  if(!inherits(x,"dist") & !is.matrix(x)){
    stop("`x' should either be a matrix, a data.frame or a `dist' object.")
  }

  if (! sideRow %in% c(2,4)){
    stop('sideRow must be either 2 or 4.')
  }
  
  if (! sideCol %in% c(1,3)){
    stop('sideCol must be either 1 or 3.')
  }
  
  ## store input
  Rowv.ori <- Rowv
  Colv.ori <- Colv
  
  
  ## check
  dendrogram <- match.arg(dendrogram)
  if ( (dendrogram %in% c("both","row")) & !inherits(Rowv,"dendrogram") ){
    warning("Discrepancy: row dendrogram is asked;  Rowv is set to `TRUE'.")
    Rowv <- TRUE
  }

  if ( (dendrogram %in% c("both","col")) & !inherits(Colv,"dendrogram") ){
    warning("Discrepancy: col dendrogram is asked;  Colv is set to `TRUE'.")
    Colv <- TRUE
  }

  
  if (identical(Rowv, FALSE) | missing(Rowv)){
    if(!identical(cluster.by.row,FALSE)){
      warning("Discrepancy: No row dendrogram is asked; cluster.by.row is set to `FALSE'.")
      cluster.by.row <- FALSE
    }
  } else {
    if(!identical(cluster.by.row,TRUE)){
      warning("Discrepancy: row dendrogram is asked; cluster.by.row is set to `TRUE'.")
      cluster.by.row <- TRUE
    }
  }

  if (identical(Colv, FALSE) | .invalid(Colv)){
    if(!identical(cluster.by.col,FALSE)){
      warning("Discrepancy: No col dendrogram is asked; cluster.by.col is set to `FALSE'.")
      cluster.by.col <- FALSE
    }
  } else {
    if(!identical(cluster.by.col,TRUE)){
      warning("Discrepancy: col dendrogram is asked; cluster.by.col is set to `TRUE'.")
      cluster.by.col <- TRUE
    }
  }
  
  if (!.invalid(kr)){
    if (is.numeric(kr)){
      if(!plot.row.partition){
        warning("Discrepancy: kr is set, therefore plot.row.partition is set to `TRUE'.")
        plot.row.partition <- TRUE
      }
    }
  }
  
  if (!.invalid(kc)){
    if (is.numeric(kc)){
      if(!plot.col.partition){
        warning("Discrepancy: kc is set, therefore plot.col.partition is set to `TRUE'.")
        plot.col.partition <- TRUE
      }
    }
  }
  
  
  ## generate dist.obj - row/col ##
  if (inherits(x,"dist")){
    dist.row <- dist.col <- x ## dist.obj
    x <- as.matrix(x)
    mat.row <- mat.col <- x ## mat.obj
    symm <- TRUE
  } else if (is.matrix(x)){
    symm <- isSymmetric(x)
    if (diss){
      if (!symm){
        stop("Dissimilary matrix should be symmetric. Please set `diss' to FALSE if `x' is not dissimilary matrix.")
      } else {
        print("`x' is treated as a dissimilary matrix.")
        flush.console()
      }
      mat.row <- mat.col <- x
      dist.row <- dist.col <- as.dist(x)
    } else{
      if (cluster.by.row) {
        if (.invalid(dist.row)){
          dist.row <- .call.FUN(dist.FUN,x,MoreArgs=dist.FUN.MoreArgs)
        }
        mat.row <- as.matrix(dist.row)
      } else {
        dist.row <- NULL
        mat.row <- NULL
      }
      if (cluster.by.col) {
        if (.invalid(dist.col)){
          dist.col <- .call.FUN(dist.FUN,t(x),MoreArgs=dist.FUN.MoreArgs)
        }
        mat.col <- as.matrix(dist.col)
      } else {
        dist.col <- NULL
        mat.col <- NULL
      }
    }
  }

  ## check input - take2: di ##
  di <- dim(x)
  cat("1.dim(x):",dim(x),"\n")
  
  if(length(di)!=2 || !is.numeric(x)){
    stop("`x' should only contain `numeric' values and can be converted to a 2-D matrix.")
  }

  ## parse param ##
  scale <- if(symm && .invalid(scale)) "none" else match.arg(scale) ## no scale on symmetric matrix
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  dist.FUN <- match.fun(dist.FUN)
  hclust.FUN <- match.fun(hclust.FUN)
  color.FUN <- match.fun(color.FUN)

  
  ## NG if both breaks and scale are specified ##
  if(!.invalid(breaks) & (scale!="none")){
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.",
            "Please consider using only one or the other.")
  }

  
  ##print(adj.main)
  ##print(mgp.main)
  ##print(cex.main)
  
  ## nr and nc ##
  nr <- di[1]
  nc <- di[2]

  ## check input - take3: nr,nc ##
  if(nr <=1 || nc <=1)
    stop("`x' must have at least 2 rows and 2 columns")

  ## font size of row/col labels ##
  cexRow0 <- 0.2+1/log10(nr)
  cexCol0 <- 0.2+1/log10(nc)

  
  if (.invalid(cexRow)) {
    cexRow <- cexRow0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for cexRow.')
    cexRow <- cexRow0*cexRow
  }
  if (.invalid(cexCol)) {
    cexCol <- cexCol0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for cexCol.')
    cexCol <- cexCol0*cexCol
  }

  
  ## cellnote ##
  ## ##if(.invalid(cellnote)) cellnote <- matrix("",ncol=ncol(x),nrow=nrow(x))

  ## ------------------------------------------------------------------------
  ## parse dendrogram ##
  ## ------------------------------------------------------------------------  
  
  if (missing(Rowv)) Rowv <- FALSE
  
  if (.invalid(Colv)) Colv <- if(symm) Rowv else FALSE
  if (Colv=="Rowv") {
    if ((!isTRUE(Rowv) | !symm) ){
      Colv <- FALSE
      warning("`Colv' is specified to use \"Rowv\", but either `Rowv' is invalid or `x' is not symmetric; Colv is suppressed.")
    } else{
      Colv <- Rowv
    }
  }

  
  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  cat("Preparing `hclust'... ")
  flush.console()

  if ( (!inherits(Rowv,"dendrogram") & !identical(Rowv,FALSE)) | (cluster.by.row & .invalid(row.clusters))){
    if (.invalid(hclust.row)){
      hclust.row <- .call.FUN(hclust.FUN,dist.row,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.row$order) != nr){
        stop("`hclust.row' should have equal size as the rows.")
      }
    }
  } else{
    hclust.row <- NULL
  }

  
  if(symm){
    hclust.col <- hclust.row
  }

  if ( !inherits(Colv,"dendrogram") & !identical(Colv,FALSE) | (cluster.by.col & .invalid(col.clusters))){
    if (.invalid(hclust.col)){
      hclust.col <- .call.FUN(hclust.FUN,dist.col,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.col$order) != nc){
        stop("`hclust.col' should have equal size as the cols.")
      }
    }
  } else {
    hclust.col <- NULL
  }
  

  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  cat("Preparing `dendrogram'... ")

  ddr <- ddc <- NULL
  ## get the dendrograms and reordering row/column indices - row ##
  if(inherits(Rowv,"dendrogram")){
    if (attr(Rowv,"members") != nr){
      stop("`Rowv' should contain equal size of members as the rows.")
    }
    ddr <- Rowv ## use Rowv 'as-is',when it is dendrogram
    rowInd <- order.dendrogram(ddr)
  } else if (is.integer(Rowv)){ ## compute dendrogram and do reordering based on given vector
    ddr <- as.dendrogram(hclust.row)
    ddr <-  reorder(ddr,Rowv) 
    rowInd <- order.dendrogram(ddr)
    if(nr != length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else if (isTRUE(Rowv)){ ## if TRUE,compute dendrogram and do reordering based on rowMeans

    Rowv <- rowMeans(x,na.rm=TRUE)
    ddr <- as.dendrogram(hclust.row)
    ddr <- reorder(ddr,Rowv)
    rowInd <- order.dendrogram(ddr)
    if(nr !=length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else{
    rowInd <- nr:1 ## from bottom.
  }
  
  ## get the dendrograms and reordering row/column indices - col ##
  if(inherits(Colv,"dendrogram")){
    if (attr(Colv,"members") != nc){
      stop("`Colv' should contain equal size of members as the cols.")
    }
    ddc <- Colv ## use Colv 'as-is',when it is dendrogram
    colInd <- order.dendrogram(ddc)
  } else if(identical(Colv,"Rowv")) {
    if(exists("ddr")){
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else{
      colInd <- rowInd
    }
  } else if(is.integer(Colv)){## compute dendrogram and do reordering based on given vector
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("`colInd' is of wrong length.")
  } else if (isTRUE(Colv)){## if TRUE,compute dendrogram and do reordering based on rowMeans
    Colv <- colMeans(x,na.rm=TRUE)
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc !=length(colInd))
      stop("`colInd' is of wrong length.")
  } else{
    colInd <- 1:nc ## from left
  }

  
  
  ## ------------------------------------------------------------------------
  ## check consistency
  ## ------------------------------------------------------------------------

  
  ## Xmisc::logme(dendrogram)
  ## Xmisc::logme(Colv)
  ## Xmisc::logme(Rowv)
  
  ## dendrogram - check consistency: Rowv ##
  if ( is.null(ddr) & (dendrogram %in% c("both","row"))){
    warning("Discrepancy: Rowv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting row dendogram.")
    if (is.logical(Colv) & (Colv.ori) & dendrogram=="both")
      dendrogram <- "column"
    else
      dendrogram <- "none"
  }

  
  ## dendrogram - check consistency: Colv ##
  if ( is.null(ddc) & (dendrogram %in% c("both","column"))){
    warning("Discrepancy: Colv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting column dendogram.")
    if (is.logical(Rowv) & (identical(Rowv.ori,TRUE) | is.numeric(Rowv.ori) | inherits(Rowv.ori,"dendrogram")) & dendrogram=="both")
      dendrogram <- "row"
    else
      dendrogram <- "none"
  }

  
  ## check consistency
  if (is.null(ddr)){
    if(isTRUE(cluster.by.row) | isTRUE(plot.row.partition) | isTRUE(plot.row.clusters) | isTRUE(plot.row.clustering) ){
      warning("Using invalid `Rowv' while allowing",
              "`cluster.by.row' or `plot.row.partition' or `plot.row.clusters' or `plot.row.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }
  
  if (is.null(ddc)){
    if(isTRUE(cluster.by.col) | isTRUE(plot.col.partition) | isTRUE(plot.col.clusters) | isTRUE(plot.col.clustering) ){
      warning("Using invalid `Colv' while allowing",
              "`cluster.by.col' or `plot.col.partition' or `plot.col.clusters' or `plot.col.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }

  if (is.null(ddr)) cluster.by.row <- plot.row.partition <- plot.row.clusters <- plot.row.clustering <- FALSE
  if (is.null(ddc)) cluster.by.col <- plot.col.partition <- plot.col.clusters <- plot.col.clustering <- FALSE


  ## ------------------------------------------------------------------------
  ## Reordering
  ## ------------------------------------------------------------------------
  cat("Reordering ... ")
  flush.console()
  ## reorder x and cellnote ##
  x <- x[rowInd,colInd]
  
  cat("2.dim(x):",dim(x),"\n")
  
  if (!.invalid(cellnote)) cellnote <- cellnote[rowInd,colInd]
  
  ## reorder labels - row ##
  if(identical(labRow,TRUE)){ ## Note: x is already reorderred 
    labRow <- if (is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
  } else if(identical(labRow,FALSE) | .invalid(labRow)){
    labRow <- rep("",nrow(x))
  } else if(is.character(labRow)){
    labRow <- labRow[rowInd]
  }
  ##cat('labRow',labRow,'\n')
  ##cat('rowInd',rowInd,'\n')

  ## reorder cellnote/labels - col ##
  if (identical(labCol,TRUE)){
    labCol <- if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
  } else if(identical(labCol,FALSE) | .invalid(labCol)){
    labCol <- rep("",ncol(x))
  } else if(is.character(labCol)){
    labCol <- labCol[colInd]
  }

  
  ## ------------------------------------------------------------------------
  ## scale
  ## center to 0 and scale to 1 in row or col but not both! ##
  ## ------------------------------------------------------------------------
  cat("Scaling ... ")
  flush.console()
  x <- .scale.data(x,scale,na.rm)
  cat("3.dim(x):",dim(x),"\n")


  ## ------------------------------------------------------------------------
  ## labels for observations/clusters/
  ## ------------------------------------------------------------------------
  ## margin for labels
  
  margin.for.labRow0 <- max(nchar(labRow))*0.75+0.2
  margin.for.labCol0 <- max(nchar(labCol))*0.75+0.2
  
  if (.invalid(margin.for.labRow)){
    margin.for.labRow <- margin.for.labRow0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labRow.')
    margin.for.labRow <- margin.for.labRow0*margin.for.labRow
  }
  
  if (.invalid(margin.for.labCol)){
    margin.for.labCol <- margin.for.labCol0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labCol.')
    margin.for.labCol <- margin.for.labCol0*margin.for.labCol    
  }
  
  ## group unique labels - row ## ##??check
  if (!.invalid(labRow.by.group) & !identical(labRow.by.group,FALSE)){
    group.value <- unique(labRow)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labRow)
    labRow <- rep("",length(labRow))
    labRow[group.index] <- group.value
  }
  
  ## group unique labels - col ## ##??check
  if (!.invalid(labCol.by.group) & !identical(labCol.by.group,FALSE)){
    group.value <- unique(labCol)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labCol)
    labCol <- rep("",length(labCol))
    labCol[group.index] <- group.value
  }


  
  ## ------------------------------------------------------------------------
  ## color breaks
  ## ------------------------------------------------------------------------
  cat("Making color breaks ... ")
  flush.console()
  
  ## set breaks for binning x into colors ##
  if(.invalid(breaks)){
    breaks <- 16
    ## print(sprintf("breaks=%s",breaks))
  }

  
  ## get x.range according to the value of x.center ##
  if (!.invalid(x.center)){ ## enhanced
    if (is.numeric(x.center)){
      x.range.old <- range(x,na.rm=TRUE)
      dist.to.x.center <- max(abs(x.range.old-x.center))
      x.range <- c(x.center-dist.to.x.center,x.center+dist.to.x.center)
    } else {
      stop("`x.center' should be numeric.")
    } 
  } else{
    x.range <- range(x,na.rm=TRUE)
  }


  ## set breaks for centering colors to the value of x.center ##
  if(length(breaks)==1){
    breaks <-
      seq(min(min(x,na.rm=TRUE),x.range[1]),
          max(max(x,na.rm=TRUE),x.range[2]),
          length.out=breaks)
  }

  ##   cat("breaks:\n")
  ##   print(breaks)
  
  ## count of breaks and colors ##
  nbr <- length(breaks)
  ncolor <- length(breaks)-1

  ## generate colors ##
  colors <- color.FUN(ncolor)
  ##cat("color.FUN:\n")
  ##print(color.FUN)
  ##cat("colors:\n")
  ##print(colors)
  
  
  ## set up breaks and force values outside the range into the endmost bins ##
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[] <- ifelse(x<min.breaks,min.breaks,x)
  x[] <- ifelse(x>max.breaks,max.breaks,x)

  
  cat("4.dim(x):",dim(x),"\n")
  
  ## ------------------------------------------------------------------------
  ## check if it is sufficient to draw side plots ##
  ## ------------------------------------------------------------------------
  if (cluster.by.row){
    if (!.invalid(row.clusters)) {## suppress kr
      if(!is.numeric(row.clusters) | length(row.clusters)!=nr | !(.is.grouped(row.clusters))){
        warning("`row.clusters' is not a grouped numeric vector of length nrow(x); cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      } else{
        row.clusters <- row.clusters[rowInd]
        kr <- length(unique(row.clusters))
      }
    } else {
      if (.invalid(kr)) kr <- 2
      if (is.numeric(kr) & length(kr)==1){
        row.clusters <- cutree(hclust.row,k=kr)
        row.clusters <- row.clusters[rowInd]
      } else {
        warning("`kr' should be numeric of length one; cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      }
    }
  }

  if (cluster.by.col){
    if (!.invalid(col.clusters)) {## suppress kc
      if(!is.numeric(col.clusters) | length(col.clusters)!=nc | !(.is.grouped(col.clusters))){
        warning("`col.clusters' is not a grouped numeric vector of length ncol(x); cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      } else{
        col.clusters <- col.clusters[colInd]
        kc <- length(unique(col.clusters))
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      }
    } else {
      if (.invalid(kc)) kc <- 2
      if (is.numeric(kc) & length(kc)==1){
        col.clusters <- cutree(hclust.col,k=kc)
        col.clusters <- col.clusters[colInd]
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      } else {
        warning("`kc' should be numeric of length one; cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      }
    }
  }

##   print("revC")
##   print(revC)
  if (!.invalid(kr) & !.invalid(kc)){
    print(sprintf("kr=%s,kc=%s",kr,kc))
    flush.console()
  }
  
  cat("5.dim(x):",dim(x),"\n")
  
  ## ------------------------------------------------------------------------
  ## Plotting
  ## ------------------------------------------------------------------------
  if (if.plot){

    ir <- length(plot.row.individuals.list)
    ic <- length(plot.col.individuals.list)
    cr <- length(plot.row.clustering.list)
    cc <- length(plot.col.clustering.list)

    
    cat("Plotting ... ")
    flush.console()
    
    if(mapratio<=1){
      sr <- 12
      sc <- sr*mapratio
    } else {
      sc <- 12
      sr <- sc/mapratio
    }
    
    ## ##print(sprintf("sr=%s,sc=%s",sr,sc))
    ## calculate the plot layout ##
    
    ## 1) for heatmap
    lmat <- matrix(1,nrow=sr,ncol=sc) 
    lwid <- c(rep(mapsize/sc,sc))
    lhei <- c(rep(mapsize/mapratio/sr,sr))

    ## 2) row.clusters
    if (plot.row.partition | plot.row.clusters){ 
      lmat <- cbind(max(lmat,na.rm=TRUE)+1,lmat) 
      lwid <- c(0.3,lwid) 
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid) 

    }
    
    ## 3) col.clusters
    if (plot.col.partition | plot.col.clusters){ 
      lmat <- rbind(c(NA,rep(max(lmat,na.rm=TRUE)+1,sc)),lmat) 
      lhei <- c(0.3/mapratio,lhei) 
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei)
    }

    if(!.invalid(RowIndividualColors)) { ## 4) add middle column to layout for vertical sidebar ##??check
      if(!is.character(RowIndividualColors) || length(RowIndividualColors) !=nr)
        stop("'RowIndividualColors' must be a character vector of length nrow(x)")
      lmat <- cbind(c(rep(NA,nrow(lmat)-sr),rep(max(lmat,na.rm=TRUE)+1,sr)),lmat)
      lwid <- c(0.2,lwid) 
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid) 
    }
    
    if(!.invalid(ColIndividualColors)) { ## 5) add middle row to layout for horizontal sidebar ##??check
      if(!is.character(ColIndividualColors) || length(ColIndividualColors) !=nc){
        stop("'ColIndividualColors' must be a character vector of length ncol(x)")
      }
      lmat <- rbind(c(rep(NA,ncol(lmat)-sc),rep(max(lmat,na.rm=TRUE)+1,sc)),lmat) 
      lhei <- c(0.2/mapratio,lhei) 
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei) 
    }

    ## 6) for row-dend
    lmat <- cbind(c(rep(NA,nrow(lmat)-sr),
                    rep(max(lmat,na.rm=TRUE)+1,sr)),
                  lmat
                  ) 
    lwid <- c(keysize,lwid)
    
    ## 7) for col-dend, 8) for kay
    lmat <- rbind(c(
                    max(lmat,na.rm=TRUE)+2,
                    rep(NA,ncol(lmat)-sc-1),
                    rep(max(lmat,na.rm=TRUE)+1,sc)
                    ),
                  lmat
                  )
    lhei <- c(keysize/mapratio,lhei)

    ##   cat("lmat:\n")
    ##   print(lmat)
    ##   print(lwid)
    ##   print(lhei)
    

    ## text.box##
    ## numbered 999 ##
    ## 9) for RowPlot (from bottom)
    ##print("plot.row.individuals")
    if(.invalid(text.box)){
      text.box <- "made by\nFunction: heatmap.3\nPackage: GMD\nin R"
    }
    if(plot.row.individuals) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),nrow(lmat)-sr),# text
                      rep((ir:1)+max(lmat,na.rm=TRUE)+(1),each=floor(sr/ir)),rep(NA,sr%%ir)
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    
    ## 10) for ColPlot from right
    ##print("plot.col.individuals")
    if(plot.col.individuals) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),ncol(lmat)-sc-1),# text
                      rep((1:ic)+max(lmat,na.rm=TRUE)+(1),each=floor(sc/ic)),rep(NA,sc%%ic),
                      ##NA # change to numeric if text.box
                      999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }
    
    ## 11) for RowPlot (from bottom)
    ##print("plot.row.clusters")
    if(plot.row.clusters) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-1), # text
                      rep((kr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sr/kr)),rep(NA,sr%%kr),
                      ##NA
                      999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    
    ## 12) for ColPlot from right 
    ##print("plot.col.clusters")
    if(plot.col.clusters) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-2),# text
                      rep((1:kc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/kc)),rep(NA,sc%%kc),
                      ##NA,NA # change to numeric if text.box
                      999,999
                      )
                    ) 
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }


    ## 13) for RowPlot (from bottom)
    ##print("plot.row.clustering")
    if(plot.row.clustering) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-2), # text
                      rep(c((cr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1)),each=floor(sr/cr)),rep(NA,sr%%cr),
                      ##NA,NA
                      999,999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    
    ## 14) for ColPlot from right
    ##print("plot.col.clustering")
    if(plot.col.clustering) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-3),# text
                      rep((1:cc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/cc)),rep(NA,sc%%cc),
                      ##NA,NA,NA # change to numeric if text.box
                      999,999,999
                      )
                    ) 
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }

    lmat[is.na(lmat)] <- 0
    if (any(lmat==999)) flag.text <- TRUE else flag.text <- FALSE
    lmat[lmat==999] <- max(lmat[lmat!=999])+1
    
    ## Graphics `output' ##
    ##op <- par(no.readonly=TRUE)
    ##suppressWarnings(on.exit(par(op)))

    ## layout
    layout(lmat,widths=lwid,heights=lhei,respect=FALSE)
    
    ## cat("lmat:\n")
    ## print(lmat)
    ## print(lwid)
    ## print(lhei)

    ## layout.show(n=max(lmat,na.rm=TRUE))
    ## readline("Touch any key to continue ..")

    
    ## reverse columns
    cat("revC:",revC,", revR:",revR,"\n")

    cat("5b.dim(x):",dim(x),"\n")
    
    if(revC){ ## x columns reversed
      iy <- nr:1
      ddc <- rev(ddc)
      x <- x[iy,]
      if (!.invalid(cellnote)) cellnote <- cellnote[iy,]
    } else {
      iy <- 1:nr
    }

    ## reverse rows
    if(revR){ ## x columns reversed
      ix <- nc:1
      ddr <- rev(ddr)
      x <- x[,ix]
      if (!.invalid(cellnote)) cellnote <- cellnote[,ix]
    } else {
      ix <- 1:nc
    }

    
    cat("5c.dim(x):",dim(x),"\n")
    
    ## 1) draw the main carpet/heatmap
    margins <- c(margin.for.labCol,0,0,margin.for.labRow)
    mgp <- c(3,1,0)
    par(mar=margins,mgp=mgp);outer=FALSE
    ##par(oma=margins,mar=c(0,0,0,0),mgp=c(0,0,0));outer=TRUE

    
    cat("scale:",scale,"\n")
    x.save <- x
    if(!symm || scale !="none"){ ##??
      x <- t(x)
      if (!.invalid(cellnote)) cellnote <- t(cellnote)
    }

    cat("5d.dim(x):",dim(x),"\n")
    
    
    ## cat(sprintf("heatmap:mar=c(%s)\n",paste(par("mai"),sep="",collapse=",")))
    image(1:nc,1:nr,
          x,
          xlim=0.5+c(0,nc),ylim=0.5+c(0,nr),
          axes=FALSE,xlab="",ylab="",col=colors,breaks=breaks,
          ...)

    
    ## plot/color NAs
    if(!.invalid(na.color) & any(is.na(x))){
      mmat <- ifelse(is.na(x),1,NA)
      image(1:nc,1:nr,mmat,axes=FALSE,xlab="",ylab="",
            col=na.color,add=TRUE)
    }

    ##
    ## labCol (?)
    if ((dendrogram %in% c("both","col")) & sideCol==3) {
      warning("Discrepancy: col dendrogram is asked; srtCol is set to 1.")
      sideCol <- 1
    }
    if (!length(srtCol)) {
      axis(sideCol,1:nc,labels=labCol,las=2,line=-0.5,tick=0,cex.axis=cexCol,outer=outer)
    } else {
      if (sideCol==1){
        if (sideCol==1) .sideCol <- par("usr")[3]-0.5*srtCol/90 else .sideCol <- par("usr")[4]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }
    
    if(!.invalid(xlab)) mtext(xlab,side=1,line=margins[1]-1.25)

    ## labRow (?)
    if ((dendrogram %in% c("both","row")) & sideRow==2) {
      warning("Discrepancy: row dendrogram is asked; sideRow is set to 4.")
      sideRow <- 4
    }
    if (!length(srtRow)) {
      axis(sideRow,iy,labels=labRow,las=2,line=-0.5,tick=0,cex.axis=cexRow,outer=outer)
    } else {
      if (sideRow==4){
        if (sideRow==4) .sideRow <- par("usr")[2]+0.5*srtRow/90 else .sideRow <- par("usr")[1]-0.5*srtRow/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)
      }
    }
    
    if(!.invalid(ylab)) mtext(ylab,side=4,line=margins[4]-1.25)

    if (!.invalid(add.expr))
      eval(substitute(add.expr))
    
    ## Enhanced: add 'sep.color' colored spaces to visually separate sections
    if (plot.row.partition | plot.row.clusters){ ##??
      plot.row.partitionList <- get.sep(clusters=row.clusters,type="row")
    } else {
      plot.row.partitionList <- NULL
    }
    if (plot.col.partition | plot.col.clusters){ ##??
      plot.col.partitionList <- get.sep(clusters=col.clusters,type="column")
    } else {
      plot.col.partitionList <- NULL
    }


    row.sepList <- sepList[[1]]
    if (!.invalid(row.sepList)){
      for (i in 1:length(row.sepList)){
        i.sep <- row.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[1]
             )
      }
    }

    col.sepList <- sepList[[2]]
    if (!.invalid(col.sepList)){
      for (i in 1:length(col.sepList)){
        i.sep <- col.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[2]
             )
      }
    }
    
        
    ## show traces
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    
    cat("6.dim(x):",dim(x),"\n")
    x.scaled  <- .scale.x(t(x),min.scale,max.scale)
    cat("7.dim(x):",dim(x),"\n")
  
    if(.invalid(hline)) hline=median(breaks)
    if(.invalid(vline)) vline=median(breaks)
    
    if(trace %in% c("both","column")){
      for( i in colInd ){
        if(!.invalid(vline)){
          vline.vals <- .scale.x(vline,min.scale,max.scale)
          abline(v=i-0.5+vline.vals,col=linecol,lty=2)
        }
        xv <- rep(i,nrow(x.scaled))+x.scaled[,i]-0.5
        xv <- c(xv[1],xv)
        yv <- 1:length(xv)-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }


    if(trace %in% c("both","row")){
      for( i in rowInd ){
        if(!.invalid(hline)){
          hline.vals <- .scale.x(hline,min.scale,max.scale)
          abline(h=i+hline,col=linecol,lty=2)
        }
        yv <- rep(i,ncol(x.scaled))+x.scaled[i,]-0.5
        yv <- rev(c(yv[1],yv))
        xv <- length(yv):1-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }

    ## cellnote
    if(!.invalid(cellnote)){
      text(x=c(row(cellnote)),
           y=c(col(cellnote)),
           labels=c(cellnote),
           col=notecol,
           cex=cex.note)
    }
    

    
    ## 2) plot.row.partition
    if(plot.row.partition |plot.row.clusters) { ##row.clusters
      par(mar=c(margins[1],0.5,0,0.1))

      row.clusters.unique <- unique(row.clusters)
      row.clusters.unique <- row.clusters.unique[!is.na(row.clusters.unique)]
      
      image(rbind(1:nr),
            xlim=0.5+c(0,1),ylim=0.5+c(0,nr),
            col=par("bg"),
            axes=FALSE)

      if (!.invalid(plot.row.partitionList)){
        for (i in 1:length(plot.row.partitionList)){
          i.sep <- plot.row.partitionList[[i]]
          rect(
               xleft=0+0.5,
               ybottom=i.sep[2]+0.5,
               xright=1+0.5,
               ytop=i.sep[4]+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- row.clusters.unique[i]
          ## ##s <- sprintf("%s (n=%s)",g,sum(g==row.clusters.unique))
          s <- g
          text(x=1,y=(i.sep[2]+0.5+i.sep[4]+0.5)/2,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=90
               )
        }
      }
      
    } 
    
    ## 3) plot.col.partition
    if(plot.col.partition | plot.col.clusters) {
      par(mar=c(0.1,0,0,margins[4]))
      col.clusters.unique <- unique(col.clusters)
      col.clusters.unique <- col.clusters.unique[!is.na(col.clusters.unique)]
      
      image(cbind(1:nc),
            xlim=0.5+c(0,nc),ylim=0.5+c(0,1),
            col=par("bg"),
            axes=FALSE)
      
      if (!.invalid(plot.col.partitionList)){
        for (i in 1:length(plot.col.partitionList)){
          i.sep <- plot.col.partitionList[[i]]
          rect(
               xleft=i.sep[1]+0.5,
               ybottom=0+0.5,
               xright=i.sep[3]+0.5,
               ytop=1+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- col.clusters.unique[i]
          ## ##s <- sprintf("%s (n=%s)",g,sum(g==col.clusters.unique))
          s <- g
          text(x=(i.sep[1]+0.5+i.sep[3]+0.5)/2,y=1,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=0
               )
        }
      }

    }
    
    
    ## 4) draw the side color bars - for row
    if(!.invalid(RowIndividualColors)) {    
      par(mar=c(margins[1],0,0,0.5))
      ##cat(sprintf("side bars - for row:mar=c(%s)\n",paste(par("mai"),sep="",collapse=",")))
      image(rbind(1:nr),col=RowIndividualColors[rowInd],axes=FALSE)
    } 
    
    ## 5) draw the side color bars - for col
    if(!.invalid(ColIndividualColors)) {
      par(mar=c(0.5,0,0,margins[4]))
      ##cat(sprintf("side bars - for col:mar=c(%s)\n",paste(par("mai"),sep="",collapse=",")))
      image(cbind(1:nc),col=ColIndividualColors[colInd],axes=FALSE)
    }

    
    ## 6) row-dend
    par(mar=c(margins[1],0,0,0))
    if(dendrogram %in% c("both","row")){
      plot(ddr,horiz=TRUE,axes=FALSE,yaxs="i",leaflab="none")
    }else{
      .plot.text(ylim=range(iy))
      if (sideRow==2){
        .sideRow <- par("usr")[2]-0.5*srtCol/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)        
      }
    }

    
    ## 7) col-dend and title
    mar3 <- (if(!is.null(main)) mar.main else 0) +
      (if(!is.null(sub)) mar.sub else 0)
    par(mar=c(0,0,mar3,margins[4]))
    ##print(par()$mar)
    
    if(dendrogram %in% c("both","column")){
      plot(ddc,axes=FALSE,xaxs="i",leaflab="none")
    } else{
      .plot.text(xlim=range(1:nc))
      if (sideCol==3){
        .sideCol <- par("usr")[3]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }

    
    ## title
    if (is.null(sub)) main.line <- 1 else main.line <- 3
    if(!is.null(main)) title(main,cex.main=cex.main,adj=adj.main,mgp=mgp.main,font.main=font.main,line=main.line)
    if(!is.null(sub)) title(sub,cex.main=cex.sub,adj=adj.main,mgp=mgp.main,font.main=font.sub,line=0)
    ##if(!is.null(main)) title(main,cex.main=1.5*op[["cex.main"]])

    
    ## 8) plot the color-key
    if(key){
      cex.key <- 0.75
      op.ori <- par()

      par(mar=c(2,1.5,0.75,1)*keysize,cex=cex.key,mgp=c(0.75,0,0),tcl=-0.05)
      z <- seq(x.range[1],x.range[2],length=length(colors))
      ##
      ## cat("z:\n");print(z)
      ## cat("col:\n");print(colors)
      ## cat("breaks:\n");print(breaks)
      
      image(z=matrix(z,ncol=1),
            col=colors,
            breaks=breaks,
            xaxt="n",
            yaxt="n",
            xlab=key.xlab,
            ylab="",
            main=""
            )
      par(usr=c(0,1,0,1))
      lv <- pretty(breaks)
      xv <- .scale.x(as.numeric(lv),x.range[1],x.range[2])
      axis(1,at=xv,labels=lv,cex.axis=cex.key*1)
      
      if(density.info=="density"){
        ## Experimental : also plot density of data
        dens <- density(x,adjust=densadj,na.rm=TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- .scale.x(dens$x,x.range[1],x.range[2])
        lines(dens$x,dens$y / max(dens$y) * 0.95,col=denscol,lwd=1)
        axis(2,at=pretty(dens$y)/max(dens$y) * 0.95,pretty(dens$y),cex.axis=cex.key*1)
        ##title("Color Key and Density",cex.lab=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,"Density",line=0.75,cex=cex.key)
      } else if(density.info=="histogram"){
        h <- hist(x,plot=FALSE,breaks=breaks)
        hx <- .scale.x(breaks,x.range[1],x.range[2])
        hy <- c(h$counts,h$counts[length(h$counts)])
        lines(hx,hy/max(hy)*0.95,lwd=1,type="s",col=denscol)
        axis(2,at=pretty(hy)/max(hy)*0.95,pretty(hy),cex.axis=cex.key*1)
        ##title("Color Key and Histogram",cex.main=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,key.ylab,line=0.75,cex=cex.key)
      } else{
        title(key.title,cex.main=cex.key,font.main=1)
      }
    } else{
      .plot.text()
    }


    ## 9)
    ##print(colnames(x))
    ##print(rownames(x))
    ##print(rownames(x.ori)[rowInd])
    if(plot.row.individuals) {
      .plot.text("Row\nIndividuals",cex=cex.text,bg="white")
      for (i in 1:ir) {
        ##.plot.text()
        tmp <- plot.row.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 10)
    if(plot.col.individuals) {
      .plot.text("Column\nIndividuals",cex=cex.text,bg="white",srt=90)
      for (i in 1:ic) {
        ##.plot.text()
        tmp <- plot.col.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }
    
    ## 11) for RowPlot from bottom
    if (plot.row.clusters){
      .plot.text("Row\nClusters",cex=cex.text,bg="white")
      
      tmp <- plot.row.clusters.list[[1]]
      row.data <- row.data[rowInd]
      for (i in unique(row.clusters)){
        i.x <- row.data[row.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Row group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }
    
    ## 12) for ColPlot from left
    if (plot.col.clusters){
      .plot.text("Col\nClusters",cex=cex.text,bg="white",srt=90)
      
      tmp <- plot.col.clusters.list[[1]]
      ##print(tmp)
      print(col.clusters)
      col.data <- if(revC) col.data[rev(colInd)] else col.data[colInd]
      for (i in unique(col.clusters)){
        i.x <- col.data[col.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Col group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }

    
    ## 13)
    if(plot.row.clustering) {
      .plot.text("Row\nClustering",cex=cex.text,bg="white")
      
      for (i in 1:cr) {
        ##.plot.text()
        tmp <- plot.row.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    
    ## 14)
    if(plot.col.clustering) {
      .plot.text("Column\nClustering",cex=cex.text,bg="white",srt=90)
      
      for (i in 1:cc) {
        ##.plot.text()
        tmp <- plot.col.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }
    
    ## 15) text
    if (!.invalid(text.box) & if.plot.info){
      .plot.text(text.box,cex=cex.text,bg="gray75")
    } else{
      if (flag.text){
        .plot.text()
      }
    }
 
  }
  cat("8.dim(x):",dim(x),"\n")
  
  ret <-
    list(x.ori=x.ori,
         x=x.save,
         rowInd=rowInd,colInd=colInd,
         row.clusters=row.clusters,col.clusters=col.clusters,
         dist.row=dist.row,dist.col=dist.col,
         hclust.row=hclust.row,hclust.col=hclust.col,
         kr=kr,kc=kc
         )
  class(ret) <- c("hclustering",class(ret))

  cat("DONE!\n\n")
  invisible(ret)
}




##' Get row or column lines of separation for \code{heatmap.3} according to clusters
##'
##' Get row or column lines of separation for \code{heatmap.3} according to clusters
##' @title Get row or column lines of separation for heatmap.3
##' @param clusters a numerical vector, indicating the cluster labels of observations.
##' @param type string, one of the following: \code{c("row","column","both")}
get.sep <-
  function(clusters,type=c("row","column","both"))
{
  ##   if(!is.numeric(clusters) | !(.is.grouped(clusters))){
  ## stop("`clusters' should be a grouped numeric vector.")
  ##   }
  tmp.whichna <- which(is.na(clusters))
  tmp.which <- which(!duplicated(clusters))
  
  tmp.sep <- data.frame(start=tmp.which,end=c(tmp.which[-1],length(clusters)+1)-1)
  tmp.sep2 <- tmp.sep[tmp.sep$start<=tmp.sep$end,]

  ## lines of separation 
  sepList <- list()
  if (type=="row"){
    xmax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(0,tmp.sep2[i.s,1]-1,xmax,tmp.sep2[i.s,2])
    }
  } else if (type=="column"){
    ymax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,0,tmp.sep2[i.s,2],ymax)
    }
  } else if (type=="both"){
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,tmp.sep2[i.s,1]-1,tmp.sep2[i.s,2],tmp.sep2[i.s,2])
    }
  }
  sepList
}


