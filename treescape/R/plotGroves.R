#'
#' Scatterplot of groups of trees
#'
#' This function displays the scatterplot of the Multidimensional
#' Scaling (MDS) output by treescape, superimposing group information
#' (derived by \code{\link{findGroves}}) using colors.
#'
#' This function relies on \code{\link[adegraphics]{s.class}}
#' from the \code{adegraphics} package.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @importFrom adegraphics s.class
#' @importFrom adegraphics s.label
#' @importFrom adegraphics s1d.barchart
#' @importFrom adegraphics insert
#' @importFrom adegenet funky
#' @importFrom adegenet bluepal
#' @importFrom adegenet transp
#'
#' @param x a list returned by \code{\link{findGroves}} or a MDS with class \code{dudi}
#' @param groups a factor defining groups of trees
#' @param xax a number indicating which principal component to be used as 'x' axis
#' @param yax a number indicating which principal component to be used as 'y' axis
#' @param type a character string indicating which type of graph to use
#' @param col.pal a color palette to be used for the groups
#' @param bg the background color
#' @param lab.show a logical indicating whether labels should be displayed
#' @param lab.col a color for the labels
#' @param lab.cex the size of the labels
#' @param lab.optim a logical indicating whether label positions should be optimized to avoid overlap; better display but time-consuming for large datasets
#' @param point.cex the size of the points
#' @param scree.pal a color palette for the screeplot
#' @param scree.size a size factor for the screeplot, between 0 and 1
#' @param scree.posi either a character string or xy coordinates indicating the position of the screeplot.
#' @param ... further arguments passed to \code{\link{s.class}}
#'
#' @return
#' An \code{adegraphics} object (class: \code{ADEgS})
#'
#' @seealso
#' \code{\link{findGroves}} to find any clusters in the tree landscape
#' \code{\link[adegraphics]{s.class}}
#' 
#'
#' @examples
#'
#' \dontrun{
#' if(require("adegenet") && require("adegraphics")){
#' ## load data
#' data(woodmiceTrees)
#'
#' ## run findGroves: treescape+clustering
#' res <- findGroves(woodmiceTrees, nf=5, nclust=6)
#'
#' ## basic plot
#' plotGroves(res)
#'
#' ## adding labels
#' plotGroves(res, lab.show=TRUE)
#'
#' ## customizing
#' plotGroves(res, lab.show=TRUE,
#' bg="black", lab.col="white", scree.size=.35)
#'
#' ## customizing
#' plotGroves(res, type="ellipse", lab.show=TRUE,
#' lab.optim=FALSE, scree.size=.35)
#'
#' ## example with no group information
#' plotGroves(res$treescape$pco)
#'
#' ## adding labels
#' plotGroves(res$treescape$pco, lab.show=TRUE, lab.cex=2)
#'
#' }
#' }
#'
#' @export
plotGroves <- function(x, groups=NULL, xax=1, yax=2,
                        type=c("chull","ellipse"), col.pal=funky, bg="white",
                        lab.show=FALSE, lab.col="black", lab.cex=1, lab.optim=TRUE,
                        point.cex=1, scree.pal=NULL, scree.size=.2,
                        scree.posi=c(.02,.02), ...){
    ## HANDLE ARGUMENTS ##
    ## checks
    type <- match.arg(type)
    if(is.null(scree.pal)) scree.pal <- function(n) rev(bluepal(n))

    ## x is a list returned by findGroves
    if(is.list(x) && !is.data.frame(x) && !inherits(x,"dudi")){
        if(is.null(x$treescape)) stop("if x is a list, it should contain a slot $treescape")
        groups <- x$groups
        x <- x$treescape$pco
    }

    ## x is a dudi object
    if(inherits(x,"dudi")){
        eig <- x$eig
        x <- x$li
    }

    ## groups missing - just s.label
    if(is.null(groups)) {
        ## with labels
        if(lab.show){
            out <- s.label(x, xax=xax, yax=yax,
                           plabels=list(optim=lab.optim, col=lab.col, cex=lab.cex),
                           ppoints=list(cex=point.cex),
                           pbackground.col=bg,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        } else {
            ## just points
            out <- s.label(x, xax=xax, yax=yax,
                           plabels=list(optim=FALSE,cex=0),
                           ppoints=list(cex=point.cex, col=lab.col),
                           pbackground.col=bg,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        }
    } else {
        ## if groups are provided
        if(!is.factor(groups)) groups <- factor(groups)
        n.lev <- length(levels(groups))


        ## MAKE GRAPH ##
        ## base scatterplot
        if(type=="chull"){
            out <- s.class(x, xax=xax, yax=yax, fac=groups, col=col.pal(n.lev),
                           ellipseSize=0, chullSize=1,
                           pbackground.col=bg,
                           ppoints.cex=point.cex,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        }
        if(type=="ellipse"){
            out <- s.class(x, xax=xax, yax=yax, fac=groups, col=col.pal(n.lev),
                           ellipseSize=1,
                           pbackground.col=bg,
                           ppoints.cex=point.cex,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        }

        ## add labels
        if(lab.show){
            out <- out + s.label(x, plabel.optim=lab.optim, plabel.col=lab.col,
                                 ppoints.cex=0, plabels.cex=lab.cex)
        }
    }
    ## add inset
    if(!is.null(scree.posi[1]) && !is.na(scree.posi[1]) && scree.posi[1]!="none"){
        screeplot <- s1d.barchart(c(0,eig), p1d.horizontal=FALSE, ppolygons.col=scree.pal(length(eig)+1),
                                  pbackground=list(col=transp("white"), box=TRUE),
                                  layout.width=list(left.padding=2),
                                  pgrid.draw=FALSE, plot=FALSE)
        out <- insert(screeplot, out, posi=scree.posi, ratio=scree.size, plot=FALSE)

    }


    ## RETURN ##
    return(out)
} # end plotGroves






######################
### ScatterD3 version
######################

#'
#' Scatterplot of groups of trees using \code{scatterD3}
#'
#' This function displays the scatterplot of the Multidimensional
#' Scaling (MDS) output by treescape, superimposing group information
#' (derived by \code{\link{findGroves}}) using colors.
#' \code{scatterD3} enables interactive plotting based on d3.js, including zooming, panning and fading effects in the legend.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @import scatterD3
#' @importFrom adegenet funky
#' @importFrom adegenet bluepal
#' @importFrom adegenet transp
#'
#' @param x a list returned by \code{\link{findGroves}} or a MDS with class \code{dudi}
#' @param groups a factor defining groups of trees. If x is a list returned by \code{\link{findGroves}} these will be detected automatically.
#' @param xax a number indicating which principal component to be used as 'x' axis
#' @param yax a number indicating which principal component to be used as 'y' axis
#' @param treeNames if a list of tree names or labels are given, these will be plotted alongside the points. Their size can be altered using \code{labels_size} - see \code{?scatterD3} for more information.
#' @param xlab the label for the 'x' axis. Defaults to use the value of 'xax'
#' @param ylab the label for the 'y' axis. Defaults to use the value of 'yax'
#' @param symbol_var a factor by which to vary the symbols in the plot
#' @param ... further arguments passed to \code{\link{scatterD3}}
#'
#' @return
#' A \code{scatterD3} plot
#' 
#' 
#' @seealso
#' \code{\link{findGroves}} to find any clusters in the tree landscape
#' 
#'
#' @examples
#'
#' \dontrun{
#' if(require("adegenet") && require("scatterD3")){
#' ## load data
#' data(woodmiceTrees)
#'
#' ## run findGroves: treescape+clustering
#' res <- findGroves(woodmiceTrees, nf=5, nclust=6)
#'
#' ## basic plot
#' plotGrovesD3(res)
#'
#' ## adding tree labels
#' plotGrovesD3(res, treeNames=1:201)
#'
#' ## customizing: vary the colour and the symbol by group
#' plotGrovesD3(res, symbol_var=res$groups)
#'
#' ## example with no group information
#' plotGrovesD3(res$treescape$pco)
#' }
#' }
#'
#' @export
plotGrovesD3 <- function(x, groups=NULL, xax=1, yax=2, treeNames=NULL, symbol_var=NULL,
                         xlab=paste0("Axis ",xax), ylab=paste0("Axis ",yax), ...){
  ## HANDLE ARGUMENTS ##
  ## checks

  ## x is a list returned by findGroves
  if(is.list(x) && !is.data.frame(x) && !inherits(x,"dudi")){
    if(is.null(x$treescape)) stop("if x is a list, it should contain a slot $treescape")
    groups <- x$groups
    x <- x$treescape$pco$li
  }
  
  ## x is a dudi object
  if(inherits(x,"dudi")){
    eig <- x$eig
    x <- x$li
  }

  scatterD3(x[,xax],x[,yax], lab=treeNames, col_var=groups, symbol_var=symbol_var,
            xlab=xlab, ylab=ylab, ...)
} # end plotGrovesD3