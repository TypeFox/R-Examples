##' Computing Generalized Minimum Distance Matrix
##'
##' Computing Generalized Minimum Distance Matrix
##' @name gmdm
##' @aliases gmdm gmdm2dist gmdm_dist print.gmdm
##' @title Generalized Minimum Distance Matrix
##' @usage
##'
##' gmdm(data,labels,pseudocount=0,sliding=TRUE,resolution=1)
##' 
##' ## S3 method for class `gmdm'
##' \method{print}{gmdm}(x, ...)
##'
##' ## convert a `gmdm' object into a `dist' object
##' gmdm2dist(m, diag=FALSE, upper=FALSE)
##'
##' ## compute GMDM and convert into a `dist' object
##' gmdm_dist(data, diag=FALSE, upper=FALSE, ...)
##' 
##' @param data a list of numeric vectors, a numeric matrix or data.frame
##' @param x a \code{gmdm} object.
##' @param m a \code{gmdm} object.
##' @param labels a character vector of the same length of x, giving the names of the numeric vectors.
##' @param pseudocount a numeric value to be allocated for each position to reduce bias;
##' by default \code{pseudocount = 0}.
##' @param sliding logical, indicating whether sliding is allowed or not for an optimal solution; 
##' by default \code{sliding = TRUE}.
##' @param resolution relative resolution, numeric (>=1), changing the size of the bin by
##' multiplying the value. A larger value (lower resolution) is more computational efficiet but
##' missing details.
##' @param diag logical value indicating whether the diagonal of the distance matrix should be 
##' printed by \code{print.dist}.
##' @param upper logical value indicating whether the upper triangle of the distance matrix should be
##' printed by \code{print.dist}.
##' @param ... arguments to be passed to method
##' @return \code{gmdm} returns an object of class \code{gmdm}, a list with components
##'
##' labels: a string vector, giving the names of distributions
##'
##' data.ori: a list of numeric vectors, giving the original input
##'
##' data: a list of numeric vectors, giving the normalized version of the original input
##' 
##' dm: a numeric numeric, the pairwise distance matrix of \emph{GM-Distances}
##' 
##' gap.pair: a numeric matrix, giving the gap pair of each alignment per row:
##' i.e. relative shifts between distributions of the optimal hit
##' 
##' sliding: logical, indicating whether sliding is performed
##' 
##' pseudocount: a numeric value that is allocated at each position in addition to original values
##' @references See \code{citation("GMD")}
##' @seealso \code{\link{plot.gmdm}}, \code{\link{gmdp}}
##' @keywords classes
gmdm <- 
  function(data,
           labels,
           pseudocount=0,
           sliding=TRUE,
           resolution=1
           )
{
  x <- data
  
  if (!is.list(x) & !is.matrix(x) & !is.data.frame(x)){
    stop("`x' should be an object of class `list', `matrix' or `data.frame'.")
  }

  if (is.matrix(x)){
    x <- data.frame(x)
  }

  if (is.data.frame(x)){
    x <- as.list(data.frame(t(x)))
  }

  
  if (!all(unlist(lapply(x,is.numeric)))){
    stop("`x' should be a list of `numeric' vectors.")
  }

  ## labels ##
  if (.invalid(labels)){
    labels <- NULL
  }

  if (is.null(labels)){
    if (!is.null(names(x)) & all(names(x)!="")){
      labels <- names(x)
    } else{
      labels <- as.character(1:length(x))
    }
  }
  
  ##
  meta <- list()
  N <- length(x)
  dm <- matrix(0, ncol=N, nrow=N)
  colnames(dm) <- rownames(dm) <- labels
  gap.pair <- matrix(0, ncol=4, nrow=N*(N-1)/2)
  shiftedLength.pair <- matrix(0, ncol=4, nrow=N*(N-1)/2)
  colnames(gap.pair) <- c("V1","V2","Gap_V1","Gap_V2")
  colnames(shiftedLength.pair) <- c("V1","V2","L_V1","L_V2")
  n <- 0
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      n <- n+1
      ij.x <- gmdp(x[[i]],x[[j]],pseudocount=pseudocount,sliding=sliding,resolution=resolution)
      dm[i,j] <- dm[j,i] <- ij.x
      ij.x.meta <- attr(ij.x,"meta")
      gap.pair[n,] <- c(labels[i],labels[j],ij.x.meta$gap.pair[1,])
      shiftedLength.pair[n,] <- c(labels[i],labels[j],ij.x.meta$gap.pair[1,]+c(length(x[[i]]),length(x[[j]])))
      if (nrow(ij.x.meta$gap.pair)>1){
        warning(sprintf("There are multiple (equally good) optimal hits/alignment
between %s and %s;
the first one will be saved.",labels[i],labels[j]))
      }
    }
  }
  gap.pair <- data.frame(gap.pair,stringsAsFactors=FALSE)
  gap.pair$Gap_V1 <- as.numeric(gap.pair$Gap_V1)
  gap.pair$Gap_V2 <- as.numeric(gap.pair$Gap_V2)
  
  length.pair <- data.frame(shiftedLength.pair,stringsAsFactors=FALSE)
  length.pair$L_V1 <- as.numeric(length.pair$L_V1)
  length.pair$L_V2 <- as.numeric(length.pair$L_V2)
  
  meta$labels <- labels
  meta$data.ori <- x
  meta$resolution <- resolution
  meta$data <- lapply(x,function(e) {ee=.resolveHist(e,resolution);ee/sum(ee)})
  meta$gap.pair <- gap.pair
  meta$length.pair <- length.pair
  meta$sliding <- sliding
  meta$pseudocount <- pseudocount

  ## return
  ret <- dm
  class(ret) <- c("gmdm",class(ret))
  attr(ret, "meta") <- meta
  return(ret)
}


print.gmdm <-
  function(x, ...)
{
  attr(x,"meta") <- c()
  class(x) <- class(x)[class(x)!="gmdm"]
  print(x, ...)
}


gmdm2dist <- 
  function(m, diag=FALSE, upper=FALSE) 
{
  if (!.is.gmdm(m)){
    stop("`m' should be an object of class `gmdm'.")
  }
  class(m) <- class(m)[class(m)!="gmdm"]
  d <- as.dist(m)
  attr(d, "Size") <- nrow(m)
  attr(d, "Labels") <- dimnames(m)[[1L]]
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- "gmdm"
  attr(d, "call") <- match.call()
  return(d)
}


gmdm_dist <- 
  function(data, diag=FALSE, upper=FALSE, ...) 
{
  gmdm2dist(gmdm(data,...),diag=diag,upper=upper)
}




##' S3 method for class \code{gmdm}
##'
##' S3 method for class \code{gmdm}
##' @title S3 method for class `gmdm'
##' 
##' @export plot.gmdm
##' 
##' @param x an object of class \code{gmdm}.
##' @param labels a string vector of the same length as \code{x$data},
##' giving the names of the numeric vectors in \code{x$data}.
##' @param colors the colors of the discrete distributions; the default is \emph{"Dark2" colors
##' in ColorBrewer palettes} if not specified.
##' @param type type of plot, as in \code{help("plot", package="graphics")}.
##' @param main an overall title for the plot. See \code{help("title", package="graphics")};
##' the default title is used if not specified.
##' @param ylab a title for the y axis. See \code{help("title", package="graphics")}.
##' @param xlab a title for the x axis. See \code{help("title", package="graphics")}.
##' @param label.length.max numeric, giving the maximum string width allowed in diagonal labels.
##' @param label.line.max numeric, giving the maximum number of lines allowed in diagonal labels.
##' @param cex.text a numerical value giving the amount by which plot text should be magnified
##' relative to the default.
##' @param cex.tickmark a numerical value giving the amount by which tickmarks should be magnified
##' relative to the default.
##' @param if.plot.new logical, indicating whether to start a new plot device.
##' @param ... arguments to be passed to methods, see \code{gmdp}.
##' @references See \code{help(GMD)}
##' @seealso \code{\link{gmdm}}, \code{\link{gmdp}}
##' @keywords methods hplot
##' @examples
##' ## ------------------------------------------------------------------------
##' ## Example1: CAGE
##' ## ------------------------------------------------------------------------
##' require("GMD") # load library
##' data(cage)     # load data
##' 
##' ## construct a distance matrix and visualize it
##' short.labels <- gsub("(.+) \\(.+","\\1",names(cage)) # get short labels
##' x <- gmdm(cage[1:6],labels=short.labels[1:6])
##' plot(x)
##'
##' \dontrun{
##' ## ------------------------------------------------------------------------
##' ## Example2: ChIP-seq
##' ## ------------------------------------------------------------------------
##' data(chipseq_mES)   # load data
##' data(chipseq_hCD4T) # load data
##' 
##' ## pairwise distance and alignment based on GMD metric
##' plot(gmdm(chipseq_mES,sliding=FALSE))
##' 
##' ## clustering on spatial distributions of histone modifications
##' x <- gmdm(chipseq_hCD4T,sliding=FALSE,resolution=10)
##' heatmap.3(x,revC=TRUE)
##' }
plot.gmdm <- 
  function(x,
           ##
           labels,
           colors,
           type=NULL,
           ##
           main,
           ylab="Fraction",
           xlab="Position",
           ## maximum size of diagonal labels
           label.length.max=8,
           label.line.max=3,
           ## style
           cex.text=1,
           cex.tickmark=0.75,
           ##
           if.plot.new=TRUE,
           ...
           ) 
{
  if (!.is.gmdm(x)){
    stop("`x' should be an object of class `gmdm'.")
  }

  ## meta ##
  x.meta <- attr(x,"meta")
  
  ## params ##
  if (.invalid(colors)){
    colors <- .brewer.pal.Dark2(8)
  }

  if (.invalid(labels)){
    labels <- x.meta$labels
  }
  
  ## dm ##
  dm <- x

  N <- nrow(dm)
  enlarge <- 0.75
  cex.text <- cex.text*N/4*enlarge

  
  ## labels for plot ##
  nchar.max <- max(nchar(labels))
  label.size.max <- label.length.max*label.line.max
  if(nchar.max>label.size.max){
    warning(sprintf("The number of characters in each `label` should not exceed %s.
The first %s characters are kept.",label.size.max,label.size.max)
            )
  }
  labels.plot <- sapply(labels,substr,start=1,stop=label.size.max)
  labels.plot <- sapply(labels.plot,.wordwrap,len=label.length.max)

  ## colors ##
  if(length(colors) < N){
    warning("The length of `colors' is shorter than that of `x' and the `colors' are recycled.")
    ##tmp <- cbind(colors,1:N)
    colors <- rep(colors,ceiling(N/length(colors)))
  }

  ## Main title ##
  if (.invalid(main)){
    main <-
      sprintf("Optimal alignments among distributions (%s sliding)",
              ifelse(x.meta$sliding,"with","without")
              )
  }
  
  ## plot ## 
  ylim.max <- max(unlist(lapply(x.meta$data,max)))
  xlim.max <- max(x.meta$length.pair[,c(-1,-2)])

  ## print(sprintf("ylim.max=%s",ylim.max))
  ## print(sprintf("xlim.max=%s",xlim.max))
  
  ## colors for distance boxs
  dm.range <- NULL # currently only an internal option
  if(is.null(dm.range)){
    dm.range <- range(dm)
  }
  
  ##colors_for_value <- .grays16
  ##colors_for_value <- gray.colors(min(c(N*(N-1)/2,16)),start=0.9,end=0.1)
  colors_for_value <- gray.colors(16,start=0.9,end=0.1)
  breaks <- length(colors_for_value) + 1

  if(length(breaks)==1){
    ##breaks <- seq(min(dm, na.rm=na.rm), max(dm,na.rm=na.rm), length=breaks)
    breaks <-
      seq(min(min(dm,na.rm=TRUE),dm.range[1]),
          max(max(dm,na.rm=TRUE),dm.range[2]),
          length=breaks
          )
  }
  
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)

  dm <- ifelse(dm<min.breaks, min.breaks, dm)
  dm <- ifelse(dm>max.breaks, max.breaks, dm)

  dm.colors.index <- apply(dm, c(1,2), function(e,breaks){ min(max(which(breaks<=e)),length(breaks)-1)}, breaks=breaks)
  dm.colors <- apply(dm.colors.index, c(1,2), function(i){colors_for_value[i]})
  
  ## 
  if (if.plot.new) {
    dev.new(width=N*2*enlarge,height=N*2*enlarge)
  }

  ##
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  par(mfrow=c(N,N),
      mgp=c(0,0.2,0),
      mar=c(1,1,1,1),
      oma=c(3,3,5,1))

  n <- 0
  for (i in 1:(N-1)){
    par(mfg=c(i,i))
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
    text(1,1,labels.plot[i],cex=cex.text/N*3.6*1.2,col=colors[i])
    for (j in (i+1):N) {
      ##cat("i=",i,", j=",j,"\n",sep="")
      n <- n+1
      par(mfg=c(i,j))
      ij.x.meta <- list()
      ij.x.meta$labels <- x.meta$labels[c(i,j)]
      ij.x.meta$v1.ori <- x.meta$data.ori[[i]]
      ij.x.meta$v2.ori <- x.meta$data.ori[[j]]
      ij.x.meta$resolution <- x.meta$resolution
      ij.x.meta$v1 <- x.meta$data[[i]]
      ij.x.meta$v2 <- x.meta$data[[j]]
      ij.x.meta$sliding <- x.meta$sliding
      ij.x.meta$pseudocount <- x.meta$pseudocount
      ij.x.meta$gap.pair <- matrix(as.numeric(matrix(x.meta$gap.pair[n,3:4],nrow=1)),nrow=1)
      ij.x.meta$n.hit <- 1
      ij.x <- dm[i,j]
      class(ij.x) <- c("gmdp",class(ij.x))
      attr(ij.x,"meta") <- ij.x.meta
      plot.gmdp(x=ij.x,
                labels=ij.x.meta$labels,
                colors=colors[c(i,j)],
                main="",
                if.text.gmd=FALSE,
                if.text.gap=FALSE,
                legend.lab=NA,
                xlab="",
                ylab="",
                ylim=c(0,ylim.max),
                xlim=c(1,xlim.max),                
                type=type,
                cex.tickmark=cex.tickmark,
                tcl=-0.1,
                ...
                )
      ##
      par(mfg=c(j,i))
      ij.text <-
        sprintf("GMD=%.3f\nGap=%s",ij.x,sprintf("(%s,%s)",ij.x.meta$gap.pair[1],ij.x.meta$gap.pair[2]))
      plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n")
      rect(xleft=0, ybottom=0, xright=1, ytop=1, col=dm.colors[i,j], border=NA)
      text(0.5,0.5,ij.text,cex=cex.text/N*3.6*1.05,col=.setTextContrastColor(dm.colors[i,j]))
    }
  }
  par(mfg=c(N,N))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  text(1,1,labels.plot[N],cex=cex.text/N*3.6*1.2,col=colors[N])
  mtext(ylab,line=0.5,side=2,outer=TRUE,cex=cex.text*0.75)
  mtext(xlab,line=0.5,side=1,outer=TRUE,cex=cex.text*0.75)
  mtext(text=main,line=1.0,side=3,outer=TRUE,cex=cex.text*0.75,adj=0.5,font=2)
  
  invisible()
  
}
