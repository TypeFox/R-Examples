#' Create 2D PCA Plots with Bootstrapping
#'
#' pcaBootPlot draws a 2D PCA plot using the first 2 principal components using
#'   the original and bootstrapped data to give some sense of variability.
#'
#' @param data A data.frame where the first column is named "ID" and contains IDs
#'   for each item measured. Measurements for each sample are in subsequent columns.
#' @param groups The default value is \strong{\code{NULL}}.\cr\cr
#'   If you want use different colors
#'   and shapes to deliniate the samples into groups, you can specify the grouping
#'   with this argument. Currently there is a limit of 9 different groups.\cr\cr
#'   For example, if you have three consecutive columns of
#'   "untreated" samples followed by three consecutive columns of "treated" samples,
#'   you can set this argument to c(1,1,1,2,2,2), and the untreated samples will be
#'   red circles and the treated samples will be blue triangles.
#' @param min.value The default value is \strong{1}.\cr\cr
#'   This allows you to filter out rows (entries) that will not conribute
#'   to the PCA. For example, if you are performing PCA on RNA-seq data, you
#'   may wish to filter out genes with less than 1 read per sample, 1 read per
#'   group or 1 read overall. If you set \code{all.min.value} to \code{TRUE},
#'   it will filter entries where at least one sample has less than
#'   \code{min.value}. If you do not set \code{all.min.value} to \code{TRUE},
#'   then filtering will be performed by group if \code{groups} are specified.
#'   In this case, an entry will be filtered out if one or more groups have
#'   less than \code{min.value}.\cr\cr
#'   If \code{groups} are not specified, then only entries where all samples
#'   have less than \code{min.value} will be removed from the analysis.\cr\cr
#'   \code{groups} will also effect filtering based on \code{min.value}.
#'   See that part of the documentation for details.
#' @param all.min.value This parameter, set to either \code{TRUE} or \code{FALSE},
#'   affects
#'   \code{min.value}. See the documentation for \code{min.value} for more
#'   details.
#' @param num.boot.samples The default value is \strong{100}. The number of bootstrap
#'   iterations to be performed.
#' @param log2.transform The default value is \strong{\code{TRUE}}. Should the data be log2
#'   transformed or not?
#' @param pdf.filename If you wish to save the the graph as a PDF, you may use
#'   this argument to specify the filename.
#' @param pdf.width If you specify a value for \code{pdf.filename}, you can specify
#'   a width for the saved graph. The default value is \strong{6} inches.
#' @param pdf.height If you specify a value for \code{pdf.filename}, you can
#'   specify a height for the saved graph. The default value is \strong{6} inches.
#' @param draw.legend The default value is \strong{\code{FALSE}}. Should there be a legend
#'   in the graph?
#' @param legend.names If \code{draw.legend} is \code{TRUE}, you can specify the
#'   names of the groups listed in the legend.
#' @param legend.x,legend.y If \code{draw.legend} is \code{TRUE}, you can specify
#'  the x and y axis coordinate for its top left corner.
#' @param transparency The default value is \strong{77}. This allow you to set how
#'   transparent the bootstrapped symbols are in the graph. Values range from 00
#'   to FF.
#' @param min.x,min.y,max.x,max.y By default, pcaBootPlot automatically
#'   determines limits for the
#'   x and y axes. Use this option to override this behavior.
#' @param correct.inversions The default value is \strong{\code{TRUE}}. Some of the
#' boostrapped PCAs may have their axes inverted. pcaBootPlot can try to correct
#' for this by ensuring that the PCA loading values are positively correlated
#' with the orginal dataset.
#' @param confidence.regions The default value is \strong{\code{FALSE}}. This option
#' will draw circles that contain confidence.size of the bootstrapped values.
#' @param confidence.size The default value is \strong{0.95}. A value betweeo 0 and 1 - the
#' proportion of bootstrapped points that need to be within the confidence regions.
#' @param step.size The default value is \strong{0.1}. This option determines how
#' the radii
#' for confidence regions are increased each iteration when trying to contain
#' confidence.size of the bootstrapped samples.
#' @param trim.proportion The default value is \strong{0.0}. This is the proportion
#' of entries that should be removed from the plot based on the size of the
#' the confidence regions. This should be a value between 0 and 1. For example,
#' if you set it to 0.1, then the top 10% of entries with the lagest confidence
#' regions will be removed from the plot.
#' @param return.samples The default value is \strong{\code{FALSE}}. If this is set
#' to \code{TRUE} then the program will return the names of the samples that were
#' included in the plot. This can be useful if \code{trim.proportion} > 0.
#' @param use.prcomp The default value is \strong{\code{FALSE}}. Usually,
#' pcaBootPlot uses FactoMineR to process samples. However, this can be
#' unnecessarily slow if there are less than 50 samples. By setting use.prcomp to
#' \strong{\code{TRUE}}, it will use prcomp() to process samples and will,
#' most likely, run much faster.
#'
#' @examples
#'
#' sample1=rnorm(n=100, mean=100, sd=10)
#' sample2=jitter(sample1, factor=10, amount=10)
#' sample3=rnorm(n=100, mean=100, sd=10)
#'
#' data <- data.frame(ID=c(1:100), sample1, sample2, sample3)
#'
#' pcaBootPlot(data, log2.transform = FALSE)
#'

#' @export
pcaBootPlot <- function(data=NULL, groups=NULL,
                        min.value=1, all.min.value=FALSE,
                        num.boot.samples=100, log2.transform=TRUE,
                        pdf.filename=NULL,
                        pdf.width=6,
                        pdf.height=6,
                        draw.legend=FALSE, legend.names=NULL,
                        legend.x=NULL, legend.y=NULL,
                        transparency=77,
                        min.x=NULL, max.x=NULL, min.y=NULL, max.y=NULL,
                        correct.inversions=TRUE,
                        confidence.regions=FALSE,
                        confidence.size=0.95,
                        step.size=0.1,
                        trim.proportion=0,
                        return.samples=FALSE,
                        use.prcomp=FALSE) {

  #library(RColorBrewer)


  if(is.null(data)) {
    return("You must provide a data.frame for the data parameter")
  }
  num.samples <- (ncol(data)-1)
  cat("Performing PCA on", num.samples, "samples\n")

  use.facto <- TRUE
  if (use.prcomp == TRUE) {
    use.facto <- FALSE
  } else if (num.samples < 50) {
    cat("\nUsing FactoMineR for analysis. However you may be able to speed up\n")
    cat("  computation by setting use.prcomp to TRUE\n\n")
  }

  ##
  ## first, we need to find duplicate entries and average the values for them.
  ##
  dup.indices <- duplicated(data$ID)
  dup.IDs <- data[dup.indices,]$ID
  avg.fpkms <- data.frame()
  cat(paste("There are ", length(dup.IDs), " duplicated entries", sep=""), "\n")
  if (length(dup.IDs) > 0) {
    cat("Averaging duplicated entries...\n")
    for (ID in levels(as.ordered(dup.IDs))) {
      avg.ID <- data.frame(ID=ID,
                           t(data.frame(colMeans(data[data$ID == ID,2:ncol(data)]))))
      row.names(avg.ID) <- 1
      avg.fpkms <- rbind(avg.fpkms, avg.ID)
      data <- data[data$ID != ID,]
    }
    data <- rbind(data, avg.fpkms)
  }

  ## convert data to a matrix
  IDs <- data[,1]
  data <- as.matrix(data[,2:ncol(data)])
  row.names(data) <- IDs

  ## DEBUG
  ##print(data[1:4,1:4])

  ##
  ## Only keep entries with a minimum value (per all samples, or per group)
  ##
  cat("Filtering entries based on min.val and groups...\n")
  if (!is.null(groups)) {
    data.factors <- as.data.frame(table(factor(groups)))
  }

  if (is.null(groups) || all.min.value) {
    if (all.min.value) {
      keep <- (apply(data, 1, min) > min.value)
    } else {
      keep <- (rowSums(data[,1:ncol(data)]) > min.value)
    }
  } else {
    all.keeps <- list()
    list.index <- 1
    for(group.id in data.factors[,1]) {
      group.keep <- which(rowSums(data[,which(groups==group.id)]) > min.value)
      all.keeps[[list.index]] <- group.keep
      list.index = list.index+1
    }
    keep <- Reduce(intersect, all.keeps)
  }
  data <- data[keep,]

  num.genes <- nrow(data)

  if (is.null(groups)) {
    if (all.min.value) {
      cat("   ", paste(num.genes, "entries had all samples with values >", min.value), "\n")
    } else {
      cat("  ", paste(num.genes, "entries had at least one sample with value >", min.value), "\n")
    }
  } else {
    cat("  ", paste(num.genes, "entries had at least one sample per group with values >", min.value), "\n")
  }
  #return(data)

  if(log2.transform) {
    cat("Adding pseudo-counts and log2 transforming the data...", "\n")
    cat("   You can turn this off by setting log2.transform to FALSE.", "\n")
    data <- log2(data+1)
  }
  #return(data)

  ## There are a bunch of ways to do PCA in R, however, I've found that
  ## FactoMineR is very fast with very large datasets (800+ samples).
  ##
  ## originally I just used the built in program, "prcomp" and that worked great
  ## until sample size > 100.

  pca.data <- data.frame()
  pc1 <- vector()
  pc2 <- vector()
  pc1.names <- vector()
  pc2.names <- vector()
  pca.var.per <- vector()
  if (use.facto) {
    cat("Using FactoMineR for analysis\n")
    pca <- FactoMineR::PCA(t(data), ncp=5, graph=FALSE, scale.unit=TRUE)
    ## arguments:
    ## t(data)    - transposed data so that samples are rows, genes are columns
    ## ncp        - the maximum number of principal components calculated
    ## graph      - draw PCA plots?
    ## scale.unit - should the variables be scaled to unit variance (yes!)
    ##
    ## return values:
    ## pca$ind$coord - the values of the "rotated data" (coordinates for each
    ##                dimension in the PCA plot). these are the scores for the
    ##                samples for each PC. That is to say:
    ##                loadings * measurements = score = rotated data = coordinates
    ## pca$eig       - a matrix with three columns:
    ##                 column 1: eigenvalues
    ##                 column 2: percentage of variance per eigenvalue
    ##                 column 3: cumulative percentage of variance
    ##print(pca$ind$coord[,c(1,2)])
    pca.data <- list(x=pca$ind$coord[,c(1,2)])
    #print(head(pca.data))
    #print(head(pca.data$x))

    pc1 <- pca$var$coord[,1]/sqrt(pca$eig[1,1])
    pc2 <- pca$var$coord[,2]/sqrt(pca$eig[2,1])

    pc1.names <- names(pc1)
    pc2.names <- names(pc2)

    pca.var.per <- round(pca$eig[,2], digits=1)
  } else {
    cat("Using prcomp for analysis\n")
    ##### This is all my old "prcomp" code
    pca <- prcomp(t(data), center=TRUE, scale. = TRUE, retx=TRUE)
    pca.data <- list(x=pca$x[,c(1,2)])
    ## arguments:
    ## t(data) - transposed data so that samples are rows, genes are columns
    ## center  - should the values be shifted to be zero centered?
    ## scale.  - should the variables be scaled to unit variance (yes!)
    ## retx    - should the rotated variables should be returned?
    ##
    ## return values:
    ## pca$rotation = the eigenvectors (the loading vectors)
    ## pca$sdev = the standard deviation of each eigen vector
    ##    - note, I'm not sure why the standard deviation is returned, since
    ##      you always want the variance. See below for code to get variances
    ##      and the proportion of variance and the cumulative proportion
    ## pca$x = the values of the "rotated data", these are the scores for the samples
    ##      for each PC. That is to say, loadings * measurements = score = rotated data
    ## pca$center = the value used for centering, if used.
    ## pca$scale = the value used for scaling, if used.
    ##
    pc1 <- pca$rotation[,1]
    pc1.names <- names(pc1)
    pc2 <- pca$rotation[,2]
    pc2.names <- names(pc2)
    ##
    ## draw a scree plot:
    ##plot(pca, las=1)
    ## las=1 make the values on the y-axis easier to read (by making them
    ## perpendicular to the y-axis)
    ##
    ## print out what proportion of the variance each PC accounts for:
    ##summary(pca)
    ##
    ## Now, let's calculate the proportion of the variance for each PC and the
    ## cumulative proportion.
    pca.var <- pca$sdev^2
    pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
    pca.var.cum <- cumsum(pca.var.per)
  }

  boot.points <- data.frame()

  if (num.boot.samples > 0) {
    ## Now let's try to bootstrap the PCA...
    cat("Bootstrapping the PCA at the entry level...\n")

    gene.names <- rownames(data)

    for (i in 1:num.boot.samples) {
      cat("Bootstrap iteration:", i, "\n")

      boot.indices <- sample(x=c(1:num.genes), size=num.genes, replace=TRUE)
      boot.data <- data[boot.indices,]

      boot.pc1 <- vector()
      boot.pc2 <- vector()
      if (use.facto) {

        ## NOTE: FactoMineR::PCA requires all rownames to be unique
        rownames(boot.data) <- c(1:nrow(boot.data))
        pca.boot <- FactoMineR::PCA(t(boot.data), ncp=5, graph=FALSE, scale.unit=TRUE)
        rownames(pca.boot$var$coord) <- gene.names[boot.indices]

        pca.boot.data <- list(x=pca.boot$ind$coord[,c(1,2)])


        boot.pc1 <- pca.boot$var$coord[,1]/sqrt(pca.boot$eig[1,1])
        boot.pc2 <- pca.boot$var$coord[,2]/sqrt(pca.boot$eig[2,1])

      } else {
        pca.boot <- prcomp(t(boot.data), center=TRUE, scale. = TRUE, retx=TRUE)

        pca.boot.data <- list(x=pca.boot$x[,c(1,2)])

        boot.pc1 <- pca.boot$rotation[,1]
        boot.pc2 <- pca.boot$rotation[,2]
      }

      if (correct.inversions) {
        ##
        ## make sure the loadings correlate with the loadings in the
        ## non-bootstrapped PCA
        ##boot.pc1 <- pca.boot$rotation[,1]
        ##boot.pc1 <- pca.boot$var$coord[,1]/sqrt(pca.boot$eig[1,1])
        boot.pc1.names <- levels(factor(names(boot.pc1)))

        #print(head(boot.pc1.names))

        pc1.cor <- cor(pc1[boot.pc1.names], boot.pc1[boot.pc1.names])
        if (pc1.cor < 0) {
          pca.boot.data$x[,1] <- pca.boot.data$x[,1] * -1
          #pca.boot$ind$coord[,1] <- pca.boot$ind$coord[,1] * -1
        }

        ##boot.pc2 <- pca.boot$rotation[,2]
        ##boot.pc2 <- pca.boot$var$coord[,2]/sqrt(pca.boot$eig[2,1])
        boot.pc2.names <- levels(factor(names(boot.pc2)))
        pc2.cor <- cor(pc2[boot.pc1.names], boot.pc2[boot.pc1.names])
        if (pc2.cor < 0) {
          pca.boot.data$x[,2] <- pca.boot.data$x[,2] * -1
          #pca.boot$ind$coord[,2] <- pca.boot$ind$coord[,2] * -1
        }
      }
      boot.points <- rbind(boot.points, pca.boot.data$x[,c(1,2)])
      #print(pca.boot$ind$coord[,c(1,2)])
      #boot.points <- rbind(boot.points, pca.boot$ind$coord[,c(1,2)])
    }
  }

  #print(head(boot.points))

  ##cat("Drawing the 2-D PCA plot with the first two PCs...\n")
  if (exists("data.factors")) {
    if (length(data.factors[,1]) == 2) {
      hex.colors <- RColorBrewer::brewer.pal(n=3, name="Set1")[1:2]
      plot.col <- paste(rep(hex.colors, data.factors[,2]), transparency, sep="")
    } else {
      hex.colors <- RColorBrewer::brewer.pal(n=length(data.factors[,1]), name="Set1")
      plot.col <- paste(rep(hex.colors, data.factors[,2]), transparency, sep="")
    }
  } else {
    hex.colors <- RColorBrewer::brewer.pal(n=3, name="Set1")[1]
    plot.col <- paste(hex.colors, transparency, sep="")
  }

  plot.pch <- 1
  if (!is.null(groups)) {
    plot.pch <- (groups)[1:ncol(data)]
  }

  radii <- NULL
  return.samples.vector <- NULL
  #num.cells <- nrow(pca$ind$coord)
  num.cells <- nrow(pca.data$x)
  if (num.boot.samples > 0) {
    if (confidence.regions | (trim.proportion > 0)) {
      cat("Calculating ", round(confidence.size * 100, digits=2), "% confidence regions\n", sep="")

      boot.points.by.cell <- cbind(boot.points, cell=c(1:num.cells))
      min.points <- round(num.boot.samples * confidence.size)
      #cat("min.points: ", min.points, "\n")

      radii <- vector(length=num.cells)

      for (i in 1:num.cells) {
        circle.center.x <- pca.data$x[i,1]
        circle.center.y <- pca.data$x[i,2]

        radius <- step.size

        ## now collect all of the bootstrapped samples that correspond to this
        ## cell...
        cell.points <- boot.points.by.cell[boot.points.by.cell$cell == i,c(1,2)]

        ## now figure out how big this circle needs to be to contain 95% of the
        ## bootstrapped points.
        done <- FALSE
        while(!done) {
          contained.points <- sum(((cell.points[,1] - circle.center.x)^2 +
                                     (cell.points[,2] - circle.center.y)^2) <=
                                    (radius^2))

          if (contained.points >= min.points) {
            done <- TRUE
          } else {
            radius <- radius + step.size
          }
        }
        radii[i] <- radius
      }

      if (trim.proportion > 0) {
        cat("Trimming samples with the most variation\n")
        cat("\tThe top ", round(trim.proportion * 100, digits=2), "% will be removed\n", sep="")
        radius.cutoff <- quantile(radii, (1-trim.proportion))
        cutoff.indices <- which(radii >= radius.cutoff)
        cat("\t", length(cutoff.indices), " samples were removed\n", sep="")
        ## just in case we need to return the genes that were not trimmed
        return.samples.vector <- c(1:num.samples)
        return.samples.vector <- return.samples.vector[-cutoff.indices]

        ## now we need to delete the original samples and the
        ## bootstrapped versions of them...
        ##pca$x <- pca$x[-cutoff.indices,]
        pca.data$x <- pca.data$x[-cutoff.indices,]

        rownames(boot.points) <- c(1:nrow(boot.points))

        for (i in 0:(num.boot.samples-1)) {
          offset <- i * num.samples
          boot.points[cutoff.indices + offset,] <- NA
        }
        na.indices <- which(is.na(boot.points[,1]))
        boot.points <- boot.points[-na.indices,]
      }
    }
  }

  if (!confidence.regions) {
    radii <- NULL
  }

  if (nrow(pca.data$x) > 0) {
    draw.pcaBootPlot(pca=pca.data, boot.points=boot.points,
                     pca.var.per=pca.var.per,
                     num.boot.samples=num.boot.samples,
                     plot.col=plot.col,
                     plot.pch=plot.pch,
                     data.factors=data.factors,
                     draw.legend=draw.legend,
                     legend.names=legend.names,
                     legend.x=legend.x,
                     legend.y=legend.y,
                     hex.colors=hex.colors,
                     transparency=transparency,
                     min.x=min.x, max.x=max.x, min.y=min.y, max.y=max.y,
                     radii=radii)

    if (!is.null(pdf.filename)) {
      pdf(file=pdf.filename, width=pdf.width, height=pdf.height)

      draw.pcaBootPlot(pca=pca.data, boot.points=boot.points,
                       pca.var.per=pca.var.per,
                       num.boot.samples=num.boot.samples,
                       plot.col=plot.col,
                       plot.pch=plot.pch,
                       data.factors=data.factors,
                       draw.legend=draw.legend,
                       legend.names=legend.names,
                       legend.x=legend.x,
                       legend.y=legend.y,
                       hex.colors=hex.colors,
                       transparency=transparency,
                       min.x=min.x, max.x=max.x, min.y=min.y, max.y=max.y,
                       radii=radii)

      dev.off()
    }
  } else {
    cat("\nThere were no points to plot! Bummer\n")
    return()
  }

  cat("\nDone! Hooray!")
  #return("Done! Hooray!")
  if (return.samples) {
    return(return.samples.vector)
  }
}

draw.pcaBootPlot <- function(pca=NULL, boot.points=NULL, pca.var.per=NULL,
                             num.boot.samples=100,
                             plot.col=NULL,
                             plot.pch=NULL,
                             data.factors=NULL,
                             draw.legend=FALSE, legend.names=NULL,
                             legend.x=NULL, legend.y=NULL,
                             hex.colors=NULL,
                             transparency=77,
                             min.x=NULL, max.x=NULL, min.y=NULL, max.y=NULL,
                             radii=NULL) {

  #print("about to draw the plot!")
  #print(head(pca$x))
  if (num.boot.samples > 0 && (nrow(boot.points) > 0)) {
    if (is.null(max.x)) {
      max.x <- max(boot.points[,1], pca$x[,1], na.rm=TRUE)
    }
    if (is.null(min.x)) {
      min.x <- min(boot.points[,1], pca$x[,1], na.rm=TRUE)
    }
    if (is.null(max.y)) {
      max.y <- max(boot.points[,2], pca$x[,2], na.rm=TRUE)
    }
    if (is.null(min.y)) {
      min.y <- min(boot.points[,2], pca$x[,2], na.rm=TRUE)
    }
    plot(boot.points, type="n", xlim=c(min.x, max.x), ylim=c(min.y, max.y), xlab=paste("PC1 (", pca.var.per[1], "%)", sep=""),
         ylab=paste("PC2 (", pca.var.per[2], "%)", sep=""))
    grid()
    points(boot.points, pch=plot.pch, col=plot.col)
  } else {
    if (is.null(legend.x)) {
      legend.x <- 0
    }
    if (is.null(legend.y)) {
      legend.y <- 0
    }
    if ((!is.null(min.x)) & (!is.null(max.x))) {
      x.lims <- c(min.x, max.x)
    }
    if ((!is.null(min.y)) & (!is.null(max.y))) {
      y.lims <- c(min.y, max.y)
    }
    if (exists("x.lims") & exists("y.lims")) {
      plot(pca$x[,c(1,2)], type="n",
           xlab=paste("PC1 (", pca.var.per[1], "%)", sep=""),
           ylab=paste("PC2 (", pca.var.per[2], "%)", sep=""),
           xlim=x.lims, ylim=y.lims)
    } else {
      plot(pca$x[,c(1,2)], type="n",
           xlab=paste("PC1 (", pca.var.per[1], "%)", sep=""),
           ylab=paste("PC2 (", pca.var.per[2], "%)", sep=""))
    }
    grid()
  }

  if(!is.null(radii)) {
    symbols(pca$x[,c(1,2)], circles=radii, add=TRUE, inches=FALSE, lty=2, col="#bdbdbd77")
  }

  points(pca$x[,c(1,2)], pch=plot.pch)

  if (draw.legend) {
    if (is.null(legend.names)) {
      legend.names <-  c(1:length(data.factors[,1]))
    }
    legend(legend.x, legend.y, legend.names, pch=1:length(data.factors[,1]), pt.cex=1, cex=0.8, col=hex.colors[1:nrow(data.factors)])
  }
}
