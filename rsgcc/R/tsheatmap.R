

#################################################################################################################
##This function plots HeatMap with different similarity measures (e.g., 1-CorCoef) for tissue/condition specific genes
##Here method could be GCC (Gini correlation coefficient), 
##PCC (Pearson product-moment correlation coefficient), 
##SCC (Spearman's rank correlation coefficient), 
##KCC (Kendall tau correlation coefficient) 
##BiWt(correlation estimates based on Tukey's biweight M-estimator)
##MI (mutual information)
##MINE
##ED
#################################################################################################################

 
gcc.tsheatmap <- function(x,
                           
                        cpus = 1,
                        
                        ## correlation method
                        method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"),
                          
                        distancemethod = c("Raw", "Abs", "Sqr"),
                           
                        #hclustfun = hclust,
                        clustermethod = c("complete", "average", "median", "centroid", "mcquitty", "single", "ward"),
                         
                        #hcdata by output gcc.tsheatmap 
                        rowhcdata = NULL,
                        colhcdata = NULL,    
                          
                          
                        keynote = "FPKM",
                        
                        ## dendrogram control
                        symm = FALSE,
                      
                        ## data scaling
                        scale = c("none","row", "column"),
                        na.rm=TRUE,

                        ## image plot
                        revC = identical(Colv, "Rowv"),
                        add.expr,

                        ## mapping data to colors
                        breaks,
                        symbreaks=min(x < 0, na.rm=TRUE) || scale!="none",

                       ## colors
                       colrange = c("yellow", "red"), 
                          
                       tissuecol= "heat.colors",

                       ## block sepration
                       colsep = 0.15,
                       rowsep,
                       sepcolor="white",
                       sepwidth=c(0.05,0.05),
          
                       ## level trace
                       trace=c("none","column","row","both"),
                       tracecol="cyan",
                       hline=median(breaks),
                       vline=median(breaks),
                       linecol=tracecol,

                       ## Row/Column Labeling
                       margins = c(5, 5),

                       ## plot labels
                       main = NULL,
                       xlab = NULL,
                       ylab = NULL,

                       ## plot layout
                       lmat = NULL,
                       lhei = NULL,
                       lwid = NULL,

                       ## extras
                       ...
                      )
{  
  
  cat("Dimension information for clustered GE matrix:", dim(x), "\n")
  
  #get axis position for tissue bar
  getAtpos <- function(genenumVec) {
    
     posVec <- rep(0, length(genenumVec))
     sum <- 0
     for( i in 1:length(genenumVec)) {
       posVec[i] <- sum + genenumVec[i]/2
       sum <- sum + genenumVec[i]
     }
     return( posVec/sum(genenumVec) )
  }
      
  #color key
  scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
  ###################################################################################
  dendrogram <- "none"  #no dendrogram, as only tissue specific genes are considered
  Rowv <- "none"  #row will be ordered by script showning below
  Colv <- "none"  #col will be ordered by script showning below
  retval <- list()
  scale <- if (symm & missing(scale)) {"none"} else { match.arg(scale) }
  trace <- match.arg(trace)
#   if (length(col) == 1 && is.character(col)) 
#       col <- get(col, mode = "function")

  if (!missing(breaks) & (scale != "none")) {
      warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  }
  if (is.null(Rowv) | is.na(Rowv)) {   Rowv <- FALSE }
  
  if (is.null(Colv) | is.na(Colv)) { Colv <- FALSE }
  else if (Colv == "Rowv" & !isTRUE(Rowv)) {   Colv <- FALSE   }
    
  di <- dim(x)
  if(length(di) != 2 ) { stop("x is not matrix?") }
  if(!is.numeric(x)) {  stop("x must be a numeric matrix") }
  nr <- di[1]
  nc <- di[2]
   
  if (nr <= 1 | nc <= 1) { stop("`x' must have at least 2 rows and 2 columns") }
  if (!is.numeric(margins) | length(margins) != 2) { stop("`margins' must be a numeric vector of length 2") }
  
  ###################################################################################
  #get cluster information and , ordered geneID
  if( is.null(rowhcdata) ) { #no hcdata, ok, caclulate it
    hcr <- gcc.hclust( x, cpus = cpus, method = method, distancemethod = distancemethod, clustermethod = clustermethod)  #clustered for rows
  }else {
    hcr <- rowhcdata
  }
  ddr <- as.dendrogram(hcr$hc)
  ddr <- reorder(ddr, TRUE)
  rowInd <- order.dendrogram(ddr)
  
  ##for hcc data
  if( is.null(colhcdata ) ) { #no colhcdata, ok, cacluate it
    hcc <- gcc.hclust( t(x), cpus = cpus, method = method, distancemethod = distancemethod, clustermethod = clustermethod)  #clustered for columns
  }else {
    hcc <- colhcdata
  }
  
  ##get ordered sampleID
  tsMatrix <- uniqueTissues(x)
  newx <- x[rowInd,]
  geneTSVec <- rep(0, nrow = nrow(newx))
  for( i in 1:nrow(newx)) {
    meanVec <- rep(0, nrow(tsMatrix))
    for( j in 1:nrow(tsMatrix)) {
      sampleIndex <- which(tsMatrix[j,] > 0)
      meanVec[j] <- mean(newx[i,sampleIndex])
    }    
    geneTSVec[i] <- which(meanVec == max(meanVec))[1]
  }

  orderMatrix <- matrix(0, nrow = nrow(tsMatrix), ncol = 4 )
  colnames(orderMatrix) <- c("TissueIndex", "startGeneIndex", "EndGeneIndex", "GeneNum" )
  rownames(orderMatrix) <- rownames(tsMatrix)
  for( i in 1:nrow(tsMatrix)) {
    GeneIndex <- which(geneTSVec == i)
    if(length(GeneIndex) > 0) { 
        if( length(GeneIndex) != (max(GeneIndex) - min(GeneIndex) + 1) ) {  #gene index not continuous, find gene index with max continuous length
          conlenmatrix <- matrix(0, nrow = length(GeneIndex), ncol = length(GeneIndex))
          for( x in 1:length(GeneIndex)-1) {
            for( y in (x+1):length(GeneIndex) ) {
              GeneIndexSub <- GeneIndex[x:y]
              if(length(GeneIndexSub) == max(GeneIndexSub) - min(GeneIndexSub) + 1) {
                conlenmatrix[x,y] <- length(GeneIndexSub)
              }
          }#end for y
        }#end for x
        
        maxconlen <- max(conlenmatrix)
        if( maxconlen > 0) {
          maxrow <- apply(conlenmatrix, 1, max)
          rowIndex <- which(maxrow == maxconlen)
          colIndex <- which(conlenmatrix[rowIndex[1],] == maxconlen)[1]
          GeneIndex <- seq(GeneIndex[rowIndex], GeneIndex[colIndex], by=1)
        }else {
          GeneIndex <- c(GeneIndex[2])
        }
      }
    }else { #no ts-gene for this tissue
      GeneIndex = 0
    }
    
    orderMatrix[i,1] <- i
    orderMatrix[i,2] <- min(GeneIndex)
    orderMatrix[i,3] <- max(GeneIndex)
    orderMatrix[i,4] <- length(GeneIndex)
    if( length(GeneIndex) == 1 & GeneIndex[1] == 0 ) 
      orderMatrix[i,4] <- 0
  } 

  colInd <- NULL
  orderMatrix <- orderMatrix[sort(orderMatrix[,2], decreasing=FALSE, index.return = TRUE)$ix,]
  print(orderMatrix)
  for( i in 1:nrow(orderMatrix) ) {
      curTissueIndex <- orderMatrix[i,1]
      sampleIndex <- which(tsMatrix[curTissueIndex,] > 0)
      colInd <- c(colInd, sampleIndex)   
  }#end for i 
  newx <- newx[, colInd]
  
  #get col for unique tissues
  TissueCol <- NULL
  if (missing(tissuecol) ) {
    TissueCol <- rainbow(nrow(orderMatrix))
  }else if( is.function(tissuecol) ) {
    TissueCol <- tissuecol(nrow(orderMatrix))
  }else {
    if( is.null(names(tissuecol)) ) {
      stop("Error: tissuecol is a vector parameters whose names are unique tissue names")
    }else {
      if( length(tissuecol) != nrow(tsMatrix) | 
        length( which( (sort(rownames(tsMatrix)) == sort(names(tissuecol))) == FALSE)) > 0 ) {
        cat("Error: the input and detected unique tissue name is not the same, pls change the input ones to detected ones.\n")
        cat("input tissue names:\n")
        print( names(tissuecol))
        cat("detected tissue namers:\n")
        print( rownames(tsMatrix))
        stop("========")
      }
      TissueCol <- rep("NA", nrow(tsMatrix))
      for( i in 1:length(TissueCol) ) {
        TissueCol[i] <- tissuecol[ which( tissuecol == rownames(TissueCol)[i])]
      }
    }
  }
  
  
  RowSideCol <- rep("NA", length(geneTSVec))
  for( i in 1:length(RowSideCol) ) {
    RowSideCol[i] <- TissueCol[geneTSVec[i]]
  }
  
  #################################################################################################
  call <- match.call()  #add 20130514
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- newx
  x.unscaled <- x  
   if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
  
  #for breaks
  if( missing(breaks) ) {
    breaks <- quantile( unique(c(x)), probs = seq(0, 1, 0.01), na.rm = TRUE)
  }else if( is.null(breaks) | length(breaks) < 1) {
         breaks <- 16
  }
  
  if (length(breaks) == 1) {
     if (!symbreaks) {
          breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), length = breaks)
     } else {
          extreme <- max(abs(x), na.rm = TRUE)
          breaks <- seq(-extreme, extreme, length = breaks)
        }
  }
  
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(colrange) == "function") {
     col <- colrange(ncol)
  }else if( is.null(colrange) ) {
    if( !require(grDevices) ) install.packages("grDevices")
    require(grDevices)
    if( !require(grDevices) )
    col <- colorRampPalette(c("yellow", "red"))(nbr - 1) #for RNA-Seq
  }else if( is.vector(colrange) & length(colrange) == 2) {
    if( !require(grDevices) ) install.packages("grDevices")
    require(grDevices)
    cat("color range is:", colrange, "\n" )
    col <- colorRampPalette(colrange)(nbr - 1) #for RNA-Seq
  }else {
    stop("Error: colrange should be a character vector containing two colors")
  }
  
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  

  #########################################################################
  ##below, plot
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
 
  if( is.null( lmat) ) {
    lmat <- t(matrix(rep(c(4, 1, 6, 2, 7, 0, 5), 3), ncol = 3))
    lmat[2,6] <- 3
  }
 
  if( is.null( lwid)) {
    lwid <- c(0.5, 0.1, 0.01, 3.5, 0.01, 0.1, 0.5)
  }
  lhei <- rep(5, 3)
  
  if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
 
  

  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  nr = nrow(newx)
  nc = ncol(newx)
  
  #tissue col bar
   par(mar = c(margins[1], 0, 2, 0.1))
   image(rbind(1:nr), col = RowSideCol, xaxt = "n", yaxt = "n")
   xv <- getAtpos(orderMatrix[which(orderMatrix[,4] > 0),4])
   axis(2, at=xv, labels=rownames(orderMatrix[which(orderMatrix[,4] > 0),]), las= HORIZONTAL<-1, cex.axis=1.2, adj = 1, xpd = TRUE )

  #heat map
  image(1:nc, 1:nr, t(newx), 
        xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), 
        axes = FALSE, xlab = "", ylab = "", 
        col = col, breaks = breaks, ...)



  min.raw <- min(newx)
  max.raw <- max(newx)
  z <- seq(min.raw, max.raw, length = length(col))
  image( t(matrix(z, ncol = 1)), col = col, xaxt = "n", yaxt = "n")  #remove breaks from orginal function
  lv <- pretty(breaks)
  xv <- scale01(as.numeric(lv), min.raw, max.raw)
  axis(4, at = xv, labels = lv) 
  mtext(side = 4, keynote, line = 3)
  
  retval$breaks <- breaks
  retval$col <- col
  retval$rownames <- rownames(newx)
  retval$colnames <- colnames(newx)
  
  return(list(retval = retval, hcr = hcr, hcc = hcc) )
  
  
}
 
