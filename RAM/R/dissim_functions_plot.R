#########################
# dissim.plot functions #
#########################

### A slightly ambitious task would be to abstract almost all the dissim.X.plot methods
### into one; the only (minor) difficulty is adjusting the plot settings for the different
### outputs

# small internal function to store the default distance methods
.get.dist.methods <- function() {
  return(c("morisita", "bray", "jaccard", "chao", "euclidean"))
}

# small internal function to store the default clustering methods
.get.clust.methods <- function() {
  return(c("average", "centroid", "complete", "mcquitty"))
}

.get.ord.methods <- function() {
  return(c("euclidean"))
}

.get.stand.methods <- function() {
  return(c("total", "max", "frequency", "normalize", "range", "standardize", "pa", "chi.square", "hellinger", "log"))
}

dissim.clust.plot <- function(data, is.OTU=TRUE, stand.method=NULL,
                          dist.methods=NULL,clust.methods=NULL,      
                          file=NULL) {
  
  .valid.data(data, is.OTU=is.OTU)
  num.data <- length(data)
  labels <- names(data)

  given.dist.methods <- !is.null(dist.methods)
  given.clust.methods <- !is.null(clust.methods)
  
  # if the user does not supply any methods, set the default ones
  if ( !given.dist.methods ) {
    dist.methods <- .get.dist.methods()
  }
  
  if ( !given.clust.methods ) {
    clust.methods <- .get.clust.methods()
  }
  
  ncols <- length(dist.methods) * length(clust.methods)
  
  # if too many methods supplied, warn the user
  if (ncols >  10) {
    warning("you specified many distance and cluster methods; some may not fit on the plot. If none are displayed, try again with fewer.")
  }
  
  save <- !is.null(file)
  if (save) { 
    # the default size is too small; so adjust size manually
    .get.dev(file, "pdf", height = 5 * num.data, width = 5 * ncols)
  }
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # determine the number of rows needed for plotting 
  # (one row for otu1, one for otu2)
  par(mfcol=c(num.data, ncols))
  
  # plot each distance/clustering method, for each dataset
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]

    for (dist.met in dist.methods) {
      for (clust.met in clust.methods) {
        title <- paste(dist.met, "/", clust.met, sep="")
      
        # calculate the cluster data (also does data validation)
        dist.clust <- dissim.clust(elem=elem, is.OTU=is.OTU, 
                            stand.method=stand.method, 
                           dist.method=dist.met, clust.met)
        # plot cluster data; hang=-1 aligns the labels
        plot(dist.clust, hang=-1, main=title)
      }
    }
  }
  
  if (save) { dev.off() }
  
  invisible()
}

dissim.eig.plot <- function(data, is.OTU=TRUE, stand.method=NULL,
                            dist.methods=NULL, file=NULL) {
  
  .valid.data(data, is.OTU=is.OTU)
  num.data <- length(data)
  labels <- names(data)

  save <- !is.null(file)  
  given.dist.methods <- !is.null(dist.methods)  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  # if the user does not supply any methods, set the default ones
  if ( !given.dist.methods ) {
    dist.methods <- .get.dist.methods()
  }
  
  if (save) { .get.dev(file, "pdf") }
  
  par(mfcol=c(num.data, length(dist.methods)))  
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]
    for (met in dist.methods) {
      # calculate the eigenvalues (also does data validation)
      eigenvals <- dissim.eig(elem=elem, is.OTU=is.OTU, 
                             stand.method=stand.method,
                             dist.method=met)
      barplot(eigenvals, ylab="Eigenvalue", xlab="Axis Number", 
              main=paste(met, label, sep=": "))
    }
  }  
  if (save) { dev.off() }  
  invisible()
}


dissim.ord.plot <- function(data, is.OTU=TRUE, stand.method=NULL,
                         dist.methods=NULL, k=NULL, file=NULL) {
  .valid.data(data, is.OTU=is.OTU)
  num.data <- length(data)
  labels <- names(data)
  save <- !is.null(file)
  
  samples <- list()
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]

    if ( is.OTU ) {
      samples[[label]] <- ncol(elem) -1
    } else {
      samples[[label]] <- nrow(elem)
    }
  } 
  # number of dimension
  k.max <- min(unlist(samples)-1)
  if (is.null(k)) {
    k <- k.max
  }
  if (!is.numeric(k) || length(k) != 1L) {
    stop("k must be a numeric vector of length 1")
  }
  
  given.dist.methods <- !is.null(dist.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
   
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  ncols <- length(dist.methods)
  
  # if too many methods supplied, warn the user
  if (ncols >  10) {
    warning("you specified many distance and ordination methods; some may not fit on the plot. If none are displayed, try again with fewer.")
  }
  
  if (save) { .get.dev(file, "tiff") }
  
  par(mfcol=c(num.data, ncols))
  
  # plot the data, recall that output of dist.ord is a list where the first 
  # item is the ordination distances, the second is the given method distances
  samples <- list()
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]    
    for (dist.met in dist.methods) {
      title <- paste(label,": dist-", dist.met, sep="")
      
      # calculate the distances (also does data validation)
      distances <- dissim.ord(elem=elem, is.OTU=is.OTU,
                              stand.method=stand.method,
                              dist.method=dist.met, k=k)
      
      plot(distances[[1]], distances[[2]], 
           xlab="Ordination Distance", ylab="Method Distance", 
           main=title)
    }
  }
 
  if (save) { dev.off() }
  invisible()
}

dissim.GOF.plot <- function(data, is.OTU=TRUE, stand.method=NULL, 
                            dist.methods=NULL, file=NULL) {
  .valid.data(data, is.OTU=is.OTU)
  num.data <- length(data)
  labels <- names(data)
  save <- !is.null(file)
  
  given.dist.methods <- !is.null(dist.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (save) { .get.dev(file, "pdf") }
  
  par(mfcol=c(num.data, length(dist.methods)))
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]

    for (met in dist.methods) {
      # calculate the goodness of fit values 
      # (also does data validation)
      GOF.vals <- dissim.GOF(elem=elem, is.OTU=is.OTU, 
                             stand.method=stand.method,
                             dist.method=met)
      plot(GOF.vals, xlab="Dimensions", 
           ylab="Goodness of Fit", 
           main=paste(met, label, sep=": "))

    }
  }
  if (save) { dev.off() }  
  invisible()
}

dissim.tree.plot <- function(data, is.OTU=TRUE, 
                             stand.method=NULL,
                             dist.methods=NULL, 
                             clust.methods=NULL, file=NULL) {
  .valid.data(data, is.OTU=is.OTU)
  num.data <- length(data)
  labels <- names(data)

  save <- !is.null(file)
  
  given.dist.methods <- !is.null(dist.methods)
  given.clust.methods <- !is.null(clust.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (!given.clust.methods) {
    clust.methods <- .get.clust.methods()
  }
  
  ncols <- length(dist.methods) * length(clust.methods)
  
  # if too many methods supplied, warn the user
  if (ncols >  10) {
    warning("you specified many distance and cluster methods; some may not fit on the plot. If none are displayed, try again with fewer.")
  }
  
  if (save) { .get.dev(file, "pdf", width=3*ncols) }
  
  par(mfcol=c(num.data, ncols))
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]

    for (dist.met in dist.methods) {
      for (clust.met in clust.methods) {
        title = paste(label, ": ", clust.met, "/", dist.met, sep="")
      
        # calculate the distances (also does data validation)
        distances <- dissim.tree(elem=elem, is.OTU=is.OTU, 
                                 stand.method=stand.method, 
                                 dist.method=dist.met,
                                 clust.method=clust.met)
        # plot the data, recall that output of dist.tree is a 
        # list where the first 
        # item is the tree distances, the second is the given 
        # method distances
        plot(distances[[1]], distances[[2]], 
              xlab="Tree Distances", 
             ylab="Method Distances", main=title)
      }
    }
  }  
  if (save) { dev.off() } 
  invisible()
}


dissim.pvar.plot <- function(data, is.OTU=TRUE, stand.method=NULL, 
                             dist.methods=NULL, file=NULL) {

  .valid.data(data, is.OTU=is.OTU)
  num.data <- length(data)
  labels <- names(data)

  save <- !is.null(file)
  
  given.dist.methods <- !is.null(dist.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (save) { .get.dev(file, "pdf") }
  
  par(mfcol=c(num.data, length(dist.methods)))
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]

    for (met in dist.methods) {
      # calculate the percent variation (also does data validation)
      pct.var <- dissim.pvar(elem=elem, is.OTU=is.OTU, 
                             stand.method=stand.method,
                             dist.method=met)
      barplot(pct.var, xlab="Axis Number", ylab="% of Variation", 
              main=paste(met, label, sep=": "))
    }
  }
  if (save) { dev.off() }
  invisible()
}

dissim.alleig.plot <- function(data, is.OTU=TRUE, 
                        stand.method=NULL, 
                        dist.methods=NULL, file=NULL) {

  .valid.data(data, is.OTU=is.OTU)
  num.data <- length(data)
  labels <- names(data)

  save <- !is.null(file)
  given.dist.methods <- !is.null(dist.methods)
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  # parameter for pco calculations; 
  # can be at most n-1 where n is the # of samples
  samples <- list()
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]
    
    if ( is.OTU ) {
      samples[[label]] <- ncol(elem) -1
    } else {
      samples[[label]] <- nrow(elem)
    }
  } 

  # number of dimension
  k.max <- min(unlist(samples)-1)  

  df.rows <- list()
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]
    
    # set region for rows
    ### can this be improved?

    for (met in dist.methods) {
      # get the eigenvalules
      eigenvals <- dissim.eig(elem=elem, is.OTU=is.OTU,
                              stand.method=stand.method,
                              dist.method=met)
      # we want to save the fraction of the sum for each eigenvalue
      lab <- paste(label, met, sep="_")
      df.rows[[lab]] <- c(met, label, eigenvals / sum(eigenvals))
    }
  }
 
    # make our data frame; set stringsAsFactors must be false, 
    # otherwise our numeric data becomes factors 
    # (this way the numeric data becomes character
    # data which we can easily convert back)
    eigs <- as.data.frame(do.call(rbind, df.rows), 
                         stringsAsFactors=FALSE)

    # the number of columns with eigenvalue data
    eig.cols <- ncol(eigs)-2
  
    # convert our eigenvalues measures back to numeric data
    for (i in 3:ncol(eigs)) {
      eigs[ ,i] <- as.numeric(eigs[ ,i])
    }
  
    # make the column names "Method", "Region", "1", ..., "n"
    colnames(eigs) <- c("Method", "Data", paste(1:eig.cols))

    #if ( !require("reshape2") ) {
    #  stop("package 'reshape2' required for this function")
   # }  
    #if ( !require("scales") ) {
   #   stop("package 'scales' required for this function")
   # }  
    eigs.m <- melt(eigs, id.vars=c("Method", "Data"),
                 variable.name="eigenval", value.name="pct_of_sum")
  
    names(eigs.m)[ncol(eigs.m)-1] <- "eigenval"
    names(eigs.m)[ncol(eigs.m)] <- "pct_of_sum"

    eigs.m$pct_of_sum <- as.numeric(eigs.m$pct_of_sum)
  
    # we need to use aes_string to pass CRAN check; see 
    # http://goo.gl/JxgZ9u
    p <- ggplot(eigs.m, aes_string(x="eigenval", 
                 y="pct_of_sum", group="Method", 
                               colour="Method")) +
       scale_y_continuous(labels = percent_format()) + 
       ylab("Fraction of Sum") +
       xlab("Axis Number") +
       facet_wrap(~Data) +
       geom_line(size=1) +
       geom_point(size=2)
  
  if (save) {
    file <- .ensure.filepath(file, "pdf")
    ggsave(filename=file)
  }  
  p
}

###############################################################
# These functions all deal with calculating measures related 
# to dissimilarity matrices
###############################################################
dissim.clust <- function(elem, is.OTU=TRUE, stand.method=NULL,
                         dist.method="morisita", 
                         clust.method="average") {
  
  dist <- .dissim.dist(elem=elem, is.OTU=is.OTU, 
                       stand.method=stand.method, 
                       dist.method=dist.method)
  dist.clust <- hclust(dist, method=clust.method)
  
  return(dist.clust)
}

dissim.eig <- function(elem, is.OTU=TRUE, stand.method=NULL, 
                       dist.method="morisita") {
  if ( is.OTU ) {
     valid.OTU(elem)
     data <- transpose.OTU(elem)
  } else {
     data <- elem
  }
  k.max <- nrow(data) - 1
  dist <- .dissim.dist(elem=elem, is.OTU=is.OTU, 
                       stand.method=stand.method, 
                       dist.method=dist.method)
  # note: we suppress warnings because we will often have 
  # negative eigenvalues
  # due to numeric error (which raises a warning in pco) 
  dist.pco <- suppressWarnings(pco(dist, k=k.max))  
  return(dist.pco$eig)
}

dissim.ord <- function(elem, is.OTU=TRUE, stand.method=NULL, 
                       dist.method="morisita", k=NULL) {
  
  if ( is.OTU ) {
     valid.OTU(elem)
     data <- transpose.OTU(elem)
  } else {
     data <- elem
  }

  dist <- .dissim.dist(elem=elem, is.OTU=is.OTU, 
                       stand.method=stand.method, 
                       dist.method=dist.method)
  k.max <- nrow(data) - 1 
  if ( is.null(k) ) {
    k <- k.max
  }  
  if (!is.numeric(k) || length(k) != 1L) {
    stop("k must be a numeric vector of length 1")
  }
  
  if (k <= k.max) {
      # calculate the ordination distances
      ord <- pco(dist, k=k)
      ord.dist <- vegdist(ord$points, method="euclidean")
      met.dist <- dist
      return(list(ord.dist, met.dist))
  } else { 
    stop(paste("specified value for k was", k, 
               "which is larger than the maximum value", k.max)) 
  } 
}

dissim.GOF <- function(elem, is.OTU=TRUE, stand.method=NULL, 
                       dist.method="morisita") {
  if ( is.OTU ) {
     valid.OTU(elem)
     data <- transpose.OTU(elem)
  } else {
     data <- elem
  }
  k.max <- nrow(data) - 1
  
  met.dist <- .dissim.dist(elem=elem, is.OTU=is.OTU, 
                       stand.method=stand.method, 
                       dist.method=dist.method)

  GOF.vals <- numeric(k.max - 1)
  
  # calculate goodness of fit values for all possible values of k 
  # note: we suppress warnings because we will often have negative eigenvalues
  # due to numeric error (which raises a warning in cmdscale) 
  for (k in 2:k.max) {
    coord <- suppressWarnings(cmdscale(met.dist, k, eig=TRUE))
    GOF.vals[k] <- coord$GOF[2]
  }
  
  return(GOF.vals)
}

dissim.tree <- function(elem, is.OTU=TRUE, stand.method=NULL, 
                dist.method="morisita", clust.method="average") {
  
  met.dist <- .dissim.dist(elem=elem, is.OTU=is.OTU, 
                       stand.method=stand.method, 
                       dist.method=dist.method)
 
  # calculate the clustering and tree distances
  clust <- hclust(met.dist, method=clust.method)
  tree.dist <- cophenetic(clust)
  
  return(list(tree.dist, met.dist))
}

dissim.pvar <- function(elem, is.OTU=TRUE, stand.method=NULL, 
                        dist.method="morisita") {
  if ( is.OTU ) {
     valid.OTU(elem)
     data <- transpose.OTU(elem)
  } else {
     data <- elem
  }
  
  met.dist <- .dissim.dist(elem=elem, is.OTU=is.OTU, 
                       stand.method=stand.method, 
                       dist.method=dist.method)
  k.max <- nrow(data) - 1
  
  # note: we suppress warnings because we will often have negative eigenvalues
  # due to numeric error (which raises a warning in pco) 
  met.pcoa <- suppressWarnings(pco(met.dist, k=k.max))
  
  # get the number of positive eigenvalues
  axes <- dim(met.pcoa$points)[2]  
  return(met.pcoa$eig / sum(met.pcoa$eig))
}

.dissim.dist <- function(elem, is.OTU, stand.method, 
                         dist.method) {
  if ( is.OTU ) {
     valid.OTU(elem)
     data <- transpose.OTU(elem)
  } else {
     data <- elem
  }

  if ( is.null(stand.method) ) {
    data.new <- data
  } else {
    if ( any(grepl(stand.method, .get.stand.methods())) &&
         (length(stand.method) == 1L) ) {
      data.new <- decostand(data, stand.method)
    } else {
       warning("not a valid standardization method, see ?dissim.clust and ?decostand for details")
    }
  }   
  # calculate the distances, then create the hclust
  dist <- vegdist(data.new, method=dist.method)
}
