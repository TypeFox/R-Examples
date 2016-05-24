# These functions all deal with creating plots of given data.

######################
# plot.top functions #
######################

### NOTE: since the change in tax.abund 
# (commit 2a420ad061ff9478fc7765d10bd63c2717f614fb), the group.top.[number|percent]
# functions have gotten very messy, and would benefit from a cleanup.

# the current control flow is as follows:
# user calls group.top.[number|percent] -> .top.samples.plot -> 
# top.samples.[ggplot2|base]

# internal function (hidden from user; they call plot.top.percent or plot.top.number)
.top.samples.plot  <-  function(data, top=10, ranks=c("p","c","o","f","g"),
                              drop.unclassified, cex.x=NULL,
                              main=NULL,file=NULL, ext=NULL, 
                              height=8, width=16, 
                              bw=FALSE, ggplot2=TRUE, mode) {

 .valid.data(data)

  save <- !is.null(file)
  .valid.plot.settings(file, ext)
  
  labels <- names(data)
  num.otus <- length(data)

  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    otu <- data[[i]]
    if ( is.null(otu) ) { break }
    valid.OTU(otu)
   
    # set the function to select the top samples
    # mode, top, and drop.unclassified inherit from this function 
    # (.top.samples.plot) call, data and rank are determined from when top.func 
    # is called later
    top.func  <-  function(otu, rank){
                       tax.abund(otu, top=top, rank=rank, 
                                 count=FALSE, mode=mode,
                            drop.unclassified=drop.unclassified)}
  }
  
  # set the appropriate titles
  if ( is.null(main) ) {
    if (mode == "number") {
    
      top.title  <-  paste("Relative Abundance of Top", top, 
                       "Taxon Groups at Five Taxonomic Ranks", sep=" ")
    
    } else if (mode == "percent") {
      top.title  <-  paste("Relative Abundance of Taxon Groups Above ", top, 
                       "% Relative Abundance at Five Taxonomic Ranks", sep="")
    }
  } else {
    top.title <- main
  }
  
  # call the appropriate plotting function
  if (ggplot2) {
    .top.samples.ggplot2(data=data, top=top, mode=mode, ranks=ranks, 
                        cex.x=cex.x, file=file, ext=ext, 
                        height=height, width=width, bw=bw, 
                        top.title=top.title, top.func=top.func)
  } else {
    .top.samples.base(data=data, top=top, mode=mode, ranks=ranks,
               cex.x=cex.x, file=file, ext=ext, height=height, 
                      width=width, bw=bw, 
                      top.title=top.title, top.func=top.func)
  }
}

.top.samples.ggplot2  <-  function(data, top, mode, 
                       ranks=c("p","c","o","f","g"), cex.x=NULL,
                       file, ext, height, width, bw, 
                       top.title, top.func) {
  
  save  <-  !is.null(file)
  .valid.plot.settings(file, ext)
  #ranks  <-  c("phylum", "class", "order", "family", "genus", "species")
  ranks <- ranks
  
  .valid.data(data)

  labels <- names(data)
  num.otus <- length(data)

  rows <- list()

  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    OTU <- data[[i]]
    if ( is.null(OTU) ) { break }
    valid.OTU(OTU)
 
    rank.list <- list()
    # get the number of samples
    num.samples  <-  ncol(OTU) - 1

    tax.names <- vector()
    
    for (j in ranks) {
      # get the data, add region/rank information
      .valid.rank(j)
      rank <- .get.rank.name(j)
      tax.names <- c(tax.names, .capitalize(rank))
      top.otus <- top.func(OTU, rank=rank)
      top.otus <- cbind(top.otus, Region=rep(label, times=num.samples),
                        Rank=rep(.capitalize(rank), times=num.samples))
      
      #top.otus[["Region"]] <- label
      #top.otus[["Rank"]] <- rank
          
      # add to our list after melting
      rank.list[[rank]]  <-  reshape2::melt(top.otus, id.vars=c("Region", "Rank"), 
                            value.name="RA", variable.name="OTU") 
  
    }
   rows[[label]]  <-  do.call(rbind, rank.list)
  }
 
  data  <-  do.call(rbind, rows)
  
  data$Rank <- factor(data$Rank,levels=tax.names, ordered=TRUE)
 
  # get a colour palette with colours for each OTU
  cols.needed  <-  length(unique(data$OTU))
  if (bw) {
    colours  <-  rep("white", times=cols.needed)
  } else {
    if (cols.needed > 12) {
      colours  <-  grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(cols.needed)
    } else {
      colours  <-  RColorBrewer::brewer.pal(12, "Set3")
    }
  }
  
  p  <-  ggplot(data, aes_string(x="OTU", y="RA", fill="OTU")) + 
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=colours) +
    scale_y_continuous(labels = percent_format()) +
    xlab("Taxon") + ylab("Relative Abundance") + 
    ggtitle(top.title) +
    theme(legend.position="none",
          panel.background = element_rect(fill="grey90"),
          panel.grid.major.x = element_line(colour="white"),
          panel.grid.major.y = element_line(colour="white")) +
    facet_grid(Region ~ Rank, scales="free_x") 
  
  if ( is.null(cex.x) ) {
    p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  } else {
    p <- p + theme(axis.text.x = element_text(size=cex.x, angle = 45, vjust = 1, hjust=1))
  }

  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  }
  
  p
}

.top.samples.base  <-  function(data, top, mode, 
                     ranks=c("p","c","o","f","g"), cex.x=NULL,
                     file, ext, height, width, bw, top.title, top.func) {
  
  save  <-  !is.null(file)
  .valid.plot.settings(file, ext)
  if (save) { .get.dev(file, ext, height=height, width=width) }

  .valid.data(data)

  labels <- names(data)
  num.otus <- length(data)

  # backup plot settings and restore afterwards
  opar  <-  par(no.readonly = TRUE)
  on.exit(par(opar))
  
  # configure plot settings
  par(oma=c(4, 0, 2, 0), mar=c(8, 2, 1, 1), las=2)
  
  # setup our layout 
  layout(matrix(1:(5 * num.otus), ncol=5, byrow=TRUE))
  

  rows <- list()

  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    valid.OTU(elem)

    elem.ranks  <-  list()
    
    # call get.rank for each taxonomic classification
    tax.names <- vector()
    for (j in ranks) {
      pretty.rank <- .get.rank(.get.rank.ind(j), pretty=TRUE)
      tax.names <- c(tax.names, pretty.rank)
      elem.ranks[[pretty.rank]] <- top.func(elem, rank=j)
    }
    
    # get color palette
    if (bw) {
      cols  <-  grey.colors(top)
    } else {
      # the max pallette size is 12; if length > 12 we get the palette with 12
      # colors (which will be recycled by boxplot, which is fine)
      cols  <-  suppressWarnings(brewer.pal(top, "Set3"))
    }
    
    # plot each rank
    for (x in 1:length(tax.names)) {
      # plot the data
      if ( is.null(cex.x) ) {
        boxplot(elem.ranks[[x]], col=cols, notch=FALSE, cex.axis=0.8)
      } else {
        boxplot(elem.ranks[[x]], col=cols, notch=FALSE, cex.axis=cex.x)
      } 
      # add a legend
      region  <-  label
      legend("topright", legend=paste0(region, " - ", tax.names[x]))
    }
  }
  
  # print title
  title(main=top.title, outer=TRUE)
  
  if (save) { dev.off() }
  
  invisible()
}

# group.top.percent and group.top.number are what the user call; the other top.samples.X
# functions are used interally to produce the graph.
group.top.percent  <-  function(data, top=10, ranks=c("p","c","o","f","g"), 
                              drop.unclassified=FALSE,  cex.x=NULL,
                              main=NULL,file=NULL, ext=NULL, height=8, width=16, 
                              bw=FALSE, ggplot2=TRUE) {

  .valid.data(data)
  labels <- names(data)
  num.otus <- length(data)
    
  .top.samples.plot(data=data, top=top, ranks=ranks,
                    drop.unclassified=drop.unclassified, main=main, cex.x=cex.x,
                    file=file, ext=ext, height=height, width=width,
                    bw=bw, ggplot2=ggplot2, mode="percent")
}

group.top.number  <-  function(data, top=10, ranks=c("p","c","o","f","g"),
                             drop.unclassified=FALSE, cex.x=NULL,
                             main=NULL, file=NULL, ext=NULL, height=8, 
                             width=16, bw=FALSE, ggplot2=TRUE) {
  .valid.data(data)
  labels <- names(data)
  num.otus <- length(data)

  .top.samples.plot(data, top=top, ranks=ranks, cex.x=cex.x,
                    drop.unclassified=drop.unclassified, main=main,
                    file=file, ext=ext, height=height, width=width, 
                    bw=bw, ggplot2=ggplot2, mode="number")
}

group.abundance  <-  function(data, rank, 
                            top=NULL, count=FALSE, drop.unclassified=FALSE,
                            cex.x=NULL, main=NULL, file=NULL, ext=NULL, 
                            height=8, width=16, bw=FALSE, ggplot2=TRUE) {
  
  .valid.rank(rank)
  .valid.data(data)
  
  save  <-  !is.null(file)
  .valid.plot.settings(file, ext)
  if (save) { .get.dev(file, ext, height=height, width=width) }

  
  labels <- names(data)
  num.otus <- length(data)
  
  if (ggplot2) {
    .abundance.ggplot2(data, rank, top, count, drop.unclassified, cex.x, main, file, ext,
                       height, width, bw)
  } else {
    .abundance.base2(data, rank, top, count, drop.unclassified, cex.x, main, file, ext,
                    height, width, bw)
  }
}


.abundance.base2 <- function(data, rank, 
                            top=NULL, count=FALSE, drop.unclassified=FALSE,
                            cex.x=NULL, main=NULL, file=NULL, ext=NULL, labels, 
                            height=8, width=9, bw=FALSE) {
  .valid.data(data)
  .valid.rank(rank)
  .valid.plot.settings(file, ext)
  
  save <- !is.null(file)
  
  #single.otu <- is.null(otu2)
  
  num.otus <- length(data)
  
  
  # set up the appropriate device if saving
  if (save) {
    .get.dev(file, ext, height, width)
  }
  
  # backup plot settings and restore afterwards
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  # set up the layout: one shared pane up top (for the title),
  # one pane for each graph and one for its legend (if applicable)
  #mat <- matrix(c(rep(1, times=num.otus), 2:((2 * num.otus) + 1)), ncol=num.otus,
  #              byrow=TRUE)
  
  mat <- matrix(c(1,1,2,4,3,5), ncol=2, byrow=TRUE)
  
  if (num.otus == 1) {
    lmat <- matrix(1:3, ncol=1)
  } else if ( num.otus == 2) {
    lmat <- matrix(c(1:3, c(1, 4, 5)), ncol=2)
  } else {
    stop("base plotting in function group.abundance only supports a maximum of TWO OTUs")
  } 
  
  layout(lmat, heights=c(1, 8, 2))
  
  pretty.rank <- .get.rank(.get.rank.ind(rank), pretty=TRUE)
  
  if ( is.null(main) ) {
    if (count) {
      top.title <- paste("Counts of Taxonomic Groups at", pretty.rank, "Level") 
    } else {
      top.title <- paste("Relative Abundance of Taxonomic Groups at", pretty.rank,
                       "Level")
    }
  } else {
    top.title <- main
  }
  
  par(mar=c(0, 0, 4, 0))
  plot.new()
  title(main=top.title, cex.main=1.5)
  
  #index <- 1
  
  for (i in 1:length(data) ) {
    elem <- data[[i]]
    valid.OTU(elem)
    # stop if not given otu2
    if (is.null(elem)) { break }
    label <- names(data)[i]
      
    abund <- tax.abund(elem, rank=rank, drop.unclassified=drop.unclassified,
                         top=top, count=count)
    
    sample.totals <- c(rowSums(abund), 0)
    breaks <- pretty(sample.totals)
    
    num.taxa <- dim(abund)[2]
    
    # set the appropriate label
    if (count) {
      y.label <- "Count"
    } else {
      y.label <- "Relative Abundance"
    }
    
    if (bw) {
      cols <- grey.colors(n=num.taxa, end=0.7)
    } else {
      # generate palette of colours, warn the user if too many OTUs
      cols <- suppressWarnings(RColorBrewer::brewer.pal(num.taxa, "Set3"))
      
      if (num.taxa > 12) {
        warning("this colour palette only has 12 distinct colours, and more than 12 taxon groups have been provided. Some colours are being recycled in your graph (be careful!)")
      }
    }
    
    # the barplot function is expecting the transpose
    abund <- t(abund)
    
    # we add the axes manually with our breaks variable later
    barplot.args <- list(abund, beside=FALSE, names.arg=colnames(abund),
                         cex.names=0.7, col=cols, las=2, axes=FALSE,
                         ylim=c(0, max(breaks)))
   
    if ( is.null(cex.x) ) {
      barplot.args$cex.names=0.7
    } else {
      barplot.args$cex.names=cex.x
    }
    
    leg.args <- list(x="center", legend=rownames(abund), cex=0.7, horiz=FALSE, 
                     ncol=3, fill=cols)
    
    # add shading to the plot
    if (bw) {
      dens <- rep(c(18, 9), length.out=num.taxa)
      ang <- c(30, 60, 90, 120, 150)
      
      barplot.args <- c(barplot.args, density=list(dens), angle=list(ang))
      leg.args <- c(leg.args, density=list(3 * dens), angle=list(ang))
    }
    
    par(mar = c(5.1, 5.1, 4.1, 2.1))
    do.call(barplot, barplot.args)
    
    title(ylab=y.label, line=3.75)
    title(xlab="Samples", line=4.2)
    title(main=label)
    
    
    axis(side=2, at=breaks, las=2)
    
    
    # set margins for legend
    par(mar=c(0, 0, 0, 0))
    plot.new()
    do.call(legend, leg.args)
  }
  
  if (save) {
    dev.off()
  }
  
  invisible()
}


.abundance.ggplot2  <-  function(data, rank, 
                               top=NULL, count=FALSE, drop.unclassified=FALSE,
                               cex.x=NULL, main=NULL, file=NULL, ext=NULL,  
                               height=8, width=16, bw=FALSE) {
  # validate inputs
  save  <-  !is.null(file)
  if (save) {
    .valid.plot.settings(file, ext)
  }  

  .valid.rank(rank)

  if ( class(data) != "list" ) {
    stop("please provide data as list, see ?RAM.input.formatting")
  }

  labels <- names(data)
  num.otus <- length(data)

  taxa  <-  list()
  for (i in 1:length(data) ) {
    elem <- data[[i]]
    if (is.null(elem)) { break }
    label <- names(data)[i]

    # get the groups
    elem.tax  <-  tax.abund(elem, rank=rank, drop.unclassified=drop.unclassified,
                          count=count, top=top)
    # melt by Sample
    elem.tax  <-  reshape2::melt(cbind(elem.tax, Sample=rownames(elem.tax)), id.vars="Sample",
                     value.name="RA", variable.name="Taxon")
    # bind together with appropriate label
    elem.tax  <-  cbind(elem.tax, Region=label)
    
    taxa[[label]]  <-  elem.tax
  }
  
  all.taxa  <-  do.call("rbind", taxa)
  
  # reorder levels based on order of appearance in table
  all.taxa$Sample  <-  ordered(all.taxa$Sample, levels=unique(all.taxa$Sample))
  
  if ( is.null(main) ) {
    if (count) {
      title  <-  paste("Counts of Taxonomic Groups at", 
                   .get.rank(.get.rank.ind(rank), pretty=TRUE),
                   "Level")
    } else {
      title  <-  paste("Relative Abundance of Taxonomic Groups at", 
                   .get.rank(.get.rank.ind(rank), pretty=TRUE),
                   "Level")
    }
  } else {
    title <- main
  }
  
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  p  <-  ggplot(all.taxa, aes_string(x="Sample", y="RA", fill="Taxon")) + 
    geom_bar(position="stack", stat="identity") +
    #coord_flip() +
    xlab("Samples") +
    ggtitle(title) 

  if ( is.null(cex.x) ) {
     p <- p + theme(legend.position="bottom", 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major.x = element_blank())
  } else {
     p <- p + theme(legend.position="bottom", 
      axis.text.x = element_text(size=cex.x, angle = 45, vjust = 1, hjust=1),
          panel.grid.major.x = element_blank())
  }
 
  if (!count) {
    p  <-  p +  scale_y_continuous(labels = percent_format()) +
      ylab("Relative Abundance")
  } else {
    p  <-  p + ylab("Count")
  }
  
  p  <-  p + facet_wrap(~Region, scales="free_x")
 
  
  if (bw) { # for black/white plots
    warning("the ggplot2 package used to create this graph cannot handle patterning for bar plots; greyscale shading is being used.")
    p  <-  p + scale_fill_grey(guide = guide_legend(direction="horizontal", ncol=5)) +
      theme_bw() + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
  } else { # use actual colours
    
    cols.needed  <-  length(unique(all.taxa$Taxon))
    
    if (cols.needed > 12) {
      # palette max is 12; so if we have more than 6 entries for otu1/2, we 
      # need to manually construct a palette to use
      
      col.func  <-  grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      p  <-  p + scale_fill_manual(values=col.func(cols.needed), 
                                 guide=guide_legend(direction="horizontal", ncol=5))
    } else {
      # otherwise, we can just use the Set3 palette from RColorBrewer
      p  <-  p + scale_fill_brewer(palette="Set3",
                                 guide=guide_legend(direction="horizontal", ncol=5))
    }
  }
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  }
  
  p
}


factor.abundance <- function(data, rank, top=NULL, count=FALSE, 
                             meta=meta, meta.factor="", drop.unclassified=FALSE, 
                             file=NULL, ext=NULL, height=8, width=16, main="") {
  # validate inputs
 # valid.OTU(otu1, otu2)
 # .valid.plot.settings(file, ext)

  .valid.rank(rank)
  .valid.data(data)
  
  save  <-  !is.null(file)
  .valid.plot.settings(file, ext)
  if (save) { .get.dev(file, ext, height=height, width=width) }

  labels <- names(data)
  num.otus <- length(data)
  
  taxa <- list()
  for (i in 1:length(data)) {
    elem <- data[[i]]
    label <- names(data)[i]
    if (is.null(elem)) {
      break
    }
    valid.OTU(elem)
    #return(elem)
    # get the groups
    if (length(meta.factor) > 1L || length(meta.factor)==0L || meta.factor == "" ) {
      stop("Please provide ONE category variable in the metadata")
    } else {
      elem.tax<-.group.rank(otu=elem, meta=meta, meta.factor=meta.factor, relative.abund=TRUE, top=top, drop.unclassified=drop.unclassified, rank=rank)
    }
  #return(elem.tax)
    # melt by Sample
  if (!requireNamespace("reshape2")) {
    stop("package 'reshape2' is required to use this function")
  }
  
    elem.tax.m <- reshape2::melt(cbind(elem.tax, Sample=rownames(elem.tax)), id.vars="Sample", value.name="RA", variable.name="Taxon")
    # bind together with appropriate label
    elem.tax.m <- cbind(elem.tax.m, Region=label)
    
    taxa[[label]] <- elem.tax.m
  }
  
  all.taxa  <-  do.call("rbind", taxa)
  
  # reorder levels based on order of appearance in table
  all.taxa$Sample <- ordered(all.taxa$Sample, levels=unique(all.taxa$Sample))
  
    title <- main
  
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  p <- ggplot(all.taxa, aes_string(x="Sample", y="RA", fill="Taxon")) + 
    geom_bar(position="stack", stat="identity") +
    #coord_flip() +
    xlab(meta.factor) +
    ggtitle(title) +
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major.x = element_blank())
    if (!requireNamespace("scales") ||   !requireNamespace("RColorBrewer")) {
       stop("packages scales and RColorBrewer are required for this function")
    }
    p <- p +  scale_y_continuous(labels = percent_format()) +
      ylab("Relative Abundance")
  
  #if (!single.otu) { # for multiple OTUs, wrap based on region
  #  p <- p + facet_wrap(~Region, scales="free_x")
  #}
  
   # use actual colours
    
    cols.needed <- length(unique(all.taxa$Taxon))
    
    if (cols.needed > 12) {
      # palette max is 12; so if we have more than 6 entries for otu1/2, we 
      # need to manually construct a palette to use
      
      col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      p <- p + scale_fill_manual(values=col.func(cols.needed), 
                                 guide=guide_legend(direction="horizontal", ncol=5))
    } else {
      # otherwise, we can just use the Set3 palette from RColorBrewer
      p <- p + scale_fill_brewer(palette="Set3",
                                 guide=guide_legend(direction="horizontal", ncol=5))
    }
    p <- p + facet_wrap(~Region, scales="free_x")
    
  print(p)
  if (save) { dev.off() }
  
        invisible()
  
}


group.temporal <- function(data, meta, date.col, factors, rank, group, 
                            file=NULL, ext=NULL, height=8, width=12) {
  
  valid.OTU(data)
  .valid.rank(rank)
  .valid.meta(otu1=data, meta=meta)
  save <- !is.null(file)
  
  # check all sample names (ignore the last column; that's taxonomy)
  if (!all(colnames(data)[-dim(data)[2]] %in% rownames(meta))) {
    stop("sample names do not match for OTU table and metadata.")
  }
  
  if (!date.col %in% colnames(meta)) {
    stop("'date.col' was not found in 'meta'.")
  }
  
  if (class(meta[ ,date.col]) != "Date") {
    warning("the date column in metadata must be of type Date; coercing it now. See ?RAM.dates.")
    meta[ ,date.col] <- as.Date(meta[ ,date.col])
  }
  
  meta.factors <- .valid.factors(meta, factors)
  
  # regular apply will coerce the data frame to a matrix, so all columns become 
  # character vectors (infuriating!); so instead lapply and unlist after (since
  # data frames are lists in disguise)
  not.numeric <- unlist(lapply(meta.factors,
                       FUN=function(column) {
                         !is.numeric(column)}))
  
  if (any(not.numeric)) {
    stop((paste("this function only supports numeric variables, the following variables are not numeric:\n",
                paste(names(meta.factors)[not.numeric], collapse=" "))))
  }
  
  dates <- meta[ ,date.col]
  
  meta.factors <- aggregate(meta.factors, by=list(dates), FUN=mean)
  meta.factors <- rename(meta.factors, c("Group.1"="Date"))
  meta.factors <- reshape2::melt(meta.factors, id.vars="Date", variable.name="Measure",
                       value.name="Value")
  names(meta.factors)[ncol(meta.factors)] <- "Value"
  names(meta.factors)[ncol(meta.factors)-1] <- "Measure"
  
  xlims <- c(min(meta.factors$Date), max(meta.factors$Date))
  
  meta.plot <- ggplot(meta.factors, aes_string(x="Date", y="Value")) + geom_line() +
               facet_wrap(~Measure, ncol=1, scales="free_y") + 
               ylab("Average Value") + 
               theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     plot.margin=unit(c(1, 1, 0, 0.5), "lines")) +
               xlim(xlims)
     
  
  groups <- .get.tax.group(data, rep(rank, times=length(group)), group)
  abund <- tax.abund(groups, rank=rank)
  
  ordering <- match(rownames(abund), rownames(meta))
  meta.sort <- meta[ordering, ]

  # tax.abund reorders the samples, so be careful when aggregating
  abund.agg <- aggregate(abund, by=list(dates[ordering]), FUN=sum)
  abund.agg <- rename(abund.agg, c("Group.1"="Date"))
  
  abund.agg <- reshape2::melt(abund.agg, id.vars="Date", variable.name="Group", 
                     value.name="Count")
  names(abund.agg)[ncol(abund.agg)] <- "Count"
  names(abund.agg)[ncol(abund.agg)-1] <- "Group"
  
  main.plot <- ggplot(abund.agg, aes_string(x="Date", y="Count", fill="Group")) + 
               geom_bar(stat="identity", position="stack") +
               theme(legend.position="bottom",
                     plot.margin=unit(c(0, 1, 1, 0.5), "lines")) +
               xlim(xlims) + 
               scale_fill_brewer(palette = "Set3")
  
  # much credit due to baptiste (on SO)
  # taken from http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
  gA <- ggplot_gtable(ggplot_build(meta.plot))
  gB <- ggplot_gtable(ggplot_build(main.plot))

  maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
  gA$widths[2:3] <- as.list(maxWidth)
  gB$widths[2:3] <- as.list(maxWidth)
  
  if (save) {
    file <- .ensure.filepath(file, ext)
    .get.dev(file, ext, height=height, width=width)
  }
  
  grid.arrange(gA, gB, ncol=1, heights=c(0.2, 0.8))
  
  if (save) { dev.off() }
  invisible()
}

group.spatial <- function(data, meta, date.col, province.col, rank, group, 
                          breaks="year", file=NULL, ext=NULL, height=8, width=10) {
  valid.OTU(data)
  .valid.rank(rank)
  .valid.meta(otu1=data, meta=meta)
  save <- !is.null(file)
  
  if (!date.col %in% colnames(meta)) {
    stop("'date.col' was not found in 'meta'.")
  }
  
  if (!province.col %in% colnames(meta)) {
    stop("province.col' was not found in 'meta'.")
  }
  
  if (!is.character(meta[ ,province.col])) {
    warning("the 'Province' column of meta is not a character vector; coercing it to one.")
    meta$Province <- as.character(meta[ ,province.col])
  }
  
  meta$Province <- .get.province.id(meta[ ,province.col])
  
  group.len <- length(group)
  groups <- .get.tax.group(data, rep(rank, times=group.len), group)
  abund <- tax.abund(groups, rank=rank)
  
  ordering <- match(rownames(abund), rownames(meta))
  meta.sort <- meta[ordering, ]
  
  meta.sort[ ,date.col] <- as.Date(meta.sort[ ,date.col])
  # note that setting ordered_result = TRUE in cut does NOT order the 
  # levels by actual date value, so we do that manually
  date.factor <- cut(meta.sort[ ,date.col], breaks=breaks)
  date.factor <- ordered(date.factor, sort(levels(date.factor)))
  
  abund.loc <- aggregate(abund, by=list(meta.sort$Province, date.factor), 
                         FUN=sum)
  
  names(abund.loc)[1:2] <- c("id", "Date")
  
  # this is only here to silence a CRAN check note; otherwise it complains 
  # that we reference map.fortify without defining it (even though loading
  # the file below will create an object map.fortify)
  map.fortify.list <- NULL
  # load file map.fortify, which contains a Canadian province map already 
  # fortified (for plotting with ggplot2)
  load(system.file("extdata", "map.fortify.list.RData", package="RAM"))
  
  # we get all missing provinces (which have count 0), and add them to 
  # abund.loc so they show up in our map later (for each time segment)
  missing <- vector(length=length(levels(abund.loc$Date)), mode="list")
  
  dates <- levels(abund.loc$Date)
  for (i in 1:length(dates)) {
    # get all provinces without counts
    missing.provinces <- setdiff(unique(map.fortify.list$map$id),
                                 unique(abund.loc[abund.loc$Date == dates[i], ]$id))
    
    # add a NA count for each group
    group.count <- vector(length=group.len, mode="list")
    group.count[1:group.len] <- NA
    
    # create and rename data frame
    missing[[i]] <- data.frame(missing.provinces, dates[i], group.count)
    # all names after the first two are groups
    names(missing[[i]]) <- names(abund.loc)
  }
  
  abund.loc <- rbind(abund.loc, do.call(rbind, missing))
  abund.loc.melt <- reshape2::melt(abund.loc, id.vars=c("id", "Date"), variable.name="Group",
                         value.name="Count")
  names(abund.loc.melt)[ncol(abund.loc.melt)] <- "Count"
  names(abund.loc.melt)[ncol(abund.loc.melt)-1] <- "Group"
  
  loc.map <- base::merge(map.fortify.list$map, abund.loc.melt, by="id")
  
  # set the fill to the given taxon group
  tax.group <- names(abund.loc)[3]

  p <- ggplot(loc.map, aes_string(x="long", y="lat", map_id="id", fill="Count")) +
       scale_fill_gradient(na.value = "white") +
       geom_map(map=loc.map, colour="black") +
       facet_grid(Group ~ Date) +
       theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
             axis.text.y = element_blank(), axis.title.y = element_blank(),
             axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
             panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
}

# this function simply maps the province codes to the names used in the GeoBase map
.get.province.id <- function(province) {
  
  valid.abbrev <- c("AB", "BC", "MB", "NB", "NL", "NT", "NS", "NU", "ON", "PE",
                    "QC", "SK", "YT")
  
  provinces <- c("ALBERTA", "BRITISH COLUMBIA", "MANITOBA", "NEW BRUNSWICK", 
                 "NEWFOUNDLAND AND LABRADOR", "NORTHWEST TERRITORIES", "NOVA SCOTIA",
                 "NUNAVUT", "ONTARIO", "PRINCE EDWARD ISLAND", "QUEBEC", 
                 "SASKATCHEWAN", "YUKON" )
  
  if (!is.character(province)) {
    warning("'province' should be a character vector; coercing it to one now.")
    province <- as.character(province)
  }
  
  if (!all(province %in% valid.abbrev)) {
    stop("invalid province code. See ?location.formatting for details.")
  }

  # get the name at the correct index
  provinces[match(province, valid.abbrev)]
}

group.indicators  <-  function(data, is.OTU=TRUE, meta, factor, rank,
                            thresholds = c(A=0.85, B=0.8, stat=0.8, p.value=0.05),
                            cex.x=NULL, file=NULL, ext=NULL,
                            height=12, width=12) {
  
  # I have seen this discussion: http://yihui.name/en/2014/07/library-vs-require/
  # but I think returning an explanatory error message is worthwhile
  
  if ( !requireNamespace("indicspecies") ) {
  #  indicspecies::multipatt
 # } else {
    stop("package 'indicspecies' is required to use this function: try 'install.packages('indicspecies')'.")
  }
  
  .valid.data(data, is.OTU=is.OTU)

  save  <-  !is.null(file)
  meta.factor  <-  .valid.factors(meta, factor, min.factors = 1, max.factors = 1)
  meta.name  <-  names(meta.factor)

  labels <- names(data)
  num.otus <- length(data)
  # to store the data we're about to generate
  rows <- list()
  sel <- list()
  indicators.df.kept <- list()

  for ( ot in 1:length(data) ) {
    elem <- data[[ot]]
    if ( is.null(elem) ) { break }
    label <- names(data)[ot]

    if ( is.OTU ) {
      otu <- elem
      valid.OTU(otu)
      if ( !is.null(rank) ) {
        .valid.rank(rank)
        .valid.meta(otu, meta=meta)
        abund  <-  tax.abund(otu, rank=rank)
      } else if (is.null(rank) ) {
        rank <- NULL
        if ( length(data) > 1 ) {
          stop("To identify otus as indicators, only provide ONE data set")
        }
        abund <- data.revamp(data=list(df=otu), is.OTU=is.OTU, 
                          stand.method=NULL, 
                          ranks=rank, top=NULL)[[1]]
      }
    } else {
      abund <- elem      
    }

    abund <- abund[match(rownames(meta), rownames(abund)),]
    if ( !identical(rownames(meta), rownames(abund)) ) {
      stop("data and metadata do not have same subjects")
    }

    if (!is.numeric(thresholds) || length(thresholds) != 4L) {
      stop("thresholds must be a numeric vector of length four (see ?indicators.plot for details).")
    }
   
    abund.stand  <-  decostand(abund, method="total")
    
    mp  <-  indicspecies::multipatt(abund.stand, meta.factor[[1]], control=how(nperm=999))
    
    #return(mp)
    indicators <- capture.output(indicspecies::summary.multipatt(mp, indvalcomp=TRUE))

    # the summary contains some human-friendly non-data lines we need to strip
    # this regex kepts only lines that have {number number number number},
    # since there are four numerical columns in the data (and none elsewhere)
    # this selects only the numerical data we want
    matches  <-  grepl("([[:digit:]\\.]+[[:space:]]+){4}", indicators)
    
    if (!any(matches)) {
      stop("no taxon groups were significant with the given parameters; try again with a lower taxon group and/or a different meta.factor.")
    }
    
    indicators  <-  indicators[matches]
    
    # split the rows at all whitespace, convert to dataframe
    indicators.df  <-  data.frame(Group=character(), A=numeric(), B=numeric(),
                                stat=numeric(), p.value=numeric(),
                                stringsAsFactors=FALSE)
    
    pieces  <-  strsplit(indicators, "[[:space:]]+")
    
    for (i in 1:length(pieces)) {
      # the [-6] removes the trailing significance code (which is just a human-friendly
      # display of the p-values)
      indicators.df[i, ]  <-  pieces[[i]][-6]
    }
    
    for (i in 2:5) {
      indicators.df[ , i]  <-  as.numeric(indicators.df[ , i])
    }
   
    # get names of all taxon groups above the given thresholds    
    df<-indicators.df
    df$Select <- NA
    for ( i in 1:nrow(df)) {
      if ( df[i, "A"] >= thresholds[1] & df[i, "B"] >= thresholds[2] & df[i, "stat"] >= thresholds[3] & df[i, "p.value"] <= thresholds[4] ) {
      df[i, "Select"] <- "Y"
      } else {
         df[i, "Select"] <- ""
     }
   }

    indicators.df.kept[[label]]<-df
    kept.all <- df$Group[df$Select=="Y"]

    if (length(kept.all) == 0L) {
      stop("no taxon groups met all thresholds. Either relax your thresholds, or try again with each otu individually.")
    } else {
      kept <- kept.all
    }

    # set up our data frame with all the data ggplot needs
    abund.filtered  <-  data.frame(Sample=rownames(abund.stand),
                                 meta=meta.factor,
                                 Region=rep(label, times=dim(abund.stand)[1]), 
                                 abund.stand[ , kept, drop=FALSE])

    # melt it & store the result
    rows[[label]]  <-  reshape2::melt(abund.filtered, id.vars=c("Sample", meta.name, "Region"),
                                                    variable.name="Indicator",
                                                    value.name="Value")
    sel[[label]] <- kept
    
  } # end otu for loop
  
  data  <-  do.call(rbind, rows)

  # we use as.formula/paste to create the faceting formula (otherwise we would 
  # have to type Region ~ meta.name, and facet_grid won't evaluate meta.name)
  p  <-  ggplot(data, aes_string(x="Sample", y="Value", fill="Indicator")) +
       geom_bar(stat="identity") +
       facet_grid(as.formula(paste("Region", "~", meta.name)),
                  scales="free_x", space="free_x") +
       ylab("Relative Abundance")
  if ( is.null(cex.x) ) {
     p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
            legend.position="bottom")
  } else {
     p <- p + theme(axis.text.x = element_text(size=cex.x, angle = 45, 
                    vjust = 1, hjust=1), 
               legend.position="bottom")
  }
 
  cols.needed  <-  length(kept)

  # color palette
  col.pal <- RAM.pal(cols.needed)
    
  if ( cols.needed > 12 ) {
    p  <-  p + scale_fill_manual(values=col.pal, 
                                 guide=guide_legend(direction="horizontal", ncol=5)) 
  } else {
    # otherwise, we can just use the Set3 palette from RColorBrewer
    p  <-  p + scale_fill_brewer(palette="Set3",
                                 guide=guide_legend(direction="horizontal", ncol=5))
  }

  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    print(p)
  }

  return(list(indicators.df.kept, sel))
}

sample.locations <- function(otu1, otu2=NULL, meta, factor=NULL, zoom=5,
                             source="google", labels=c("ITS1", "ITS2"),
                             lat.col="Latitude", long.col="Longitude",
                             file=NULL, ext=NULL, height=10, width=12) {
  
  valid.OTU(otu1, otu2)
  .valid.meta(otu1, otu2, meta=meta)
  
  if (is.null(otu2)) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  .valid.labels(num.otus, labels)
  
  save <- !is.null(file)
  
  if (!source %in% c("google", "osm")) {
    stop("'source' must be one of 'google' or 'osm'.")
  }
  # .valid.factors will check this later, but we do it here to give descriptive
  # error message
  if (!lat.col %in% names(meta)) {
    stop("'lat.col' was not found in meta (see ?sample.locations for help).")
  }
  
  if (!long.col %in% names(meta)) {
    stop("'long.col' was not found in meta (see ?sample.locations for help).")
  }
  
  if (source == "google") {
    if (zoom < 3 || zoom > 21) {
      stop("when source == 'google', zoom must be between 3 and 21 inclusive.")
    }
  } else if (source == "osm") {
    if (zoom < 3 || zoom > 18) {
      stop("when source == 'osm', zoom must be between 3 and 18 inclusive.")
    }
  }
  
  columns <- c(Latitude=lat.col, Longitude=long.col)
  
  # to be used in the ggmap call later
  points.aes <- aes_string(x="Longitude", y="Latitude", size="Counts")
  
  if (!is.null(factor)) {
    
    if (is.null(attr(factor, which="names"))) {
      stop("'factor' must be a named character vector.")
    }
    
    columns <- c(columns, factor)
    
    points.aes <- c(points.aes, aes_string(colour=names(factor)))
    class(points.aes) <- "uneval" # the call to c above strips the class info
  }
  
  if (num.otus == 2) {
    points.aes <- c(points.aes, aes_string(shape="Region"))
    class(points.aes) <- "uneval" # the call to c above strips the class info
  }
  
  meta.data <- .valid.factors(meta, columns, 
                              min.factors=length(columns), max.factors=3)
  
  index <- 1
  for (elem in list(otu1, otu2)) {
    if (is.null(elem)) { break }
    
    meta.data <- cbind(meta.data, colSums(elem[ ,-dim(elem)[2]]))
    # the first two columns are lat/long data, then possibly a metadata factor
    names(meta.data)[index + 2 + !is.null(factor)] <- labels[index]
                       
    index <- index + 1
  }
  
  meta.data <- reshape2::melt(meta.data, id.vars=c(lat.col, long.col, names(factor)), 
                    variable.name="Region", value.name="Counts")
  names(meta.data)[ncol(meta.data)] <- "Counts"
  names(meta.data)[ncol(meta.data)-1] <- "Region"
  
  meta.data.agg <- aggregate(Counts ~ ., data=meta.data, FUN=sum)
  
  buffer <- 0.5
  
  left <- min(meta.data$Longitude) - buffer
  bottom <- min(meta.data$Latitude) - buffer
  right <- max(meta.data$Longitude) + buffer
  top <- max(meta.data$Latitude) + buffer
  
  # create bounding box for map, with a little extra room
  bounding <- c(left=left, bottom=bottom, right=right, top=top)
  
  get_map.args <- list(location=bounding, zoom=zoom, maptype="roadmap", 
                       color="bw", source=source)
  
  if (source == "osm") {
    get_map.args$scale <- OSM_scale_lookup(zoom)
  }
  
  map <- do.call(get_map, get_map.args)
  
  #p <- ggmap(map, extent="device") +
  p <- ggmap(map) +
       geom_point(points.aes, data=meta.data.agg, alpha=0.7) +
       scale_color_brewer(palette = "Set1") +
       scale_size(range=c(2,10))
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p 
  }
}

group.abundance.meta  <-  function(data, rank, top=NULL, count=FALSE, 
                            drop.unclassified=FALSE, cex.x=NULL, main=NULL, 
                            file=NULL, ext=NULL, height=8, width=16, bw=FALSE, 
                            meta=NULL, meta.factor=NULL) {
  
  .valid.rank(rank)
  .valid.data(data)
  
  save  <-  !is.null(file)
  .valid.plot.settings(file, ext)
  if (save) { .get.dev(file, ext, height=height, width=width) }

  
  labels <- names(data)
  num.otus <- length(data)
  
  .abundance.ggplot2.meta(data, rank, top, count, drop.unclassified, cex.x, main, file, ext,
                       height, width, bw, meta=meta, meta.factor=meta.factor)
}



.abundance.ggplot2.meta  <-  function(data, rank, 
                               top=NULL, count=FALSE, drop.unclassified=FALSE,
                               cex.x=NULL, main=NULL, file=NULL, ext=NULL,  
                               height=8, width=16, bw=FALSE, meta=NULL, meta.factor=NULL) {
  # validate inputs
  save  <-  !is.null(file)
  if (save) {
    .valid.plot.settings(file, ext)
  }  

  .valid.rank(rank)

  if ( class(data) != "list" ) {
    stop("please provide data as list, see ?RAM.input.formatting")
  }

  labels <- names(data)
  num.otus <- length(data)

  taxa  <-  list()
  for (i in 1:length(data) ) {
    elem <- data[[i]]
    if (is.null(elem)) { break }
    label <- names(data)[i]
    valid.OTU(elem)

    if ( !is.null(meta) ) {
       .valid.meta(otu1=elem, meta=meta)
    }
    
    if ( is.null(meta) && !is.null(meta.factor) ) {
      stop("Metadata was not provided")
    } else if ( !is.null(meta) && is.null(meta.factor)) {
      meta <- meta
    } else if ( !is.null(meta) && ! is.null(meta.factor)) {
      fac.len <- length(meta.factor)
      vec.fac <- vector()
      for ( i in 1:fac.len) {
        if ( ! any(meta.factor[i] %in% names(meta)) ) {
           vec.fac <- c(vec.fac, meta.factor[i])
        }
      }
      if (length(vec.fac) !=0L ) {
        stop (paste(paste(vec.fac, collapse=", "), " are not in metadata", sep=""))
      } 
      meta <- meta[, meta.factor, drop=FALSE]  
    } else {
      meta <- meta
    }
      
    # get the groups
    elem.tax  <-  tax.abund(elem, rank=rank, drop.unclassified=drop.unclassified,
                          count=count, top=top)
    #return(elem.tax)
    if ( !is.null(meta) && ! is.null(meta.factor)) {
      if ( identical(rownames(elem.tax), rownames(meta)) ) {
        elem.tax1 <- cbind(elem.tax, meta[, meta.factor])
        names(elem.tax1)[(ncol(elem.tax)+1):ncol(elem.tax1)] <- meta.factor
      } else {
        stop("otu table and metadata do not have same sample or sampleIDs are in different order")
      }
    } else {
        elem.tax1 <- elem.tax
    }
    #return(elem.tax1)
    # melt by Sample

    elem.tax2  <-  reshape2::melt(cbind(elem.tax1, Sample=rownames(elem.tax)), id=meta.factor, id.vars=c("Sample", meta.factor),
                     value.name="RA", variable.name="Taxon")
    
    # bind together with appropriate label
    elem.tax2  <-  cbind(elem.tax2, Region=label)
    #return(elem.tax2)
 
    taxa[[label]]  <-  elem.tax2
  }
  
  all.taxa  <-  do.call("rbind", taxa)
  
  # reorder levels based on order of appearance in table
  all.taxa$Sample  <-  ordered(all.taxa$Sample, levels=unique(all.taxa$Sample))
  
  if ( is.null(main) ) {
    if (count) {
      title  <-  paste("Counts of Taxonomic Groups at", 
                   .get.rank(.get.rank.ind(rank), pretty=TRUE),
                   "Level")
    } else {
      title  <-  paste("Relative Abundance of Taxonomic Groups at", 
                   .get.rank(.get.rank.ind(rank), pretty=TRUE),
                   "Level")
    }
  } else {
    title <- main
  }
  
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  p  <-  ggplot(all.taxa, aes_string(x="Sample", y="RA", fill="Taxon")) + 
    geom_bar(position="stack", stat="identity") +
    #coord_flip() +
    xlab("Samples") +
    ggtitle(title) 

  if ( is.null(cex.x) ) {
     p <- p + theme(legend.position="bottom", 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major.x = element_blank())
  } else {
     p <- p + theme(legend.position="bottom", 
      axis.text.x = element_text(size=cex.x, angle = 45, vjust = 1, hjust=1),
          panel.grid.major.x = element_blank())
  }
 
  if (!count) {
    p  <-  p +  scale_y_continuous(labels = percent_format()) +
      ylab("Relative Abundance")
  } else {
    p  <-  p + ylab("Count")
  }
  

  if ( !is.null(meta) && ! is.null(meta.factor)) {
    if (length(meta.factor) ==1L) {
       formula <- paste("Region", " ~ ", meta.factor, sep="") 
        p  <-  p + facet_grid(as.formula(formula), scales="free_x")
    } else {
       formula <- paste("Region", " ~ ", paste(meta.factor, collapse=" + "), sep="") 
        p  <-  p + facet_grid(as.formula(formula), scales="free_x")
    }
  } else {    
    p  <-  p + facet_wrap(~Region, scales="free_x")
  }
  
  if (bw) { # for black/white plots
    warning("the ggplot2 package used to create this graph cannot handle patterning for bar plots; greyscale shading is being used.")
    p  <-  p + scale_fill_grey(guide = guide_legend(direction="horizontal", ncol=5)) +
      theme_bw() + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
  } else { # use actual colours
    
    cols.needed  <-  length(unique(all.taxa$Taxon))
    
    if (cols.needed > 12) {
      # palette max is 12; so if we have more than 6 entries for otu1/2, we 
      # need to manually construct a palette to use
      
      col.func  <-  grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      p  <-  p + scale_fill_manual(values=col.func(cols.needed), 
                                 guide=guide_legend(direction="horizontal", ncol=5))
    } else {
      # otherwise, we can just use the Set3 palette from RColorBrewer
      p  <-  p + scale_fill_brewer(palette="Set3",
                                 guide=guide_legend(direction="horizontal", ncol=5))
    }
  }
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  }
  
  p
}


