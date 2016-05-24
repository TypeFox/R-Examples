correlation<-function(data=NULL, is.OTU=TRUE, meta=NULL, rank="g", 
                      sel=NULL, sel.OTU=TRUE, data.trans=NULL, 
                      method="pearson", main=NULL, file=NULL, 
                      ext=NULL, width=8, height=8) {
  
  save <- !is.null(file)
  if (save) { .get.dev(file, ext, height=height, width=width) }
  
  # make sure either data or metadata was provided
  if ( is.null(data) & is.null(meta) ) {
    stop("Error: provide at least one of the following: OTU table, taxonomy abundance matrix, or metadata")
  } 
  
  if ( !is.null(meta) ) {
    # metadata of numeric variables only
    col <- sapply(meta, is.numeric)
    meta.sel <- meta[, col, drop=FALSE]
    if ( ncol(meta) == 0 ) {
      warning("no columns in metadata is numeric, will only plot data")
      meta.sel <- data.frame(matrix(vector(), nrow(meta.sel), 0))
    } else {
      meta.sel <- meta.sel
    }
  } else {
    if ( is.OTU ) {
      meta.sel <- data.frame(matrix(vector(), ncol(data)-1, 0))
    } else {
      meta.sel <- data.frame(matrix(vector(), nrow(data), 0))
    }
  }
  
  # validate data and metadata
  if ( !is.null(data) & !is.null(meta) ) {
    if ( is.OTU ) {
      valid.OTU(data)
      .valid.meta(otu1=data, meta=meta)         
    } else {
      data <- data[match(rownames(meta), rownames(data)), , drop=FALSE]
      if ( !identical(rownames(data), rownames(meta)) ) {
        stop("Error: not same subjects in data and metadata")
      } 
    }
  } else {
    data <- data
  }
  
  # require count data
  if ( !is.null(data) ) {
    # check whether integer data being provided
    if ( is.OTU ) {
      int <- unique(sapply(data[, -ncol(data), drop=FALSE], is.integer))
    } else {
      int <- unique(sapply(data, is.integer))
    }
    if ( !isTRUE(int) & !is.null(data.trans) ) {
      warning("data contain non-integer records, are you sure to transform the data?")
    }
  }
  
  if ( !is.null(data) ) { 
    # plot OTU
    if ( sel.OTU ) {
      if ( !is.OTU ) {
        stop("Error: provided outIDs but data is not an OTU table")
      } else {
        if ( is.null(data.trans) ) {
          data.sel <- data
        } else {
          data.sel <- cbind(vegan::decostand(data[, -ncol(data), drop=FALSE], 
                                      MARGIN=2, data.trans), data$taxonomy)
          names(data.sel)[ncol(data.sel)] <- "taxonomy"
        }
      }
      if ( !is.null(sel) ) {
        sel.otu<- which(rownames(data.sel) %in% sel)
        if ( length(sel.otu) != 0 ) {
          data.sel <- data.sel[sel.otu, ]
        } else {
          warning("did not find provided otuIDs in data,
                  use all data")
          data.sel <- data.sel
        }
      } else {
        data.sel <- data.sel
      }
      if ( class(data.sel$taxonomy) != "character" ) {
        data.sel$taxonomy<-as.character(data.sel$taxonomy)
      }
      data.tax <- LCA.OTU(data.sel, strip.format=FALSE, drop=TRUE)
      rownames(data.tax) <- paste( rownames(data.tax), 
                                   data.tax$LCA, sep="_")
      data.tax <- data.tax[, -ncol(data.tax), drop=FALSE]
      data.tax <- as.data.frame(t(data.tax))
      }
    
    # plot taxa at given rank    
    if ( !sel.OTU ) {
      if ( is.null(rank) ) {
        stop("please provide rank of the selected taxa")
      }
      #get rank names
      rank.name <- .get.rank.name(rank, plural=TRUE)
      rank_name <- .get.rank.name(rank)
      rank_pat <- .get.rank.pat(rank)
      
      if ( is.OTU ) {
        data.tax <- tax.abund(data, count=TRUE, rank=rank, 
                              drop.unclassified=FALSE)
      } else {
        data.tax <- data
      }
      # transform data
      if ( !is.null(data.trans) ) {
        data.tax <- decostand(data.tax, data.trans)
      } else {
        data.tax <- data.tax
      }
      # remove unclassified taxa
      remove.pat <- gsub(.get.rank.pat(rank), "", 
                         paste0(.blacklist(.get.rank.pat(rank)), "|no_taxonomy"))
      data.tax <- data.tax[, !grepl(remove.pat, names(data.tax), 
                                    ignore.case=TRUE), drop=FALSE]
      # select only taxa in sel
      if ( !is.null(sel) ) {
        sel.tax <- which(colnames(data.tax) %in% sel)
        if ( length(sel.tax) != 0 ) {
          data.tax <- data.tax[ , sel.tax, drop=FALSE]
        } else {
          warning("No taxa in the provided list was found in data, 
                  will plot all taxa")
          data.tax <- data.tax
        }
      } else {
        data.tax <- data.tax
      }
      } 
  } else {
    data.tax <- data.frame(matrix(vector(), nrow(meta), 0))
  }
  
  # combine data
  dat <- cbind(data.tax, meta.sel)
  
  if ( ncol(dat) == 0 ) {
    stop("Error: nothing to plot")
  } else if ( ncol(dat) > 200 ) {
    stop("Too many variables to be plotted, please use sel 
         option to reduce the number of variables")
  } else {
    
    # method <- c("pearson", "kendall", "spearman")
    cor.mat<-cor(dat, method=method, use = "pairwise.complete.obs")
    x.scale <- list(cex=0.5, alternating=1, col='black') 
    y.scale <- list(cex=0.5, alternating=1, col='black') 
    
    # levelplot
    bp1 <- lattice::levelplot(cor.mat,xlab=NULL,ylab=NULL,cex=0.08,  
                     at=do.breaks(c(-1.01,1.01),101),
                     scales=list(x=list(rot=90)),
                     colorkey=list(space="top"), 
                     col.regions=colorRampPalette(c("red","white","blue")), 
                     main=main,cex=1)
    print(bp1)
    
    if (save) { dev.off() }
    invisible()
  }
  
}

reset.META <- function(meta, factor=NULL, numeric=NULL, date=NULL) {

  if( !is.null(factor) ) {
    for ( i in factor ) {
      if ( ! (i %in% names(meta) ) ) {
        warning(paste (i," is not in in metadata!", sep=""))
      } else {
        meta[[i]] <- as.factor(meta[[i]])
      }
    }
    meta<-meta
  } else {
    meta<-meta
  }

  if( !is.null(numeric) ) {
    for ( i in numeric ) {
      if ( ! ( i %in% names(meta) ) ) {
         warning(paste (i," is not in in metadata!", sep=""))
      } else if ( !is.numeric(meta[[i]]) ) {
         warning(paste(i, " is not a numeric data type, but will
                       be converted to numeric as required!", sep=""))
         meta[[i]] <- as.numeric(as.factor(meta[[i]]))
      } else {
        meta[[i]] <- as.numeric(meta[[i]])
      }
    }
    meta<-meta
  } else {
    meta<-meta
  }

  if( !is.null(date) ) {
    for ( i in date ) {
      if( ! ( i %in% names(meta) ) ) {
        warning(paste (i," is not in in metadata!", sep=""))
      } else {
        meta[[i]] <- as.Date(meta[[i]])
      }
    }
    meta<-meta
  } else { 
    meta<-meta
  }
  
  return(meta)
}


sample.map <- function(meta, siteID="City", maptype="roadmap",
                       shape=16, colour=c("red", "blue"),
                       lat="Latitude", lon="Longitude", zoom=3, 
                       file=NULL, ext=NULL, width=10, height=10) {
  
  save <- !is.null(file)
  if (save) { .get.dev(file, ext, height=height, width=width) }
  
  if (length(siteID)==0 || length(siteID)>1 || siteID=="" ) {
    stop("Error: please provide ONE uniuqe siteID for each pair of Latitude and Longitude record")
  } 
  
  if ( !(siteID %in% names(meta)) ) {
    stop(paste(siteID, " is not in the metadata", sep=""))
  }
  
  if( length(levels(factor(meta[[siteID[1]]])))>6 ) {
    stop("Error: Too many siteID levels, please use RAM::sampling.map instead")
  }  
  
  if (! isTRUE(lat %in% names(meta)) & (lon %in% names(meta))) {
    stop ("Error: no Latitude or Longitude information, please check metadata")
  }
  
  
  sample <- .samples(meta=meta, siteID=siteID, lat=lat, lon=lon)
  
  
  #if ( !require("gtable") ) {
  if(!requireNamespace('gtable')) {
    stop("package 'gtable' is required for this function")
  }
  
  #if ( !require("ggmap") ) {
  #  stop("package 'ggmap' is required for this function")
  # }
  #if ( !require("grid") ) {
  #   stop("package 'grid' is required for this function")
  #}
  
  if ( !requireNamespace("mapproj") ) {
    stop("package 'mapproj' is required for this function")
  }
  # if ( !require("RColorBrewer") ) {
  #   stop("package 'RColorBrewer' is required for this function")
  # }
  
  
  #osmMap <- get_map(location="Canada", zoom=4, source = 'google')
  map_type = c("terrain", "satellite", "roadmap", "hybrid", 
               "toner")
  if ( maptype %in%  map_type) {
    location <- c(mean(sample$Longitude), mean(sample$Latitude))
    osmMap <- ggmap::get_map(location=location, zoom=zoom, 
                      maptype=maptype)
  } else {
    warning(paste("maptype must be one of the following: ", map_type, ". Will use roadmap as default", sep=""))
    osmMap <- ggmap::get_map(location=location, zoom=zoom, 
                      maptype="roadmap")
  }
  #return(sample)
  # points1: color gradients for number of samples collected
  #points1 <- geom_point(data = sample, aes_string(x = "Longitude", 
  #                y = "Latitude", colour = "Freq", shape="siteID"), 
  #                    alpha = 0.6,  size = 10)
  points1 <- ggplot2::geom_point(data = sample, aes_string(x = "Longitude", 
                                                  y = "Latitude", colour = "Freq"), shape=shape,
                        alpha = 0.6,  size = sample$Freq)
  # points2: shape for locations
  #points2<-geom_point(data=sample, aes(Longitude, Latitude, 
  #                    shape=rownames(sample)), size=3)
  
  p <- ggmap::ggmap(osmMap) + 
    points1 + 
    scale_colour_continuous(name="Frequency", low = colour[1], 
                            high = colour[2], space = "Lab", guide = "colorbar") + 
    scale_shape_discrete(name=siteID, breaks=rownames(sample),
                         labels=paste(rownames(sample), sep=":", sample$Freq)) + 
    labs(title="Sampling Locations") + 
    theme(plot.title=element_text(size=14, colour="black", 
                                  face="bold", family="serif")) + xlab("Longitude") + 
    ylab("Latitude") + 
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12, 
                                  face="bold", family="serif"))
  
  assign("last.warning", NULL, envir = baseenv())  
  
  print(p)
  
  # capture warnings: e.g. zoom issues
  if( length(tryCatch(print(p))) > 0) {
    warning("Change zoom setting if missing points")
  }
  
  #  return(sample)
  
  if (save) { dev.off() }
  
  invisible()
}


sample.sites <- function(meta, siteID="City", marker.size="small", 
                         lat="Latitude", lon="Longitude", maptype="hybrid", 
                         zoom=5, file=NULL, ext=NULL, width=8, height=8) {
  
  save <- !is.null(file)
  if (save) { .get.dev(file, ext, height=height, width=width) }
  
  if (length(siteID)==0 || length(siteID) > 1 || siteID=="" ) {
    stop("Error: please provide ONE uniuqe siteID for each pair of Latitude and Longitude record")
  } 
  
  if ( !(siteID %in% names(meta)) ) {
    stop(paste(siteID, " is not in the metadata", sep=""))
  }
  
  if (! isTRUE(lat %in% names(meta)) & (lon %in% names(meta))) {
    stop ("Error: no Latitude or Longitude information, please check metadata")
  }
  
  sample <- .samples(meta=meta, siteID=siteID, lat=lat, lon=lon)
  
  if ( !requireNamespace("gtable") ) {
    stop("package 'gtable' is required for this function")
  }
  
  #if ( !requireNamespace("ggmap") ) {
  #  stop("package 'ggmap' is required for this function")
  #}
  
  #if ( !requireNamespace("grid") ) {
  #  stop("package 'grid' is required for this function")
  #}
  
  if ( !requireNamespace("mapproj") ) {
    stop("package 'mapproj' is required for this function")
  }
  
  #if ( !requireNamespace("RColorBrewer") ) {
  #  stop("package 'RColorBrewer' is required for this function")
  #}    
  
  #if ( !requireNamespace("RgoogleMaps") ) {
  #  stop("package 'RgoogleMaps' is required for this function")
  #}
  
  sample.2<-sample
  sample.2$size <- marker.size  #create a column indicating size of marker
  sample.2$col <- "red"   #create a column indicating color of marker
  sample.2$char <- ""   #normal Google Maps pinpoints will be drawn
  names(sample.2)[1]<-"lat"
  names(sample.2)[2]<-"lon"
  
  range.lat <- range(sample.2$lat) #define our map's ylim
  range.lon <- range(sample.2$lon) #define our map's xlim
  center = c(mean(range.lat), mean(range.lon))  #tell what point to center 
  zoom <- min(RgoogleMaps::MaxZoom(range.lat, range.lon))
  
  # map range
  map.range <- RgoogleMaps::qbbox(lat = sample.2[,"lat"], lon = sample.2[,"lon"]);
  #terrain_close <- GetMap.bbox(lonR= range.lon, latR= range.lat, center= center, destfile= "terrclose.png", markers= sample.2[,c(1,2,4,5,6)],  zoom=6, maptype="terrain")
  
  map_close <- RgoogleMaps::GetMap.bbox(lonR= range.lon, latR= range.lat, 
                           center= center, markers= sample.2[,c(1,2,4,5,6)],  
                           zoom=zoom, maptype=maptype)
  RgoogleMaps::PlotOnStaticMap(map_close)
  
  if (save) { dev.off() }
  
  invisible()
}

.make.siteID <- function(meta, lat, lon) {
  factor<-c(lat, lon)
  .valid.factor(meta, factor)
  env <- meta[, which(names(meta) %in% factor), drop=FALSE]
  env$ll <- paste(env[[lat]], env[[lon]], sep=",")
  env$siteID <- NA
  env.uni <- unique(env)
  env.uni$ll <- paste(env.uni[[lat]], env.uni[[lon]], sep=",")
  env.uni$siteID <- paste0("sample_", 1:nrow(env.uni))
  for ( i in 1:nrow(env.uni) ) {
    n <- which(env$ll == env.uni[i, "ll"])
    env[n, "siteID"] <- env.uni[i, "siteID"]
  }    
  env <- env[, c("siteID", lat, lon)] 
  return(env)
} 

.samples <- function(meta, siteID=NULL, lat="Latitude", lon="Longitude") {
  
  if ( is.null(siteID) ) {
    env <- .make.siteID(meta=meta, lat=lat, lon=lon)
  } else {
    #drop levels if metadata is a subset of a bigger dataset.
    meta[[siteID]] <- factor(as.character(meta[[siteID]]))
    
    # only siteID, lat, lon from metadata
    sel <- which(names(meta) %in% c(siteID, lat, lon))
    env <- meta[, sel, drop=FALSE]
    names(env)[names(env)==siteID] <- "siteID"
    
    if ( nrow(unique(env[, c(lat, lon)])) != length(unique(env[,"siteID"]))) {
      warning("siteID do not match latitude and logitude records, will not use the site IDs")
      env <- .make.siteID(meta=meta, lat=lat, lon=lon)
    }
  }
  siteID <- "siteID"
  count <- as.data.frame(t(table(env[[siteID]])))
  rownames(count) <- count$Var2
  count <- count[order(rownames(count)),-1]
  
  env_sub <- unique(env)
  rownames(env_sub) <- env_sub[,1]
  env_sub <- env_sub[order(rownames(env_sub)), -1]
  names(env_sub)[1] <- "Latitude"
  names(env_sub)[2] <- "Longitude"
  env_sub <- env_sub[match(rownames(count), rownames(env_sub)), , drop=FALSE]  
  
  if( identical(rownames(count), rownames(env_sub)) ) {
    sample<-cbind(env_sub, count)
    sample<-sample[, -3]
    sample[["siteID"]] <- rownames(sample)
  } else {
    stop("error: not all sites has longitude and latitude info")
  }
  
  return(sample)
}

group.diversity <- function(data, meta, factors="", indices="", 
                            diversity.info=FALSE,  
                            x.axis=NULL, compare=NULL, 
                            facet=NULL, facet.y=TRUE, 
                            facet.x.cex=NULL, 
                            facet.y.cex=NULL, scale.free=NULL, 
                            xlab=NULL, ylab=NULL, 
                            legend.title=NULL, legend.labels=NULL,                 
                            file=NULL, ext=NULL, width=8, height=8) {
  #  x can be one of the metadata variable or "SampleID"
  save <- !is.null(file)
  
  labels <- names(data)
  for ( i in 1:length(data) ) {
    valid.OTU(otu1=data[[i]])
    if ( is.null(data[[i]]) ) { break }
    .valid.meta(data[[i]], meta=meta)
  }
  
  # calcluate diveristy 
  if ( diversity.info ) {
    warning("make sure using RAM::OTU.diversity output for this function")
    meta.diversity <- meta
  } else {
    meta.diversity <- OTU.diversity(data=data, meta=meta)
  }
  
  # valid metadata variables
  meta.fac <- unique(c(factors, x.axis, compare, facet))
  
  .valid.factor(meta, meta.fac) 
  
  Indices <- c("spec", "sim", "invsim", "shan", "sim_even", "shan_even", 
               "sim_trudiv", "shan_trudiv", "chao", "ACE")
  
  if ( length(indices) !=0 || indices != "" ) {
    if ( !any(indices %in% Indices) ) {
      stop(paste("indices has to be one or more of the following: ", paste(Indices, collapse=", "), 
                 " See ?OTU.diversity for details!", sep=""))
    } else {
      clp <- vector()
      # colnames of the indices in meta.diversity (with labels)
      for ( i in indices ) {
        clp <- c(clp, paste(i, labels, sep="_"))
      }
      div.ind <- which(names(meta.diversity) %in% clp)
      div.ind <- names(meta.diversity)[div.ind]
    }  
  } else {
    stop(paste("what diversity indices to be compared? indices has to be one or more of the following: ", 
               paste(Indices, collapse=", "), " See ?OTU.diversity for details!", sep=""))
  }
  
  #if ( !requireNamespace("reshape2") ) {
  #   stop("package 'reshape2' is required for this function")
  #}
  
  if(  length(factors) ==0 || factors=="" || !any(factors %in% names(meta)) ) {
    stop("Error: please provide metadata factor(s) for comparison ")
  } else {
    meta.m <- reshape2::melt(cbind(meta.diversity[, which(colnames(meta.diversity) 
                                                %in% c(meta.fac, div.ind)), drop=FALSE], 
                         SampleID=rownames(meta.diversity)), 
                   variable.name="Index", value.name="Value")
    names(meta.m)[ncol(meta.m)] <- "Value"
    names(meta.m)[ncol(meta.m)-1] <- "Index"
  }
  
  len <- length(indices)
  if ( len == 1 ) {
    levels(meta.m$Index) <- labels
  } else {
    levels(meta.m$Index) <- levels(meta.m$Index)
  }
  
  if ( is.null(x.axis) || x.axis=="") {
    warning(paste("you did not provide x.axis to plot, will use: ", factors[1], " as default!, See ?group.diversity for detail", sep=""))
    x.axis <- factors[1]
  } else {
    x.axis <- x.axis
  }
  
  if ( is.null(compare) || compare=="" ) {
    warning(paste("you did not provide a metadata factor to compare, will use: ", factors[1], " as default!, See ?group.diversity for detail", sep=""))
    compare <- factors[1]
  } else {
    compare <- compare
  }
  
  p <- ggplot(meta.m, aes_string(x=x.axis, y="Value", color=compare)) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, hjust=1, size=8))
  
  if ( is.null(facet.x.cex) ) { 
    p <- p + theme(strip.text.x = element_text(face="bold"))
  } else {
    p <- p + theme(strip.text.x = element_text(size=facet.x.cex, 
                                               face="bold"))
  }
  
  if ( is.null(facet.y.cex) ) { 
    p <- p + theme(strip.text.y = element_text(face="bold"))
  } else {
    p <- p + theme(strip.text.y = element_text(size=facet.y.cex,
                                               face="bold"))
  }
  
  
  if (is.null(legend.title)) {
    p <- p
  } else {
    p <- p + scale_color_discrete(name = legend.title,
                                  labels=legend.labels)
  }
  
  if ( is.null(facet) || length(facet) == 0 || facet=="" ) {
    p <- p + facet_grid(. ~ Index)
  } else {
    if( facet.y ) {
      formula <- paste("Index", facet[1], sep=" ~ ")
    } else {
      formula <- paste(facet[1], "Index", sep=" ~ ")
    }
    if( is.null(scale.free) ) {
      p <- p + facet_grid(as.formula(formula))
    } else {
      p <- p + facet_grid(as.formula(formula), scales="free", 
                          space=scale.free)
    }
  }
  
  if (is.null(xlab)) {
    p <- p
  } else {
    p <- p+xlab(xlab)
  }
  
  if (is.null(ylab)) {
    p <- p
  } else {
    p <- p + ylab(ylab)
  }
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
  
}

filter.META <- function(meta=meta, excl.na=TRUE, excl.NNF=TRUE, 
                        exclude=NULL) {
  
  # remove empty strings with NA
  #meta[ meta == "" ] <- NA
  # convert all character columns to as factor
  meta <- .CharToFac(meta)
  
  # whether provided list of variables to be excluded
  if ( !is.null(exclude) ) {
    .valid.factor(meta, exclude)
    suppressWarnings(excl <- .excl.col(meta, exclude=exclude))
  } else {
    excl <- NULL
  }
  # whether or not exclude columns with missing data
  if ( excl.na ) {
    suppressWarnings(na.col <- .na.col(meta))
  } else {
    na.col <- NULL
  }
  # columns with only one level   
  suppressWarnings(level1.col <- .level1.col(meta))
  # NNF stands for non-numeric/factor variables
  # excl.NNF=TRUE: exclude NNF; 
  suppressWarnings(NNF.col <- .NNF.col(meta))
  
  # combine all column numbers needs to be removed
  remove.all <- unique(c(excl, na.col, level1.col, NNF.col))
  if( length(remove.all) == 0 ) {
    meta.new <- meta
  } else {
    meta.new <- meta[,-c(remove.all), drop=FALSE]
  }
  
  print(paste("A total of ", length(remove.all), 
              " columns are removed: ", 
              paste(names(meta)[remove.all], collapse=", ")))
  print(paste("columns in exclude list: ", sep="", 
              paste(sort(excl), collapse=", ")))
  print(paste("column# with missing data: ", sep="", 
              paste(sort(na.col), collapse=", ")))
  print(paste("columns# with one level: ", sep="", 
              paste(sort(level1.col), collapse=", ")))
  print(paste("column# that are not numeric or factor: ", 
              sep="", paste(NNF.col, collapse=" ")))  
  return(meta.new)
  
}

# columns with missing data
.na.col <- function(df) {
 df.na <- df[sapply(df, function(x) any(is.na(x)))] 
  if ( ncol(df.na) == 0 ) {
      na.col <- NULL
  } else {
      na.col <- which(names(df) %in% names(df.na))
  }
  return(na.col)
}

# convert all character columns to as factor
.CharToFac <- function(df) {
  for ( i in 1:ncol(df) ) {
    if( is.character(na.omit(df[, i])) ) {
      df[,i] <- as.factor(df[, i])
    }
  }
  return(df)
}

# check variables only has one level
.level1.col <- function(df) {
  df <- .CharToFac(df)
  df.l1 <- df[sapply(df, function(x)
                  length(levels(factor(na.omit(x)))) == 1 )]
  if ( ncol(df.l1) == 0 ) {
      level1.col <- NULL
  } else {
      level1.col <- which(names(df) %in% names(df.l1))
  }
  return(level1.col)
}

# NNF stands for non-numeric/factor variables
  # excl.NNF=TRUE: exclude NNF; 
.NNF.col <- function(df) {
  df <- .CharToFac(df)
  df.NNF <- df[sapply(df, function(x) 
              !is.numeric(na.omit(x)) & !is.factor(na.omit(x)) ) ]
  if ( ncol(df.NNF) == 0 ) {
      NNF.col <- NULL
  } else {
      NNF.col <- which(names(df) %in% names(df.NNF))
  } 
  return(NNF.col)
}

# select only factor/character variables
.fac.col <- function(df) {
  df <- .CharToFac(df)
  fac.col <- numeric()
  df.fac <- df[sapply(df, function(x) 
              is.factor(na.omit(x)) ) ]
  if ( ncol(df.fac) == 0 ) {
      fac.col <- NULL
      #print("No variables are factor data type")
  } else {
      fac.col <- which(names(df) %in% names(df.fac))
      #print(paste("The following variables are factors: ", 
      #            paste(names(df)[fac.col], collapse=", "), sep=""))
  } 

  return(fac.col)
}

# select only factor/character variables
.num.col <- function(df) {
  num.col <- numeric()
  df.num <- df[sapply(df, function(x) 
              is.numeric(na.omit(x)) ) ]
  if ( ncol(df.num) == 0 ) {
      num.col <- NULL
  } else {
      num.col <- which(names(df) %in% names(df.num))
  } 
  return(num.col)
}
           
# exclude columns of a df
.excl.col <- function(df, exclude) {
  df <- .CharToFac(df)
  excl <- numeric()
  if ( is.null(exclude) ) {
    # print("No columns in exclude list")
    excl <- NULL 
  } else {
    # obtain column numbers of excluded variables
    excl <- numeric()
    # column numbers in exclude list
    if ( is.numeric(exclude) ) {
      if ( all(exclude %in% 1:ncol(df)) ) {
        excl <- exclude
      } else {
        excl <- exclude[ which( exclude %in% 1:ncol(df) ) ]
        excluded <- exclude[-which(1:ncol(df) %in% exclude)]
        warning( paste("data frame does not have ", sep="", 
                       paste(excluded, collapse=", "), " column numbers", 
                       ", will be ignored"))
      }
    }
    # column names in exclude list
    if ( is.character(exclude) ) {
      if ( all(exclude %in% names(df)) ) {
        excl <- which(names(df) %in% exclude)
      } else {
        excl <- which(names(df) %in% exclude)
        excluded <- exclude[!grepl(paste(names(df), 
                                         collapse="|"), exclude)]
        warning( paste("data frame does not have variables 
                 named ", sep="", paste(excluded, collapse=", "), 
                       ", will be ignored"))
      }
    }
  }
  return(excl)
}





