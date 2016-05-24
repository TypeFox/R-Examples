group.Taxa.box <- function(data, is.OTU=TRUE, rank="g", 
                      taxa="", meta, meta.factor="",     
                      cex.y=5, cex.x=5, cex.main=10, RAM.theme=NULL, 
                      col.pal=NULL, main="", 
                      file=NULL, ext=NULL, height=8, width=16) {

  # check whether data is list
  if ( !is.list(data) ) {
      stop("Please provide input data as list, see ?RAM.input.formatting")
  }

  # valide rank
  rank.name <- .get.rank.name(rank, plural=TRUE)
  rank_name <- .get.rank.name(rank)
  rank_pat <- .get.rank.pat(rank)
  
  # get the groups
  if( !length(meta.factor)==1 ) {
    stop ("Error: please provide one factor")
  } else {
    fac.levels  <-  levels(factor(meta[[meta.factor]]))
  }

  if( length(taxa)==0 || taxa=="" ) {
        stop ("Error: please provide list of taxa to plot")
    }

  taxa.m  <-  list()
  total.count  <- list()

  # labels for each dataset
  labels <- names(data)

  # process each data set  
  # process each data set  
  tax.list <- list()
  nocap.list <- list()
  # find out what in the taxa list not present in all datasets
  for (i in 1:length(data)) {
      if ( is.null(data[[i]]) ) {
         break
      }

      label <- names(data)[i]
      elem <- data[[i]]

      if ( is.OTU ) {
          valid.OTU(elem)
          sum.count <- sum(transpose.OTU(elem))
          tax <- tax.abund(elem, rank=rank, drop.unclassified=TRUE, 
                       count=TRUE)
      } else {
          sum.count <- sum(elem)
          tax <- elem
      }

      if( length(which(names(tax) %in% taxa))==0) {
           stop ("Error: no taxa present in the otu table")
      }

      tax.sel <- tax[, which(names(tax) %in% taxa), drop=FALSE]
      nocap <- setdiff(taxa, names(tax.sel))
      if ( length(nocap) != 0 ) {
          warning(paste(length(nocap), " taxa are not in the ", label, " dataset", sep=""))
      } else {
        nocap <- NULL
      }

      tax.list[[label]] <- tax.sel
      nocap.list[[label]] <- nocap
  }

  # taxa that not in all datasets

  if ( is.null(unlist(nocap.list)) ) {
     taxa <- taxa
  } else {
    all.nocap <- Reduce(intersect, nocap.list)
    if ( length(all.nocap) != 0 ) {
       warning(paste(length(nocap), " taxa are not in all datasets: ", paste(all.nocap, collapse=", "), sep=""))
       rm <- which(taxa %in% all.nocap)
       # remaining taxa list
       taxa <- taxa[-rm]
    } else {
       taxa <- taxa
    }
  }
 
  # process the tax.abund matrices
  taxa.m <- list()
  total.count <- list()
  for (i in 1:length(tax.list)) {
    
    label <- names(tax.list)[i]
    tax <- tax.list[[i]] 
    if ( is.null(tax) ) { break }
    df <- data.frame(matrix(vector(), nrow(tax), 0))
    for ( i in taxa ) {
      if ( i %in% names(tax) ) {
        df[[i]] <- tax[[i]]
      } else {
        df[[i]] <- 0
      }
      rownames(df) <- rownames(tax)
      tax.sel <- df
    }
   
    # relative abundance of each taxa among samples
    tax.sel <- tax.sel[match(rownames(meta), 
                      rownames(tax.sel)) , , drop=FALSE]
    if ( identical(rownames(tax.sel), rownames(meta)) ) {
        tax.sel.p  <-  vegan::decostand(tax.sel, "total", MARGIN=2)
        tax.sel.p.fac  <-  cbind(tax.sel.p, meta[[meta.factor]])
        names(tax.sel.p.fac)[ncol(tax.sel.p.fac)] <- meta.factor
        # for calculate total.count
        tax.sel.fac = cbind(tax.sel, meta[[meta.factor]])
        names(tax.sel.fac)[ncol(tax.sel.fac)] <- meta.factor  
    } else {
      stop(paste(label, " and metadata do not have same samples", sep=""))
    }
  
    # shorten names
    df.count  <-  tax.sel.fac
    df.ra  <-  tax.sel.p.fac
 
    # reshape2::melt by sample
    #if ( require("reshape2") ) {
    #  reshape2::melt
    #} else {
    #   stop("package 'reshape2' is required to use this function")
    #}

    df.ra.m <- reshape2::melt(cbind(df.ra, sample=rownames(df.ra)),
                variable.name=rank_name, value.name="RA")
    names(df.ra.m)[ncol(df.ra.m)-1] <- rank_name
    names(df.ra.m)[ncol(df.ra.m)] <- "RA"
    df.count.m <- reshape2::melt(cbind(df.count, sample=rownames(df.count)),
                variable.name=rank_name, value.name="count")

    names(df.count.m)[ncol(df.count.m)-1] <- rank_name
    names(df.count.m)[ncol(df.count.m)] <- "count"

    #return(list(df.count.m, df.ra.m))
        # bind together with appropriate label
    df.ra.m  <-  cbind(df.ra.m, Region=label)
    df.count.m  <-  cbind(df.count.m, Region=label)

    #calculate total count of each Taxon
    total= sapply(levels(df.count.m[[rank_name]]), 
         function(x){sum(df.count.m$count[df.count.m[[rank_name]]==x])}, 
          USE.NAMES=F)
    total.df <- data.frame(Taxon=levels(df.count.m[[rank_name]]),total)
    total.df[["Region"]] <- label
    names(total.df)[1]  <-  rank_name
    # reorder levels based on order of appearance in total.count, 
     df.ra.m[[rank_name]]  <-  factor(as.character(df.ra.m[[rank_name]]),
                              levels=total.df[[rank_name]], 
                              ordered=TRUE)
    
   
      total.count[[label]]  <-  total.df
      taxa.m[[label]]  <-  df.ra.m
    }
  
  all.taxa  <-  do.call("rbind", taxa.m)
  all.count  <-  do.call("rbind", total.count)

  #return(list(all.taxa=all.taxa, all.count=all.count))
  # reorder by total counts
  all.count <- all.count[order(all.count$total, decreasing=FALSE),]
  all.count[[rank_name]] <- factor(all.count[[rank_name]],
                              levels=rev(levels(all.count[[rank_name]])), 
                              ordered=TRUE)
  all.taxa[[rank_name]] <- factor(all.taxa[[rank_name]],
                              levels=levels(all.count[[rank_name]]), 
                              ordered=TRUE)

  # plot title, default is "":
    title  <-  main

   # use actual colours
  cols.needed <- length(unique(all.taxa[[meta.factor]]))
  cols <- .cols(cols.needed=cols.needed, col.pal)
  names(cols) <- levels(factor(all.taxa[[meta.factor]])) 
      
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  # the distribution of each taxon in each category of meta.factor
  #if ( !require("scales") ) {
  #    stop("package 'scales' is required for this function")
  #}
 # if ( !require("RColorBrewer") ) {
 #    stop("package 'RColorBrewer' is required cex.y=5, cex.x=5, for this function")
 # }
 # if ( !require("grid") ) {
 #    stop("package 'grid' is required for this function")
 # }

  # make informative labels
  lab.lst <- list()
  for ( i in taxa ) {
      lab <- list()
      for ( j in 1:length(labels) ) {
          sel <- which(total.count[[labels[j]]][[rank_name]]==i)
          #lab[[labels[j]]] <- paste(labels[j], ":", 
          #         total.count[[labels[j]]]$total[sel], sep="")
          lab[[labels[j]]] <- round(total.count[[labels[j]]]$total[sel], digits=0)

      }      
      #lab1 <- paste(unlist(lab), collapse="\n")
      #lab.lst[[i]] <- paste(i, "\n", lab1, sep="")
      lab1 <- paste(unlist(lab), collapse=",")
      lab.lst[[i]] <- paste(i, "(", lab1, ")", sep="")
  }
  lab2 <- unlist(lab.lst)
  
 # boxplot 
 p <- ggplot2::ggplot(all.taxa, aes_string(x=rank_name,
            y="RA",fill=meta.factor)) + 
           geom_boxplot() 

  # facet by region 
  p <- p + scale_x_discrete(labels = lab2) + 
           facet_wrap(~Region) 
  # flip x & y
  p <- p +  scale_fill_manual(values=cols[names(cols) %in% 
                   levels(factor(all.taxa[[meta.factor]]))])+
            scale_y_continuous(labels=scales::percent_format()) + 
            xlab(.capitalize(.get.rank.name(rank))) +
            ylab("Percentage") + 
            coord_flip() +
            ggtitle(title) 

  if ( is.null(RAM.theme) ) {     
     p <- p + theme(legend.key.size=unit(6,"mm"), 
              legend.text = element_text(size=10),
              axis.text.y=element_text(size=cex.y, face="bold"),
              axis.text.x=element_text(size=cex.x, face="bold"),
              #axis.ticks.length=grid::unit(-0.1, "cm"), 
              #axis.ticks.margin=grid::unit(0.5, "cm"),
              axis.ticks.y=element_line(size=1),
              legend.position="right",
              plot.title = element_text(face="bold", size=cex.main))
  } else {
    p <- p + RAM.theme
  } 

  .valid.plot.settings(file, ext)
  save <- !is.null(file)
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
  
}


group.Taxa.bar <- function(data, is.OTU=TRUE, rank="g", taxa="", 
                        meta, meta.factor="", func="sum", cex.y=5, 
                        cex.x=5, cex.main=10, bar.width=NULL, 
                        RAM.theme=NULL, col.pal=NULL, main="", 
                        file=NULL, ext=NULL, 
                        height=8, width=16) {

  # valide rank
  rank.name <- .get.rank.name(rank, plural=TRUE)
  rank_name <- .get.rank.name(rank)
  rank_pat <- .get.rank.pat(rank)

  # check whether data is list
  if ( !is.list(data) ) {
      stop("Please provide input data as list, see ?RAM.input.formatting")
  }

  # get the groups
  .valid.factor(meta, meta.factor)
  if(!length(meta.factor)==1) {
     stop ("Error: please provide one factor")
  }
  if(length(taxa)< 1) {
    stop ("Error: please provide a vector of selected taxa names")
  }

  # labels for each dataset
  labels <- names(data)

  # process each data set  
  tax.list <- list()
  nocap.list <- list()
  # find out what in the taxa list not present in all datasets
  for (i in 1:length(data)) {
      if ( is.null(data[[i]]) ) {
         break
      }

      label <- names(data)[i]
      elem <- data[[i]]

      if ( is.OTU ) {
          valid.OTU(elem)
          sum.count <- sum(transpose.OTU(elem))
          tax <- tax.abund(elem, rank=rank, drop.unclassified=TRUE, 
                       count=TRUE)
      } else {
          sum.count <- sum(elem)
          tax <- elem
      }

      if( length(which(names(tax) %in% taxa))==0) {
           stop ("Error: no taxa present in the otu table")
      }

      tax.sel <- tax[, which(names(tax) %in% taxa), drop=FALSE]
      nocap <- setdiff(taxa, names(tax.sel))
      if ( length(nocap) != 0 ) {
          warning(paste(length(nocap), " taxa are not in the ", label, " dataset", sep=""))
      } else {
         nocap <- NULL
      }

      tax.list[[label]] <- tax.sel
      nocap.list[[label]] <- nocap
  }
  
  # taxa that not in all datasets
  if ( is.null(unlist(nocap.list)) ) {
     taxa <- taxa
  } else {
    all.nocap <- Reduce(intersect, nocap.list)
    if ( length(all.nocap) != 0 ) {
       warning(paste(length(nocap), " taxa are not in all datasets: ", paste(all.nocap, collapse=", "), sep=""))
       rm <- which(taxa %in% all.nocap)
       taxa <- taxa[-rm]
    } else {
       rm <- NULL
       taxa <- taxa
    }
    # remaining taxa list    
  }

  # process the tax.abund matrices
  taxa.m <- list()
  total.count <- list()
  
  for (i in 1:length(tax.list)) {
    label <- names(tax.list)[i]
    tax <- tax.list[[i]] 
    if ( is.null(tax) ) { break }
    df <- data.frame(matrix(vector(), nrow(tax), 0))
    for ( i in taxa ) {
      if ( i %in% names(tax) ) {
        df[[i]] <- tax[[i]]
      } else {
        df[[i]] <- 0
      }
    } 
    rownames(df) <- rownames(tax)     
    tax.sel <- df

    tax.sel.fac <- stats::aggregate(tax.sel, by=list(meta[[meta.factor]]), FUN=func)
   
    rownames(tax.sel.fac) <- tax.sel.fac$Group.1
    tax.sel.fac<-as.data.frame(tax.sel.fac[,-1, drop=FALSE])

    # reshape2::melt by factor
    #if (!require("reshape2")) {
    #   stop("package 'reshape2' is required to use this function")
    #}
    tax.sel.fac.m<-reshape2::melt(cbind(tax.sel.fac, 
                       factor=rownames(tax.sel.fac)), 
                       is.vars=c('factor'), value.name="count", 
                       variable.name=rank_name)
    
    names(tax.sel.fac.m)[ncol(tax.sel.fac.m)-1] <- rank_name
    names(tax.sel.fac.m)[ncol(tax.sel.fac.m)] <- "count"
    names(tax.sel.fac.m)[names(tax.sel.fac.m)=="factor"] <- meta.factor
  
    # bind together with appropriate label
    tax.sel.fac.m <- cbind(tax.sel.fac.m, Region=label)
 
    #calculate total count of each Taxon
    total= sapply(levels(tax.sel.fac.m[[rank_name]]), 
           function(x){sum(tax.sel.fac.m$count[tax.sel.fac.m[[rank_name]]==x])}, 
           USE.NAMES=F)
    total.df<-data.frame(Taxon=levels(tax.sel.fac.m[[rank_name]]),total)

    names(total.df)[1] <- rank_name
    total.df[["Region"]] <- label

    # reorder levels based on order of appearance in total.count, 
    tax.sel.fac.m[[rank_name]] <- factor(as.character(tax.sel.fac.m[[rank_name]]),
                              levels=total.df[[rank_name]], 
                              ordered=TRUE)
      
    # bind together with appropriate label
    total.count[[label]] <- total.df
    taxa.m[[label]] <- tax.sel.fac.m
  }
   
  all.taxa  <-  do.call("rbind", taxa.m)
  all.count  <-  do.call("rbind", total.count)

  # reorder by total counts
  all.count$total <- round(all.count$total, digits=0)
  #all.count <- all.count[order(all.count$total, decreasing=FALSE),]
  all.count[[rank_name]] <- factor(all.count[[rank_name]],
                              levels=rev(levels(all.count[[rank_name]])), 
                              ordered=TRUE)
  all.taxa[[rank_name]] <- factor(all.taxa[[rank_name]],
                              levels=levels(all.count[[rank_name]]), 
                              ordered=TRUE)

  #return(list(all.taxa, all.count))
  title = main

  # use actual colours
  cols.needed <- length(unique(all.taxa[[meta.factor]]))
  cols <- .cols(cols.needed, col.pal)
  names(cols) <- levels(factor(all.taxa[[meta.factor]])) 
    
   
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  # the distribution of each taxon in each category of meta.factor
  #if ( !require("scales") ) {
  #   stop("package 'scales' is required for this function")
  #}
 # if ( !require("RColorBrewer") ) {
  #   stop("package 'RColorBrewer' is required for this function")
 # }
 # if ( !require("grid") ) {
 #    stop("package 'grid' is required for this function")
 # }
  
  # make informative labels
  lab.lst <- list()
  for ( i in taxa ) {
      lab <- list()
      for ( j in 1:length(labels) ) {
          sel <- which(total.count[[labels[j]]][[rank_name]]==i)
          #lab[[labels[j]]] <- paste(labels[j], ":", 
          #         total.count[[labels[j]]]$total[sel], sep="")
          lab[[labels[j]]] <- round(total.count[[labels[j]]]$total[sel], digits=0)

      }      
      #lab1 <- paste(unlist(lab), collapse="\n")
      #lab.lst[[i]] <- paste(i, "\n", lab1, sep="")
      lab1 <- paste(unlist(lab), collapse=",")
      lab.lst[[i]] <- paste(i, "(", lab1, ")", sep="")
  }
  lab2 <- unlist(lab.lst)

  
 # barplot 
 if ( is.null(bar.width) ) {
    p <- ggplot2::ggplot(all.taxa, aes_string(x=rank_name,
            y="count",fill=meta.factor)) + 
           geom_bar(colour="white",stat="identity",
             position="fill") 
 } else {
   p <- ggplot2::ggplot(all.taxa, aes_string(x=rank_name,
            y="count",fill=meta.factor)) + 
           geom_bar(width=bar.width, colour="white",stat="identity",
             position="fill")
 }
           
  # facet by region,   # flip x & y
  p <- p + scale_fill_manual(values=cols[names(cols) %in% 
                   levels(factor(all.taxa[[meta.factor]]))]) +
           scale_y_continuous(labels=scales::percent_format()) +
           scale_x_discrete(labels = lab2) + 
           facet_wrap(~Region) + 
           ggtitle(title) +
           xlab(.capitalize(rank_name)) +
           ylab("Percentage") + 
           coord_flip()
           

  if ( is.null(RAM.theme) ) {
    p <- p + theme(legend.key.size=unit(6,"mm"), 
                   axis.text.y=element_text(size=cex.y, face="bold"),
                   axis.text.x=element_text(size=cex.x, face="bold"),
                   #axis.ticks.length=grid::unit(-0.25, "cm"), 
                   #axis.ticks.margin=grid::unit(0.5, "cm"),
                   axis.ticks.y=element_line(size=0.5),
                   legend.position="right",
                   legend.text = element_text(size=10), 
                   plot.title = element_text(face="bold", size=cex.main),
                   strip.text.x = element_text(size=8, face="bold", angle=0))
  } else {
    p <- p + RAM.theme
  } 

  .valid.plot.settings(file, ext)
  save <- !is.null(file)
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
}


group.abund.Taxa <- function(data, is.OTU=TRUE, rank="g", taxa, 
                             drop.unclassified=FALSE, bar.width=NULL,
                             meta, meta.factor="", RAM.theme=NULL, 
                             col.pal=NULL, main="",
                             file=NULL, ext=NULL, height=8, width=16) {
  
  save <- !is.null(file)
  
  # valide rank
  rank.name <- .get.rank.name(rank, plural=TRUE)
  rank_name <- .get.rank.name(rank)
  rank_pat <- .get.rank.pat(rank)
  
  # get the groups
  if( !length(meta.factor)==1 || 
        !any(names(meta) %in% meta.factor) ) {
    stop ("Error: please provide one factor")
  } 
  
  if( length(taxa)==0 || taxa=="" ) {
    stop ("Error: please provide list of taxa or top# of groups to plot")
  }
  
  # check whether data is list
  if ( !is.list(data) ) {
    stop("Please provide input data as list, see ?RAM.input.formatting")
  }
  
  # labels for each dataset
  labels <- names(data)
  
  # process each data set  
  tax.list <- list()
  for (i in 1:length(data)) {
    if ( is.null(data[[i]]) ) {
      break
    }
    
    label <- names(data)[i]
    elem <- data[[i]]
    
    if ( is.OTU ) {
      valid.OTU(elem)
      # get the groups
      tax <-.group.rank2(data=elem, meta=meta, meta.factor=meta.factor,
                        relative.abund=TRUE, taxa=taxa,  
                        drop.unclassified=drop.unclassified, 
                        rank=rank)  
    } else {
      tax <- .group.rank2(data=elem, is.OTU=FALSE, meta=meta, meta.factor=meta.factor,
                          relative.abund=TRUE, taxa=taxa,  
                          drop.unclassified=drop.unclassified, 
                          rank=rank) 
    }
    #print(tax)
    # reshape2::melt by Sample
    #if (!require("reshape2")) {
    #   stop("package 'reshape2' is required to use this function")
    # }
    if (ncol(tax) == 0 ) { 
      stop("selected taxa was not found")
    }
    tax.m <- reshape2::melt(cbind(tax, Sample=rownames(tax)), 
                            id.vars="Sample", value.name="RA", 
                            variable.name=rank_name)
    names(tax.m)[ncol(tax.m)-1] <- rank_name
    names(tax.m)[ncol(tax.m)] <- "RA"
    # bind together with appropriate label
    tax.m <- cbind(tax.m, Region=label)
    
    tax.list[[label]] <- tax.m
  }
  
  all.taxa <- do.call("rbind", tax.list)
  
  # reorder levels based on order of appearance in table
  all.taxa$Sample <- ordered(all.taxa$Sample, levels=unique(all.taxa$Sample))
  
  title <- main
  
  # if ( !require("scales") ) {
  #     stop("package 'scales' is required for this function")
  # }
  # if ( !require("RColorBrewer") ) {
  #     stop("package 'RColorBrewer' is required for this function")
  # }
  # if ( !require("grid") ) {
  #      stop("package 'grid' is required for this function")
  # }
  
  
  # use actual colours
  cols.needed <- length(unique(all.taxa[[rank_name]]))
  cols <- .cols(cols.needed=cols.needed, col.pal)
  names(cols) <- levels(factor(all.taxa[[rank_name]])) 
  
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  if ( is.null(bar.width) ) {
    p <- ggplot2::ggplot(all.taxa, aes_string(x="Sample", y="RA", 
                                              fill=rank_name)) + 
      geom_bar(position="stack", stat="identity") 
  } else {
    p <- ggplot2::ggplot(all.taxa, aes_string(x="Sample", y="RA", 
                                              fill=rank_name)) + 
      geom_bar(position="stack", stat="identity", width=bar.width)
  }
  
  if ( is.null(RAM.theme) ) {
    p <- p + theme(legend.position="bottom", 
                   axis.text.x=element_text(angle = 45, vjust = 1,  
                                            hjust=1),
                   panel.grid.major.x = element_blank())
  } else {
    p <- p + RAM.theme
  } 
  #coord_flip() +
  p <- p +  scale_y_continuous(labels = scales::percent_format()) +
    facet_wrap(~Region, scales="free_x") +
    xlab(meta.factor) +
    ylab("Relative Abundance") + 
    ggtitle(title)+
    scale_fill_manual(values=cols[names(cols) %in% 
                                    levels(factor(all.taxa[[rank_name]]))], 
                      guide=ggplot2::guide_legend(direction="horizontal", ncol=5))+
    theme(legend.position="bottom") 
  
  
  .valid.plot.settings(file, ext)
  save <- !is.null(file)
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
  
}

shared.Taxa <- function(data, is.OTU=TRUE, rank="g") {
  .valid.data(data=data, is.OTU=is.OTU)
  # valide rank
  if ( is.OTU ) {
    rank.name <- .get.rank.name(rank, plural=TRUE)
    rank_name <- .get.rank.name(rank)
    rank_pat <- .get.rank.pat(rank)
  } else {
    warning("data are not otu tables, will ignore the rank provided")
    rank_name <- "taxa"
    rank.name <- "taxa"
    rank_pat <- ""
  }
  
  # check whether data is list
  if ( !is.list(data) ) {
    stop("Please provide input data as list, see ?RAM.input.formatting")
  }
  
  # labels for each dataset
  labels <- names(data)
  
  # process each data set  
  shared.list <- list()
  for (i in 1:length(data)) {
    if ( is.null(data[[i]]) ) {
      break
    }
    
    label <- names(data)[i]
    elem <- data[[i]]
    
    if ( is.OTU ) {
      valid.OTU(elem)
      total.count <- sum(transpose.OTU(elem))
      tax <- tax.abund(elem, rank=rank, drop.unclassified=TRUE, 
                       count=TRUE)
    } else {
      total.count <- sum(elem)
      tax <- elem
    }
    
    ###check for zero entries?
    if(dim(tax)[2] == 1) {
      stop(paste("Error: Only one taxon at the ", rank_name, " level! Please choose another taxonomic rank to compare", sep=""))
    } else { 
      tax.pa <- decostand(tax, "pa")
    }
    
    # for readability
    x <- dim(tax.pa)[1] # number of samples
    y <- dim(tax.pa)[2] # number of taxa
    
    num.tax.one.sample <- dim(tax.pa[ ,colSums(tax.pa) == 1, drop=FALSE])[2]
    num.tax.mult.sample <- dim(tax.pa[ ,colSums(tax.pa) > 1, drop=FALSE])[2] 
    ### test this
    num.tax.all.sample <- dim(tax.pa[ ,colSums(tax.pa) == x, drop=FALSE])[2]
    
    per.tax.one.sample <- num.tax.one.sample / y
    per.tax.all.sample <- num.tax.all.sample / y
    
    if( length(which(colSums(tax.pa) == x ))==0 ) {
      warning("NO taxa are shared by all samples")
      num.seq.shared.tax <- 0
    } else {
      num.seq.shared.tax <- sum(tax[ , colnames(tax) %in% 
                                      names(tax.pa)[which(colSums(tax.pa) == x)]])
    }
    per.seq.shared.tax <- num.seq.shared.tax / total.count
    
    # get the taxa names for all taxs present in all samples
    tax <- names(tax.pa)[which(colSums(tax.pa) == x)]
        
    val <- list(num.tax.one.sample, num.tax.mult.sample, 
                num.tax.all.sample, per.tax.one.sample, 
                per.tax.all.sample, num.seq.shared.tax,
                per.seq.shared.tax, tax)
    
    # format numerical data?
    #val <- format(round(val[1:7], 2), nsmall=2)
    
    names(val) <- c(paste("#_of_", rank.name, "_in_1_sample", sep=""), 
                    paste("#_of_", rank.name, "_in_>1_sample", sep=""), 
                    paste("#_of_", rank.name, "_in_all_samples", sep=""),
                    paste("%_of_", rank.name, "_in_one_sample", sep=""),
                    paste("%_of_", rank.name, "_in_all_samples", sep=""), 
                    paste("#_of_sequence_in_shared_", rank.name, sep=""),
                    paste("%_of_sequence_in_shared_", rank.name, sep=""), 
                    paste(rank.name, "_in_all_samples", sep=""))
    
    shared.list[[label]] <- val
  }
  return(shared.list)
}


.cols <- function(cols.needed, col.pal=NULL) {
  if ( is.null(col.pal) ) {    
    if ( cols.needed <= 2 ) {
      cols <- c("#E41A1C", "#377EB8")
    } else if ( cols.needed <= 12 && cols.needed > 2 ) {
      cols <- RColorBrewer::brewer.pal(cols.needed, "Set3")
    } else {
      
      # if palette max is 12; so if we have more than 6 entries for otu1/2, we 
      # need to manually construct a palette to use
      
      #cols <- rainbow(cols.needed) 
      col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      cols <- col.func(cols.needed)
    } 
  } else {
    cols <- col.pal
  }
}
        

.group.rank2<-function(data, is.OTU=TRUE, meta=meta, meta.factor="",
                       drop.unclassified=FALSE, relative.abund=FALSE, 
                       rank="g", func="sum", taxa=NULL){
  ## whether data an OTU or taxonomic abundance matrix
  if ( is.OTU ) {
    valid.OTU(data)
    tax <- tax.abund(data, rank=rank, drop.unclassified=drop.unclassified, 
                     count=TRUE)
  } else {
    tax <- data
  }
  
  # get the groups
  if( !length(meta.factor)==1 || 
        !any(names(meta) %in% meta.factor) ) {
    stop ("Error: please provide one factor")
  } 
  
  if(all(rownames(tax)==rownames(meta))) {
    tax.factor <- stats::aggregate(tax, by=list(meta[[meta.factor]]), FUN=func)
    rownames(tax.factor) <- tax.factor[,"Group.1"]
    tax.factor <- tax.factor[,-1]
  } else {
    stop("Error: otu and metadata have different samples")
  }
  tax.factor <- tax.factor[, order(colSums(tax.factor), decreasing=TRUE)]
  if(!isTRUE(relative.abund)) {
    tax.factor <- tax.factor
  } else {
    tax.factor <- vegan::decostand(tax.factor, "total")    
  }
  
  if (drop.unclassified) {
    # this selects all columns NOT containing in the blacklist
    # drop.unclassified using blacklist
    remove.pat <- gsub(.get.rank.pat(rank), "", 
                       paste0(.blacklist(.get.rank.pat(rank)), 
                              "|no_taxonomy"))
    tax.factor <- tax.factor[ , !grepl(remove.pat, names(tax.factor), 
                                       ignore.case=TRUE), drop=FALSE]
  } else {
    tax.factor<-tax.factor
  }
  # keep only the 'top' most abundant groups, where top is user-given 
  if ( is.null(taxa) ) {
    tax.factor.sel<-tax.factor
  } else if ( is.numeric(taxa) ) {
    if ( taxa <= ncol(tax.factor) ) {
      tax.factor.sel<-tax.factor[,1:taxa]
    } else {
      tax.factor.sel<-tax.factor
    }
  } else if ( is.character(taxa) && length(taxa) >= 1 ) {
    if( !any(names(tax.factor) %in% taxa) ) {
      stop ("Error: no taxa present in the dataset")
    } else {
      tax.factor.sel <- tax.factor[, which(names(tax.factor) %in% taxa), drop=FALSE]
      nocap <- setdiff(taxa, names(tax.factor.sel))
      if ( length(nocap) != 0 ) {
        warning(paste("Taxa not in the dataset: ", paste(nocap, collapse=" ,"), sep=""))
        # if a taxon not being found in dataset, add 0 
        for (i in nocap) {
          tax.factor.sel[[i]] <- 0
        } 
      } 
    }
  } else {
    tax.factor.sel<-tax.factor
  }
  
  return(tax.factor.sel)    
}
