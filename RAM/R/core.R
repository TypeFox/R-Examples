core.OTU <- function(data, meta, meta.factor="", percent=1){
  
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list. See ?RAM.input.formatting for details")
  }
  
  .valid.data(data=data, is.OTU=TRUE)
  
  list.core <-list()
  labels <- names(data)
  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    otu <- data[[i]]
    if ( is.null(otu) ) { break }
    
    # valid otu & meta
    valid.OTU(otu)
    
    #.valid.meta(otu1=otu, meta=meta)
    samp <- colnames(otu)[1:(ncol(otu)-1)]
    meta <- meta[match(samp, rownames(meta)), ]
    if ( !identical(rownames(meta), samp) ) {
      stop("Error: metadata and otu table do not have same samples")
    }
    
    # check meta category
    .valid.factor(meta, meta.factor)
    fac.levels <- levels(factor(meta[[meta.factor]]))
    otu <- LCA.OTU(otu, strip.format=FALSE, drop=TRUE)
    otu$LCA <- paste(otu$LCA, rownames(otu), sep="_")
    
    ###check for zero entries?
    otu.fac.list=list()
    otu.fac.id = list()
    otu.fac.rank = list()
    for (j in 1:length(fac.levels)) {
      # each fac levels
      meta.fac <- meta[grep(fac.levels[j], meta[[meta.factor]]), ]
      
      # presence/absence
      otu.pa <- decostand(otu[, -ncol(otu)], "pa") 
      
      # select only otus present in more than percent of samples 
      # in each category
      col <- which(colnames(otu.pa) %in% rownames(meta.fac))
      if ( length(col) == 1 ) {
        vec <- otu.pa[, col]
        sel <- which(vec == 1 )
      } else {
        sel <- which(rowSums(otu.pa[, col, drop=FALSE]) >= round(percent*nrow(meta.fac)))
      } 
      
      otu.fac <- otu[, c(col, ncol(otu))]
      
      if ( length(sel) != 0) {
        otu.fac.ids <- as.character(rownames(otu.fac)[sel])
        otu.fac.id[[j]] <- as.character(rownames(otu.fac)[sel])
        otu.fac.ranks <- unique(as.character(factor(otu.fac[["LCA"]][sel])))
        otu.fac.rank[[j]] <- unique(as.character(factor(otu.fac[["LCA"]][sel])))
        otu.fac.ids.num <- length(otu.fac.ids)
        
      } else {
        otu.fac.ids <- ""
        otu.fac.id[[j]] <- ""
        otu.fac.ranks <- ""
        otu.fac.rank[[j]] <- ""
        otu.fac.ids.num <- 0
      } 
      # summarize
      description <- paste(otu.fac.ids.num, "_otus_found_in_", 
                           percent*100, "%_", 
                           fac.levels[j], "_samples;",  sep="")
      otu.fac.list[[j]] <- list(summary=description, otuID=otu.fac.ids,
                                LCA=otu.fac.ranks)
      names(otu.fac.list)[[j]] <- fac.levels[j]
    }
    
    #return(otu.fac.list)
    # find otus across all categories        
    otu.fac.list.num <- length(otu.fac.list)
    
    otu.fac.id[[otu.fac.list.num+1]] <- as.character(factor(
      Reduce(intersect, otu.fac.id)))
    otu.fac.rank[[otu.fac.list.num+1]] <- as.character(factor(
      Reduce(intersect, otu.fac.rank)))
    #return(list(otu.fac.id, otu.fac.list.num))
    # otu ids in all samples
    # percent of sequences in core otus across all categories
    
    if ( length(otu.fac.id[[otu.fac.list.num+1]]) == 0 ||
           unique(unique(otu.fac.id[[otu.fac.list.num+1]]) == "")  ) {
      otu.fac.id.all.num <- 0
      otu.fac.rank.all.num <- 0
      otu.fac.perc <- 0
    } else {
      otu.fac.id.all.num <- length(unique(otu.fac.id[[otu.fac.list.num+1]]))
      otu.fac.rank.all.num <- length(unique(otu.fac.rank[[otu.fac.list.num+1]]))
      otu.sel <- otu[rownames(otu) %in% otu.fac.id[[otu.fac.list.num+1]], ]
      otu.fac.perc <- 100*sum(otu.sel[, -ncol(otu.sel)])/sum(otu[, -dim(otu)[2]])
    } 
    
    description <- paste(otu.fac.id.all.num, "_otus_present_in_core_OTU_of_each_", 
                         meta.factor, sep="")
    
    otu.fac.list[[otu.fac.list.num+1]] <- list(summary=description, 
                                               otuID=otu.fac.id[[otu.fac.list.num+1]], 
                                               LCA=otu.fac.rank[[otu.fac.list.num+1]])
    
    names(otu.fac.list)[[otu.fac.list.num+1]] <- meta.factor
    list.core[[label]] <- otu.fac.list
  }
  
  return(list.core)
}


core.OTU.rank<-function(data, rank="g", drop.unclassified=TRUE,
                        meta, meta.factor="", percent=1){
  
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list. See ?RAM.input.formatting for details")
  }
  
  .valid.data(data=data, is.OTU=TRUE)
  
  list.core <-list()
  labels <- names(data)
  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    otu <- data[[i]]
    if ( is.null(otu) ) { break }
    
    # valid otu & meta
    valid.OTU(otu)
    
    #.valid.meta(otu1=otu, meta=meta)
    samp <- colnames(otu)[1:(ncol(otu)-1)]
    meta <- meta[match(samp, rownames(meta)), ]
    if ( !identical(rownames(meta), samp) ) {
      stop("Error: metadata and otu table do not have same samples")
    }
    
    # check meta category
    .valid.factor(meta, meta.factor)
    fac.levels <- levels(factor(meta[[meta.factor]]))
    
    # get rank name & pattern
    rank.name <- .get.rank.name(rank, plural=TRUE)
    rank_name <- .get.rank.name(rank)
    rank_pat <- .get.rank.pat(rank)
    
    # tax.split
    otu.all <- tax.split(otu1=otu, rank=rank)
    otu <- otu.all[otu.all[[rank_name]] !="", ]
    rm <- paste(.blacklist(), collapse="|")
    if ( drop.unclassified ) {
      rm.row <- which(grepl(rm, otu[[rank_name]], ignore.case=TRUE))
      if ( length(rm.row) == 0 ) {
        rm.row <- NULL
      } else {
        rm.row <- rm.row
      }
    } else {
      rm.row <- NULL
    }
    
    if ( !is.null(rm.row) ) {
      otu <- otu[-rm.row, ]       
    } else {
      otu <- otu
    }
    # presence/absence
    otu.pa <- decostand(otu[, -ncol(otu)], "pa") 
    
    ###check for zero entries?
    otu.fac.list=list()
    otu.fac.id = list()
    otu.fac.rank = list()
    for (j in 1:length(fac.levels)) {
      # each fac levels
      meta.fac <- meta[grep(fac.levels[j], meta[[meta.factor]]), ]
      
      # select only otus present in more than percent of samples 
      # in each category
      col <- which(colnames(otu.pa) %in% rownames(meta.fac))
      if ( length(col) == 1 ) {
        vec <- otu.pa[, col]
        sel <- which(vec == 1 )
      } else {
        sel <- which(rowSums(otu.pa[, col]) >= round(percent*nrow(meta.fac)))
      } 
      
      otu.fac <- otu[, c(col, ncol(otu))]
      otu.all.fac <- otu.all[, c(col, ncol(otu.all))]
      if ( length(sel) != 0) {
        otu.fac.ids <- as.character(rownames(otu.fac)[sel])
        otu.fac.id[[j]] <- as.character(rownames(otu.fac)[sel])
        otu.fac.ranks <- unique(as.character(factor(otu.fac[[rank_name]][sel])))
        otu.fac.rank[[j]] <- unique(as.character(factor(otu.fac[[rank_name]][sel])))
        otu.fac.ids.num <- length(otu.fac.ids)
        otu.fac.ranks.num <- length(otu.fac.ranks)
        
        otu.fac.ids.perc <-  100*sum(otu.fac[c(otu.fac.ids), 
                                             -ncol(otu.fac)])/sum(otu.all.fac[, 
                                                                              -ncol(otu.all.fac)])
        
      } else {
        otu.fac.ids <- ""
        otu.fac.id[[j]] <- ""
        otu.fac.ranks <- ""
        otu.fac.rank[[j]] <- ""
        otu.fac.ids.num <- 0
        otu.fac.ranks.num <- 0
        otu.fac.ids.perc <- 0
      } 
      # summarize
      description <- paste(otu.fac.ids.num, "_otus_found_in_", percent*100, "%_", 
                           fac.levels[j], "_samples; ","assigned_to_", otu.fac.ranks.num, "_", 
                           rank.name, "; ", otu.fac.ids.perc, "%_of_total_sequences_in_", 
                           fac.levels[j], sep="")
      otu.fac.list[[j]]<- list(summary=description, otuID=otu.fac.ids, taxa=otu.fac.ranks)
      names(otu.fac.list)[[j]] <- fac.levels[j]
    }
    
    # find otus across all categories        
    otu.fac.list.num <- length(otu.fac.list)
    
    otu.fac.id[[otu.fac.list.num+1]] <- as.character(factor(Reduce(intersect, otu.fac.id)))
    otu.fac.rank[[otu.fac.list.num+1]] <- as.character(factor(Reduce(intersect, otu.fac.rank)))
    
    # otu ids in all samples
    # percent of sequences in core otus across all categories
    
    if ( length(otu.fac.id[[otu.fac.list.num+1]]) == 0 ||
           unique(unique(otu.fac.id[[otu.fac.list.num+1]]) == "")  ) {
      otu.fac.id.all.num <- 0
      otu.fac.rank.all.num <- 0
      otu.fac.perc <- 0
    } else {
      otu.fac.id.all.num <- length(unique(otu.fac.id[[otu.fac.list.num+1]]))
      otu.fac.rank.all.num <- length(unique(otu.fac.rank[[otu.fac.list.num+1]]))
      otu.sel <- otu[rownames(otu) %in% otu.fac.id[[otu.fac.list.num+1]], ]
      otu.fac.perc <- 100*sum(otu.sel[, -ncol(otu.sel)])/sum(otu.all[, -dim(otu.all)[2]])
    } 
    
    description <- paste(otu.fac.id.all.num, "_otus_present_in_core_OTU_of_each_", 
                         meta.factor, "; ", "assigned_to_", otu.fac.rank.all.num, "_", 
                         rank.name, "; ", otu.fac.perc, "%_of_total_sequences",sep="")
    
    otu.fac.list[[otu.fac.list.num+1]] <- list(summary=description, 
                                               otuID=otu.fac.id[[otu.fac.list.num+1]], 
                                               taxa=otu.fac.rank[[otu.fac.list.num+1]])
    
    names(otu.fac.list)[[otu.fac.list.num+1]] <- meta.factor
    list.core[[label]] <- otu.fac.list
  }
  return(list.core)
}


core.Taxa<-function(data, is.OTU=FALSE, rank="g", 
                    drop.unclassified=TRUE,
                    meta, meta.factor="", percent=1){
  
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list. See ?RAM.input.formatting for details")
  }
  
  .valid.data(data=data, is.OTU=is.OTU)
  
  # valide rank
  if ( is.OTU ) {
    rank.name <- .get.rank.name(rank, plural=TRUE)
    rank_name <-  .get.rank.name(rank)
    rank_pat <- .get.rank.pat(rank)
  } else {
    warning("data are not otu tables, will ignore the rank provided")
    rank_name <- "taxa"
    rank.name <- "taxa"
    rank_pat <- ""
  }
  
  list.core <-list()
  labels <- names(data)
  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    
    if ( is.OTU ) {
      valid.OTU(elem)
      .valid.meta(otu1=elem, meta=meta)
      
      # get rank name & pattern
      rank.name <- .get.rank.name(rank, plural=TRUE)
      rank_name <- .get.rank.name(rank)
      rank_pat <- .get.rank.pat(rank)
      
      # tax.abund
      tax.all <- tax.abund(elem, rank=rank, drop.unclassified=FALSE, 
                           count=TRUE)
    } else {
      elem <- elem[match(rownames(meta), rownames(elem)), ]
      if (identical(rownames(elem), rownames(meta))) {
        tax.all <- elem
      } else {
        stop("Error: metadata and data don't have same subjects")
      }
    }
    
    # check meta category
    if( !length(meta.factor)==1 ) {
      stop("Error: please provide one factor variable in metadata")
    } else {
      .valid.factor(meta, meta.factor)
      fac.levels <- levels(factor(meta[[meta.factor]]))
    }
    
    if ( drop.unclassified ) {
      # remove unclassified columns
      # this selects all columns NOT containing in the blacklist
      remove.pat <- paste(.blacklist(), collapse="|")
      tax <- tax.all[ , !grepl(remove.pat, names(tax.all), 
                               ignore.case=TRUE), drop=FALSE]
    } else {
      tax <- tax.all
    }
    
    ###check for zero entries?
    tax.fac.list=list()
    tax.fac.id = list()
    
    for ( j in 1:length(fac.levels) ) {
      
      # each fac levels
      meta.fac <- meta[grep(fac.levels[j], meta[[meta.factor]]), ]
      tax.fac <- tax[match(rownames(meta.fac), rownames(tax)),]
      tax.fac <- tax.fac[rowSums(tax.fac)>0, colSums(tax.fac)>0]
      # presence or absence
      tax.fac.pa <- vegan::decostand(tax.fac, "pa")
      
      tax.all.fac <- tax.all[match(rownames(meta.fac), rownames(tax.all)),]
      tax.all.fac <- tax.all.fac[rowSums(tax.all.fac)>0, colSums(tax.all.fac)>0]
      # presence or absence
      tax.all.fac.pa <- vegan::decostand(tax.all.fac, "pa")
      
      # select taxa present in more than percent of samples in each category
      sel <- which(colSums(tax.fac.pa) >= round(percent*nrow(tax.fac.pa)))
      
      if (length(sel) != 0) {
        tax.fac.ids <- as.character(colnames(tax.fac)[sel])
        tax.fac.id[[j]] <- as.character(colnames(tax.fac)[sel])
        tax.fac.ids.num <- length(tax.fac.ids)
        tax.fac.ids.perc <-  100*sum(tax.fac[, c(tax.fac.ids)])/sum(tax.all.fac)
        
      } else {
        tax.fac.ids <- ""
        tax.fac.id[[j]] <- ""
        tax.fac.ids.num <- 0
        tax.fac.ids.perc <- 0
      } 
      # summarize
      description <- paste(tax.fac.ids.num, "_", rank.name, "_found_in_", 
                           percent*100, "%_", fac.levels[j], "_samples; ", 
                           tax.fac.ids.perc, "%_of_total_sequences_in_", 
                           fac.levels[j], sep="")
      tax.fac.list[[j]]<- list(summary=description, taxa=tax.fac.ids)
      names(tax.fac.list)[[j]] <- fac.levels[j]
    }
    
    # find taxa across all categories        
    tax.fac.list.num <- length(tax.fac.list)
    tax.fac.id[[tax.fac.list.num+1]] <- as.character(factor(
      Reduce(intersect, tax.fac.id)))
    
    # taxa in all samples
    # percent of sequences in core taxa across all categories
    if( unique(tax.fac.id[[tax.fac.list.num+1]]) == "" || 
          length(unique(tax.fac.id[[tax.fac.list.num+1]])) == 0 ) {
      tax.fac.id.all.num <- 0
      tax.fac.perc <- 0
      tax.fac.id.all <- ""
    } else {
      tax.fac.id.all.num <- length(tax.fac.id[[tax.fac.list.num+1]])
      tax.fac.id.all <- tax.fac.id[[tax.fac.list.num+1]]
      tax.sel <- tax[, which(colnames(tax) %in% tax.fac.id.all)]
      tax.fac.perc <- 100*sum(tax.sel)/sum(tax.all)
    } 
    
    description <- paste(tax.fac.id.all.num, "_", rank.name, "_found_in_", 
                         percent*100, "%_of_samples_at_each_", meta.factor,
                         "; ", tax.fac.perc, "%_of_total_sequences", 
                         sep="")
    tax.fac.list[[tax.fac.list.num+1]] <- list(summary=description, taxa=tax.fac.id.all)
    names(tax.fac.list)[[tax.fac.list.num+1]] <- meta.factor
    list.core[[label]] <- tax.fac.list
  }
  return(list.core)
  
}


group.OTU <- function(otu, rank="g", otuIDs="", meta, 
                      meta.factor="", boxplot=TRUE,  
                      main="", file=NULL, ext=NULL, 
                      height=8, width=16) {
  
  if(!length(meta.factor)==1) {
    stop("Error: please provide one factor variable in metadata")
  } else {
    fac.levels <- levels(factor(meta[[meta.factor]]))
  }
  
  # get the groups
  if(!length(meta.factor)==1) {
    stop ("Error: please provide one factor")
  }
  
  if(length(otuIDs)==0) {
    stop ("Error: please provide list of otuIDs to plot")
  }
  
  # validate inputs
  valid.OTU(otu)
  .valid.meta(otu1=otu, meta=meta)
  
  # tax.split
  if ( is.null(rank) ) {
    otu <- LCA.OTU(otu, strip.format=FALSE, drop=TRUE)
    rank_name <- "LCA"
  } else {
    otu <- tax.split(otu1=otu, rank=rank)
    rank.name <- .get.rank.name(rank, plural=TRUE)
    rank_name <- .get.rank.name(rank)
    rank_pat <- .get.rank.pat(rank)
  }
  otu.p <- cbind(decostand(otu[, -ncol(otu)], "total", 
                           MARGIN=2), otu[[rank_name]])
  colnames(otu.p)[ncol(otu.p)] <- rank_name
  otu.p <- otu.p[otu.p[[rank_name]] !="", ]
  
  if( length(which(rownames(otu.p) %in% otuIDs))==0) {
    stop ("Error: no selected otuIDs present in the otu table")
  }
  
  otuIDs <- unique(otuIDs)
  otu.p.sel <- otu.p[which(rownames(otu.p) %in% c(otuIDs)), ]
  otu.sel <- otu[which(rownames(otu) %in% c(otuIDs)), ]
  
  otus.ex<-vector()
  for (i in otuIDs) {
    if(!(i %in% rownames(otu.p))) {
      otus.ex<-unique(c(otus.ex, i))
    }
  } 
  if (length(otus.ex) != 0) {
    print(paste("The following selected otus in your list were not in the otu table: ", paste(otus.ex, collapse=", "), sep=""))
  }
  
  # combine otu with taxonomy
  rownames(otu.p.sel) <- paste(paste0("OTU_ID|",rownames(otu.p.sel)), 
                               otu.p.sel[[rank_name]], sep=": ")
  rownames(otu.sel) <- paste(paste0("OTU_ID|",rownames(otu.sel)), 
                             otu.sel[[rank_name]], sep=": ")
  
  # transpose OTU table
  otu.p.sel.t <- as.data.frame(t(otu.p.sel[, -ncol(otu.p.sel)]))
  otu.p.sel.tax <- otu.p.sel[[rank_name]]
  
  otu.sel.t <- as.data.frame(t(otu.sel[, -ncol(otu.sel)]))
  otu.sel.tax <- otu.sel[[rank_name]]
  
  if(identical(rownames(otu.p.sel.t), rownames(meta))) {
    otu.p.sel.fac <- cbind(otu.p.sel.t, meta[[meta.factor]])
    names(otu.p.sel.fac)[ncol(otu.p.sel.fac)] <- meta.factor
    otu.sel.fac <- cbind(otu.sel.t, meta[[meta.factor]])
    names(otu.sel.fac)[ncol(otu.sel.fac)] <- meta.factor
  } else {
    stop("samples not identical in otu and metadata")
  }
  
  # melt by sample
  #if (!require("reshape2")) {
  # stop("package 'reshape2' is required to use this function")
  # }
  
  otu.sel.fac.m <- reshape2::melt(cbind(otu.sel.fac, Sample=rownames(otu.sel.fac)), 
                                  tax=otu.sel.fac$tax, variable.name="OTU", 
                                  value.name="Count")  
  names(otu.sel.fac.m)[ncol(otu.sel.fac.m)-1] <- "OTU"
  names(otu.sel.fac.m)[ncol(otu.sel.fac.m)] <- "Count"
  
  # boxplot or barplot
  if(boxplot) {
    # boxplot
    otu.p.sel.fac.m <- reshape2::melt(cbind(otu.p.sel.fac, 
                                            Sample=rownames(otu.p.sel.fac)), 
                                      variable.name="OTU", value.name="RA")
    names(otu.p.sel.fac.m)[ncol(otu.p.sel.fac.m)-1] <- "OTU"
    names(otu.p.sel.fac.m)[ncol(otu.p.sel.fac.m)] <- "RA"
    
    otu.p.sel.fac.m <- cbind(otu.p.sel.fac.m, do.call(rbind,
                                                    strsplit(x=as.character(otu.p.sel.fac.m$OTU),split=":")))
    names(otu.p.sel.fac.m)[ncol(otu.p.sel.fac.m)-1] <- "otu"
    names(otu.p.sel.fac.m)[ncol(otu.p.sel.fac.m)] <- "tax" 
    
    otu.p.sel.fac.m <- otu.p.sel.fac.m[order(otu.p.sel.fac.m$tax, 
                                             otu.p.sel.fac.m$otu), ]    
  } else {
    # barplot
    otu.sel.fac.agg = stats::aggregate(otu.sel.t, 
                                       by=list(meta[[meta.factor]]), FUN=sum)
    otu.sel.fac.agg <- cbind(otu.sel.fac.agg[,1], 
                            vegan::decostand(otu.sel.fac.agg[, -1], 
                                             MARGIN=2, "total"))
    names(otu.sel.fac.agg)[1] <- meta.factor
    
    # label counts as "RA" to be consistent with previous df
    otu.sel.fac.agg.m <- reshape2::melt(otu.sel.fac.agg, 
                                        is.vars=c(meta.factor), value.name="RA", 
                                        variable.name="OTU")
    names(otu.sel.fac.agg.m)[ncol(otu.sel.fac.agg.m)-1] <- "OTU"
    names(otu.sel.fac.agg.m)[ncol(otu.sel.fac.agg.m)] <- "RA"
    otu.sel.fac.agg.m <- cbind(otu.sel.fac.agg.m, 
                               do.call(rbind,strsplit(x=as.character(otu.sel.fac.agg.m$OTU), 
                                                      split=":")))
    names(otu.sel.fac.agg.m)[ncol(otu.sel.fac.agg.m)-1] <- "otu"
    names(otu.sel.fac.agg.m)[ncol(otu.sel.fac.agg.m)] <- "tax"
    
    otu.sel.fac.agg.m <- otu.sel.fac.agg.m[order(otu.sel.fac.agg.m$tax, 
                                                 otu.sel.fac.agg.m$otu), ]
    otu.p.sel.fac.m <- otu.sel.fac.agg.m
    #return(otu.p.sel.fac.m)
  }
  
  #calculate total count of each Taxon
  total= sapply(levels(otu.sel.fac.m$OTU), 
                function(x){sum(otu.sel.fac.m$Count[otu.sel.fac.m$OTU==x])}, 
                USE.NAMES=F)
  total.count <- data.frame(OTU=levels(otu.sel.fac.m$OTU),total)
  # split otu&tax
  total.count<-cbind(total.count, 
                     do.call(rbind,strsplit(x=as.character(total.count$OTU),
                                            split=":")))
  names(total.count)[ncol(total.count)-1] <- "otu"
  names(total.count)[ncol(total.count)] <- "tax"
  
  #total.count <- total.count[order(total.count$tax, 
  #total.count$OTU, total.count$total, decreasing=TRUE), ]
  
  #return(list(otu.sel.fac.m, otu.p.sel.fac.m, total.count))
  
  # reorder levels based on order of appearance in total.count, 
  otu.p.sel.fac.m$OTU <- factor(as.character(otu.p.sel.fac.m$OTU), 
                                levels=total.count$OTU, ordered=TRUE)
  
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  # the distribution of each taxon in each category of meta.factor
  if (!requireNamespace("scales")) {
    stop("package 'scales' is required to use this function")
  }

  if (!requireNamespace("RColorBrewer")) {
    stop("package 'RColorBrewer' is required to use this function")
  }

  if (!requireNamespace("grid")) {
    stop("package 'grid' is required to use this function")
  }
  
  ylab <- "Relative Abundance (%)"
  
  if(boxplot) {
    
    p <- ggplot(otu.p.sel.fac.m, aes_string(x="OTU", y="RA", 
                                                     fill=meta.factor)) + 
      geom_boxplot() 
  } else {
    p <- ggplot(otu.p.sel.fac.m, aes_string(x="OTU", y="RA", 
                                                     fill=meta.factor)) + 
      geom_bar(colour="white",stat="identity",
               position="fill")
  }
  
  cols.needed <- length(levels(factor(meta[[meta.factor]])))
  if (cols.needed <= 12 ) {
    p <- p + scale_fill_brewer(type="div", palette="Set3") 
  } else {
    p <- p + scale_fill_manual(values=RAM.pal(cols.needed))  
  } 
  
  p <- p + scale_x_discrete(labels = paste(total.count$otu, "\n", 
                                           total.count$tax,":",total.count$total, sep="")) + 
    theme(legend.key.size=unit(6,"mm"),
          legend.text = element_text(size=10), 
          axis.text.y=element_text(size=7),
          legend.position="right") + 
    xlab("OTU") +
    ylab(ylab) + 
    ggtitle(main) +
    coord_flip()
  
  p <- p + scale_y_continuous(labels=scales::percent_format())
  
  
  .valid.plot.settings(file, ext)
  save <- !is.null(file)
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
  
}

# venn for classified taxa
group.venn <- function(vectors, cat.cex=1.5, cex=1,
                       cat.pos=NULL, cat.dist=NULL, 
                       label=TRUE, lab.cex=1, 
                       lab.col= "black", fill=NULL, 
                       file=NULL, ext=NULL, width=8, height=8) {
  
  save <- !is.null(file)
  if (save) { .get.dev(file, ext, height=height, width=width) }
  
  if ( !requireNamespace("VennDiagram") ) {
    stop("package 'VennDiagram' is required for this function")
  }
  if ( !requireNamespace("RColorBrewer") ) {
    stop("package 'RColorBrewer' is required for this function")
  }

  if (!requireNamespace("grid")) {
    stop("package 'grid' is required to use this function")
  }

  
  # Generate plot
  # number of vectors to plot
  len <- length(vectors)
  if ( is.null(fill) ) {
    if ( len == 2 ) {
      fill = c("lightpink", "lightblue")
    } else {
      fill = RColorBrewer::brewer.pal(len, "Pastel1")
    }
  } else {
    if ( length(fill) == len ) {
      fill = fill
    } else if ( length(fill) > len ) {
      warning(paste("more colors being provided than required, will ignore ", length(fill)-len, " colors", sep=""))
      fill = fill[1:len]
    } else {
      warning("not enough colors being provided, will use default")
      if ( len == 2 ) {
        fill = c("lightpink", "lightblue")
      } else {
        fill = RColorBrewer::brewer.pal(len, "Pastel1")
      }
    }
  }
  
  if ( len > 2 && label )  {
    warning("currently only support 2 groups to have actual item labels; will only use numbers")
  } else if ( len > 5 || len < 2 ) {
    stop("please provide 2 to 5 vectors")
  }
  
  alpha = rep(0.5, len)
  
  if ( !is.null(cat.pos) && !is.null(cat.dist) ) { 
    v <- VennDiagram::venn.diagram(vectors, fill = fill, alpha = alpha, 
                                   cat.dist=cat.dist, cat.pos=cat.pos, 
                                   cat.fontface = "bold", cat.cex = cat.cex, cex=cex, 
                                   filename=NULL)
  } else if ( !is.null(cat.pos) && is.null(cat.dist) ) { 
    v <- VennDiagram::venn.diagram(vectors, fill = fill, alpha = alpha, 
                                   cat.pos=cat.pos, cat.fontface = "bold", 
                                   cat.cex = cat.cex, cex=cex, 
                                   filename=NULL)
  } else if ( is.null(cat.pos) && !is.null(cat.dist) ) { 
    v <- VennDiagram::venn.diagram(vectors, fill = fill, alpha = alpha, 
                                   cat.fontface = "bold", cat.dist= cat.dist, 
                                   cat.cex = cat.cex, cex=cex, 
                                   filename=NULL)
  } else {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, alpha = alpha, 
                                   cat.fontface = "bold", cat.cex = cat.cex, cex=cex, 
                                   filename=NULL)
  }
  
  if ( len > 2 && len <=5 ) {
    grid::grid.newpage()
    grid::grid.draw(v)
  }
  
  if ( len==2 ) {
    if ( !label ) {
      grid::grid.newpage()
      grid::grid.draw(v)
    } else {
      name <- lapply(v,  names)
      lapply(v, function(i) i$label)
      
      # plot with labels
      # first find out whether the vectors got rearranged
      v.labels <- lapply(v, function(i) i$label)
      v.lab <- vector()
      for ( i in 1:length(v.labels) ) {
        if ( length(v.labels[[i]] %in% names(vectors)) !=0 && 
               isTRUE(v.labels[[i]] %in% names(vectors))  ) {   
          v.lab <- c(v.lab, v.labels[[i]])
        }
      }
      v1 <- vectors[[v.lab[1]]]
      v2 <- vectors[[v.lab[2]]]
      v[[5]]$label  <- paste(c(v[[5]]$label,setdiff(v1, v2)), collapse="\n")
      v[[5]]$gp$cex <- lab.cex
      v[[5]]$gp$col <- lab.col
      # in baa only
      v[[6]]$label <- paste(c(v[[6]]$label,setdiff(v2, v1))  , collapse="\n") 
      v[[6]]$gp$cex <- lab.cex
      v[[6]]$gp$col <- lab.col
      # intesection
      v[[7]]$label <- paste(c(v[[7]]$label,intersect(v1, v2)), collapse="\n")
      v[[7]]$gp$cex <- lab.cex
      # plot with labels
      v[[7]]$gp$col <- lab.col
      
      # plot with labels
      grid::grid.newpage()
      grid::grid.draw(v)
    }
  } 
  
  if (save) { dev.off() }
  invisible()
}
