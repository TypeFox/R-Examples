META.clust <- function(meta, group=4, data.trans=NULL, dist=NULL, 
                       clust=NULL, type=NULL, main="", 
                       file=NULL, ext=NULL, width=8, height=8) {
  
  save <- !is.null(file)
  if (save) { .get.dev(file, ext, height=height, width=width) }
  
  suppressWarnings(meta.new <- filter.META(meta))
  
  # check wither mixed or contain factor/character variables
  cl <- vector()
  fac <- vector()
  num <- vector()
  for ( i in 1:ncol(meta.new) ) {
    if ( class(meta.new[,i]) == "factor" || 
           class(meta.new[,i]) == "character" ) {
      fac <- unique(c(fac, i))
    } 
    if ( class(meta.new[,i]) == "numeric" ) {
      num <- unique(c(num, i))
    }
    cl <- unique(c(cl, class(meta.new[,i]) ))
  } 
  
  # data transforation for numeric variables
  if ( is.null(data.trans) ) {
    meta.new.tr <- meta.new
  } else {
    if ( length(num) != 0 && length(fac) !=0 ) {
      meta.new.num <- vegan::decostand(meta.new[, num], data.trans)
      meta.new.tr <- cbind(meta.new[, fac], meta.new.num)
      names(meta.new.tr) <- c(names(meta.new)[fac], 
                              names(meta.new)[num])
    } else {
      warning("no numeric variables to be transformed")
      meta.new.tr <- meta.new
    }       
  }
  
  if ( any(cl %in% c("factor", "character")) ) {
    # if (!require("FD")) {
    #     stop("package 'FD' is required to use this function: 
    #     try 'install.packages('FD')'.")
    #  }
    dis <- FD::gowdis(meta.new.tr)
  } else {
    dis <- vegan::vegdist(meta.new.tr, method=dist)
  }
  
  # cluster methodes
  if ( is.null(clust) ) {
    hc <- stats::hclust(dis)
  } else {
    hc <- stats::hclust(dis, clust)
  }
  # plot(hc, hang=-1)
  # background color
  
  #if (!require("RColorBrewer")) {
  #         stop("package 'RColorBrewer' is required to use this function: 
  #         try  'install.packages('RColorBrewer')'.")
  #}
  
  
  if ( length(group) != 1L ) {
    warning(" group must be a digit or a metadata variable; 
            will cut tree in 4 groups as default")
    mem <- stats::cutree(hc, 4)
  } else {   
    if ( is.numeric(group) ) {
      # cut dendrogram in 4 clusters
      mem <- stats::cutree(hc, group)
    } else {
      if ( any(group %in% names(meta.new.tr)) ) {
        # group subjects by metadata factor
        mem <- factor(as.numeric(factor(meta.new.tr[[group]])))
        attr(mem, "names") <- rownames(meta.new.tr)
      } else {
        warning(" group is NOT a metadata variable, 
                will cut tree in 4 groups as default")
        mem <- stats::cutree(hc, 4)
      }
    }
    }
  
  .pCLUST(hc, type=type, mem=mem)
  
  if (save) { dev.off() }
  invisible()
}


data.clust <- function(data, is.OTU=TRUE, meta, rank=NULL, top=NULL, 
                       mode="number", group=4, data.trans=NULL, 
                       dist=NULL, clust=NULL, type=NULL, main=NULL, 
                       file=NULL, ext=NULL, width=8, height=8) {
  
  save <- !is.null(file)
  if (save) { .get.dev(file, ext, height=height, width=width) }
  
  # make sure only one level of data.new
  if ( !is.null(rank) && length(rank) != 1L ) {
    warning("rank should be one length, multiple ranks were      
            provided, will only process the first rank")
    rank <- rank[1]
    .valid.rank(rank)
  } else {
    rank <- rank
  }
  
  # transform the data
  elem <- data
  data.new <- data.revamp(data=list(data=elem), is.OTU=is.OTU, 
                          ranks=rank, stand.method=data.trans, 
                          top=top, mode=mode)
  
  # filter METAdata, exclude variable with only 1 level, with
  # missing data, and non numeric&&factor/charactor (NNF)
  suppressWarnings(meta.new <-filter.META(meta))
  
  # extract the first in the data.new list.
  data.new <- data.new[[1]]
  if ( !is.null(dist) ) {
    dis <- vegan::vegdist(data.new, method=dist)
  } else {
    dis <- vegan::vegdist(data.new)
  }
  
  # cluster methodes
  if ( is.null(clust) ) {
    hc <- stats::hclust(dis)
  } else {
    hc <- stats::hclust(dis, clust)
  }
  # plot(hc, hang=-1)
  # background color
  
  #if (!require("RColorBrewer")) {
  #         stop("package 'RColorBrewer' is required to use this function: 
  #         try  'install.packages('RColorBrewer')'.")
  #}
  
  op = graphics::par(bg = "#E8DDCB")
  
  if ( length(group) != 1L ) {
    warning(" group must be a digit or a metadata variable; 
            will cut tree in 4 groups as default")
    mem <- stats::cutree(hc, 4)
  } else {   
    if ( is.numeric(group) ) {
      # cut dendrogram in 4 clusters
      mem <- stats::cutree(hc, group)
    } else {
      if ( any(group %in% names(meta.new)) ) {
        # group subjects by metadata factor
        mem <- factor(as.numeric(factor(meta.new[[group]])))
        attr(mem, "names") <- rownames(meta.new)
      } else {
        warning(" group is NOT a metadata variable, 
             will cut tree in 4 groups as default")
        mem <- stats::cutree(hc, 4)
      }
    }
  }
  
  .pCLUST(hc, type=type, mem=mem)
  
  if (save) { dev.off() }
  invisible()
}

.type.phylo <-  function() {
  return(c("phylogram", "cladogram", "fan", "unrooted", "radial"))
}  
   
.type.dendro <-  function() {
  return(c("triangle", "rectangle"))
}
   

# function to plot clusters
.pCLUST <- function(hc, type, mem, main="") {
  #if (!require("RColorBrewer")) {
  #    stop("package 'RColorBrewer' is required to use 
  #this function: try  'install.packages('RColorBrewer')'.")
  # }
  
  num.mem <- length(levels(factor(mem)))
  if ( num.mem <= 2 ) {
    leaf.col <- c("#E41A1C", "#377EB8")
  } else if ( num.mem > 2 && num.mem <= 9 ) {
    leaf.col <- RColorBrewer::brewer.pal(num.mem, "Set1")
  } else if ( num.mem > 9 && num.mem <= 12 ) {
    leaf.col <- RColorBrewer::brewer.pal(num.mem, "Set3")
  } else if ( num.mem > 12 && num.mem <=20 ) {
    leaf.col <- .ram.pal(num.mem)
  } else {
    leaf.col <- grDevices::rainbow(num.mem) 
  } 
  
  # set background color
  op = graphics::par(bg = "#E8DDCB")
  
  if ( any(type %in% .type.phylo() ) ) {
    #if (!require("ape")) {
    #    stop("package 'ape' is required to use this function: 
    #    try  'install.packages('ape')'.")
    #}
    fhc <- ape::as.phylo(hc)
    plot(fhc, type = type, edge.width = runif(20, 0.5, 3), 
         tip.color = leaf.col[mem], col="red", cex = 1, 
         label.offset = 0.01, main=main)
  }
  
  # function to get color labels
  .Col.group <- function(x) {
    if ( stats::is.leaf(x) ) {
      att <- attributes(x)
      LAB_col <- leaf.col[mem[which(names(mem) == att$label)]]
      attr(x, "nodePar") <- c(att$nodePar, lab.col = LAB_col)
    }
    x
  }
  
  if ( any(type %in% .type.dendro()) || is.null(type) ) {
    # using dendrogram objects
    dhc <- stats::as.dendrogram(hc)
    # load code of A2R function
    Dendro = stats::dendrapply(dhc, .Col.group)
    plot(Dendro, type = type, main=main)
  }
}
