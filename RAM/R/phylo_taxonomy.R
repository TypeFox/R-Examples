.sumcol = function(X){
  ret = X
  X = as.matrix(X)
  if(dim(X)[1]>1)
    ret = colSums(X)
  return(ret)
}

.as.phylo.formula_SL<-function (x, data = parent.frame(), ...){
  # credit goes to Liam Revell, https://stat.ethz.ch/pipermail/r-sig-phylo/2013-August/003017.html 
  if ( !requireNamespace('phytools') ) {
    stop("function read.newick in package phytools is required for this function")
  }

  if ( !requireNamespace('ape') ) {
    stop("functions in package ape is required for this function")
  }
  
  err <- "use formula \"~column1/column2/.../columnN\"."
  if (length(x) != 2) 
    stop(err)
  if (x[[1]] != "~") 
    stop(err)
  f <- x[[2]] # left of the formula
  taxo <- list()
  while (length(f) == 3) {
    if (f[[1]] != "/") 
      stop(err)
    if (!is.factor(data[[deparse(f[[3]])]])) 
      stop(paste("Variable", deparse(f[[3]]), "must be a factor."))
    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
    if (length(f) > 1) 
      f <- f[[2]]
  }
  if (!is.factor(data[[deparse(f)]])) 
    stop(paste("Variable", deparse(f), "must be a factor."))
  taxo[[deparse(f)]] <- data[[deparse(f)]]
  taxo.data <- as.data.frame(taxo)
  leaves.names <- as.character(taxo.data[, 1])
  taxo.data[, 1] <- 1:nrow(taxo.data)
  f.rec <- function(subtaxo) {
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[, u])
    if (u == 1) {
      if (length(levels) != nrow(subtaxo)) 
        warning("Error, leaves names are not unique.")
      return(as.character(subtaxo[, 1]))
    }
    t <- character(length(levels))
    for (l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)])
      t[l] <- paste("(", paste(x, collapse = ","), ")", sep = "")
    }
    return(t)
  }
  string <- paste("(", paste(f.rec(taxo.data), collapse = ","),");", sep = "")
  phy <- phytools::read.newick(text = string) ## so that singles will be read without error
  phy$edge.length <- rep(1,nrow(phy$edge))
  phy <- ape::collapse.singles(phy)
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  return(phy)
}

phylo_taxonomy <- function(otu, rank="order", rank.sep="; ", meta, factors, plot.type="phylogram", 
                           edge.width=1, cex=0.7, font = 1, x.lim = NULL, tip.offset=0, tip.cex=0.5, 
                           thermo=FALSE, thermo.horiz=TRUE, thermo.width=0.5, thermo.height=1, 
                           node.frame="r", node.bg="white", node.col="black", node.width=0.5, 
                           node.height=0.6, node.cex=0.6, node.font=1) {
  #if (!require("ape")) {
  #  stop("R packages RAM and ape are required for this function")
  #}
  valid.OTU(otu)
  otu.tax <- gsub("Incertae sedis", "Incertae_sedis", x=otu$taxonomy, ignore.case=TRUE)
  otu$taxonomy <- otu.tax
  data <- tax.fill(otu, downstream=TRUE)
  data <- data[order(data$taxonomy), ]
  rank <- rank
  if (.get.rank(.get.rank.ind(rank)) == "species") {
    last.rank=TRUE
  } else {
    last.rank=FALSE
  }
  
  all.tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  if ( !last.rank ) {
    rank.split <- paste0(rank.sep, .get.rank.pat(.get.rank(.get.rank.ind(rank)+1)))
    tax.classes <- c("taxonomy", "exclude")
    data.split <- col.splitup(data, col="taxonomy", sep=rank.split,  drop=TRUE, names=tax.classes)
    data.split <- data.split[, -grep("exclude", names(data.split))]
  } else {
    data.split <- data
  }
  #head(data.split)
  
  data.ag <- stats::aggregate(data.split[, 1:(ncol(data.split)-1)], by=list(factor(data.split$taxonomy)), FUN=.sumcol)
  names(data.ag)[1] <- "taxonomy"
  data.ag <- data.ag[, c(names(data.ag)[2:ncol(data.ag)], "taxonomy")]
  data.ag$taxonomy <- as.character(data.ag$taxonomy)
  
  ## meta
  if ( !identical(colnames(data.ag)[1:(ncol(data.ag)-1)], rownames(meta) )) {
    stop("not same samples or not the same order in otu table and metadata")
  }
  
  pie.list <- list()
  fac.lev <- list()
  for ( i in 1:length(factors) ) {
    fac <- factors[i]
    if ( length(which(is.na(meta[, factors[i]]))) != 0 ) {
      warning(paste("missing data in ", factors[i], ", the plot may not be completed; consider remove this variable!", sep=""))
    } 
    if ( is.numeric(meta[, fac]) ) {
      stop("factors must be category data")
    } else {
      meta[, fac] <- as.character(meta[, fac])
    }
    fac.lev[[factors[i]]] <- length(levels(factor(meta[, fac])))
    temp <- RAM::transpose.OTU(data.ag)
    meta.ag <- stats::aggregate(temp, by=list(meta[, fac]), FUN=.sumcol)
    rownames(meta.ag) <- meta.ag[,1]
    meta.ag <- meta.ag[,-1]
    meta.ag.t <- as.data.frame(t(meta.ag))
    pie <- as.matrix(vegan::decostand(meta.ag.t, "total", 1))
    pie <- cbind(pie, rowSums(meta.ag.t))
    colnames(pie)[ncol(pie)] <- "Seq_count"
    pie.list[[factors[i]]] <- pie
  } 
  
  
  sel <- which(all.tax.classes == .get.rank(.get.rank.ind(rank)))
  sel.tax.classes <- all.tax.classes[1:sel]
  data.ag.split <- col.splitup(data.ag, col="taxonomy", sep=rank.sep,  drop=TRUE, names=sel.tax.classes)
  #head(data.ag.split)
  
  
  t1 <- .as.phylo.formula_SL(as.formula(paste0("~", paste(sel.tax.classes, collapse="/"))), data=data.ag.split)
  t1$tip.label
  t1$node.label
  t1$edge
  t1$tip.label <- paste(t1$tip.label, pie[,ncol(pie)], sep=":")
  
  # plot the tree
  plot(t1, type=plot.type, edge.width=edge.width, cex=cex, font = font, label.offset = tip.offset, x.lim = x.lim, no.margin = FALSE)
  
  if (!requireNamespace("RColorBrewer")) {
    stop("R package RColorBrewer is required for this function")
  }
  
  col.list <-list()
  col.all <- RAM.pal(50)
  remove <- vector()
  remove[1] <- 1
  for ( i in 1:length(factors) ) {
    #col.list[[factors[i]]] <- sample(brewer.pal(12, "Set3"), fac.lev[[i]])
    #if ( fac.lev[[i]] <= 12 ) {
    #  col.list[[factors[i]]] <- brewer.pal(fac.lev[[i]], "Set3")
    #} else {   
    remove <- c(remove, (remove[i] + fac.lev[[i]]))   
    col.list[[factors[i]]] <- col.all[remove[i]:(remove[i+1]-1)]
    #}
    names(col.list[[factors[i]]]) <- levels(factor(meta[, factors[i]]))
  }
  
  
  #tiplabels(pch = 21, bg = gray(pie[,3]/sum(pie[,3])), cex = 2, adj = 1.4)
  #tiplabels(pch = 19, col = c("yellow", "red", "blue"), adj = 2.5, cex = 2)
  for ( i in 1:length(factors) ) {
    if ( !thermo ) {
      ape::tiplabels(piecol=col.list[[i]], pie=pie.list[[i]][,1:(ncol(pie.list[[i]])-1)], bg = "white", 
                     cex = tip.cex, adj = 0+ (i-1)*0.15,  label.offset = tip.offset, no.margin=TRUE)
    } else {
      #     tiplabels(piecol=col.list[[i]], thermo = pie.list[[i]][,1:(ncol(pie.list[[i]])-1)], 
      #     bg = "white", cex = tip.cex, adj = 0+ (i-1)*0.5, width=thermo.width, height=thermo.height, 
      #     horiz=thermo.horiz, label.offset = tip.offset, x.lim = x.lim, no.margin=TRUE)
      ape::tiplabels(piecol=col.list[[i]], thermo = pie.list[[i]][,1:(ncol(pie.list[[i]])-1)], bg = "white", 
                     cex = tip.cex, adj = 1, width=thermo.width, height=thermo.height, horiz=thermo.horiz, 
                     label.offset = tip.offset, x.lim = x.lim, no.margin=TRUE)
    }
  }
  
  #nodelabels()
  
  for (i in 1:length(t1$tip.label)) { 
    #tip.v<-c(tip.v, i); 
    tax.tip <- unlist(strsplit(t1$tip.label[i], ":"))[1]
    if (tax.tip != "") { 
      tax <- data.ag[grep(tax.tip, data.ag$taxonomy), "taxonomy"]
      tax.split <- unlist(strsplit(tax, split=rank.sep))
      if ( !last.rank ) {
        tax.split <- rev(tax.split[1:(length(tax.split)-1)])
      } else {
        tax.split <- rev(tax.split[1:(length(tax.split))])
      }
      a <- phangorn::Ancestors(t1, i, type=c("all")); 
      dis <- ape::dist.nodes(t1)[i, c(i, a)]
      #if (length(a) == length(tax.split)) {
      for ( j in 1:length(a) ) {
        #if ( (dis[j+1]-dis[j]) >1 ) {
        #  inter <- dis[j+1] - dis[j] -1
        #  inter.tax <- tax.split[(dis[j]+1):(dis[j+1]-1)]
        #edgelabels ( tax.split[dis[j+1]]
        
        #print(paste(tax.split[dis[j+1]], a[j], sep="   "))
        ape::nodelabels(tax.split[dis[j+1]], a[j], frame=node.frame, bg=node.bg, col=node.col, 
                        width=node.width, height=node.height, cex=node.cex, font=node.font)        
      }
    }
  }
  #return(pie.list)
  return(list(col.list=col.list, tree=t1, pie.list=pie.list))
}



phylog_taxonomy <- function(otu, rank="order", rank.sep="; ", meta, factors=NULL, sel.taxon=NULL, sel.rank=NULL, root="root", cleaves=1, cnodes=1, clabel.leaves = 0.5, clabel.nodes = 0.5, f.phylog = 0.5, sub = TRUE, csub = 1.25, possub = "bottomleft", draw.box = TRUE) {

  if (!requireNamespace("ade4") ) {
  #if (!require("ade4") ) {
    stop("R packages ade4 is required for this function")
  }
 
  valid.OTU(otu)
  ## meta
  if ( !identical(colnames(otu[, -ncol(otu)]), rownames(meta) )) {
    stop("not same samples or not the same order in otu table and metadata")
  }


  # extract subdata for each level of factor
  data.list <- list()
  meta.list <- list()
  fac.lev <- list()

  if ( is.null(factors) ) {
    data.list[["ALL"]] <- otu
    meta.list[["ALL"]] <- meta
  } else {
    for ( i in 1:length(factors) ) {
      fac <- factors[i]
      if ( is.numeric(meta[, fac]) ) {
        stop("factors must be category data")
      } else {
         meta[, fac] <- as.character(meta[, fac])
      }
      fac.lev[[factors[i]]] <- length(levels(factor(meta[, fac])))
      for (j in levels(factor(meta[,fac]))) {
        name <- paste(fac, j, sep=": ")
        meta1 <- meta[meta[fac]==j, ]
        otu1 <- otu[, match(c(rownames(meta1), "taxonomy"), colnames(otu))]
        if ( is.null(sel.taxon) ) {
          meta.list[[name]] <- meta1
          data.list[[name]] <- otu1
        } else {
          if ( is.null(sel.rank) ) {
             stop("what's the taxonomic rank of the selected taxon? eg., 'order' or 'o' ")
          } else {
             pat <- .get.rank.pat(sel.rank)
             pat1 <- paste(pat, sel.taxon, ";", sep="")
             otu2 <- otu1[grep(pat1, otu1$taxonomy), ]
             # exclude samples don't have the taxon
             otu2 <- otu2[rowSums(otu2[, -ncol(otu2)])>0, colSums(otu2[, -ncol(otu2)])>0]
             if ( ncol(otu2) ==0 || nrow(otu2) == 0 ) {
                stop(paste("no ", sel.taxon, " was found in samples from ", name, sep=""))
             } else {
               meta2 <- meta1[match(names(otu2)[1:(ncol(otu2)-1)], rownames(meta1)), ]
             }
             meta.list[[name]] <- meta2
             data.list[[name]] <- otu2
           }
         }
       }
     }
  }
  
     
  # deal with each subdata
  data.ag.split.list <- list()
  tax.phy.list <- list()
  for ( k in 1:length(data.list) ) {
    data<-data.list[[k]]
    data.name <- names(data.list)[k]
    #return(data)
    data.tax<-gsub("Incertae sedis", "Incertae_sedis",x=data$taxonomy, ignore.case=TRUE)
    data$taxonomy <- data.tax
    data <- tax.fill(data, downstream=TRUE)
   
    #return(data)
    data <- data[order(data$taxonomy), ]
    rank <- rank
    if (.get.rank(.get.rank.ind(rank)) == "species") {
      last.rank=TRUE
    } else {
      last.rank=FALSE
    }

    if ( !last.rank ) {
      rank.split <- paste0(rank.sep, .get.rank.pat(.get.rank(.get.rank.ind(rank)+1)))
      tax.classes <- c("taxonomy", "exclude")
      data.split <- col.splitup(data, col="taxonomy", sep=rank.split,  drop=TRUE, names=tax.classes)
     data.split <- data.split[, -grep("exclude", names(data.split))]
    } else {
      data.split <- data
    }
    #head(data.split)

    data.ag <- aggregate(data.split[, 1:(ncol(data.split)-1)], by=list(factor(data.split$taxonomy)), FUN=.sumcol)
    names(data.ag)[1] <- "taxonomy"
    data.ag <- data.ag[, c(names(data.ag)[2:ncol(data.ag)], "taxonomy")]
    data.ag$taxonomy <- as.character(data.ag$taxonomy)
    #return(data.ag)
    # split 
    all.tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    sel <- which(all.tax.classes == .get.rank(.get.rank.ind(rank)))
    sel.tax.classes <- all.tax.classes[1:sel]
    data.ag.split <- col.splitup(data.ag, col="taxonomy", sep=rank.sep,  drop=TRUE, names=sel.tax.classes)
    dup <- which(duplicated(data.ag.split[,ncol(data.ag.split)]) )  
    if ( length(dup) != 0L) {
      for (w in dup) {
        data.ag.split[, ncol(data.ag.split)] <- as.character(data.ag.split[, ncol(data.ag.split)])
        data.ag.split[w, ncol(data.ag.split)] <- paste(data.ag.split[w, ncol(data.ag.split)], w, sep="_")
      }
    }
    data.ag.split[,ncol(data.ag.split)] <- paste(data.ag.split[, ncol(data.ag.split)], "|", rowSums(data.ag.split[, 1:(ncol(data.ag)-1)]), sep="")
   
    #return(data.ag.split)
 
    data.ag.split.list[[data.name]] <- data.ag.split
    #return(data.ag.split)
    tax <- data.ag.split[, rev(names(data.ag.split)[ncol(data):ncol(data.ag.split)])]
    
    
    rownames(tax) <- tax[,1]
    #head(tax)
    tax<-tax[,-1]
    #head(tax)
    vec <- vector()
    # remove ranks with ONLY ONE level, i.e. all from k__Fungi, then kingdom level will be removed.
    for ( z in 1:ncol(tax) ) {
      if ( length(levels(factor(tax[,z]))) != 1 ) {
        vec <- c(vec, names(tax)[z])
      }
    }
    tax1 <- tax[, which(names(tax) %in% vec), drop=FALSE]
    if ( is.null(root) ) {
      root1 <- "root"
    } else {
      root1 <- root
    }
    #ade4
    #return(tax1)
    tax.phy <- ade4::taxo2phylog(ade4::as.taxo(tax1), abbrev=FALSE, root=root1)
    tax.phy.list[[data.name]] <- tax.phy
  }
  

  # plot
  if ( is.null(factors) ) {
    plot(tax.phy.list[[1]], cleaves=cleaves, cnodes=cnodes, clabel.leaves = clabel.leaves, clabel.nodes = clabel.nodes, f.phylog = f.phylog, sub = rank)
  } else {
    par(mfrow = c(length(factors), max(unlist(fac.lev))))
  # plot.phylo
    
    
    for ( x in 1:length(factors) ) {
      fac1 <- factors[x]
      lev.len <- length(levels(factor(meta[, fac1])))
      n.blank <- max(unlist(fac.lev)) - lev.len
      tax.phy.list1 <- tax.phy.list[grep(fac1, names(tax.phy.list))]
      for ( l in 1:length(tax.phy.list1) ) {
        tax.phy1 <- tax.phy.list1[[l]]
        n <- names(tax.phy.list1)[l]
        if ( !sub ) {
          sub1 <- ""
        } else {
          sub1 <- paste(names(tax.phy.list1)[l], " | ", rank, sep="")      
        }
        plot(tax.phy1, cleaves=cleaves, cnodes=cnodes, clabel.leaves = clabel.leaves, clabel.nodes = clabel.nodes, f.phylog = f.phylog, sub = sub1, csub = csub, possub = possub, draw.box = draw.box)
      }
      if ( n.blank != 0 ) {
        for ( p in 1:n.blank) {
          plot(1, type="n", axes=F, xlab="", ylab="")
        }
      } 
    }
  }
#detach(package:ade4)
par()
}

#detach(package:ade4)

