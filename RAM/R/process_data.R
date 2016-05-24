transpose.OTU <- function(data) {
  valid.OTU(data)
  # return the transpose (without the taxonomy column, 
  # which should be the last column)
  return(as.data.frame(t(data[, -dim(data)[2]])))
}

filter.OTU <- function(data, percent=NULL, number=NULL) {

  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  

  otu_sel <- list()
  otu.sel <- list()
  for ( i in 1:length(data) ){
    otu <- data[[i]]
    if (is.null(otu)) { break }
    label <- names(data)[i]    
    valid.OTU(otu)

    # both filter are NULL
    if (is.null(number) & is.null(percent)) {
       warning("No filter was provided, all otus will be used. Consider filtering otu by total number of sequence (e.g. 50) or maximum relative abundance (e.g. 0.01)")
      sel<-1:nrow(otu)
    } else if (!is.null(number) & !is.null(percent)) {
      stop("Error: can only provide one type of filter by total number of sequence (e.g. 50) or maximum relative abundance (e.g. 0.01)")
    } else if (!is.null(percent)) {
        # select OTUs more than percent in at least on sample (otu columns)
        otu.p <- cbind(decostand(otu[, -ncol(otu)], 
                       MARGIN=2, "total"), otu$taxonomy)
        names(otu.p)[ncol(otu.p)]<-"taxonomy"
        sel <- rownames(otu.p)[which(apply(otu.p[, -ncol(otu.p)], 
                     MARGIN=1, FUN=max) > percent)]
    } else if (!is.null(number)) {
        # select OTUs more than number sequences in total
        sel<-rownames(otu)[which(rowSums(otu[, -ncol(otu)]) > number)]
    } else {
        warning("no filtering requirment provided, will return the original otus")
    } 
    otu_sel[[label]] <- sel
    #return(otu_sel)  

    if ( length(sel) == 0 ) {
      warning("no otus in ", label, " met the filter requirement, original OTU returned")
      sel <- ""
    } else {
      print(paste(length(sel), " otus in ", label,
                  " met the filter requirment.", sep=""))
        otu.sel[[label]] <- otu[sel, ]
    }
  }
  return(otu.sel)
}


data.subset <- function(data, meta, factors="", values="", and=TRUE){

 if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  

  otu_sel <- list()
  for ( i in 1:length(data) ){
    otu <- data[[i]]
    if (is.null(otu)) { break }
    label <- names(data)[i]    
    valid.OTU(otu)
  
    # meta factor and value for selection
    factors <- .valid.factor(meta=meta, meta.factor=factors)
    factors1 <- parse(text = paste("meta$",factors,sep=""))
    values1 <- paste(values, collapse="|")

    if( length(factors1) > 0 ) {
      if ( !and ) {
          # and=FALSE; any subjects belong to the values should be kept
          meta_sel <- list()
          for ( i in 1:length(factors1) ) {
              meta_sel[[i]] <- meta[grepl(values1, eval(factors1[i])), ]
          }
          meta_sel <- meta_sel[lapply(meta_sel,length)>0]
          meta_sub <- unique(do.call("rbind", meta_sel))

      } else {
          # and=TRUE: subjects should belong to all values to be kept
          arg <- vector()
          for ( i in 1:length(factors1) ) {
              new <- paste("grepl( \"", values1, "\", ", factors1[i], " )", sep="")
              arg <- c(arg, new)
          } 
          arg.all <- paste(arg, collapse = " & ")
          sel <- eval(parse(text = arg.all))        
          meta_sub <- meta[sel, ]          
      }
    }

    if ( nrow(meta_sub) == 0 ) {
      stop("Error: no subject met the filter requirement")
    } else {
      for (i in factors) {
          if (is.factor(meta[[i]])) {
              # drop levels
              meta_sub[[i]] <- factor(meta_sub[[i]])
          }
      }
    }
  }

  # match samples among OTUs and metadata
  result=list()
  for (i in 1:length(data)) {
      otu <- data[[i]]
      label <- names(data)[i]
      if ( is.null(otu) ) { break }
      label <- names(data)[i]
      otu_sub = otu[, match(c(rownames(meta_sub),"taxonomy"), colnames(otu))]  
      ### remove empty rows and columns and order the data.frame ###
      otu_sub = otu_sub[rowSums(otu_sub[,-ncol(otu_sub)])>0, colSums(otu_sub[,-ncol(otu_sub)])>0]
      otu_sub = otu_sub[order(rowSums(otu_sub[,-ncol(otu_sub)]),decreasing=TRUE),]  
      ### Verify if the data.frame is non empty ###
      if(nrow(otu_sub) == 0){
          return ("Length of factor data is 0, please check if there were errors in your factors and factor.values")
      }
      if ( identical(rownames(meta_sub), colnames(otu_sub)[-ncol(otu_sub)]) ) {
          result[[label]]<-otu_sub
      } else {
          stop("meta and `otu` do not have same samples")      
      }
  }
  result[["meta"]] <- meta_sub
  return(result)
}


data.revamp <- function(data, is.OTU=TRUE, ranks=NULL,
                        stand.method=NULL, top=NULL, 
                        mode="number") {

  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=is.OTU)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  

  data.new <- list()
  for ( i in 1:length(data) ){
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]   
 
    if ( is.OTU ) {
      valid.OTU(elem)
      if ( is.null(ranks) || ranks == "" ) {
         # LCA of each otu
         lca <- LCA.OTU(elem)
         rownames(lca) <- paste(lca$LCA,
                                    rownames(lca), 
                                    sep="_")
         lca <- lca[, -ncol(lca)]
         lca <- as.data.frame(t(lca))
         lca.sort <- lca[, order(colSums(lca), decreasing=TRUE)]
         data.new[[label]]  <- lca.sort     
      } else {
         num.rank <- length(ranks)
         rank.list <- list()
         for ( i in 1:num.rank ) {
             .valid.rank(ranks[i])
             rank_name <- .get.rank.name(ranks[i])
             tax_rank <- tax.abund(elem, rank=ranks[i],     
                                    drop.unclassified=FALSE, 
                                    count=TRUE)
             tax_rank.sort <- tax_rank[, order(colSums(tax_rank),
                                          decreasing=TRUE)]
             lab <- paste(label, rank_name, sep="_")
             data.new[[lab]] <- tax_rank.sort
             
        }
     }    
    } else {
     data.new[[label]] <- elem[, order(colSums(elem), decreasing=TRUE)]
    }
   
  }

  # order by total counts or relative abund
  for ( i in 1:length(data.new) ) {
    elm <- data.new[[i]]
    if ( is.null(elm) ) { break }
    label <- names(data.new)[i]

    elm.percent <- decostand(elm, "total")
    ord.num <- order(colSums(elm), decreasing=TRUE)
    ord.per <- order(apply(elm.percent, MARGIN=2, 
                             FUN=max), decreasing=TRUE)
    if ( !is.null(top) ) {
      # check whether have so many top groups
      if ( top <= ncol(elm) ) {
        top <- top
      } else {
        warning(paste(label, " doesn't have ", top, " taxa groups, will retrun all", sep=""))
        top <- ncol(elm)
      }
      if ( mode == "number" ) {
        elm <- elm[, ord.num[1:top]]
      } else if ( mode == "percent" ) {
        elm <- elm[, ord.per[1:top]]
      } else {
        warning("seletion mode must be \"number\" or \"percent\", 
        the default is \"number\".  Will keep the taxa with 
        highest total count") 
        elm <- elm[, ord.num[1:top]]
      }
    } else {
       elm <- elm
    }
   data.new[[label]] <- elm
 }
  
  # data transformation & remove unclassified
  remove.pat <- paste(.blacklist(), collapse="|")
  for ( i in 1:length(data.new) ) {
    elm <- data.new[[i]]
    if ( is.null(elm) ) { break }
    label <- names(data.new)[i]
      
    keep <- !grepl(remove.pat, names(elm),ignore.case=TRUE)
    if ( length(keep) != 0 ) {
        if ( !is.null(stand.method) ) {
          elm <- decostand(elm, stand.method)
          elm <- elm[, keep, drop=FALSE]  
        } else { 
          elm <- elm[, keep, drop=FALSE]
        }
    } else { 
         warning("No identified groups met requirement, will retrun 
         all groups without filtering")
         elm <- elm
    }
    data.new[[label]] <- elm
  }

  return(data.new)
}


filter.Taxa <- function(taxa, drop.unclassified=TRUE,
                        percent=NULL, number=NULL) {

 # whether or not keep unclassified groups
 remove.pat <- paste(.blacklist(), collapse="|")
 if ( drop.unclassified ) {
    taxa <- taxa[, !grepl(remove.pat, names(taxa))]
 } else {
    taxa <- taxa
 }
  # both filter are NULL
  if (is.null(number) & is.null(percent)) {

        warning("No filter was provided, all taxa will be used. Consider filtering taxonomic abundance matrix by total number of sequence (e.g. 50) or maximum relative abundance (e.g. 0.01)")
    sel<-1:nrow(taxa)

  } 
  if ( !is.null(number) & !is.null(percent) ) {
    warning("should only provide one type of filter by total number of sequence (e.g. 50) or maximum relative abundance (e.g. 0.01), will filter by percent as default")
    taxa.p <- decostand(taxa, MARGIN=1, "total")
    sel <- which(apply(taxa.p, MARGIN=2, FUN=max) > percent)
  } 
  if ( !is.null(percent) ) {
    # select taxa more than percent in at least on sample 
    taxa.p <- decostand(taxa, MARGIN=1, "total")
    sel <- names(taxa.p)[which(apply(taxa.p, MARGIN=2, FUN=max) > percent)]
  } 
  if ( !is.null(number) ) {
    # select taxa more than number sequences in total
    sel <- names(taxa)[which(colSums(taxa) > number)]
 } 

 if ( length(sel) == 0) {
   warning("no taxa met the filter requirement, original taxa abundance matrix returned")
   taxa.new <- taxa    
 } else {
   taxa.new <- taxa[, names(taxa) %in% sel]
   print(paste(length(names(taxa.new)), " taxa met the filter requirment", sep=""))
 }  
 return(taxa.new)
} 


OTU.rarefy <- function(data, sample=NULL) {
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  

  data.new <- list()
  for ( i in 1:length(data) ){
    otu <- data[[i]]
    if ( is.null(otu) ) { break }
    label <- names(data)[i]   
    valid.OTU(otu)

    # transpose OTU
    otu.t<-transpose.OTU(otu)
    # OTU taxonomy
    otu_tax <- otu$taxonomy
  
    # rarefy
    if (is.null(sample)) {
       otu.t.rf<-rrarefy(otu.t, sample=min(apply(otu.t, 1, sum)))
    } else if (is.numeric(sample)) {
      if (sample <= min(apply(otu.t, 1, sum))) {
         otu.t.rf<-rrarefy(otu.t, sample=sample)
      } else {
        warning("sampling number is larger than the minimum sequence of a sample, using minimum number instead ")
        otu.t.rf<-rrarefy(otu.t, sample=sample)
      }
    } else {
       stop("Error: what's the sampling number?")
    }

    otu.new <-as.data.frame(t(otu.t.rf))
    otu.new$taxonomy<-otu_tax
    data.new[[label]] <- otu.new     
  }
  return(data.new)
}


combine.OTU <- function(data, meta) {
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  

  matched <-match.data(data=data, is.OTU=TRUE, meta=meta)
  meta <- matched$meta
  # check sampleID orders, 
  # make the samples IDs in otu tables match metadata
  order <- rownames(meta)
  otu.list <- list()
  samples <- list()
  for ( i in 1:(length(matched)-1) ){
    otu <- matched[[i]]
    if ( is.null(otu) ) { break }
    label <- names(matched)[i]   
    valid.OTU(otu)
    .valid.meta(otu1=otu, meta=meta)
    # give a prefix 
    sampleID <- colnames(otu)[-ncol(otu)]
      if ( ! identical(sampleID, order) ) {
          otu <- otu[, match(c(order, "taxonomy"), colnames(otu))]
          if ( ! identical(c(order, "taxonomy"), colnames(otu)) ) {
              stop(paste(label, " doesn't have same samples as metadata", 
                       sep=""))
          }
      }
      #rownames(elm) <- paste(lab, rownames(elm), sep="_")
      otu.list[[label]] <- otu
  }
  df <- do.call("rbind", otu.list)
  return(df)
}


match.data <- function(data, is.OTU=TRUE, meta) {

  samples <- list()
  samples[["meta"]] <- rownames(meta)

  elems <- list()
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]
    elems[[label]] <- elem
    # otus
    if ( is.OTU ) {
       samples.elem <- colnames(elem)[-ncol(elem)]
    } else {
       samples.elem <- rownames(elem)
    }
    samples[[label]] <- samples.elem
  }
  
  # samples exist in all data and metadata
  samples.all <- Reduce(intersect, samples)
  
  # new metadata
  sel <- which(samples[["meta"]] %in% samples.all)
  meta.new <- meta[sel, ]
       
  # new data
  elems.new <- list()
  for ( i in 1:length(elems) ) {
    label <- names(elems)[i]
    elem <- elems[[i]]
    if ( is.null(elem) ) { break }
    if ( is.OTU ) {
      to_match <- c(rownames(meta.new), "taxonomy")
      elem.new <- elem[, match(to_match, 
                           colnames(elem))] 
      elem.new <- elem.new[rowSums(elem.new[, 
                                     -ncol(elem.new)])>0,    
                          colSums(elem.new[, 
                                     -ncol(elem.new)]) >0]
      if ( identical(to_match, colnames(elem.new)) ) {
        elems.new[[label]] <- elem.new
      } else {
        stop("couldn't match samples in data and metadata, please check your datasets")
      } 
    } else {
      to_match <- rownames(meta.new)
      elem.new <- elem[match(to_match, 
                           rownames(elem)), ]
      elem.new <- elem.new[rowSums(elem.new)>0,    
                           colSums(elem.new) >0]
      if ( identical(to_match, rownames(elem.new)) ) {
        elems.new[[label]] <- elem.new
      } else {
        stop("couldn't match samples in data and metadata, please check your datasets")
      } 
    }
  }
  elems.new[["meta"]] <- meta.new
  return(elems.new)
}
