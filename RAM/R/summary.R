shared.OTU <- function(data) {
  
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  
  
  shared <- list()
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]
    
    otu <- transpose.OTU(elem)
    otu.pa <- decostand(otu, "pa")
    otu.lca <- LCA.OTU(otu=elem)
    
    # for readability
    x <- nrow(otu.pa)
    y <- ncol(otu.pa)
    
    if(length(which(colSums(otu.pa) == x))==0) {
      warning(paste("NO OTUs are shared by all samples in ", label, sep=""))
    } 
    
    num.otu.one.sample <- length(which(colSums(otu.pa) == 1))
    num.otu.mult.sample <- length(which(colSums(otu.pa) > 1))
    num.otu.all.sample <- length(which(colSums(otu.pa) == x))
    
    per.otu.one.sample <- num.otu.one.sample / y
    per.otu.all.sample <- num.otu.all.sample / y
    
    if(length(colnames(otu.pa)[which(colSums(otu.pa) == x)]) == 0){
      num.seq.shared.otu <- 0
    } else {
      num.seq.shared.otu <- sum(otu[ ,colnames(otu) %in% 
                                      colnames(otu.pa)[which(colSums(otu.pa) == x)]])
    }
    per.seq.shared.otu <- num.seq.shared.otu / sum(otu)
    
    # get the OTU #s for all OTUs present in all samples
    otu.IDs <- colnames(otu.pa[ ,colSums(otu.pa) == x])
    # make a list of the form "OTU-taxonomic_information"
    sel <- which(rownames(otu.lca) %in% otu.IDs)
    tax <- paste(otu.IDs, otu.lca[sel, "LCA"], sep=": ")
    
    val <- list(num.otu.one.sample, num.otu.mult.sample, 
                num.otu.all.sample, per.otu.one.sample, 
                per.otu.all.sample, num.seq.shared.otu,
                per.seq.shared.otu, otu.IDs, tax)
    
    # format numerical data?
    #val <- format(round(val[1:7], 2), nsmall=2)
    
    names(val) <- c("#_of_OTUs_in_1_sample", 
                    "#_of_OTUs_in_>1_sample", 
                    "#_of_OTUs_in_all_samples", 
                    "%_of_OTUs_in_one_sample",
                    "%_of_OTUs_in_all_samples", 
                    "#_of_sequence_in_shared_OTUs",
                    "%_of_sequence_in_shared_OTUs", 
                    "OTUs_in_all_samples", 
                    "LCA_of_OTUs_in_all_samples")
    shared[[label]] <- val
  }
  if (length(shared) == length(data) ) {
    return(shared)
  } else {
    warning(paste("not all data being processed: only ", names(shared), " being processed", sep=""))
  }
}

    
percent.classified <- function(data, ranks=c("f","g")) {

  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  

  classified <- list()
  for ( i in 1:length(data) ) {
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    label <- names(data)[i]
    valid.OTU(elem)

    rk.classified <- list()
    for ( rk in ranks ) {
      .valid.rank(rk)
      rank_name <- .get.rank.name(rk)
      rank_pat <- .get.rank.pat(rk)
      blacklist <- paste(.blacklist(rk), "|", rank_pat, "$", sep="")
  
      # number of rows in get.rank gives us the number of OTUs classified at that
      # rank; divide by total number of OTUs and multiply to get %
      # return(100 * dim(get.rank(otu1=data, rank=rank))[1] / dim(data)[1])
    
      elem_rank <- get.rank(otu1=elem, rank=rk)
          
      if(length(grep(blacklist, elem_rank$taxonomy)) == 0) {
          elem_ided <- elem_rank
      } else {      
          elem_ided <- elem_rank[!grepl(blacklist, elem_rank$taxonomy), , drop=FALSE]
      }
      name <- paste("%_OTU_classified_at_", rank_name, "_level:", sep="")
      rk.classified[[name]] <- (100 * nrow(elem_ided) / nrow(elem))
    }   
    classified[[label]] <- rk.classified
  }
  return(classified)
}


OTU.recap <- function(data, ranks=c("p", "c", "o", "f", "g"), 
                      brewer.pal="Pastel1", file=NULL, ext="pdf", 
                      width=12, height=8) {

  save <- !is.null(file)
 
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  
    
  if (length(ranks)==0) {
    stop("Error: please provide taxonomic ranks to summarize, e.g. ranks=c('p', 'o', 'c', 'f', 'g', 's')")
  }

  names<-c("total_otu", "total_seq", "singleton", "percent_singleton", 
"percent_singleton_seq", "otu", "seq", "identified_otu", "identified_seq", 
"percent_identified_otu_of_total", "percent_identified_seq_of_total", 
"percent_identified_otu_at_rank", "percent_identified_seq_at_rank", 
"identified_tax")

  len <- length(names)
  len.ranks<-length(ranks)  
  
  df<-list()
  otu_tax<-list()
  df.total <- list()
  
  level <- vector()
  for ( j in ranks ) {
     rk <- .get.rank.name(j)
     level <- unique(c(level,  rk))
  }
 
 
  for (i in 1:length(data) ){
    label <- names(data)[i]
    otu <- data[[i]]

    # exit if no more otus 
    if ( is.null(otu) ) { break }
  
    valid.OTU(otu)
    total_OTUs_num <- nrow(otu)
    total_seqs <-sum(otu[, -ncol(otu)])

    # singletons
    singleton <- which(rowSums(otu[, -ncol(otu)]) == 1)

    if( length(singleton) ==0 ) {
        warning("NO singletons (otus that being observed only once)")
	singleton_OTUs_num <- 0
	singleton_OTUs_percent <- 0
	singleton_seqs_percent <- 0
    } else {

      singleton_OTUs_num <-  nrow(otu[singleton, ])
      singleton_OTUs_percent <- round(100*nrow(otu[rowSums(otu[, 
                           -ncol(otu)])==1,])/nrow(otu), digits=3)
      singleton_seqs_percent <- round(100*nrow(otu[rowSums(otu[, 
                               -ncol(otu)])==1,])/total_seqs, digits=3)
    }

      df.names <- names 

      df[[label]]<-data.frame(matrix(vector(), 0, len, dimnames=list(c(), 
                          df.names)), stringsAsFactors=F)    
    
    tax<-list()
    for (j in ranks) {
      .valid.rank(j)
      rank<- .get.rank.name(j)
      rank_pat<-.get.rank.pat(j)
      blacklist <- paste(.blacklist(j), "|", rank_pat, "$", sep="")
      
      # All otus with rank_pat, e.g. f__      
      rank_otu <- otu[grep(rank_pat, otu$taxonomy),]
      if ( nrow(rank_otu) == 0 ) {
        rank_otu_num <- 0
        rank_otu_seq <- 0
      } else {
        rank_otu_num <- nrow(rank_otu)
        rank_otu_seq <- sum(rank_otu[, -ncol(rank_otu)])
      }


      # only identified otus at the rank
      rank_otu_ided <- get.rank(otu1=otu, rank=j)
      col.ided <- ncol(rank_otu_ided)

      # split the taxonomy column
      tax.classes <- c("kingdom", "phylum", "class", "order", 
                          "family", "genus", "species") 
      rank_otu_ided_splitup <- col.splitup(rank_otu_ided, col="taxonomy", 
                                    sep="; ", drop=TRUE, names=tax.classes)
      
      
      if ( (unique(rank_otu_ided_splitup[[rank]])) != "" ) {
          # identified taxa at the rank  
          rank_otu_ided_splitup.rank <- cbind(rank_otu_ided[, 1:(col.ided-1)], 
                                           rank_otu_ided_splitup[[rank]])
          names(rank_otu_ided_splitup.rank)[ncol(rank_otu_ided_splitup.rank)] <- rank
          tax[[j]] = unique(rank_otu_ided_splitup.rank[[rank]])
          rank_tax_ided_num = length(tax[[j]])
     
          rank_otu_ided_num <- nrow(rank_otu_ided)
          rank_otu_ided_seq <- sum(rank_otu_ided[, -ncol(rank_otu_ided)])
          rank_otu_ided_num_percent <- round(100*nrow(rank_otu_ided)/nrow(rank_otu), digits=3)       
          rank_otu_ided_seq_percent <- round(100*sum(rank_otu_ided[, 
                  -ncol(rank_otu_ided)]) / sum(rank_otu[, -ncol(rank_otu)]),  digits=3)
          percent_identified_otu_of_total <- round(100*nrow(rank_otu_ided)/nrow(otu), digits=3)
          percent_identified_seq_of_total <- round(100*sum(rank_otu_ided[, 
                  -ncol(rank_otu_ided)]) / sum(otu[, -ncol(otu)]),  digits=3)

         
      } else {
	  warning(paste("no otus being identified at ", rank, " level"))
	  rank_otu_ided_num <- 0
          rank_otu_ided_seq <- 0
          rank_otu_ided_num_percent <- 0       
          rank_otu_ided_seq_percent <- 0
          percent_identified_otu_of_total <- 0
          percent_identified_seq_of_total <- 0	
          tax[[j]] = ""
          rank_tax_ided_num = 0
      }

      df[[label]][rank,] <- c(total_OTUs_num, total_seqs, 
                  singleton_OTUs_num, singleton_OTUs_percent, 
                   singleton_seqs_percent, rank_otu_num, 
                 rank_otu_seq, rank_otu_ided_num, rank_otu_ided_seq,  
                 percent_identified_otu_of_total, percent_identified_seq_of_total,  
                 rank_otu_ided_num_percent, rank_otu_ided_seq_percent, rank_tax_ided_num)
    }
    otu_tax[[label]] <- tax
  }
  #return(otu_tax)

 # multiple data sets to be compared for taxa lists  
  if ( length(data) > 1 ) {
    df_tax_unique <- list()
    for ( i in 1:length(data) ) {
      label <- names(data)[i]
      df_tax_unique_rank <- list()
      for (j in ranks) {
        rank<- .get.rank.name(j)
         # find what taxa in data[[i]] not in others
        target <- as.character(otu_tax[[i]][[j]])
        other <- otu_tax[-i]
        vec <- vector()
        for ( o in 1:length(other) ) { 
            vec <- c(vec, other[[o]][j])
        }
        others <- vector()
        for (v in 1:length(vec) ) {
           others <- unique(c(others, as.character(vec[[v]])))
        }
        
        tax_num <- length(target)
        tax_num_unique <- length(setdiff(target, others))
        tax_unique <- paste(setdiff(target, others), collapse=",")
        name <- paste(rank, ": ", "classified_to_", 
                       tax_num, "; ", tax_num_unique, 
                       "_only_in_", label, sep="")
        df_tax_unique_rank[[name]] <- tax_unique
      }
     df_tax_unique[[label]] <- df_tax_unique_rank
    }
    df[["tax_unique"]] <- df_tax_unique
  }
   

   # plot percent classified
  #if ( !require("reshape2") ) {
  #   stop("package 'reshape2' is required for this function")
  #}
 # if ( !require("RColorBrewer") ) {
 #    stop("package 'RColorBrewer' is required for this function")
 # }
  
  df.m <- list()
  
  for ( i in 1:length(data) ) {
     label <- names(data)[i]
     df.m[[label]] <- melt(cbind(df[[label]][,10:13], 
                        ind=rownames(df[[label]])))
     df.m[[label]]$Region <- label
  }
  
  if ( length(data) == 1 ) { 
     df.m.total <- df.m[[1]]
  }

  if ( length(data) > 1 ) { 
     df.m.total <- do.call("rbind", df.m)
  }
 
  # make sure rank order as it was provided
  df.m.total$ind <- factor(df.m.total$ind, levels=level)
  len <- length(df.m.total$variable)
  breaks <- levels(df.m.total$variable)
  lab <- c("Total_OTUs", "Total_Seqs", "OTUs_at_rank", 
               "Seqs_at_rank")
  p <- ggplot(df.m.total, aes_string(x="ind", y="value", fill="variable")) + 
           geom_bar(stat="identity", position="dodge")  + 
           scale_fill_manual(values=brewer.pal(len, brewer.pal), 
                    name="Percent Classified", 
                    breaks=breaks, labels=lab) + 
                    xlab("Taxonomic Ranks") + 
                    ylab("Percent ( % )")

  if ( length(data) == 1 ) { 
     p <- p
  }   
  if ( length(data) > 1 ) {
   
    p <- p +
         facet_wrap(~Region)
  }

  
  if (save) {
      .ggsave.helper(file, ext, width, height, plot=p) 
  } else {
       print(p)
  }
  
   return(df)
}


.true.calc <- function(data, index, mode) {
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  

  sampleids <- .valid.data(data, is.OTU=TRUE, 
                           export.sampleIDs=TRUE)
  # empty dataframe to store the output
  output <- data.frame(matrix(vector(), 0, length(sampleids), 
                   dimnames=list(c(), sampleids)), 
                    stringsAsFactors=F)

  for ( i in 1:length(data) ){
    otu <- data[[i]]
    if (is.null(otu)) { break }
    label <- names(data)[i]    
    valid.OTU(otu)
  
    otu.t <- transpose.OTU(otu)
  
    if (!is.character(index)) {
      stop("argument 'index' must be either 'simpson' or 'shannon'.")
    }
 
    if (mode == "diversity") {
      func <- .true.div
    } else if (mode == "evenness") {
      func <- .true.even
    }
  
    methods <- c("simpson", "shannon")
    ### is this error message descriptive enough?
    index <- match.arg(index, methods)
    res <- func(otu.t, index)
    output[label,] <- res
  } 
  return(output)
}

true.diversity <- function(data, index="simpson") {
  .true.calc(data=data, index=index, mode="diversity")
}

evenness <- function(data, index="simpson") {
  .true.calc(data=data, index=index, mode="evenness")
}

.true.div <- function(OTU, index) {
  
  if (index == "simpson") {
    return(vegan::diversity(OTU, index="invsimpson"))
    
  } else if (index == "shannon") {
    return(exp(vegan::diversity(OTU, index="shannon", base = exp(1))))
  }
}

.true.even <- function(OTU, index) {
  
  if (index == "simpson") {
    return(vegan::diversity(OTU, index="invsimpson") / specnumber(OTU))
    
  } else if (index == "shannon") {
    return(vegan::diversity(OTU, index="shannon", base = exp(1)) / 
             log(vegan::specnumber(OTU), base=exp(1)))
  }
}


OTU.diversity<-function(data, meta) {
  
  if ( class(data) != "list" ) {
    stop("please provide otu tables as list; see ?RAM.input.formatting")
  } 
  .valid.data(data, is.OTU=TRUE)
  ###check for zero entries?
  num.data <- length(data)
  labels <- names(data)  
  
  for (i in 1:length(data) ){
    label <- names(data)[i]
    otu <- data[[i]]
    if (is.null(otu)) { break }
    # name of the otu
    
    valid.OTU(otu)
    .valid.meta(otu, meta=meta)  
    
    otu.t<-transpose.OTU(otu)
    
    # add columns for different diversity indices
    meta[[paste("spec", label, sep="_")]] <- vegan::specnumber(otu.t)
    meta[[paste("sim", label, sep="_")]] <- vegan::diversity(otu.t, index="simpson", MARGIN=1)
    meta[[paste("invsim", label, sep="_")]] <- vegan::diversity(otu.t, index="invsimpson", MARGIN=1)
    meta[[paste("shan", label, sep="_")]] <- vegan::diversity(otu.t, index="shannon", MARGIN=1)
    data2<-list(otu)
    names(data2)[1] <- label
    meta[[paste("sim_even", label, sep="_")]] <- as.vector(t(
      evenness(data=data2)[label,]))
    meta[[paste("shan_even", label, sep="_")]] <- as.vector(t(
      evenness(data=data2, index="shannon")[label,]))
    meta[[paste("sim_trudiv", label, sep="_")]] <- as.vector(t(
      true.diversity(data=data2)[label,]))
    meta[[paste("shan_trudiv", label, sep="_")]] <- as.vector(t(
      true.diversity(data=data2, index="shannon")[label,]))
    meta[[paste("chao", label, sep="_")]] <- vegan::estimateR(otu.t)["S.chao1",]
    meta[[paste("ACE", label, sep="_")]] <- vegan::estimateR(otu.t)["S.ACE",] 
    
  }
  
  return (meta)
}

group.spec  <-  function(otu, meta, factor, file=NULL, ext=NULL,
                         width=8, height=8) {
 
  save  <-  !is.null(file)
  
  .valid.meta(otu1=otu, meta=meta)
  m <- meta
  m$rn <- rownames(m)
  m <- m[, c("rn", factor)]
  names(m) <-c("rownames", factor)
  m <- m[complete.cases(m),]
  not <- setdiff(rownames(meta), rownames(m))
  if ( length(not) !=0 ) {
    warning(paste("The following samples contain missing data in ", factor, ": ", paste(not,collapse=", "), "; which are excluded in the analysis", sep=""))
  }
  meta <- m

  ex <- vector()
  for ( i in levels(factor(meta[[factor]])) ) {
   n <- rownames(meta)[which(meta[[factor]]==i)]
   if ( length(n) == 1 ) {
     warning(paste("Only 1 sample, ", n, ", is from ", i, "; some diversity indices cannot be calculated, this sample is removed from the analysis", sep=""))
     ex <- c(ex, n)
    }
  }
  if ( length(ex) == 0 ) {
    meta <- meta
  } else {
    rm <- paste(ex, collapse="|")
    meta <- meta[!grepl(rm, rownames(meta)), ]
  }
  meta[[factor]] <- factor(meta[[factor]])
  otu <- otu[, match(c(rownames(meta),"taxonomy"), colnames(otu))]
  
  # calculate specnumber, specpool & specpool2vect
  # specpool for pooled samples (e.g. metadata categories): 
  # Function 'specpool' uses presence data. The incidence-based, 
  # i.e. presence/absence, estimates in 'specpool' use the 
  # frequencies of species in a collection of sites.

  pool  <-  specpool(transpose.OTU(otu), meta[[factor]])
  pool.vect <- specpool2vect(pool)

  # specpool2vect: 'specpool2vect' maps the pooled values into a 
  # vector giving the value of selected 'index' for each original 
  # site. specpool2vect(pool, index = c("jack1","jack2", "chao", 
  # "boot","Species")) 
  
  index = c("jack1","jack2", "chao", "boot")
  df<-data.frame(matrix(vector(), 0, nrow(meta), 
         dimnames=list(c(), rownames(meta))), stringsAsFactors=F) 

  OTU_observed<-specnumber(transpose.OTU(otu))
 
  for ( i in index ) {
      rn <- paste("Percent_OTU_observed", i, sep=":")
      df[rn,] <- 100*OTU_observed/specpool2vect(pool, index=i)
  }
  
  df["OTU_observed",]  <-  specnumber(transpose.OTU(otu))
  df.t <- as.data.frame(t(df))
  if ( identical(rownames(df.t), rownames(meta)) ) {
      df.t[[factor]] <- meta[[factor]]
  }
  
  # melt data frame
  #if (!require("reshape2")) {
  #  stop("package 'reshape2' is required to use this function: try 'install.packages('reshape2')'.")
  #}
  
  df.m  <-  melt(cbind(df.t, rownames(df.t)))
  names(df.m)[2]  <-  "Sample"
  names(df.m)[3]  <-  "Index"

  # plot
  p <- ggplot(data=df.m, aes_string(x=factor, y="value", 
            fill=factor)) + 
            geom_boxplot() + 
            facet_wrap(~Index, ncol=2, scales="free") + 
            theme(axis.text.x = element_text(angle = 45, 
                 vjust = 1, hjust=1)) + ylab("")

#  return(list(pool, spec.obs, percent.spec.obs, fac, df, df.m))
  if (save) {
    file <- .ensure.filepath(file, ext)
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
  
}


group.rich <- function(otu, meta, factor, file=NULL, ext=NULL,
                       width=8, height=8) {
  
  #library("reshape2") #melt()
  #  library("vegan") # specpool()
  #  library("reshape") # rename(). colsplit()
  
  save  <-  !is.null(file)
  # vegan
  V1 = se = NULL
  .valid.meta(otu1=otu, meta=meta)
  .valid.factor(meta, factor)
  m <- meta
  m$rn <- rownames(m)
  m <- m[, c("rn", factor)]
  names(m) <-c("rownames", factor)
  m <- m[complete.cases(m),]
  not <- setdiff(rownames(meta), rownames(m))
  if ( length(not) !=0 ) {
    warning(paste("The following samples contain missing data in ", factor, ": ", paste(not,collapse=", "), "; which are excluded in the analysis", sep=""))
  }
  meta <- m
  
  ex <- vector()
  for ( i in levels(factor(meta[[factor]])) ) {
    n <- rownames(meta)[which(meta[[factor]]==i)]
    if ( length(n) == 1 ) {
      warning(paste("Only 1 sample, ", n, ", is from ", i, "; some diversity indices cannot be calculated, this sample is removed from the analysis", sep=""))
      ex <- c(ex, n)
    }
  }
  if ( length(ex) == 0 ) {
    meta <- meta
  } else {
    rm <- paste(ex, collapse="|")
    meta <- meta[!grepl(rm, rownames(meta)), ,drop=FALSE]
  }
  meta[[factor]] <- factor(meta[[factor]])
  otu <- otu[, match(c(rownames(meta),"taxonomy"), colnames(otu)), drop=FALSE]
  pool <- vegan::specpool(transpose.OTU(otu), meta[[factor]])
  pool <- pool[, -dim(pool)[2]]
  
  #if (!require("reshape2")) {
  #  stop("package 'reshape2' is required to use this function: try 'install.packages('reshape2')'.")
  #}
  
  #if (!require("reshape")) {
  #  stop("package 'reshape' is required to use this function: try 'install.packages('reshape')'.")
  #}
  
  pool.melt <- reshape2::melt(cbind(pool,ind=rownames(pool)), is.vars=c('ind'))
  names(pool.melt)[names(pool.melt)=="ind"] <- factor
  names(pool.melt)[names(pool.melt)=="variable"] <- "Index"
  #pool.melt <- rename(pool.melt,c(ind=paste0(factor), variable="Index"))
  pool.melt.split <- cbind(pool.melt, 
                           Index = reshape2::colsplit(pool.melt$Index, 
                                                      "\\.",  c('type', 'se')))
  pool.melt.split.cast = reshape::cast(pool.melt.split, 
                                       paste0(factor, "+ Index.type ~ Index.se"))
  
  # plot
  # bar posittion. bit overlap
  dodge  <- ggplot2::position_dodge(width = 0.8)
  
  x <- colnames(pool.melt.split.cast)[1]
  y <- colnames(pool.melt.split.cast)[3]
  fill <- colnames(pool.melt.split.cast)[2]
  legend <- levels(pool.melt.split.cast$Index.type)
  xlab <- colnames(pool.melt.split.cast)[1]
  ylab <- "richness"
  

  specpool.plot <- ggplot2::ggplot(pool.melt.split.cast, 
                          aes_string(x=x, y=y, fill=fill)) + 
    geom_bar(stat="identity",position=dodge) +  
    geom_errorbar(aes(ymin=V1-se, ymax=V1+se), 
                  position = dodge, size = 0.5, 
                  linetype = 1, width = 0.1) + 
    scale_fill_brewer(legend, type = "div" , 
                      palette = "Set3" ) +  
    xlab(xlab) + ylab(ylab)
  #  return(list(pool, pool.melt, pool.melt.split.cast))
  
  if (save) {
    file <- .ensure.filepath(file, ext)
    .ggsave.helper(file, ext, width, height, plot=specpool.plot)
  } else {
    specpool.plot
  }
  
}

