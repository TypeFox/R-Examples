.analysis.helper <- function(otu1, otu2=NULL, meta, full, exclude, mode, rank, na.action=na.exclude) {
  
  valid.OTU(otu1, otu2)
  .valid.meta(otu1, otu2, meta)
  .valid.rank(rank)
  
  single.otu <- is.null(otu2)
  
  if (!single.otu) {
    return(list(otu1=.analysis.helper(otu1, meta=meta, full=full,
                                      exclude=exclude, mode=mode, rank=rank, na.action=na.exclude),
                otu2=.analysis.helper(otu2, meta=meta, full=full,
                                      exclude=exclude, mode=mode, rank=rank, na.action=na.exclude)))
  }
  
  tax<-tax.abund(otu1, otu2=NULL, rank=rank, drop.unclassified=TRUE,
                 top=NULL, count=TRUE, mode="number")
  
  meta <- filter.META(meta=meta, excl.na=FALSE)
  
  if (is.null(exclude)) {
    meta.remain <- meta
  } else {
    meta.remain <- meta[ , -exclude]
  }
  
  if (full) {
    form <- stats::formula(tax ~ .)
  } else {
    form <- stats::formula(tax ~ 1)
  }
  
  if (mode == "cca") {
    mod <- vegan::cca(form, data=meta.remain, na.action=na.action)
  } else if (mode == "rda") {
    mod <- vegan::rda(form, data=meta.remain, na.action=na.action)
  }
  
  good <- vegan::goodness(mod,summ = TRUE)
  vif.scores <- vegan::vif.cca(mod)
  #pct.var <- vegan:::summary.cca(mod)$concont$importance
  pct.var <- utils::getFromNamespace("summary.cca", "vegan")(mod)$concont$importance
  CCA.eig <- mod$CCA$eig
  CA.eig <- mod$CA$eig
  anv <- vegan::anova.cca(mod)
  
  return(list(GOF=good, VIF=vif.scores, percent_variation=pct.var, CCA_eig=CCA.eig, 
              CA_eig=CA.eig, anova=anv))
}

assist.cca <- function(otu1, otu2=NULL, meta, full=TRUE, exclude=NULL, rank, na.action=na.exclude) {
  .analysis.helper(otu1, otu2, meta, full, exclude, mode="cca", rank=rank, na.action=na.action)
}


assist.rda <- function(otu1, otu2=NULL, meta, full=TRUE, exclude=NULL, rank, na.action=na.exclude) {
  .analysis.helper(otu1, otu2, meta, full, exclude, mode="rda", rank=rank, na.action=na.action)
}


# plot rda analysis for OTUs
OTU.ord <- function(otu, meta=meta, factors=NULL, group=NULL, 
                 na.action=c("na.fail", "na.omit", "na.exclude"),
                 rank="g", taxa=NULL, data.trans="total", 
                 plot.species=TRUE, plot.scaling=-1, 
                 biplot.scale=NULL, biplot.sig=NULL, biplot.label= TRUE, 
                 mode=c("rda", "cca"), choice=c(1,2), 
                 main="", cex.point=3, cex.leg=12, cex.bp=3,
                 file=NULL, ext=NULL, width=10, height=10) {

  .valid.meta(otu, meta=meta)
  # plot otus on ordination space.
  # get identified ranks
  rank_name <- .get.rank.name(rank)
  rank_pat <-.get.rank.pat(rank)
  
  if ( is.null(factors) ) {
     factors1 <- names(meta)
  } else {
     factors1 <- .valid.factor(meta, factors)
  }
  meta.new <- meta[, which(names(meta) %in% factors1)]
  meta.new <- filter.META(meta=meta.new, excl.na=FALSE)
  
  # return(meta.new)
  # problem passing dt to the formula
  # error: Error in parse(text = x) : 
  #<text>:2:0: unexpected end of input 1: dt ~ 
  #formula <- paste("dt", paste(factors, collapse=" + "), 
  #                 sep=" ~ ")
  #form <- as.formula(formula)

  if ( length(mode) != 1L || 
    !any(grepl(mode, c("rda", "cca"))) ) {
    warning("Only either 'rda' or 'cca' can be used as mode,
    will use 'rda' as default for analysis")
    mode <- "rda"
  } else {
    mode <- mode
  }
  
 # handling missing data in metadata 

  na.sel=c("na.fail", "na.omit", "na.exclude")
  if ( length(na.action) != 1L || !any(na.action %in% na.sel)) {
    na.action <- "na.fail"
    warning("na.action was not properly selected, will used default 'na.fail'")
    suppressWarnings(meta.new <- filter.META(meta.new, excl.na=TRUE))  
  } else { 
    if ( na.action == "na.fail" ) {
      na.action <- "na.fail"
      suppressWarnings(meta.new <- filter.META(meta.new, excl.na=TRUE))  
    } else if (na.action == "na.exclude") {
      na.action <- "na.exclude"
      suppressWarnings(meta.new <- filter.META(meta.new, excl.na=FALSE)) 
    } else if (na.action == "na.omit") {
      na.action <- "na.omit"
      suppressWarnings(meta.new <- filter.META(meta.new, excl.na=FALSE)) 
      meta.new <- na.omit(meta.new)
      if ( nrow(meta.new) == 0 ) {
        stop("no records left after removing missing data")
      }
    } else {
       stop("how do you want to handle the missing data?")
    }      
  }
  
# prepare dataset after removing missing data in metadata table
  otu <- otu[, match(c(rownames(meta.new), "taxonomy"), colnames(otu))]
  OTU <- .subset.otu(otu=otu, rank=rank, taxa=taxa)
  # transpose to sample x taxa 
  OTU_otu <- OTU[, -ncol(OTU)]
  OTU_otu.t <- as.data.frame(t(OTU_otu))

  # taxa levels at the selected rank
  OTU_tax <- OTU[[rank_name]]
  tax <- OTU_tax
  tax.lvl <- levels(factor(OTU_tax))
  #return(OTU)
  # RDA / cca
  # problem in passing output to anova testing.
  # so have to use the original otu table
  if ( is.null(data.trans) ) {
    OTU_otu.t <- OTU_otu.t
  } else {
    OTU_otu.t <- decostand(OTU_otu.t, data.trans)
  }

  # check sample # & order 
  if( !identical(rownames(OTU_otu.t), rownames(meta.new))) {
     stop ("Error: samples doesn't match between data and metadata after handling missing data")
  }
  
  dt <- OTU_otu.t
  if (mode=="rda") {
    # scale=TRUE would scale species to unit variance, but 
    # since we have standadized data, so this being set FALSE
    ord_model = vegan::rda(dt ~., data=meta.new, scale=FALSE,
                           na.action=na.action)
  } else {
    ord_model = vegan::cca(dt ~., data=meta.new, scale=FALSE,
                           na.action=na.action)
  }

  #return(list(meta.new, OTU_otu.t, ord_model))
  # anova test for each factor
  #test.anova <- anova(ord_model, by = "terms", perm = 1000)
  
  # res <- .perform_ord(OTU, count=count, meta=meta, factors=factors, mode=mode)

  choice <- choice
  res1 <- .eig(ord_model, choice, mode=mode)
  n.x <- choice[1]
  n.y <- choice[2]
  eig <- res1[[1]]
  xlab <- res1[[2]][1]
  ylab <- res1[[2]][2]

  
 # site, species and cn scores
 score <- scores(ord_model, display = c("sites", "species", "bp", 
               "cn"), scaling = plot.scaling)
 ## now get the bp scores 
 bp.score <- as.data.frame(score$biplot)

 ## biplot scale 
 #mul <- vegan:::ordiArrowMul(bp.score, fill=0.75) 

 # bp.score.sel <- .prune.biplot(dt, meta.new, ord_model, bp.score, 
 #                              biplot_sig=biplot.sig)
 # I have problem using anova to test significance of factors
 # so I'll just plot all.

 bp.score.sel <- bp.score
 # species/sites/centroids scores
 species.score <- as.data.frame(score$species)
 sites.score <- as.data.frame(score$sites)
 cn.score <- as.data.frame(score$centroids)
       
 if( plot.species ) {
   xlim <- range(species.score[, n.x], cn.score[, n.x], 
                      bp.score.sel[, n.x] )
   ylim <- range(species.score[, n.y], cn.score[, n.y], 
                      bp.score.sel[, n.y] )
 } else {
   xlim <- range(sites.score[, n.x], cn.score[, n.x], 
                      bp.score.sel[, n.x] )
   ylim <- range(sites.score[, n.y], cn.score[, n.y], 
                      bp.score.sel[, n.y] )
 }
 
 ### Plot OTUs ###
 meta <- meta[match(rownames(meta.new), rownames(meta)), ]

 # plot levels
 if( plot.species ) {
    df <- as.data.frame(species.score)
    df <- df[match(rownames(OTU), rownames(df)),]
    if ( identical(rownames(df), rownames(OTU)) ) {
       df[[rank_name]]<-OTU[[rank_name]]
    } else {
       stop(" not same otuIDs in rda output and OTU table provided")
    }
 } else if ( !plot.species & !(is.null(group)) ) {
    .valid.factor(meta, group)
    df<-as.data.frame(sites.score)
    df <- df[match(rownames(meta), rownames(df)),]
    if ( identical(rownames(df), rownames(meta)) ) {
       if ( length(group) > 2 ) {
         warning(" currently only support 2 groups for plotting, will only use the first 2 factor groups") 
         group <- group[1,2]
       } else {
         group <- group
       }
       labels <- names(group)
       for ( i in 1:length(group) ) {
         gr <- group[[i]]
         if ( is.null(gr) ) { break }
         label <- names(group)[i]
         df[[label]] <- meta[[gr]]
       }
    } else {
       stop(" not same subjects in rda output and metadata table provided")
    }
    
 } else {
    warning("No levels to plot, provide one group")
 }
 
 #if ( require("grid") ) {
 #  grid::arrow
 #} else {
 #  stop("package 'grid' is required for this function; try install.packages('grid')") # arrow()
# }

 ## ggplot
 # set up vairables for plotting
 if ( mode=="rda" ) {
    x.string <- paste("RDA", n.x, sep="")
    y.string <- paste("RDA", n.y, sep="")
 } else {
    x.string <- paste("CCA", n.x, sep="")
    y.string <- paste("CCA", n.y, sep="")
 }

 shape <- c(16, 17, 15, 3, 7, 8, 1:2, 4:6, 9:15, 18, 48:57) 
 if ( plot.species ) {
   p <- ggplot(data=df) + 
        geom_point(data=df, aes_string(x=x.string, y=y.string, 
                    col=rank_name, shape=rank_name), 
                    size=cex.point)
   len <- length(levels(factor(df[[rank_name]])))
   if (len <= 6) {
     p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=rank_name) + 
             scale_color_discrete(name=rank_name) 
   } else {
     p <- p + scale_shape_identity() + 
              # # define shapes
              scale_shape_manual(name=rank_name, 
                            values=shape[0:len]) + 
              scale_color_discrete(name=rank_name)
   }
 } else {
   if ( length(group) == 1 ) {
     len <- length(levels(factor(df[[labels[1]]])))
     p <- ggplot(data=df) + 
         geom_point(data=df, aes_string(x=x.string, y=y.string, 
                    col=labels[1], shape=labels[1]), 
                    size=cex.point)
     if (len <= 6) {
       p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=labels[1]) + 
             scale_color_discrete(name=labels[1]) 
     } else {
       p <- p + scale_shape_identity() + 
             scale_shape_manual(name=labels[1], 
                            values=shape[0:len]) + 
             scale_color_discrete(name=labels[1])
     }
    } else {
     len1 <- length(levels(factor(df[[labels[1]]])))
     len2 <- length(levels(factor(df[[labels[2]]])))
     len <- max(len1, len2)
     p <- ggplot(data=df) + 
         geom_point(data=df, aes_string(x=x.string, y=y.string, 
                    col=labels[1], shape=labels[2]), 
                    size=cex.point)
     if (len <= 6) {
       p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=labels[2]) + 
             scale_color_discrete(name=labels[1]) 
     } else {
       p <- p + scale_shape_identity() + 
             scale_shape_manual(name=labels[2], 
                            values=shape[0:len]) + 
             scale_color_discrete(name=labels[1])
     }
    }
 }

 if ( len <= 12 ) {
    p <- p + scale_color_brewer(palette="Set1")
 } else {
    col.func  <-  grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    p <- p + scale_color_manual(values=col.func(len))
 }
    
 # biplot arrows
 if ( is.null(biplot.scale) ) {
   bp.score.sel <- as.data.frame(bp.score.sel)
 } else {
   bp.score.sel <- as.data.frame(bp.score.sel*biplot.scale)
 }
    
 if ( biplot.label ) {
    bp.score.sel$label <- rownames(bp.score.sel)
 } else {
    bp.score.sel$label <- 1:nrow(bp.score.sel)
    lab.biplot <- paste(1:nrow(bp.score.sel), 
                        rownames(bp.score.sel), sep=": ")
    print(paste("each number of biplot label represents: ", 
                paste(lab.biplot, collapse="; "), sep=""))
 }
 
 p<-p+ geom_segment(data=bp.score.sel, 
       aes_string(x=0,xend=x.string,y=0,yend=y.string), 
           arrow = arrow(length = unit(0.3, "cm")), 
              colour="grey",inherit_aes=FALSE) 
 p <- p + geom_text(data = bp.score.sel, aes_string(x= x.string, 
          y = y.string, label ="label"), size = cex.bp,  hjust = 0.5) + 
          theme_bw() +
          theme(legend.title = element_text(size = cex.leg, 
                face = "bold", hjust = 0))
               
 p <- p + xlab(xlab) + 
          ylab(ylab) +
          ggtitle(main)

 save <- !is.null(file)
 if (save) {
   .ggsave.helper(file, ext, width, height, plot=p)
 } else {
    print(p)
 }
 return(list(plot=p, ord_model=ord_model, dt=dt, meta.new=meta.new))
}


.subset.otu <- function (otu, rank, taxa) {
  valid.OTU(otu)
  .valid.rank(rank)

  # get identified ranks
  rank_name<- .get.rank.name(rank)
  rank_pat<-.get.rank.pat(rank)
  remove.pat <- .blacklist()
                
  # split taxonomy at rank
  OTU <- tax.split(otu, rank=rank)
  OTU <- OTU[ OTU[[rank_name]] != "", ]
  move <- which(OTU[[rank_name]] %in% remove.pat)
  if ( length(move) == 0 ) {
    OTU <- OTU
  } else {
    OTU <- OTU[-which(OTU[[rank_name]] %in% remove.pat),]
  }
  OTU.tax <- OTU[[rank_name]]

  if ( is.null(taxa) ) {
    sel.taxa <- 1:length(OTU.tax)
  } else {
    if (is.character(taxa)) {
      sel.taxa <- which(OTU.tax %in% taxa)
        
      if (length(sel.taxa) ==0) {
         warning ("no taxa are found in OTUs")
         sel.taxa <- NULL
      } else {
         sel.taxa <- sel.taxa
      }  
    }      
    if (is.numeric(taxa)) {
      # select the OTUs from the most abundant taxa
      OTU.agg <- aggregate(OTU[, -ncol(OTU)], 
                          by=list(OTU[[rank_name]]), FUN=sum)
      OTU.agg <- OTU.agg[order(rowSums(OTU.agg[, -1]), 
                             decreasing=TRUE),]
      #return(OTU.agg)
      if (length(levels(factor(OTU.agg[,1]))) <= taxa) {
         sel<-OTU.agg[, 1]
      } else {
         sel<-OTU.agg[1:taxa, 1]
      }
      sel.taxa <- which(OTU.tax %in% sel)        
    }
  }   
  OTU <-OTU[sel.taxa, ]
  return(OTU)
}


.perform_ord <- function(OTU, meta, factors, 
                          mode=c("rda", "cca")) {
  
  # transpose to sample x taxa 
  OTU_otu <- OTU[, -ncol(OTU)]
  OTU_otu.t <- as.data.frame(t(OTU_otu))
  #if ( count ) {
  #   OTU_otu.t <- OTU_otu.t
  #} else {
  #   OTU_otu.t <- decostand(OTU_otu.t, "total")
  #}
  # check sample # & order 
  if( !identical(rownames(OTU_otu.t), rownames(meta))) {
     stop ("Error: samples doesn't match between data and metadata")
  }
   
  # prepare metadata
  suppressWarnings(meta <- filter.META(meta))  

  if ( is.null(factors) ) {
     factors1 <- names(meta)
  } else {
     factors1 <- .valid.factor(meta, factors)
  }
  meta.new <- meta[, which(names(meta) %in% factors1)]

  # return(meta.new)
  # problem passing dt to the formula
  # error: Error in parse(text = x) : 
  #<text>:2:0: unexpected end of input 1: dt ~ 
  #formula <- paste("dt", paste(factors, collapse=" + "), 
  #                 sep=" ~ ")
  #form <- as.formula(formula)

  if ( length(mode) != 1L || 
    !any(grepl(mode, c("rda", "cca"))) ) {
    warning("Only either 'rda' or 'cca' can be used as mode,
    will use 'rda' as default for analysis")
    mode <- "rda"
  } else {
    mode <- mode
  }
  
  if (mode=="rda") {
    # scale=TRUE would scale species to unit variance, but 
    # since we have standadized data, so this being set FALSE
    dt_ord = vegan::rda(OTU_otu.t ~., data=meta.new, scale=FALSE)
  } else {
    dt_ord = vegan::cca(OTU_otu.t ~., data=meta.new, scale=FALSE)
  }
  #return(list(meta.new, OTU_otu.t, dt_ord))
  # anova test for each factor
  test.anova <- anova(dt_ord, by = "terms", perm = 1000)

  return(list(dt_ord, test.anova))
}              


.prune.biplot <- function(data, meta, ord_model, bp_score, biplot_sig=NULL) {
        
  #prune ordination biplot arrows: include only the environment 
  # constraints that are significant at conf=0.05, 
  # doesn't actually work because for factors with multi-levels, 
  # the rownames of test.anova do not match rownames of biplot 
  # scores, biplot scores are for each level of each factor, but 

  # anova test for each factor
  dt <- data
  meta.new <- meta
  ord_model <- ord_model
  test.anova <- anova(ord_model, by = "terms", perm = 1000)

  ## get ones with p values less than or equal to some threshold 
  if ( is.null(biplot_sig) ) {
    biplot_sig <- 100
  } else {
    biplot_sig <- biplot_sig
  }
  anova.sig <- rownames(test.anova[with(test.anova, 
                        head(`Pr(>F)` <= biplot_sig, -1)),])
  #return(anova.sig)    
  ## now get the bp scores 
  bp_score <- as.data.frame(bp_score)
  ## and take the ones that are signif 
  if( length(anova.sig) == 0 ) {
     warning(paste("No factor had significant impact at significance level of ", biplot_sig, ", will plot all", sep=""))
     sig <- 1:nrow(bp_score)
  } else {
     sig.pat <- paste(anova.sig, collapse="|")
     sig <- which(grepl(sig.pat, rownames(bp_score)))
  }

  bp.score.sel <- bp_score[sig, ] 
  return(bp.score.sel)
}

.eig <- function(ord_model, choice=c(1,2), 
                   mode=c("rda", "cca")) {
    
  ### Choices of axis to plot ###
  choices <- choice
  n.x <- choice[1]
  n.y <- choice[2]

  #extract eigenvalues
  eig<-eigenvals(ord_model)
  #make nice xlab and ylab with "variance explained"
  if (mode=="rda") {
    xlab <- paste("RDA", n.x, " ", 
          sprintf("(%.2f%%)",100*(eig[n.x]/(sum(eig)))), sep="") 
    # sprintf("%f", 
    ylab <- paste("RDA", n.y, " ", 
          sprintf("(%.2f%%)",100*(eig[n.y]/(sum(eig)))), sep="")
 } else {
    xlab <- paste("CCA", n.x, " ", 
          sprintf("(%.2f%%)",100*(eig[n.x]/(sum(eig)))), sep="") 
    ylab <- paste("CCA", n.y, " ", 
          sprintf("(%.2f%%)",100*(eig[n.y]/(sum(eig)))), sep="")
 }
 return(list(eig, labs=c(xlab, ylab)))
}
 

# plot rda analysis for sites and taxa at a given rank
Taxa.ord <- function(data, is.OTU=TRUE, meta=meta, factors=NULL, 
                 group=NULL, rank="g", taxa=10, data.trans="total", 
                 plot.species=TRUE, plot.scaling=-1, 
                 biplot.scale=NULL, biplot.sig=NULL, biplot.label= TRUE, 
                 mode=c("rda", "cca"), choice=c(1,2), main="", 
                 cex.point=3, cex.label=1, cex.leg=12, cex.bp=3,
                 cex.text=3, file=NULL, ext=NULL, width=10, height=10) {


  # plot otus on ordination space.
  # get identified ranks
  rank_name <- .get.rank.name(rank)
  rank_pat <-.get.rank.pat(rank)

  if ( is.OTU ) {
    abund  <-  tax.abund(data, rank=rank, drop.unclassified=TRUE)
  } else {
    abund <- data
  }
 
  abund <- abund[, order(colSums(abund), decreasing=TRUE)]
  
  if (!is.null(data.trans)) {
    abund  <-  decostand(abund, data.trans)
  } else {
    abund <- abund
  }
       
  if ( is.numeric(taxa) && length(taxa) == 1L ) {
    if ( ncol(abund) <= taxa ) {
       abund <- abund
    } else {
       abund <- abund[, 1:taxa]
    }
  } else {
    sel <- which(names(abund) %in% taxa)
    if ( length(sel)== 0 ) {
       warning("selected taxa are not in the data set, will plot the top taxa as default")
       if ( ncol(abund) <=10 ) {
         abund <- abund[, 1:ncol(abund)]
       } else {
         abund <- abund[, 1:10]
       }
    } else {
       abund <- abund[, sel]
    }
  }
  
  tax <- names(abund)
  tax.lvl <- levels(factor(tax))

  # check sample # & order 
  if( !identical(rownames(abund), rownames(meta))) {
     stop ("Error: samples doesn't match between data and metadata")
  }
   
  # prepare metadata
  suppressWarnings(meta <- filter.META(meta))  

  if ( is.null(factors) ) {
     factors1 <- names(meta)
  } else {
     factors1 <- .valid.factor(meta, factors)
  }
  meta.new <- meta[, which(names(meta) %in% factors1)]

  # return(meta.new)
  # problem passing dt to the formula
  # error: Error in parse(text = x) : 
  #<text>:2:0: unexpected end of input 1: dt ~ 
  #formula <- paste("dt", paste(factors, collapse=" + "), 
  #                 sep=" ~ ")
  #form <- as.formula(formula)

  if ( length(mode) != 1L || 
    !any(grepl(mode, c("rda", "cca"))) ) {
    warning("Only either 'rda' or 'cca' can be used as mode,
    will use 'rda' as default for analysis")
    mode <- "rda"
  } else {
    mode <- mode
  }
  
  dt <- abund
  if (mode=="rda") {
    # scale=TRUE would scale species to unit variance, but 
    # since we have standadized data, so this being set FALSE
    ord_model = vegan::rda(dt ~., data=meta.new, scale=FALSE)
  } else {
    ord_model = vegan::cca(dt ~., data=meta.new, scale=FALSE)
  }
  
  #return(list(meta.new, OTU_otu.t, ord_model))
  # anova test for each factor
  #abund <- abund
  #meta.new <- meta.new
  #ord_model <- ord_model
  #test.anova <- anova(ord_model, by = "terms", perm = 1000)
  #return(test.anova)
  # res <- .perform_ord(OTU, count=count, meta=meta, factors=factors, mode=mode)

  choice <- choice
  res1 <- .eig(ord_model, choice, mode=mode)
  n.x <- choice[1]
  n.y <- choice[2]
  eig <- res1[[1]]
  xlab <- res1[[2]][1]
  ylab <- res1[[2]][2]

  
 # site, species and cn scores
 score <- scores(ord_model, display = c("sites", "species", "bp", 
               "cn"), scaling = plot.scaling)
 ## now get the bp scores 
 bp.score <- as.data.frame(score$biplot)

 ## biplot scale 
 #mul <- vegan:::ordiArrowMul(bp.score, fill=0.75) 

 #bp.score.sel <- .prune.biplot(dt, meta.new, ord_model, bp.score, 
 #                             biplot_sig=biplot.sig)
 #return(bp.score.sel)
 # I have problem using anova to test significance of factors
 # so I'll just plot all.

 bp.score.sel <- bp.score
 # species/sites/centroids scores
 species.score <- as.data.frame(score$species)
 sites.score <- as.data.frame(score$sites)
 cn.score <- as.data.frame(score$centroids)
       
 if( plot.species ) {
   xlim <- range(sites.score[, n.x], species.score[, n.x])           
   ylim <- range(sites.score[, n.y], species.score[, n.y]) 
 } else {
   xlim <- range(sites.score[, n.x])           
   ylim <- range(sites.score[, n.y])
 }
 
 ### Plot sites and/or taxa ###

 # plot levels
 if ( !is.null(group) ) {
   .valid.factors(meta, group)
   df<-as.data.frame(sites.score)
   df <- df[match(rownames(meta), rownames(df)),]
   if ( identical(rownames(df), rownames(meta)) ) {
     if ( length(group) > 2 ) {
        warning(" currently only support 2 groups for plotting, will only use the first 2 factor groups") 
        group <- group[1,2]
     } else {
        group <- group
     }
     labels <- names(group)
     len <- list()
     legend.name <- list()
     for ( i in 1:length(group) ) {
       gr <- group[[i]]
       if ( is.null(gr) ) { break }
       label <- names(group)[i]
       df[[label]] <- meta[[gr]]
     }
   } else {
     stop(paste("not same subjects in ", mode, " output and metadata table provided", sep=""))
   }    
 } else {
    warning("No levels to plot, provide one group")
 }
 
 ## ggplot
 # set up vairables for plotting
 if ( mode=="rda" ) {
    x.string <- paste("RDA", n.x, sep="")
    y.string <- paste("RDA", n.y, sep="")
 } else {
    x.string <- paste("CCA", n.x, sep="")
    y.string <- paste("CCA", n.y, sep="")
 }
 
 # plot sites.score
 shape <- c(16, 17, 15, 3, 7, 8, 1:2, 4:6, 9:15, 18, 48:57) 
 if ( length(group) == 1 ) {
   len <- length(levels(factor(df[[labels[1]]])))
   p <- ggplot(data=df) + 
         geom_point(data=df, aes_string(x=x.string, y=y.string, 
                    col=labels[1], shape=labels[1]), 
                    size=cex.point)
   if (len <= 6) {
     p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=labels[1]) + 
             scale_color_discrete(name=labels[1]) 
   } else {
       p <- p + scale_shape_identity() + 
             scale_shape_manual(name=labels[1], 
                            values=shape[0:len]) + 
             scale_color_discrete(name=labels[1])
   }
 } else {
   len1 <- length(levels(factor(df[[labels[1]]])))
   len2 <- length(levels(factor(df[[labels[2]]])))
   len <- max(len1, len2)
   p <- ggplot(data=df) + 
        geom_point(data=df, aes_string(x=x.string, y=y.string, 
                    col=labels[1], shape=labels[2]), 
                    size=cex.point)
   if (len <= 6) {
     p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=labels[2]) + 
             scale_color_discrete(name=labels[1]) 
   } else {
       p <- p + scale_shape_identity() + 
             scale_shape_manual(name=labels[2], 
                            values=shape[0:len]) + 
             scale_color_discrete(name=labels[1])
   }
 }
 
 if ( len <= 12 ) {
    p <- p + scale_color_brewer(palette="Set1")
 } else {
    col.func  <-  grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    p <- p + scale_color_manual(values=col.func(len))
 }
    
 if ( plot.species ) {
   num.species  <-  length(rownames(species.score))
   df.tax <- as.data.frame(species.score)
  
   # create a vector for the vertical position of labels; make half the labels 
   # justified up/down on y-axis, the others constant on y-axis
   v.jitter  <-  sample(c(2.8, 0.5), size=num.species, 
                              replace=TRUE)
   v.jitter  <-  sapply(v.jitter, 
                     FUN=function(x){if (x != 0.5 && runif(1) < 0.5) {x * -1} else {x}})
  
   # for the labels with no y-axis jitter, add x-axis jitter
   h.jitter  <-  sapply(v.jitter, 
                     FUN=function(x){if (x != 0.5) {0.5} else {sample(c(-0.4, 1.4), size=1)}})
  
   # jitter the sample labels (using vectors from above)
   df.tax[["name"]] <- paste0(rownames(df.tax), "")
   p  <-  p + geom_point(data=df.tax, aes_string(x=x.string, y=y.string)) +
              geom_text(data=df.tax, aes_string(x=x.string, y=y.string, 
                        label="name"), size=cex.text, colour="grey50", 
                        vjust=v.jitter, hjust=h.jitter)
 }
 
 # biplot arrows
 if ( is.null(biplot.scale) ) {
   bp.score.sel <- as.data.frame(bp.score.sel)
 } else {
   bp.score.sel <- as.data.frame(bp.score.sel*biplot.scale)
 }
    
 if ( biplot.label ) {
    bp.score.sel$label <- rownames(bp.score.sel)
 } else {
    bp.score.sel$label <- 1:nrow(bp.score.sel)
    lab.biplot <- paste(1:nrow(bp.score.sel), 
                        rownames(bp.score.sel), sep=": ")
    print(paste("each number of biplot label represents: ", 
          paste(lab.biplot, collapse="; "), sep=""))
 }
 
 p<-p+ geom_segment(data=bp.score.sel, 
       aes_string(x=0,xend=x.string,y=0,yend=y.string), 
           arrow = arrow(length = unit(0.3, "cm")), 
              colour="grey",inherit_aes=FALSE) 
 p <- p + geom_text(data = bp.score.sel, aes_string(x= x.string, 
          y = y.string, label ="label"), size = cex.bp,  hjust = 0.5) + 
          theme_bw() +
          theme(legend.title = element_text(size = cex.leg, 
                face = "bold", hjust = 0))
               
 p <- p + xlab(xlab) + 
          ylab(ylab) +
          ggtitle(main) +
          xlim(xlim) +
          ylim(ylim)

 save <- !is.null(file)
 if (save) {
   .ggsave.helper(file, ext, width, height, plot=p)
 } else {
    print(p)
 }
 return(list(plot=p, ord_model=ord_model, dt=dt, meta.new=meta.new))
}


pcoa.plot <- function(data, is.OTU=TRUE, meta, factors, rank, stand.method=NULL,
                      dist.method="morisita", sample.labels=TRUE,
                      top=20, ellipse=FALSE, main=NULL, file=NULL, ext=NULL,
                      height=8, width=10, ggplot2=TRUE, bw=FALSE) {
  
  
  num.facs <- length(factors)
  
  if (!is.numeric(top) || length(top) != 1L || top < 0) {
    stop("'top' must be a numeric vector of length 1 and have a value > 0.")
  }
  
  if (!is.logical(sample.labels) || length(sample.labels) != 1L || 
        is.na(sample.labels)) {
    stop("'sample.labels' must be a logical vector of length 1.")
  }
  
  if (!is.logical(bw) || length(bw) != 1L || is.na(bw)) {
    stop("'bw' must be a logical vector of length 1.")
  }
  
  # we have to check with identical, because in R, TRUE == 1 but 
  # identical(TRUE, 1) is FALSE
  if (!any(identical(ellipse, 1), identical(ellipse, 2), 
           identical(ellipse, FALSE))) {
    stop("'ellipse' should be either 1, 2, or FALSE; see ?pcoa.plot for details.")
  }
  
  if (ellipse > num.facs) {
    stop("argument 'ellipse' cannot be greater than the number of factors.")
  }
  
  if (ellipse && ggplot2) {
    warning("drawing the ellipses for groups is not currently supported when 'ggplot2=TRUE'.")
  }
  
  # get the metadata factors
  meta.factors <- .valid.factors(meta, factors, min.factors=1, max.factors=2)
  
  # if we're plotting 2 factors in b&w, we need to 'cross' the two factors and 
  # assign a unique symbol to each
  if (bw && num.facs >= 2) {
    
    # ggplot complains if we include /, -, " ", etc. in names
    if (ggplot2) {separator <- "_"} else {separator <- "/"}
    
    cross.name <- paste(names(meta.factors), collapse=separator)
    meta.factors <- cbind(paste(meta.factors[ ,1], meta.factors[ ,2], 
                                sep=separator),
                          meta.factors)
    
    names(meta.factors)[1] <- cross.name
  }
  
  for (column in colnames(meta.factors)) {
    if (!is.factor(meta.factors[ ,column])) {
      warning(paste("column", column,
                    "from 'factors' is not a factor; attempting to coerce now", 
                    "(see ?RAM.factors for help)."))
      
      meta.factors[ ,column] <- as.factor(meta.factors[ ,column])
    }
    
    if (length(levels(meta.factors[ ,column])) > 9) {
      warning(paste("there are more than 9 levels in column", column,
              "from 'factor'; only 9 will be shown."))
    }
    
    # to plot ellipses, we need more than 2 counts for each level
    if (any(summary(meta.factors[ ,column]) <= 2) && ellipse) {
      warning(paste("column", column, "from 'factor' has less than two",
                    "observations for a level, this prevents ellipses from",
                    "being plotted."))
    }
  }
  
  if ( is.OTU ) {
    valid.OTU(data)
    if ( is.null(rank) ) {
      rank <- NULL
      abund <- transpose.OTU(data)
    } else {
      rank <- rank
      .valid.rank(rank)
      .valid.meta(otu1=data, meta=meta)
      #otu.t <- transpose.OTU(data)
      abund <- tax.abund(data, rank=rank, drop.unclassified=TRUE)
    }
  } else if ( !is.OTU) {
    rank <- NULL
    abund <- data
    abund <- abund[match(rownames(meta), rownames(abund)),]
    if ( !identical(rownames(meta), rownames(abund)) ) {
      stop("data and metadata do not have same samples")
    }
  }    
  
  if (!is.null(stand.method)) {
    abund <- decostand(abund, stand.method)
  } else {
    abund <- abund
  }
  
  dists <- vegdist(abund, method=dist.method)
  
  k.max <- nrow(abund)-1
  
 # if (!require("labdsv")) {
 #     stop("package 'labdsv' is required to use this function")
 # }
  
  pcoa <- suppressWarnings(pco(sqrt(dists), k.max))
  # sqrt the distance to avoid negative eigenvalues or reduce them
  
  # if we don't get at least two axes of ordination, throw error
  if (dim(pcoa$points)[2] < 2) {
    stop("less than two axes of ordination were calculated. Please try again with different values for 'stand.method' and/or 'dist.method'.")
  }
  
  sp.scores <- wascores(pcoa$points, abund)
  
  if (ggplot2) {
    .pcoa.ggplot2(abund, pcoa, rank, sp.scores, meta.factors, sample.labels, 
                  top, ellipse, main, file, ext, height, width, bw)
  } else {
    .pcoa.base(abund, pcoa, rank, sp.scores, meta.factors, sample.labels, top, 
               ellipse, main, file, ext, height, width, bw)
  }
}

.pcoa.ggplot2 <- function(abund, pcoa, rank, sp.scores, meta.factors, 
                          sample.labels, top, ellipse, main, file, ext, height, width,
                          bw) {
  
  num.facs <- length(meta.factors)
  save <- !is.null(file)
  
  # set up the data frame with pcoa data
  samples <- as.data.frame(cbind(Sample=rownames(pcoa$points), 
                                 X=pcoa$points[ , 1], Y=pcoa$points[ , 2]))
  
  for (i in 1:num.facs) {
    samples <- cbind(samples, meta.factors[[i]])
    names(samples)[i + 3] <- names(meta.factors)[i]
  }

  samples$X <- as.numeric(as.character(samples$X))
  samples$Y <- as.numeric(as.character(samples$Y))
  samples$Sample <- as.character(samples$Sample)
  
  # set up the data frame with taxonomic data
  otus <- as.data.frame(cbind(OTU=rownames(sp.scores), X=sp.scores[ ,1], 
                              Y=sp.scores[ ,2]))
  otus$X <- as.numeric(as.character(otus$X))
  otus$Y <- as.numeric(as.character(otus$Y))
  
  # filter out the top samples
  if (dim(otus)[1] < top) {
    if ( !is.null(rank) ) {
      warning(paste("there are less than", top, 
                  "taxon groups at the given rank; plotting them all."))
    } else {
      warning(paste("there are less than", top, 
                  "taxon groups in the input data; plotting them all.")) 
    }   
    top <- dim(sp.scores)[1]
  }
  
  otus <- otus[1:top, ]
  
  # determine aes based on number of meta factors and bw setting
  # (recall the special case when bw=T and num.facs >=2 since we 'crossed' the
  # factors)
  if (num.facs >= 2) {
    if (bw) {
      samples.aes <- aes_string(x="X", y="Y", label="Sample",
                                shape=names(samples)[4])
    } else {
      samples.aes <- aes_string(x="X", y="Y", label="Sample", 
                                shape=names(samples)[4],
                                colour=names(samples)[5])
    }
  } else if (num.facs == 1) {
    samples.aes <- aes_string(x="X", y="Y", label="Sample", shape=names(samples)[4])
  }

   # add points for samples (with labels)
  p <- ggplot(samples, samples.aes) + geom_point(alpha=0.65, size=7)
  
  shape <- c(16, 17, 15, 3, 7, 8, 1:2, 4:6, 9:15, 18, 48:57) 
  # default shapes in ggplot only 6, need to add more
  if ( num.facs == 1 ) {
    len <- length(levels(factor(samples[[names(samples)[4]]])))
    if (len <= 6) {
      p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=names(samples)[4]) 
      } else {
       p <- p + scale_shape_identity() + 
             scale_shape_manual(name=names(samples)[4], 
                            values=shape[0:len])  
   }
 } else if (num.facs >=2 ) {
   len1 <- length(levels(factor(samples[[names(samples)[4]]])))
   len2 <- length(levels(factor(samples[[names(samples)[5]]])))
   len <- max(len1, len2)
   if (bw) {
     if (len <= 6) {
       p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=names(samples)[4])    
     } else {
       p <- p + scale_shape_identity() + 
             scale_shape_manual(name=names(samples)[4], 
                            values=shape[0:len])  
     }
   } else {
     if (len <= 6) {
       p <- p + scale_shape_identity() + 
             scale_shape_discrete(name=names(samples)[4]) + 
             scale_color_discrete(name=names(samples)[5]) 
     } else {
       p <- p + scale_shape_identity() + 
             scale_shape_manual(name=names(samples)[4], 
                            values=shape[0:len]) + 
             scale_color_discrete(name=names(samples)[5])
     }
   }
 }
 
  num.samples <- length(unique(samples$Sample))
  
  # create a vector for the vertical position of labels; make half the labels 
  # justified up/down on y-axis, the others constant on y-axis
  v.jitter <- sample(c(2.8, 0.5), size=num.samples, replace=TRUE)
  v.jitter <- sapply(v.jitter, 
                     FUN=function(x){if (x != 0.5 && runif(1) < 0.5) {x * -1} else {x}})
  
  # for the labels with no y-axis jitter, add x-axis jitter
  h.jitter <- sapply(v.jitter, 
                     FUN=function(x){if (x != 0.5) {0.5} else {sample(c(-0.4, 1.4), size=1)}})
  
  
 
  
  if (sample.labels) {
    # jitter the sample labels (using vectors from above)
    p <- p + geom_text(size=2, colour="black", vjust=v.jitter, hjust=h.jitter)
  }
    
  if (top != 0) {
    # add taxon groups
    p <- p +
      geom_text(aes_string(x="X", y="Y", label="OTU", colour=NULL, shape=NULL),
                data=otus, size=3, colour="darkgrey", alpha=0.8)
  }

  x.lab <- paste0("Axis I (", round(100 * pcoa$eig[1] / sum(pcoa$eig), digits=2), "%)")
  y.lab <- paste0("Axis II (", round(100 * pcoa$eig[2] / sum(pcoa$eig), digits=2), "%)")
  
  #main.title <- paste("PCoA for Top", top, "Taxon Groups at",
  #                    .get.rank(.get.rank.ind(rank), pretty=TRUE),
  #                     "Level")
  if (is.null(main)) {
        main.title = ""
  } else {
      main.title = main
  }
  
  # add titles
  p <- p + ggtitle(main.title) + xlab(x.lab) + ylab(y.lab)
  
  # add theme
  p <- p + theme(panel.background = element_rect(fill="white"),
                 panel.border = element_rect(colour="black", fill="transparent"))

  # add colour
  if ( len <= 12 ) {
    p <- p + scale_color_brewer(palette="Set1")
  } else {
    col.func  <-  grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    p <- p + scale_color_manual(values=col.func(len))
  }
    
    
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  }
  
  p
}

.pcoa.base <- function(abund, pcoa, rank, sp.scores, meta.factors, 
                       sample.labels, top, ellipse, main, file, ext, height, width,
                       bw) {
  
  .valid.plot.settings(file, ext)
  save <- !is.null(file)
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  num.facs <- min(length(meta.factors), 2)
  
  # ugly hack, if bw=T and we have >=2 factor, we've "crossed" two other factors
  # so manually set num.facs <- 1
  
  if (bw && num.facs >=2) { num.facs <- 1 }
  
  if (save) {
    .get.dev(file, ext, height=height, width=width)
  }
  
  lmat <- matrix(c(rep(1, times=num.facs), 2:(num.facs + 1)), ncol=2)
  layout(lmat, widths=c(9, 1))
  
  # get label names
  x.lab <- paste0("Axis I (", round(100 * pcoa$eig[1] / sum(pcoa$eig), digits=2), "%)")
  y.lab <- paste0("Axis II (", round(100 * pcoa$eig[2] / sum(pcoa$eig), digits=2), "%)")
  # main.title <- paste("PCoA for Top", top, "Taxon Groups at",
  #                    .get.rank(.get.rank.ind(rank), pretty=TRUE),
  #                    "Level")
  
  if (is.null(main)) {
      main.title = ""
  } else {
      main.title = main
  }
  
  plot.args <- list(pcoa$points[ , 1:2], xlab=x.lab, ylab=y.lab, main=main.title,
                    lwd=4, cex=2.5)
  
  par(mar=c(5.1, 4.1, 4.1, 0))
  
  # set up colours for meta factors 
  meta.cols <- vector(length=num.facs, mode="list")
  palettes <- list(RColorBrewer::brewer.pal(9, "Pastel1"), RColorBrewer::brewer.pal(9, "Set1"))
  
  for (i in 1:num.facs) {
    
    meta.cols[[i]] <- palettes[[i]][as.numeric(meta.factors[[i]])]
    names(meta.cols[[i]]) <- as.character(meta.factors[[i]])
  }
  
  # the shapes to be used for bw plotting (we pass these numbers to pch)
  shapes <- c(0:2, 5:6, 15:18)
  
  if (num.facs == 2) {
    
      plot.args <- c(plot.args, pch=21, bg=list(meta.cols[[1]]), 
                     col=list(meta.cols[[2]]))
    
  } else { # only one metadata factor
    if (bw) {
      plot.args <- c(plot.args, pch=list(shapes[as.numeric(meta.factors[[1]])]))
      
    } else {
      plot.args <- c(plot.args, pch=16, col=list(meta.cols[[1]]))
    }
  }
  
  do.call(plot, plot.args)
  
  if (ellipse) {
    
    # if bw=TRUE & num.facs >=2, the first column is the "crossed" factor, which
    # we correct for by adding 1
    if (bw && num.facs >= 2) { ellipse <- ellipse + 1}
    
    # add ellipses with group names as centroids
    for (level in unique(meta.factors[[ellipse]])) {
      
      if (bw) {
        ordiellipse(pcoa, meta.factors[[ellipse]], show.groups = level,
                    col="black", label=TRUE)
        
      } else { # plot in colour
        ordiellipse(pcoa, meta.factors[[ellipse]], show.groups = level,
                    # this call gets the appropriate colour palette, then extracts
                    # the colour corresponding to the value of the factor level
                    col=meta.cols[[ellipse]][names(meta.cols[[ellipse]]) == level],
                    label=TRUE)
      }
    }
  }
  
  if (sample.labels) {
  # plot the sample labels
  text(pcoa$points[ , 1:2], labels=rownames(abund), 
       # randomly place the labels; call the function until this looks nice
       pos=sample(4, size=length(rownames(abund)), replace=TRUE), 
       cex=0.6, offset=1)
  }
  
  # filter out the top samples
  if (dim(sp.scores)[1] < top) {
    if ( !is.null(rank) ) {
      warning(paste("there are less than", top, 
                  "taxon groups at the given rank; plotting them all."))
    } else {
      warning(paste("there are less than", top, 
                  "taxon groups in the input data; plotting them all.")) 
    }   
    top <- dim(sp.scores)[1]
  }
  
  if (top != 0) {
    # plot the taxonomic information
    text(sp.scores[1:top, 1:2, drop=FALSE], labels=rownames(sp.scores)[1:top], 
         col="darkgrey", cex=0.8)
  }
  
  for (i in 1:num.facs) {
    
    # align the legend boxes with the edges of the plot
    if (i == 1) {
      par(mar=c(0, 0, 4.1, 0))
      plot.new()
      dir <- "topleft"
      # pch value for filled dot
      leg.shape <- 16
      
    } else if (i == 2) {
      par(mar=c(5.1, 0, 0, 0))
      plot.new()
      dir <- "bottomleft"
      # pch value for filled dot
      leg.shape <- 21
    }
    
    leg.args <- list(dir, legend=levels(meta.factors[[i]]), 
                     title=names(meta.factors)[i], xpd=NA)
    
    # if black and white, use shapes for legend, otherwise use colour
    if (bw) {
      leg.args <- c(leg.args, pch=list(shapes[1:length(levels(meta.factors[[i]]))]))
      
    } else {
      leg.args <- c(leg.args, pch=leg.shape, bg="white", pt.lwd=3, pt.cex=1.5,
                    col=list(unique(meta.cols[[i]])))
    }
    
    do.call(legend, leg.args)
  }
  
  if (save) {
    dev.off()
  }
  
  invisible()
}

