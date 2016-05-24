.valid.meta.NB <- function (meta) {
  suppressWarnings(meta.sel<-filter.META(meta))
  suppressWarnings(col.fac <- .fac.col(meta.sel))
  error<-vector()
  
  i=0
  if( length(col.fac) > 1 ) {
    for (j in 1:length(col.fac)) {
      for (m in 2:length(col.fac)) {
        if(names(meta.sel)[col.fac[j]]==names(meta.sel)[col.fac[m]]) {
          next
        } else {
          tab<-table(meta.sel[,col.fac[j]], meta.sel[,col.fac[m]])
      
          if ( is.infinite(sum(log(tab))) ) {
            i=i+1
            error[i]<-paste(names(meta.sel)[col.fac[j]], " & ", names(meta.sel)[col.fac[m]], sep="")
          }
        }
      }
      if (j==length(col.fac)) {
        break
      }
    }
  }
  if(length(error)==0) {
    error<-NULL
    return(meta.new=meta.sel)
  } else {
    warning("missing records in the following meta factor pairs, 
    suggest to remove one in each pair and then re-run: ")
    warning(paste("\n", unique(na.omit(error)), collapse="; ")  )
    return(meta.new=meta.sel)
  }
}


assist.NB <- function(data, meta, is.OTU=TRUE, rank=NULL, 
                      meta.factors=NULL, anov.fac=NULL, 
                      taxon="") {

  #if (!require("MASS")) {
  #  stop("package 'MASS' is required to use this function: try 'install.packages('MASS')'.")
  #}
 
  meta.new <- .valid.meta.NB(meta)
  # make sure only one level of data.new
  if ( !is.null(rank) && length(rank) != 1L ) {
     warning("rank should be one length, multiple ranks were      
     provided, will only process the first rank")
     rank <- rank[1]
     .valid.rank(rank)
  } else {
     rank <- rank
     .valid.rank(rank)
  }
        
  if ( length(taxon) != 1L || !is.character(taxon) ) {
     stop("Error: provide either an otuID or taxon name to test")
  }

  if ( is.OTU ) {
     valid.OTU(data)
     rownames(data) <- factor(rownames(data))
     if ( any(taxon %in% rownames(data)) && !is.null(rank) ) { 
         warning("Test an otuID, will ignore the ranks")
         data.new <- data.revamp(list(otu=data), is.OTU=is.OTU, ranks=NULL, 
                       stand.method=NULL, top=NULL)
     } else {
         data.new <- data.revamp(list(otu=data), is.OTU=is.OTU, ranks=rank, 
                       stand.method=NULL, top=NULL)
     }
  } else {
     data.new <- data.revamp(list(taxa=data), is.OTU=is.OTU, 
                             ranks=rank, 
                             stand.method=NULL, top=NULL)
  }
   
  data.new <- data.new[[1]]
  #names(data.new) <- gsub("[ _]+", "", names(data.new))
  data.new <- data.new[match(rownames(meta.new), 
                       rownames(data.new)) ,]

  # if taxon is an otuID, since data.revamp renames each otuID
  # with otuID_LCA, need to find out the changed name
  if ( is.null(taxon) || taxon=="" ) {
      stop("must provide an otuID or taxon name for test")
  } else {
      taxon <- as.character(taxon)
      pat <- paste("_", taxon, sep="")
      if ( any(taxon %in% rownames(data) ) ) {              
          sel <- which(grepl(pat, names(data.new)))
          taxon.n <- names(data.new)[sel]
      } else if ( any(taxon %in% names(data.new)) ) {
          taxon.n <- taxon
      } else {
          stop("the taxon is not in the data")
      }
  }
  # combine the taxon and metadata
  if ( identical(rownames(meta.new), rownames(data.new)) ) {
      sel <- which(names(data.new) == taxon.n)
      tax.met <- cbind(data.new[, sel], meta.new)
      names(tax.met) <- c(taxon.n, names(meta.new))
  } else {
      stop("Error: data and metadata do not have same subjects")
  }

  #return(tax.met)
  if ( is.null(meta.factors) ) {
    factors <- names(meta.new)
    formula<-paste(taxon.n, paste(factors,collapse=" + "), sep=" ~ ")
    #return(list(taxon.n, formula))
    m_nb<-MASS::glm.nb(as.formula(formula), data=tax.met)
#   m_nb<-glm.nb(noquote(taxon) ~ noquote(paste(names(meta.sel),
#             collapse=" + ")), data=tax.met)
  } else {
    factors <- names(meta.new)[which(names(meta.new) %in% 
                                 meta.factors) ]
    warnings(paste("meta.factors has ", length(meta.factors), 
    " variables; among which, ", setdiff(meta.factors, factors), 
    " were not in metadata table"))
    if ( length(factors) != 0 )      
        formula<-paste(noquote(taxon.n), paste(factors, collapse=" + "), sep="~")
        m_nb<-MASS::glm.nb(as.formula(formula), data=tax.met)
  }
  plot(m_nb)
  m_nb_summary<-summary(m_nb)

  # factors that has significant impact
  pvalue=NULL
  pv <- as.data.frame(summary(m_nb)$coef[, "Pr(>|z|)"])
  names(pv)[1] <- "pvalue"
  pv.sel <- subset (pv, pvalue<=0.05 )
  # test output show all levels of each variable
  # but we only need the matching metadata variable names
  factor.sel <- vector()
  for ( i in names(meta.new) ) {
      find <- grep(i, rownames(pv.sel))
      if ( length(find) == 0 ) {
          next 
      } else {
          factor.sel <- unique(c(factor.sel, i))
      }
  }
                                 
  if ( length(factor.sel) == 0 ) {
      warning("No metadata variable was tested significant in 
      NB model")
      factor <- factor.sel
  } else {
      factor <- factor.sel
  }
#  factor.sel<-paste(rownames(pv.sel),collapse="',' ")
#  factor.sel<-paste(rownames(pv.sel),collapse=" + ")
  
  if( !is.null(anov.fac) && is.character(anov.fac) ) {
    anov <- list()
    for ( i in anov.fac ) {
        m2 <- paste(".", "-", i, sep=" ")
        formula<-paste(".", m2, sep="~")
        m_nb_update <- update(m_nb, as.formula(formula))
        m_nb_anova <- anova(m_nb, m_nb_update)
        anov[[i]] <- m_nb_anova
    }
    return(m_NB <- list(NB.model=m_nb, tax.met=tax.met, 
                       taxon=taxon.n, factors=factor, anova=anov))
  } else {
    return(m_NB <- list(NB.model=m_nb, tax.meta=tax.met, 
                       taxon=taxon.n, factors=factor))
  }
}


.newdf.NB <- function(data, num.col=NULL, length.out=100) {
  # This function is to create new dataframe for NB plotting
  # data is the 3rd element of the output list from function
  # assist.NB 
  # e.g 
  # data(ITS1, meta)
  # data <- assist.NB(ITS1, is.OTU=TRUE, taxon="141562", 
  #              rank="g", meta=meta)[[3]]
  
  # check how many metadata variables
  var <- names(data)[2:ncol(data)]
  if ( length(var) != 0 ) {
      meta <- data[, 2:ncol(data)]
      cols.fac <- var[.fac.col(meta)]
      cols.num <- var[.num.col(meta)]
  } else {
      stop("The data does not contain metadata variables for plotting") 
  }
  
  len.num <- length(cols.num)
  len.fac <- length(cols.fac)
  len.all <- len.num + len.fac

  # type of variables for plotting
  if ( len.all == len.fac ) {
      # all factor variables
      mode <- "all_fac"
  } else if ( len.all == len.num ) {
      # all numeric varibles
      mode <- "all_num"
  } else {
      # mixed variables
      mode <- "mix"
  }
  #return(list(len.fac, len.num, len.all, mode))
  # set initial replicate number for each variable as 1
  rep.number <- 1

  if ( len.fac != 0 ) {

      # first, find out the max replicates levels required 
      # for each factor variables
      # we will use this to create an empty dataframe for newdata
      # of factor variables
      for ( i in cols.fac ) {
          level.num <- length(levels(factor(data[[i]])))
          rep.number <- rep.number*level.num
      }
      #print(rep.number)
      
      # create a empty dataframe for factor variables
      df.fac = data.frame(matrix(vector(), 
                             rep.number*length.out, 
                             length(cols.fac)))

      # determine how many times each factor variable needs to 
      # be replicated, so that we will have a combination of each 
      # factor variable in the new df.
      levels<-numeric()
      levels[1] <- 1
      
      for (i in 1:length(cols.fac) ) {
          # number of levels of a factor
          level.num <- length(levels(factor(data[[cols.fac[i]]])))
          # numbers to replicate for next factor will 
          # be levels[i] * level.num
          levels[i+1] <- level.num*levels[i]
          # names of each replicate
          labels <- levels(factor(data[[cols.fac[i]]]))
          # numer of each level for currentl factor 
          each <- (length.out*rep.number)/(level.num*levels[i])
          warning(paste("factor: ", cols.fac[i] ,"; level.num: ", 
          level.num, "; labels: ", paste(labels, collapse=" "), 
          "; each: ", each, sep=""))
          df.fac[,i] <- factor(rep(1:level.num, 
                                   each = each), 
                                   levels = 1:level.num, 
                                   labels = labels)  
         if ( i == length(cols.fac) ) {
             break  
         }
         names(df.fac)<-c(cols.fac) 
      }
  } else {
      df.fac = df.fac
      #names(df.fac)<-c(cols.fac)
  } 
  #return(df.fac)
  # numeric variables
  if ( len.num != 0) {
      
      # replicate number should be inherited from after construct 
      # df.fac
      df.num = data.frame(matrix(vector(), 
                        rep.number*length.out, length(cols.num)))

      # determine whith numeric factor to be plotted
      if ( is.null(num.col) ) {
          # don't define which numeric variable to be ploted
          # will plot the first numeric variable in line of 
          # cols.num, and hold other cols.num    
          # set a range of first cols.num
          sel.rep <- cols.num[1] 
          sel.const <- setdiff(cols.num, sel.rep)          
      } else if ( !is.null(num.col) && 
                  any(num.col %in% cols.num) ) {
          sel.rep <- num.col
          sel.const <- setdiff(cols.num, sel.rep)
      } else {
          sel.rep <- cols.num[1] 
          sel.const <- setdiff(cols.num, sel.rep)  
          warning(paste(num.col, " is not in the data; will 
          ignore and use ", sel.rep, " for plotting", sep=""))
      }
      
      for ( i in 1:length(cols.num) ) {
          df.num[,1] <- rep(seq(from = min(data[[sel.rep]]), 
                             to = max(data[[sel.rep]]), 
                             length.out = length.out), 
                             rep.number)
          names(df.num) <- sel.rep
          if ( length(sel.const) != 0 ) {
              for ( i in 1:length(sel.const) ) {
                  df.num[, i+1] <- mean(data[[sel.const[i]]], na.rm=T)
                  names(df.num)[i+1] <- sel.const[i]
              }
          }
      } 
  } 
  #else {
  #   df.num = df.num
  #}
  
  if ( mode == "all_fac" ) {
     df <- df.fac
  } else if ( mode == "all_num" ) {
     df <- df.num
  } else if ( mode == "mix" ) {
    df<-cbind(df.fac, df.num)
  }
  # write.csv(df, file="test.csv", quote=FALSE)
  return(df)
  
}  

envis.NB <- function(NB.model="", tax.meta,  taxon="", 
                x="", num.col=NULL, group=NULL, group.order=NULL, 
                xlab=NULL, ylab=NULL, fill=NULL, facet=NULL, 
                file=NULL, ext=NULL, width=8, height=8) {
  
  save <- !is.null(file)
  
  # check fac & num variables
  var <- names(tax.meta)[2:ncol(tax.meta)]
  if ( length(var) != 0 ) {
      cols.fac <- var[.fac.col(tax.meta[, -1])]
      cols.num <- var[.num.col(tax.meta[, -1])]
  } else {
      stop("no data for plotting") 
  }
  
  len.num <- length(cols.num)
  len.fac <- length(cols.fac)
  len.all <- len.num + len.fac

  # type of variables for plotting
  if ( len.all == len.fac ) {
      # all factor variables
      mode <- "all_fac"
  } else if ( len.all == len.num ) {
      # all numeric varibles
      mode <- "all_num"
  } else {
      # mixed variables
      mode <- "mix"
  }

  #if ( len.num != 0 ) {
  #    if ( num.col == "" || is.null(num.col) ) { stop("Must choose one numeric variable as one predictor, it should be set as num.col, can also be as x axis")
  #    } else {
  if ( is.numeric(tax.meta[, x]) && !identical(x, num.col)) {
warning(paste("num.col=", num.col, " has been set as the predictor, 
will be used as x axis for plotting; ", "will ignore x=", x, sep=""))
    x <- num.col
  } else {
    x <- x
  }
     
  # create the new data for NB plotting
  suppressWarnings(newdf_NB <- .newdf.NB(tax.meta, num.col=num.col, length.out=100))
  
  # calculate predicted value
  df2 <- cbind(newdf_NB, predict(NB.model, newdf_NB, 
               type = "link", se.fit=TRUE)) 
         # type="link" gives logged scale of type="response"
  df2[[taxon]] <- exp(df2$fit)
  df2$LL <- exp(df2$fit - 1.96 * df2$se.fit)
  df2$UL <- exp(df2$fit + 1.96 * df2$se.fit)
  # return(df2)
  # reorder factor levels for plot
  if( !is.null(group) ) {
      if ( any(group %in% names(df2) ) ) {
          if ( !is.null(group.order) && 
                identical(sort(group.order), 
          sort( levels( factor( df2[[group]] ) ) ) ) ) {
              
              df2[[group]] <- factor(df2[[group]], c(group.order))
              warning(paste("Changed the ", group, " order as: " , 
              paste(group.order, collapse=", "), sep=""))
          } else {
              df2[[group]] <- factor(df2[[group]])
          }
      } else {
         stop(paste(group, " is not in data for plotting; will 
          ignore", sep=""))
      }        
  } 
  #return(list(cols.fac, cols.num, newdf_NB,df2))

  ## Plot
  # x & y labels
  if ( !is.null(xlab) ) {
    xlab <- xlab
  } else {
    xlab <- x
  }

  if ( !is.null(ylab) ) {
    ylab <- ylab
  } else {
    ylab <-paste("Predicted ", taxon, " (count)", sep="")
  }

  if ( is.null(x) || x =="" || !any(x %in% names(df2)) ) {
     stop(paste("must provide a predictor as x axis for plotting, ", 
      "it should be one of the following: ", 
      paste(names(tax.meta)[2:ncol(tax.meta)], collapse=", "), 
      sep=""))
  } else {
     if ( is.numeric(df2[,x]) ) {
          plot <- "x_num"
      } 
      if ( is.factor(df2[, x]) || is.character(df2[, x]) ) {
          plot <- "x_fac"
      }
  } 


  if ( plot == "x_num" ) {
    if ( is.null(fill) ) {
      p <- ggplot(df2, aes_string(x=x, y=taxon)) +
           geom_line(size = 2) + 
           labs(x = xlab, y = paste("Predicted ", taxon, sep=""))         
    } else {
      p <- ggplot(df2, aes_string(x=x, y=taxon, fill=fill)) +
           geom_line(aes_string(colour = fill), size = 2) +
           geom_ribbon(aes_string(ymin = "LL", ymax = "UL", 
                                   #fill = fill), alpha = .25) +
                                   fill = fill), alpha = 0.25) +
           labs(x = xlab, y = paste("Predicted ", taxon, sep="")) 
    }     

  }
  
  if ( plot == "x_fac" ) {
    if ( is.null(fill) ) {
         p <- ggplot(df2, aes_string(x=x, y=taxon)) + 
              geom_boxplot(postition="dodge") 
    } else {
         p <- ggplot(df2, aes_string(x=x, y=taxon, col=fill)) +
              geom_boxplot(postition="dodge", 
                           outlier.colour = "black", 
                           outlier.shape = 16, 
                           outlier.size = 6, 
                           notch = TRUE, 
                           notchwidth = 0.5,) 
            #   geom_jitter()
    }
  }
  
  if( !is.null(facet) ) {
      facet.new <- vector()
      for ( i in facet ) {
          if ( is.factor(df2[[i]]) || 
               is.character(df2[[i]]) ) {
             df2[[i]] <- factor(df2[[i]])
             facet.new <- unique(c(facet.new, i))
          }
       }
       
      if( length(facet.new) == 1) {
          ncol <- length(levels(tax.meta[[facet.new]]))
      } else {
          gg.ncol <- numeric()
          gg.ncol[1] <- 1
          for (i in 1:length(facet.new) ) {
            gg.ncol[i+1] <- gg.ncol[i]*
                   length(levels(factor(df2[[facet.new[i+1]]])))
          }
          ncol <- max(gg.ncol)
      }
  } else {
      facet.new <- NULL
      ncol <- 1
  }
  warning(paste("Divide plots to ", ncol, " columns", sep=""))

  #return(ncol)

  if ( length(facet.new) == 0 ) {
      p <- p
  } else {
      if ( length(facet.new) == 2  ) {
          formula<-paste(facet.new, collapse=" ~ ")
          p <- p + facet_grid(as.formula(formula))
      } else { 
          formula<-paste("", paste(facet.new, collapse="+"), sep="~")
          p <- p + facet_wrap(as.formula(formula), ncol=ncol)
      }
  }

  # add xlab, ylab
  p <- p + labs(x = xlab, y = ylab)
  # align xlab 
  p <- p + theme(axis.text.x = element_text(angle = 45, 
                                vjust = 1, hjust=1)) 

  if (save) {
    file <- .ensure.filepath(file, ext)
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
  
}
