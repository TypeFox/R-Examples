valid.taxonomy <- function(data) {
  prefix <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

  if ( class(data) != "list" ) {
    stop("please provide input data as a list, see ?RAM.input.formatting")
  }

  labels <- names(data)
  alert.list <- list()
  for ( i in 1:length(data) ) {
    otu <- data[[i]]
    if ( is.null(otu) ) { break }
    label <- names(data)[i]
    # check whether missing prefix for each rank
    valid.OTU(otu)
    # need to convert all taxonomy to character, NOT factor, 
    # otherwise, gsub wouldn't work
    otu$taxonomy <- as.character(otu$taxonomy)
    alert <- vector()
    for ( i in 1:nrow(otu) ) {
       # remove the ; at the end of the lineage
       if ( grepl(";[ ]+$", otu[i,"taxonomy"], perl=TRUE) ) {         
           otu[i,"taxonomy"] <- gsub(";[ ]+$", "", otu[i,"taxonomy"], 
                                   perl=TRUE)
       } else if ( grepl("[;]+$", otu[i,"taxonomy"], perl=TRUE) ) { 
          otu[i,"taxonomy"] <- gsub("[;]+$", "", otu[i,"taxonomy"], 
                                       perl=TRUE)
       } else {
           otu[i,"taxonomy"] <- otu[i,"taxonomy"]
       }
        # now check # of '; ' and ';'
      count1 <- sapply(regmatches("; ", gregexpr("; ", 
                          otu[i,"taxonomy"])), length) 
      count2 <- sapply(regmatches(";", gregexpr(";", 
                          otu[i,"taxonomy"])), length)
      if ( count1 == 0 && count2 ==0 ) {
         # only kingdom information avaiable 
         # otuID  k__fungi
         next
      } 
      if (count1 == 0 && count2 > 0 ) {
         # otu was separated by ';', not '; '
         alert <- c(alert, "taxonomy lineages are not properly formatted, please check ?RAM.rank.formatting; taxonomic ranks should be separated by '; ', i.e. ';' and a white space")
      } 
      if ( count2 >=7 ) {
         alert <- c(alert, "RAM accept 7 taxonomic ranks, see ?RAM.rank.formatting; the otu table has more than 7 taxonomic ranks") 
      }
      miss_pre <- vector()
      miss_rk <- vector()
      for ( j in prefix ) {
         if ( !grepl(j, otu[i,"taxonomy"]) ) {
              miss_pre <- c(miss_pre, j)
              miss_rk <- c(miss_rk, .get.rank.name(j))
           }
       }
       if ( length(miss_pre) != 0 ) {
         alert<- c(alert, paste("missing prefix or missing taxonomy ranks in the lineages; 1) if missing prefix, consider reformatting the taxonomy; 2) if missing ranks only, it's probably ok" ))
       }
     }
     len <- unique(alert)
     alert <- paste(unique(alert), collapse="\n")
     if ( len != 0 ) {
        y <- "Consider reformatting the OTU table's taxonomy, check ?RAM::reformat.taxonomy"
        warning(paste0(paste("For ", label, ": ", sep=""), alert, "\n", y, sep="; "))
     }
 
   }
}


.capitalize.vec <- function(vec) {
  output <- vector()
  for ( i in vec ) {
    output <- c(output, .capitalize(i))
  }
  return(output)
}


reformat.taxonomy <- function(data, input.ranks=NULL, sep="; ") {
  
  if ( class(data) != "list" ) { 
    stop("provide all otu tables as list. See ?RAM.input.formatting")
  }
  new.otu <- list()
  labels <- names(data)
  # process each otu
  for (i in 1:length(data) ) {
    label <- names(data)[i]
    otu <- data[[i]]
    if (is.null(otu)) {
      break
    }
    valid.OTU(otu)
    ls <- list()
    ls[[label]] <- otu
    valid.taxonomy(data=ls)
    
    # default outputs for RAM taxonomy
    prefix <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    names <- c("kingdom", "phylum", "class", "order", "family",
               "genus", "species")
    otu$taxonomy <- as.character(otu$taxonomy)
    
    # reformat 
    # check whether the 'input.ranks include names
    input.ranks <- tolower(input.ranks)
    if ( is.null(input.ranks) || !any(input.ranks %in% names) ) {
      stop("please provide ALL taxonomic ranks in the input OTU tables, see ?reformat.taxonomy")
    } 
    
    # split the current taxonomy column
    # remove the last ";" at the end of line
    otu$taxonomy <- gsub("[; ]+$", "", otu$taxonomy, perl=TRUE)
    suppressWarnings(otu.split <- col.splitup(otu, col="taxonomy", 
                                              sep=sep, drop=TRUE,
                                              names=input.ranks))
    
    # select columns only belong to 
    n <- which(colnames(otu.split) %in% names)
    n.names <- intersect(colnames(otu.split), names)
    n.otu <- ncol(otu)
    otu.split <- otu.split[, c(1:(n.otu-1), n)]
    
    # replace the 'unclassified' taxonomy annotation to proper format
    for (i in n.names) {
      rank_pat <- .get.rank.pat(i)
      rank_pat2 <- paste0(rank_pat, rank_pat)
      otu.split[[i]] <- paste0(substring(i,1,1), "__", 
                               otu.split[[i]])
      otu.split[[i]] <- gsub(rank_pat2, rank_pat, otu.split[[i]])
    }
    # combine to 'taxonomy' column
    n.split <- ncol(otu.split)
    for (i in 1:nrow(otu.split) ) {
      row<-vector();
      n <- n.otu:n.split
      for (j in n) {
        row<-c(row, otu.split[i,j])
      }
      otu.split[i, "taxonomy"] <- paste(row, collapse="; ");
      otu.split[i, "taxonomy"] <- gsub(";[ ]+$", "", 
                                       otu.split[i,"taxonomy"], 
                                       perl=TRUE)
    }  
    otu.split <- otu.split[, c(1:(n.otu-1), ncol(otu.split))]
    new.otu[[label]] <- otu.split
  }      
  return(new.otu)
}
  

get.rank <- function(otu1, otu2=NULL, rank=NULL) {
  single.otu <- is.null(otu2)
  single.rank <- !is.null(rank)
  valid.OTU(otu1, otu2)
  
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", 
                   "genus", "species")
  
  # if given both otu and otu2, call get.rank for both
  if (!single.otu) {
    
    output <- list()
    output$otu1 <- get.rank(otu1=otu1, rank=rank)
    # this looks ugly, but we are just calling get.rank with a single
    # OTU argument (which is named otu2)
    output$otu2 <- get.rank(otu1=otu2, rank=rank)
    
   # if (single.rank) {
   #   names(output)$otu1 <- .get.rank(.get.rank.ind(rank))
   #   names(output$otu2) <- .get.rank(.get.rank.ind(rank))
   # } else {
   #   names(output$otu1) <- tax.classes
   #   names(output$otu2) <- tax.classes
   # }
    
    return(output)
  }
  
  if (single.rank) {
    .valid.rank(rank)
    
    # select only the rows with the given taxonomic prefix in their taxonomy
    pat <- .get.rank.pat(rank)
    output <- otu1[grep(pat, otu1$taxonomy), ]
    
    # now only select the ones that are NOT on the blacklist
    remove.pat <- paste0(.blacklist(rank), "|no_taxonomy", paste0("|", pat,"$"))
    
    output <- output[!grepl(remove.pat, output$taxonomy, ignore.case = TRUE, perl=TRUE), ]
    
    if (dim(output)[1] == 0) {
      warning(paste("no OTUs classified at the", .get.rank(.get.rank.ind(rank)),
        "level."))
    }
    return(output)
    
  } else { # call get.rank for each taxonomic classification
    output <- vector(mode="list", length=length(tax.classes))
    
    for (i in tax.classes) {
      output[.get.rank.ind(i)] <- list(get.rank(otu1=otu1, rank=i))
    }
    
    names(output) <- tax.classes
    return(output)
  }
}

tax.split <- function(otu1, otu2=NULL, rank=NULL) {
  valid.OTU(otu1, otu2)
  single.otu <- is.null(otu2)
  single.rank <- !is.null(rank)
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # if given otu and otu2, call tax.split for both
  if (!single.otu) {
    
    output <- list()
    
    output$otu1 <- tax.split(otu1=otu1, rank=rank)
    output$otu2 <- tax.split(otu1=otu2, rank=rank)
    
    return(output)
  }
  
  if (single.rank) {
    # get the index for rank (also does data validation for rank)
    .valid.rank(rank)
    tax.ind <- .get.rank.ind(rank)
  } 
  
  # split OTU table taxonomy information into individual columns
  suppressWarnings(otu1.split <- col.splitup(otu1, col="taxonomy", sep="; ", 
                   max=length(tax.classes), names=tax.classes, drop=TRUE))
  
  # columns from 1 to max contain the numeric data, the other taxonomic
  max <- dim(otu1)[2] - 1
  
  # we need seven taxonomy columns; if one is missing (because nothing classified
  # at that level), add empty column
  
  # while we have less columns than we should...
  #while (dim(otu1.split)[2] < max + length(tax.classes)) {
    # ...add a column filled with empty strings
  #  otu1.split <- cbind(otu1.split, rep("", times=dim(otu1.split)[1]))
  #}
  
  # strip the 'formatting' bits of taxonomic info
  otu1.split[ ,-(1:max)] <- gsub("k__|p__|c__|o__|f__|g__|s__|;", "", 
        as.matrix(otu1.split[ ,-(1:max)]))
  
  if (single.rank) {
    # add the single taxonomic column to the numeric data, return that data frame
    return(otu1.split[ , names(otu1.split) %in% c(names(otu1)[1:max], tax.classes[tax.ind])])
    
  } else {
    # set up list for output
    otu1.taxa <- vector("list", length(tax.classes))
    names(otu1.taxa) <- tax.classes 
    
    for (i in 1:length(tax.classes)) {
      # creating a list of data frames is surprisingly difficult in R;
      # if you do not use the list(df1, df2, ...) syntax, you get a list composed
      # of the first data frame in its entirety, followed be the individual columns
      # of the other data frames. Instead we wrap each data frame in a list itself 
      # before we add it, then we can access them with [[]]
      otu1.taxa[i] <- list(otu1.split[ , names(otu1.split) %in% c(names(otu1)[1:max],
                   tax.classes[i])])
    }
    
    return(otu1.taxa)
  }
}

tax.abund <- function(otu1, otu2=NULL, rank=NULL,            
                      drop.unclassified=FALSE, 
                      top=NULL, count=TRUE, mode="number") {
  
  single.otu <- is.null(otu2)
  valid.OTU(otu1, otu2)
  single.rank <- !is.null(rank)
  # data validation for top is done later in the function (when the dimensions of 
  # the taxonomy matrix are known)
  filter <- !is.null(top)
  
  if (!is.logical(drop.unclassified) || length(drop.unclassified) != 1L) {
    stop("argument 'drop.unclassified' must be a logical vector of length one.")
  }
  
  if (!is.logical(count) || length(count) != 1L) {
    stop("argument 'count' must be a logical vector of length one.")
  }
  
  if (!is.character(mode) || !any(mode %in% c("number", "percent"))) {
    stop("argument 'mode' must be one of 'number' or 'percent'.")
  }
  
  # if given otu and otu2, call tax.abund for both
  if (!single.otu) {
    
    drop <- drop.unclassified
    output <- list()
    
    output$otu1 <- tax.abund(otu1, rank=rank, drop.unclassified=drop, top=top,
         count=count, mode=mode)
    output$otu2 <- tax.abund(otu2, rank=rank, drop.unclassified=drop, top=top,
         count=count, mode=mode)
    
    return(output)
  }
  
  # get the OTU table in proper format
  if (single.rank) {
    .valid.rank(rank)
    tax.list <- list(tax.split(otu1, rank=rank))
  } else {
    tax.list <- tax.split(otu1)
  }
  
  for (i in seq(along=tax.list)) { 
    
    # update taxonomy label to "taxonomy"
    names(tax.list[[i]])[dim(tax.list[[i]])[2]] <- "taxonomy" 
    # aggregate the otu table by taxonomy rank names 
    tax.list[[i]] = stats::aggregate(tax.list[[i]][ , -dim(tax.list[[i]])[2]],  
          by = list(tax.list[[i]]$taxonomy), FUN = .sumcol) 
    # change row names to header provided by aggregate
    rownames(tax.list[[i]]) <- tax.list[[i]][ , 1] 
    # remove first column (that information is now in the row names)
    tax.list[[i]] <- tax.list[[i]][ , -1]
    # transpose table (the as.data.frame generates the V1 heading) 
    tax.list[[i]] <- as.data.frame(t(tax.list[[i]])) 
    
    # if count is false, return relative abundance
    if (!count) {
      tax.list[[i]] <- vegan::decostand(tax.list[[i]], method="total")
    }
    
    # order the table by column sums
    tax.list[[i]] <- tax.list[[i]][ , order(colSums(tax.list[[i]]), 
            decreasing = TRUE), drop=FALSE] 
    
    # remove all zero entries
    tax.list[[i]] <- tax.list[[i]][ , colSums(tax.list[[i]]) > 0, drop=FALSE]
    
    # rename the "V1" column
    if (!single.rank) {
      names(tax.list[[i]])[names(tax.list[[i]]) == "V1"] <- 
    paste("unclassified_below_", .get.rank(i - 1), sep="")
    } else {
      # if we are only processing one element, we cannot use the i index 
      # to get the correct rank
      names(tax.list[[i]])[names(tax.list[[i]]) == "V1"] <- 
    paste("unclassified_below_", .get.rank(.get.rank.ind(rank) - 1), sep="")
    }
    
    # remove unclassified columns
    if (drop.unclassified) {
      # this selects all columns NOT containing in the blacklist
      remove.pat <- gsub(.get.rank.pat(.get.rank(i)), "", paste0(.blacklist(.get.rank(i)), "|no_taxonomy"))
      
      tax.list[[i]] <- tax.list[[i]][ , !grepl(remove.pat, names(tax.list[[i]]),
               ignore.case=TRUE), 
             drop=FALSE]
    }
    
    # keep only the 'top' most abundant groups, where top is user-given 
    if (filter) {
      tax.list[[i]] <- .select.top(tax.list[[i]], top, count, mode)
    }
  }

  tax.out <- list()
  index <- 1
  for ( i in 1:length(tax.list) ) {
    tax <- tax.list[[i]]
    if ( is.null(tax) ) { break }
    lab <- names(tax.list)[i]
    if (is.null(lab)) {
      lab1 <-index
    } else {
       lab1 <- lab
    }
    names(tax) <- gsub(" ", "_", names(tax))
    tax.out[[lab1]] <- tax
    index <- index + 1
  }
  
  if (single.rank) {
    return(tax.out[[1]])
    
  } else {
    return(tax.out)
  }
}

.select.top <- function(abund, top, count, mode) {
  if (!is.numeric(top) || length(top) != 1) {
    stop("argument 'top' must be a numeric vector of length one.")
  }
  
  if (!is.character(mode) || !any(mode %in% c("number", "percent"))) {
    stop("argument 'mode' must be one of 'number' or 'percent'.")
  }
  
  if (top <= 0) {
    stop("argument 'top' must be greater than zero.")
  }
  
  # take top X samples when mode == "number"
  if (mode == "number") {
    num.groups <- dim(abund)[2]
    
    if (top > num.groups) {
      warning("argument 'top' is greater than the total number of groups, returning all groups.")
      top <- num.groups
    }
    
    abund <- abund[ ,1:top, drop=FALSE]
    
    
  # find all samples with RA > X% when mode == "percent"
  } else if (mode == "percent") {
    
    if (top > 100) {
      stop("argument 'top' must be in the range (0, 100].")
    }
    
    percent <- top / 100
    
    # we need the 'count' parameter to determine whether we need to standardize
    # the data ourselves or not
    if (count) {
      abund.stand <- vegan::decostand(abund, "total")
    } else {
      abund.stand <- abund
    }
    
    abund.max <- apply(abund.stand, MARGIN=2, FUN=max)
    # find samples w/ max relative read abundance < 'top'% and remove
    exclude <- names(which(abund.max < percent))
    
    # if all groups have been excluded
    if (length(exclude) == dim(abund)[2]) {
      stop("no taxon groups with relative abundance above 'top' percent.")
    }
    
    # select the samples above 'top'%
    abund <- abund[ ,-(which(names(abund.max) %in% exclude)), drop=FALSE]
  }
  
  abund
}

tax.fill <- function(data, downstream=TRUE) {
  
  valid.OTU(data)
  
  if (length(downstream) != 1L || !is.logical(downstream) || is.na(downstream)) {
    stop("'downstream' must be a character vector of length one.")
  }
  
  if (downstream) {
    range <- 1:7
    belong <- "belongs_to_"
  } else {
    range <- 7:1
    belong <- "classified_to_"
  }
  
  taxonomy <- strsplit(data$taxonomy, "; ")
  taxonomy.length <- length(taxonomy)
  
  # initialize the vector of classified names
  classified.names <- rep("no_taxonomy", times=taxonomy.length)
  # initialize list for fixed taxonomy
  new.taxonomy <- vector(mode="list", length=taxonomy.length)
  
  for (i in range) {
    # this selects the i-th element of every list and returns a character vector
    current.groups <- unlist(lapply(taxonomy, FUN="[", i))

    # any entries that are TRUE (i.e. match blacklist) need to be replaced
    # we add that pattern to match any group prefix with no name 
    # (since we split on semicolons the rank argument in .blacklist won't work)
    blacklist <- paste0(paste(.blacklist(), collapse="|"), "|[kpcofgs]__$", sep="")
    
    replace <- is.na(current.groups) | grepl(blacklist, current.groups,
             ignore.case = TRUE)
    
    num.groups <- length(current.groups)

    for (j in 1:num.groups) {
      
      # if we need to replace the classification, use last good name
      if (replace[j]) {
    
          new.name <- paste0(.get.rank.pat(.get.rank(i)), belong,
           gsub("__", "_", classified.names[j]))
    
          new.taxonomy[[j]][i] <- new.name
      } else { # if the classification is good, store it
          classified.names[j] <- current.groups[j]
          new.taxonomy[[j]][i] <- classified.names[j]
      }
    }

  }
  
  new.taxonomy <- lapply(new.taxonomy, FUN=paste, collapse="; ")
  
  # return a new data frame with the updated taxonomy
  # (do not mutate the old data frame)
  data.frame(data[ ,-dim(data)[2], drop=FALSE], taxonomy=unlist(new.taxonomy),
     stringsAsFactors = FALSE)
}


LCA.OTU <- function(otu, strip.format=FALSE, drop=TRUE) {

  valid.OTU(otu)
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", 
                   "genus", "species")

  # split taxonomy column by ;
  suppressWarnings(otu.split <- col.splitup(otu, col="taxonomy", sep="; ", 
                 max=length(tax.classes), names=tax.classes, drop=TRUE))
  # otu table column#  
  max <- ncol(otu) - 1
  
 # strip the 'formatting' bits of taxonomic info
 blacklist<-vector()
 for ( i in c("k__", "p__", "c__", "o__", "f__", "g__", "s__") ) {
     blacklist<-c(blacklist,paste(i, .blacklist(), sep=""))
 }
 blacklist <- paste(blacklist, collapse="|")
 blacklist <- paste(blacklist, "|k__|p__|c__|o__|f__|g__|s__|;", sep="")

# remove entries matched blacklist (unclassified taxa)
  otu.split[ ,-(1:max)] <- gsub(blacklist, "", 
        as.matrix(otu.split[ ,-(1:max)]))

  otu.split[ ,-(1:max)] <- gsub("^[ ]+", "", 
       as.matrix(otu.split[ ,-(1:max)]), perl=TRUE)

  # obtain LCA of each OTU
  # replace "" to NA
  otu.split[ otu.split == "" ] <- NA
  

  # whether or not strip formatting: e.g. g__
  otu.split$LCA <- apply(otu.split, 1, function(x) tail(unlist(x [ 
                        which(!is.na(x)) ]), n=1))

  if ( strip.format ) {
      otu.split$LCA <- otu.split$LCA
  } else {
      n.r <- nrow(otu.split)
      n.c <- ncol(otu.split)
      for ( i in 1:n.r ) {
           num <- which(otu.split[i, 1:(n.c-1)]==otu.split$LCA[i])
           num <- num[length(num)]
           pat <- .get.rank.pat(names(otu.split)[num])
           otu.split$LCA[i] <- paste(pat, otu.split$LCA[i], 
                                  sep="")
      }
  }
  
  # whether or not remove taxonomy split ranks
  output <- otu.split[, setdiff(colnames(otu.split),"taxonomy")]
  if ( drop ) { 
      # keep only LCA for taxonomy info
      output <- cbind(otu[, 1:max], output$LCA)
      names(output)[ncol(output)] <- "LCA"
  } else { 
      output <- cbind(output, otu$taxonomy)
      names(output)[ncol(output)] <- "taxonomy"
  }
  return(output)
}


col.splitup <- function(df, col="", sep="", max=NULL, names=NULL, drop=TRUE) { 
    # validate all inputs
    if ( sep == "" ) {
        stop(paste("\n", "    Error: separator ?  check ?col.splitup  ", "\n", sep=""))
    }
    if ( col == "" )  {
       stop(paste("\n", "    Error: column to split? check ?col.splitup  ", "\n", sep=""))
    }
    if ( !(col %in% names(df)) ) {
        stop(paste("\n", "    Error: column to be split is not in the dataframe", "\n", sep=""))
    }  
    
    # col position in df
    if ( col %in% names(df) ) {
        num.col <- which( names(df) %in% col )
    }
    
    # split the column to list;
    list <- strsplit(df[[col]], sep, fixed=FALSE); 
    vec <- vector();

    # determine max number of split columns to be remained in output
    for (i in 1:length(list) ) {
      vec<-c( vec, length(list[[i]]) );
    }

    def.max <- max(c(max, length(names)))
    maximum <- c( max, max(vec), length(names) )
    max <- max(maximum)
    
    if ( max(vec) > def.max ) {
#        warning(paste("\n", "    ", col, " can be split to: ", max(vec), " coloumns; ", 
#                 "\n", "    required to be split to ", def.max, " columns; ", "\n", 
#                  "    Will KEEP all ", max(vec), " and IGNORE user defined ",  def.max, 
#                 " columns", sep=""))
     } else if ( max(vec) < def.max )  {
            #warning(paste("\n", "    ", col, " can be split to: ", max(vec), 
#                " columns; ", "\n", "    required to be split to: ", max, 
 #               " columns; ", "\n", "    column names provided: ", length(names), "\n", 
 #               "    Will fill empty strings in the additional ", def.max-max(vec), 
 #              " column(s)", sep=""))
      } else {
        #warning(paste("\n", "    ", col, " can be split to: ", max(vec), " coloumns; ", 
#                 "\n", "    required to be split to ", def.max, " columns; ", "\n", 
#                  "    ", col, " will be split to: ", max(vec), " columns", sep=""))
     }
   

    for ( i in 1:length(list) ) {
      # since max is equal to or larger than length(list[[i]]) as defined by section above
      if ( length(list[[i]]) < max ) {
        x=rep( "", max-length(list[[i]]) );
        list[[i]] <- c(list[[i]], x)
      } else {
        list[[i]] <- list[[i]]
      }
    } 
    
    # rbind to form a data frame of split columns    
    new<-as.data.frame( do.call("rbind", list) )
    
    # names of new columns
    if ( is.null(names) ) {
       new.name <- colnames(new)
    } else {
      if ( length(names) == max ) {
        new.name <- names
      } else if ( length(names) < max ) {
      #warning(paste("\n", "    ", col, " being split to: ", max, 
      #          " columns;", "\n", "    column names provided: ", length(names), "; ", "\n",
      #          "    will only change the first ", max, " of split columns", sep=""))
        new.name <- c(names, colnames(new)[(length(names)+1):max])
      # } 
      # since max is the maximum of length(names), pre-defined max and max(vec), so 
      # following condition is not possible
      #  else if ( length(names) > max ) {
      #     warning(paste("\n", "    remained number of split columns from ", col, " is: ", max, 
      #          ";", "\n", " number of names provided: ", length(names), "; ", "\n",
      #         "     will ignore the last ", length(names) -max, " of names", sep=""))
      #  new.name <- names[1:max]
      } else {
        new.name <- colnames(new)
      }
    }
    
    colnames(new) <- new.name    
    # for data.table dt[, setdiff(colnames(dt),col), with=FALSE] 
    if ( ! drop ) {
        #warning(paste("    Keep ", col, " column in output!", sep=""))
        df.new <- cbind(df, new)
    } else {
        #warning(paste("    Drop ", col, " column in output!", sep=""))
        df.new <- cbind(df[, setdiff(colnames(df),col)], new)
    }
    return(df.new)
}
    
transpose.LCA <- function(data) {
  .valid.data(data)
  labels <- names(data)

  lca.t.list <- list()  
  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    otu <- data[[i]]
    if ( is.null(otu) ) { break }
    
    lca <- LCA.OTU(otu, strip.format=FALSE, drop=TRUE)
    lca <- lca[order(rowSums(lca[,-ncol(lca)]), decreasing=TRUE),]
    rownames(lca) <- paste(lca$LCA, rownames(lca), sep="_")
    lca <- lca[, -ncol(lca)]
    lca.t <- as.data.frame(t(lca))
    lca.t.list[[label]] <- lca.t
  }
  return(lca.t.list)
}

