# all functions with names beginning with a period are NOT exported for the
# end user's use, but are available internally. Removing the leading period
# and rebuilding will make them available, or you can edit the NAMESPACE file
# in the root folder (which determines what functions are exported).

# open the appropriate device, given an extension
.get.dev <- function(file, ext, height=8, width=16, dpi=1000) {
  .valid.ext(ext)
  file <- .ensure.filepath(file, ext)
  
  args <- list(file=file, height=height, width=width, res=dpi, bg="white",
               units="in")
  
  if (ext == "tiff") {
    args <- c(args, compression="lzw")
  }
  
  if (ext == "pdf" || ext == "svg") {
    args$res <- NULL
    args$units <- NULL
  }
  
  # this looks strange, but since the functions to open the devices are just called
  # pdf, jpeg, ... themselves, we simply parse our ext value to get the appropriate 
  # function.
  expr <- parse(text=ext)
  
  # evaluate the expression (which gives a function); call with other arguments
  do.call(eval(expr), args)
}

.get.rank.ind <- function(rank) {
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax.classes.plural <- c("kingdoms", "phyla", "classes", "orders", "families", "genera", "species")
  tax.classes.short <- c("k", "p", "c", "o", "f", "g", "s")
  tax.classes.pattern <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax.classes.all <- c(tax.classes, tax.classes.plural, tax.classes.short, tax.classes.pattern)
  
  # we do not validate the rank here, as some methods call this 
  # in a for loop when rank=NULL
  
  # convert to upper case as ignore.case in grepl
  # the length call below should return 7 (barring massive taxonomical discoveries...)

   val <- unique(which(toupper(rank)==toupper(tax.classes.all)) %% length(tax.classes))
  # since we took the index mod 7, all species values were given index 0
  # however, we want them to have index 7:
   if (val == 0) {
      val <- 7;
    }
 
  
  return(val)
}


.get.rank <- function(ind, pretty=FALSE) {
  # we do not validate the rank here, as some methods call this 
  # in a for loop when rank=NULL
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  if (pretty) {
    pretty.rank <- tax.classes[ind]
    # capitalize the first letter (for plots)
    pretty.rank <- .capitalize(pretty.rank)
    pretty.rank
    
  } else {
    tax.classes[ind]
  }
}

.get.rank.pat <- function(rank) {
    # we do not validate the rank here, as some methods call this
    # in a for loop when rank=NULL
    tax.classes.pattern <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    
    # get the rank index and return the appropriate pattern
    return(tax.classes.pattern[.get.rank.ind(rank)])
}


.get.rank.name <- function(rank, plural=FALSE, pretty=FALSE) {
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax.classes.plural <- c("kingdoms", "phyla", "classes", "orders", "families", "genera", "species")
  tax.classes.short <- c("k", "p", "c", "o", "f", "g", "s")
  tax.classes.pattern <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax.classes.all <- c(tax.classes, tax.classes.short, tax.classes.pattern, tax.classes.plural)
  
    # we do not validate the rank here, as some methods call this 
    # in a for loop when rank=NULL
  
    # the length call below should return 7 (barring massive taxonomical discoveries...)
    val <- which(rank==tax.classes.all) %% length(tax.classes)
   
    # since we took the index mod 7, all species values were given index 0
    # however, we want them to have index 7:
    if (val == 0) {
      val <- 7;
    }
  
  if (!isTRUE(plural)) {
    val.name<-tax.classes[val]
  } else {
    val.name<-tax.classes.plural[val]
  }

  if (!is.na(val.name)) {
      if (pretty) {
          # capitalize the first letter (for plots)
          val.name <- .capitalize(val.name)       
      } else {
          val.name <- val.name
      }
      return(val.name)
  }
}


# capitalize the first letter of a string; works for vectors as well
.capitalize <- function(string) {
  gsub("(.)(.*)", "\\U\\1\\E\\2", string, perl=TRUE)
}

.sumcol <- function(x){
  x <- as.matrix(x)
  
  if (dim(x)[1] > 1) {
    return(colSums(x))
  } else {
    return(x)
  }
}

.ensure.filepath <- function(file, extension) {
  if (!is.character(file) || length(file) != 1L) {
    stop("file must be a character vector of length 1.")
  }
  
  if (!is.character(extension) || length(extension) != 1L) {
    stop("extension must be a character vector of length 1.")
  }
  
  # if the given file path is just a name, (i.e. contains no slashes); 
  # prepend the working directory
  if(!grepl("[\\\\\\/]", file)) {
    file <- paste(getwd(), "/", file, sep="")
  }
  
  # split at any / or \
  dirs <- strsplit(file, "\\\\|\\/")[[1]]
  path <- paste(dirs[-length(dirs)], collapse="/")
  file.name <- dirs[length(dirs)]
  
  if (!file.exists(path)) {
    stop(paste("the given directory \'", path, "\' does not exist.", sep=""))
  }
  
  # if the path ends without a file name; e.g. test/directory/
  # (note that strsplit will not add a "" item to the end of the list in this case)
  if (grepl("[\\\\\\/]$", file)) {
    stop(paste("please provide a valid file name (after the last slash in \'", file, "\').",
               sep=""))
  }
  
  # if the given file name does not end in the right extension
  pat <- paste(extension, "$", sep="")
  if (!grepl(pat, file, ignore.case=TRUE)) {
    file <- paste(file, ".", extension, sep="")
  }
  
  return(file)
}

.factor.settings <- function(row.factor, col.factor, hmap.args, leg.args) {
  given.rfactor <- !is.null(row.factor)
  given.cfactor <- !is.null(col.factor)
  
  if (given.rfactor) {
    # create a palette with a color for each level
    # it raises an error when n is 2, which we silence; however we do need
    # to catch the case where there are > 12 levels
    rfactor.name <- names(row.factor)
    row.factor <- as.matrix(row.factor)
    
    if (!is.factor(row.factor)) {
      warning("row.factor is not a factor; coercing it to one now (see ?factor).")
      row.factor <- as.factor(row.factor)
    }
    
    rfactor.numcols <- length(levels(row.factor))
    
    if (rfactor.numcols > 9) {
      col.func  <-  grDevices::colorRampPalette(brewer.pal(9, "Pastel1"))
      rfactor.pal <- col.func(rfactor.numcols)
      # warning("the meta factor palette can only contain 8 colours, and there are more than 8 levels; dropping some levels.")
    } else {
      rfactor.pal <- suppressWarnings(RColorBrewer::brewer.pal(length(levels(row.factor)), "Pastel1"))
    }
    # get a vector with the appropriate colour at each index
    rfactor.cols <- rfactor.pal[as.numeric(row.factor)]
    # sort the colours to group based on metadata
    rfactor.cols <- rfactor.cols[order(row.factor)]
    hmap.args$RowSideColors <- rfactor.cols
    
    # if there is no column factor, make the legend now
    if (!given.cfactor) {
      leg.args$title <- rfactor.name
      leg.args$legend <- levels(row.factor)
      leg.args$fill <- rfactor.pal[1:rfactor.numcols]
    }
  }
  
  if (given.cfactor) {
    # create a palette with a color for each level
    # it raises an error when n is 2, which we silence; however we do need
    # to catch the case where there are > 12 levels
    cfactor.name <- names(col.factor)
    col.factor <- as.matrix(col.factor)
    
    if (!is.factor(col.factor)) {
      warning("col.factor is not a factor; coercing it to one now (see ?factor).")
      col.factor <- as.factor(col.factor)
    }
    
    cfactor.numcols <- length(levels(col.factor))
    
    if (cfactor.numcols > 9) {
      col.func  <-  grDevices::colorRampPalette(brewer.pal(9, "Pastel1"))
      cfactor.pal <- col.func(cfactor.numcols)
      # warning("the meta factor palette can only contain 8 colours, and there are more than 8 levels. Some colours are being recycled in your plot (be careful!).")
    } else {
      # create a palette with a color for each level
      cfactor.pal <- suppressWarnings(brewer.pal(length(levels(col.factor)), "Accent"))
    }
    # get a vector with the appropriate colour at each index
    cfactor.cols <- cfactor.pal[as.numeric(col.factor)]
    # sort the colours to group based on the row.factor, if given
    if (!given.rfactor) {
      cfactor.cols <- cfactor.cols[order(col.factor)]
    } else {
      cfactor.cols <- cfactor.cols[order(row.factor)]
    }
    
    hmap.args$ColSideColors <- cfactor.cols
    
    # if there is no row factor, make the legend now
    if (!given.rfactor) {
      leg.args$title <- cfactor.name
      leg.args$legend <- levels(col.factor)
      leg.args$fill <- cfactor.pal[1:cfactor.numcols]
    }
  }
  
  if (given.rfactor && given.cfactor) {
    leg.args$title <- paste(rfactor.name, "&", cfactor.name)
    leg.args$legend <- c(levels(row.factor), levels(col.factor))
    # the default palette size is 3, so we need to only select the elements of the 
    # palette we need
    leg.args$fill <- c(rfactor.pal[1:rfactor.numcols], cfactor.pal[1:cfactor.numcols])
  }
  
  return(list(hmap.args, leg.args))
}

.get.tax.group <- function(data, rank, group) {
  valid.OTU(data)
  group.len <- length(group)
  
  if (!is.character(rank)) {
    stop("rank must be a character vector.")
  }
  
  if (!is.character(group)) {
    stop("group must be a character vector.")
  }
  
  if (length(rank) != group.len) {
    stop("rank and group must be vectors of equal length.")
  }
  ### add validation to ensure no duplicates
  
  rows <- vector(length=group.len, mode="list")
  
  for (i in 1:group.len) {
    search.pat <- paste0(.get.rank.pat(rank[i]), group[i])
    
    hits <- grepl(search.pat, data$taxonomy, ignore.case=TRUE)
    
    if (length(which(hits)) == 0) {
      stop(paste0("no results found for ", search.pat, "."))
    }
    
    rows[[i]] <- data[hits, , drop=FALSE]
  }
  
  do.call(rbind, rows)
}

.ggsave.helper <- function(file, ext, width, height, plot) {
  
  file <- .ensure.filepath(file, ext)
  
  ggsave.args <- list(filename=file, width=width, height=height, plot=plot,
                      units="in", dpi=1000)
  
  if (ext == "tiff") {
    ggsave.args <- c(ggsave.args, compression="lzw")
  }
  
 # it does not work because ggsave wants an object of class ggplot, while you're passing a grob. arrangeGrob will sometimes trick ggsave in pretending inheritance from ggplot, but only when at least one of the grobs belongs to this class; here, however, you're only passing a gtable. Perhaps the easiest workaround is to clone ggsave and bypass the class check,

  ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
  do.call(ggsave, ggsave.args)
}

# this function returns the blacklist of group names we consider
# "missing" (it essentially acts as a global variable)
.blacklist <- function(rank = NULL) {
  
  blacklist <- c("unclassified", "unidentified", "incertae_sedis", 
                 "incertae sedis", "unassignable")
    
  # example: if we "phylum" as the rank, we want to add "p__;" to the blacklist
  if (!is.null(rank)) {
    .valid.rank(rank)
    
    # convert whatever the rank is to its pattern
    rank <- .get.rank.pat(rank)
    blacklist <- paste(paste(rank, .blacklist(), collapse="|", sep=""), 
                       paste0(rank, ";"), sep="|")
  }
  
  blacklist
}

.ram.pal <- function(n) {
  col_pal <- c("red4", "darkslategray3", "dodgerblue1", "darkcyan",
               "gray79", "black", "skyblue2", "dodgerblue4",
               "purple4", "maroon", "chocolate1", "bisque3", "bisque",
               "seagreen4", "lightgreen", "skyblue4", "mediumpurple3",
               "palevioletred1", "lightsalmon4", "darkgoldenrod1")
  ram.pal <- col_pal[1:n]
  return(ram.pal)
}


 
.varTostr <- function(var) {
  deparse(substitute(var))
}

.group.rank<-function(otu, meta=meta, meta.factor="", drop.unclassified=FALSE, relative.abund=FALSE, rank="g", top=10){
  tax <- tax.abund(otu, rank=rank, drop.unclassified=FALSE, count=TRUE)
  #return(tax)
  if(all(rownames(tax)==rownames(meta))) {
    tax.factor <- stats::aggregate(tax, by=list(meta[[meta.factor]]), FUN=sum)
    rownames(tax.factor)<-tax.factor[,"Group.1"]
    tax.factor <- tax.factor[,-1]
  } else {
    stop("Error: otu and metadata have different samples")
  }
  #return(tax.factor)
  tax.factor <- tax.factor[, order(colSums(tax.factor), decreasing=TRUE)]
  if(!isTRUE(relative.abund)) {
    tax.factor <- tax.factor
  } else {
    tax.factor <- vegan::decostand(tax.factor, "total")
    
  }

  if (drop.unclassified) {
      # this selects all columns NOT containing in the blacklist
      # drop.unclassified using blacklist
     remove.pat <- gsub(.get.rank.pat(rank), "", paste0(.blacklist(.get.rank.pat(rank)), "|no_taxonomy"))
          tax.factor <- tax.factor[ , !grepl(remove.pat, names(tax.factor), ignore.case=TRUE), drop=FALSE]
    } else {
    tax.factor<-tax.factor
  }
  # keep only the 'top' most abundant groups, where top is user-given 

  if ( top > ncol(tax.factor) ) {
    warning(paste("There are ", ncol(tax.factor), " taxa at ", .get.rank.name(rank), " will return all!", sep=""))
    top <- ncol(tax.factor)
  } else {
    top <- top
  }
  tax.factor<-tax.factor[,1:top]

    return(tax.factor)    
}


