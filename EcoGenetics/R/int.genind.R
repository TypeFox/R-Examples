
####################---INT.GENIND AND GENDATA CLASSES----#######################
# A product of the genind class for internal processing of ecogen class.
#                    ----------o----------

# The present code is a modified version of that written by T. Jombart
# for genind class of the package Adegenet. 
# The internal genind object is used only as an interface between 
# data frames <-> ecogen objects, and A slot <-> G slot,
# without S4 methods or validations. Additional slots were added for handling
# coding/missing data information and were removed the slots loc.nall, loc.names
# and ind.names. The int.genind objects are temporal (they only exist as
# intermediate stages of the information in several 
# processes): they are are unfolded 
# into a frequency matrix and a class "int.gendata" with pure information 
# about the matrix with the function "int.genind2gendata".
# These two objects are joined together in temporal operations with the function
# "int.gendata2genind".
# The information of int.gendata objects, used in internal processing,
# is stored in the slot INT, invisible at user level.
# There are also defined: an int.genind constructor, and 
# an importation - exportation mechanism to data frame (both also modified
# versions of the original ones). 
# With the exception of the first importation of the genetic data
# into ecogen information, the importation-exportation mechanism 
# is used in the internal processes of some ecogen methods (e.g., resizing). 
# The conversion into allelic frequencies preserves options of previous 
# Adegenet versions, and the ploidy of the input data must be unique. 
# Leandro Roser, July 2015.

#------------------------------------------------------------------------------#
#---------------------------CLASSES DEFINITIONS--------------------------------#
#------------------------------------------------------------------------------#
#' int.gendata
#' @name int.gendata-class
#' @slot loc.fac locus of each allele
#' @slot all.names alleles names
#' @slot  ploidy ploidy
#' @slot type type of data ("codominant" or "dominant")
#' @slot NA.char NA character
#' @slot sep separator
#' @slot ncod number of digits coding each allele (codominant data)
#' @slot  missing missing character
#' @slot  missing cells character
#' @slot removed.image removed columns (coded as 2) or rows (coded as 1)
#' @keywords internal

setClass("int.gendata", 
         representation(loc.fac = "factorORnull",
                        all.names = "characterORnull",
                        ploidy = "numeric",
                        type = "character",
                        NA.char = "character",
                        sep = "characterORnull",
                        ncod = "intORnull",
                        missing = "character",
                        missing.cells = "intORnull",
                        removed.image = "list"))


#' int.genind
#' @name int.genind-class
#' @param tab allelic frequency matrix
#' @keywords internal

setClass("int.genind", 
         representation(tab = "matrix"),
         contains =  "int.gendata")


#------------------------------------------------------------------------------#
#-TWO CONVERSORS TO UNFOLD/ JOIN A FREQUENCY MATRIX AND ITS INT.GENDATA OBJECT-#
#------------------------------------------------------------------------------#

#' int.genind2gendata
#' @param input int.genind object
#' @keywords internal

int.genind2gendata <- function(input) {
  out <- new("int.gendata")
  out@loc.fac <- input@loc.fac
  out@all.names <- input@all.names
  out@ploidy <- input@ploidy
  out@type <-  input@type
  out@NA.char <- input@NA.char
  out@sep <-   input@sep
  out@ncod <- input@ncod
  out@missing <-  input@missing
  out@missing.cells <- input@missing.cells
  out@removed.image <-  input@removed.image
  out
}

#' int.gendata2genind
#' @param tab frequency matrix
#' @param input int.gendata object
#' @keywords internal

int.gendata2genind <- function(tab, input) {
  out <- new("int.genind")
  out@tab <- as.matrix(tab)
  out@loc.fac <- input@loc.fac
  out@all.names <- input@all.names
  out@ploidy <- input@ploidy
  out@type <-  input@type
  out@NA.char <- input@NA.char
  out@sep <-   input@sep
  out@ncod <- input@ncod
  out@missing <-  input@missing
  out@missing.cells <- input@missing.cells
  out@removed.image <-  input@removed.image
  out
}

#------------------------------------------------------------------------------#
#--------------------------INT.GENIND CONSTRUCTOR------------------------------#
#------------------------------------------------------------------------------#
#' constructor
#' @keywords internal

int.genind <- function(X, 
                       ploidy = 2, 
                       type = c("codominant", "dominant"),
                       NA.char = "NA",
                       missing = c("0", "NA", "MEAN"),
                       missing.cells = integer(0),
                       sep = "",
                       ncod = NULL,
                       removed.image) {
  
  
  #----GENERAL CONFIGURATION---------------------------------------------------#
  
  type <- match.arg(type)
  missing <- match.arg(missing)
  
  if(class(X) != "data.frame" && class(X) != "matrix") {
    stop("X is not of class <data.frame> or <matrix>")
  }
  
  # X configuration
  X <- as.matrix(X, rownames.force = TRUE)
  nind <- nrow(X)
  nloc <- ncol(X)
  
  # labels configuration 
  ## column names configuration
  temp <- colnames(X)
  temp <- gsub("[.][^.]*$", "", temp)
  loc.names <- aue.rmspaces(temp)
  loc.names <- unique(loc.names)
  
  ## ind.names configuration. If duplicated or no names
  ## present, create generic labels
  X <- int.check.rownames(X)
  ind.names <- rownames(X)
  
  
  # MARKER SPECIFIC CONFIGURATION----------------------------------------------#
  
  #------------------------------
  # CODOMINANT MARKERS
  #------------------------------
  if(type == "codominant"){
    
    ## loc.fac
    loc.nall <-  table(temp)[match(loc.names, names(table(temp)))]
    loc.nall <- as.integer(loc.nall)
    loc.fac <- factor(rep(loc.names, loc.nall), levels = loc.names)
    
    ## alleles name
    temp <- colnames(X)
    temp <- gsub("^.*[.]","",temp)
    all.names <- aue.rmspaces(temp)
    names(all.names) <- loc.fac
    # all.names is now a vector. It can be splitted into a list with:
    # all.names <- split(as.vector(all.names), loc.fac)
    # all.names <- all.names[loc.names]
    
    # END CODOMINANT---------------
    
    
    #------------------------------
    # DOMINANT MARKERS
    #------------------------------
  } else { # end if type=="codominant" <=> if type=="dominant"
    loc.fac <- as.factor(loc.names)
    all.names <- NULL
  }
  # END DOMINANT----------------
  
  # OUTPUT CREATION------------------------------------------------------------#
  
  res <- new("int.genind")
  res@tab <- X
  res@loc.fac <- loc.fac
  res@all.names <- all.names
  
  if(ploidy < as.integer(1)) {
    stop("ploidy inferior to 1")
  }
  res@ploidy <- as.integer(ploidy)
  
  res@type <- as.character(type)
  res@NA.char <- NA.char
  res@sep <- sep
  res@ncod <- as.integer(ncod)
  res@missing <- missing
  res@missing.cells <- as.integer(missing.cells)
  res@removed.image <- removed.image
  
  res
  
} 

#------------------------------------------------------------------------------#
#---------------------------IMPORTER-------------------------------------------#
#------------------------------------------------------------------------------#
#' importer
#' @keywords internal

int.df2genind <- function(indata, 
                          sep = "", 
                          ncod = NULL,
                          NA.char = "NA",
                          ploidy = 2, 
                          type = c("codominant","dominant"),
                          missing = c("0", "NA", "MEAN"),
                          rm.empty.ind = FALSE,
                          poly.level = 5) {
  
  
  
  # DATA CHECK-----------------------------------------------------------------#
  
  type <- match.arg(type)
  missing <- match.arg(missing)
  
  # check ploidy
  if(ploidy < 1L) {
    stop("ploidy can not be less than 1")
  }
  
  # check format and type congruence
  X.format <- class(indata)
  
  if(!any(c("data.frame","matrix") %in% X.format)) {
    stop("data has not a valid format (<data.frame> or <matrix>).
         Check the class of your data.")
  }
  
  #if data has not rows, return an empty int.genind object
  if(dim(indata)[1] == 0) {
    return(new("int.genind"))
  } 
  
  #if data has not columns but n rows, return an int.genind object with 
  # a matrix of dimension n x 0
  if(dim(indata)[2] == 0) {
    temp <- new("int.genind")
    temp@tab <- matrix(nrow = dim(indata)[1], ncol = dim(indata)[2])
    rownames(temp@tab) <- rownames(indata)
    rownames(temp@tab) <- rownames(indata)
    return(temp)
  }
  
  
  # GENERAL INPUT CONFIGURATION------------------------------------------------#
  
  # X configuration
  X <- as.matrix(indata, rownames.force = TRUE)
  mode(X) <- "character"
  X <- aue.rmspaces(X)
  
  nind <- nrow(X)
  nloc <- ncol(X)
  
  # check row and column names
  X <- int.check.rownames(X)
  X <- int.check.colnames(X)
  loc.names <- colnames(X)
  ind.names <- rownames(X)
  
  
  # NA CONFIGURATION---------------------------------------------------------#
  
  # NA translated into character
  if(is.na(NA.char)) {
    X[is.na(X)] <- "NA"
    NA.char <- "NA"
  } 
  
  # find NA strings
  NA.list <- unlist(lapply(unique(ploidy), 
                           function(nrep) {
                             paste(rep(NA.char, nrep),  collapse = sep)
                           }))
  
  NA.list <- unique(c(NA.list, NA.char))
  
  # replace NAs
  X[X %in% NA.list] <- NA
  
  
  # HANDLING MISSING DATA----------------------------------------------------#
  # REMOTION OF DATA IS REVERSIBLE
  
  # image of the columns and rows that will be removed tagged as 1 and 2
  removed.temp <- is.na(X)
  mode(removed.temp) <- "numeric"
  removed.image <- removed.temp - removed.temp
  col.rem <- apply(removed.temp, 2, sum) == nrow(removed.temp)
  removed.image[, col.rem] <- 1 
  row.rem <- apply(removed.temp, 1, sum) == ncol(removed.temp)
  removed.image[row.rem, ] <- 2
  
  # erase entirely non-type loci
  remove.loc <- which(colSums(is.na(X)) == nrow(X))
  
  if(length(remove.loc) > 0) {
    
    ## preserve information about the original matrix
    ## reset if the cells contain a 2
    removed.image[, remove.loc] <- 1
    
    ## remove non informative loci
    X <- X[, -remove.loc, drop = FALSE]
    loc.names <- colnames(X)
    nloc <- ncol(X)
    message("Note: removed noninformative loci -pure NAs column(s)- from slots G and A")
    
  }
  
  # erase entirely non-type individuals
  remove.ind <- which(rowSums(is.na(X)) == ncol(X))
  if(length(remove.ind) > 0) {
    
    ## preserve information about the matrix, for restoring
    ## NA individuals. In case of rm.empty.ind = FALSE,
    ## reset removed image for these individuals 
    ## (because removed image is defined as is.na(X))
    if(!rm.empty.ind) {
      removed.image[remove.ind, ] <- 0
    }
    ## preserve order
    old.order.row <- seq(nrow(indata))
    new.order.ind <- old.order.row[-remove.ind]
    ind.oldnames <- ind.names
    
    ## remove individuals
    X <- X[-remove.ind, ]
    ind.names <- rownames(X)
    nind <- nrow(X)
  }
  
  # MARKER SPECIFIC CONFIGURATION--------------------------------------------#
  
  #----------------------------
  # DOMINANT MARKERS 
  #----------------------------
  
  if(type == "dominant") {
    
    out <- X
    
    # check that data values are "0", "1" and NA
    if(!all(out %in% c(NA, "1", "0"))) { 
      stop("dominant data must be binary (0 for absence, 
            1 for presence")
    }
           
    
    # restore missing individuals if required
    if(!rm.empty.ind && length(remove.ind) > 0) {
      temp <- matrix(nrow = nrow(indata), ncol = ncol(out))
      temp[new.order.ind, ] <- out
      rownames(temp) <- ind.oldnames
      colnames(temp) <- colnames(out)
      out <- temp
    } else if(rm.empty.ind && length(remove.ind) > 0) {
      message("Note: removed noninformative individuals -pure NAs row(s)- from slots G and A")
    }
    
    # remove non polymorphic data
    mode(out) <- "numeric"
    isPoly <- aue.is.poly(out, poly.level)
    out <- X[,  isPoly,  drop = FALSE]
    
    loc.names <- colnames(out)
    nloc <- ncol(out)
    
    if(ncol(removed.image) > ncol(out)) {
      if(poly.level == 0) {
        message("Note: non-polymorphic marker(s) deleted")
      } else {
        message(paste("Note: marker(s) with polymorphism level <", 
                      paste(poly.level, "%", sep = ""), "deleted"))
      }  
    }
    # save data image
    removed.image[, !isPoly] <- 1
  }
  
  # END DOMINANT----------------
  
  
  #----------------------------
  # CODOMINANT MARKERS 
  #----------------------------
  
  if(type == "codominant") {
    
    # ncod control 
    ncod <- int.check.ncod(X, ploidy = ploidy, ncod = ncod, sep = sep)
    
    # Handling separators
    
    if(sep == "" && ploidy > 1) {
      
      ## add "/" as separator
      X <- gsub(paste("([[:alnum:]]{",ncod,"})", sep = ""), "\\1/", X)
      X <- sub("/$", "", X)
      sep <- "/"
      
      ## non missing case, checking if <sep> is metacharacter
    } else  {
      X <- gsub(meta2char(sep),"/",X)
      sep <- "/"
    }
    
    
    # Translate data into allelic frequencies 
    
    ## unfold data for each cell of the table
    if (ploidy > 1) {
      allele.data <- strsplit(X, "/")
      n.items <- sapply(allele.data, length)
      locus.data <- rep(rep(loc.names, each = nind), n.items)
      ind.data <- rep(rep(ind.names,ncol(X)), n.items)
      allele.data <- unlist(allele.data)
    } else {
      n.items <- rep(1, length(X))
      locus.data <- rep(rep(loc.names, each = nind), n.items)
      ind.data <- rep(rep(ind.names, ncol(X)), n.items)
      allele.data <- unlist(X)
    }
    
    ## identify NAs
    NA.posi <- which(is.na(allele.data))
    NA.ind <- ind.data[NA.posi]
    NA.locus <- locus.data[NA.posi]
    
    ## remove NAs
    if(length(NA.posi) > 0){
      allele.data <- allele.data[-NA.posi]
      locus.data <- locus.data[-NA.posi]
      ind.data <- ind.data[-NA.posi]
    }
    
    ## get matrix with frequencies
    allele.data <- paste(locus.data, allele.data, sep = ".")
    allele.data <- factor(allele.data, levels = unique(allele.data))
    out <- table(ind.data, allele.data)
    out <- out[ind.names, , drop = FALSE] # table sorts alphabetically. This resets.
    out <- out/2
    ## force type 'matrix'
    class(out) <- NULL
    dimnames(out) <- list(rownames(out), colnames(out))
    
    ## restore NAs
    if(length(NA.posi) > 0) {
      out.colnames <- colnames(out)
      for(i in 1:length(NA.ind)) {
        loc <- paste0(NA.locus[i], "\\.")
        out[NA.ind[i], grep(loc, out.colnames)] <- NA
      }
    }
    
  }
  
  # END CODOMINANT---------------
  
  
  # GENERAL OUTPUT CONFIGURATION-----------------------------------------------#
  
  # restore missing individuals if required
  if(!rm.empty.ind && length(remove.ind) > 0) {
    temp <- matrix(nrow = nrow(indata), ncol = ncol(out))
    temp[new.order.ind, ] <- out
    rownames(temp) <- ind.oldnames
    colnames(temp) <- colnames(out)
    out <- temp
  } else if(rm.empty.ind && length(remove.ind) > 0) {
    message("Note: removed noninformative individuals -pure NAs row(s)- from slots G and A")
  }
  
  # missing data manipulation
  
  missing.cells <-  which(is.na(out))
  if(length(missing.cells) == 0) {
    missing.cells <- integer(0)
  }
  
  ## 0 case
  if (missing == "0") {
    out[is.na(out)] <- 0
  }
  
  ## mean case
  if (missing == "MEAN") {
    mode(out) <- "numeric"
    moy <- round(apply(out, 2, function(x) mean(x, na.rm = TRUE)), 3)
    
    for (j in 1:ncol(out)) {
      out[, j][is.na(out[, j])] <- moy[j]
    }
  }  
  
  # removed image configuration
  rem.cols <- which(apply(removed.image, 2, function(x) any(x == 1)))
  rem.rows <- which(apply(removed.image, 1, function(x) any(x == 2)))
  removed.image <- list(init.dim = dim(indata), 
                        rem.rows = as.vector(rem.rows), 
                        rem.cols = as.vector(rem.cols),
                        names.rows = rownames(indata), 
                        names.cols = colnames(indata))
  
  # numeric matrix
  mode(out) <- "numeric"
  
  # OUTPUT CREATION------------------------------------------------------------#
  
  #dominant configuration
  if(type == "dominant") {
    sep <- NULL
    ncod <- NULL
  }
  
  res <- int.genind(X = out,
                    ploidy = ploidy,
                    type = type,
                    NA.char = NA.char,
                    sep = sep,
                    ncod = ncod,
                    missing = as.character(missing),
                    missing.cells = as.integer(missing.cells),
                    removed.image = removed.image)
  
  res
  
  }

#------------------------------------------------------------------------------#
#---------------------------EXPORTER-------------------------------------------#
#------------------------------------------------------------------------------#

#' export
#' @keywords internal

int.genind2df <- function(x, sep = "",                   #the product is a matrix
                          NA.char = "NA", 
                          restore.removed = FALSE) {
  
  # restore missing cells
  output <- x@tab
  loc.names <- levels(x@loc.fac)
  row.names <- rownames(output)
  
  #--(1/2)---dominant case-----------------------------------------------------#
  if(x@type == "dominant"){
    return(output) 
  }
  
  #--(2/2)---codominant case---------------------------------------------------#
  # make separate tables
  kX <- list()
  loc.fac <- as.character(x@loc.fac)
  for(i in loc.names){
    kX[[i]] <- output[, i == loc.fac, drop = FALSE]
  }
  kX <- lapply(kX, function(X) round(X * x@ploidy))
  
  # function to recode a genotype in form 
  # "A1[sep]...[sep]Ak" from frequencies--------#
  recod <- function(vec, lab){
    if(any(is.na(vec))) {
      return(NA)
    }
    res <- paste(rep(lab, vec), collapse = sep)
    return(res)
  }
  #---------------------------------------------#
  
  # kGen is a list of nloc vectors of genotypes--------------------------------#
  # all.names is splitted into a list 
  all.names <- x@all.names
  all.names <- split(as.vector(all.names), loc.fac)
  all.names <- all.names[loc.names]
  
  kGen <- lapply(1:length(kX), 
                 function(i) { 
                   apply(kX[[i]], 1, recod, all.names[[i]])
                 }
  )
  
  names(kGen) <- loc.names
  
  ## build the final data.frame
  res <- do.call(cbind, kGen)
  res[res == ""] <- NA
  
  
  # RESTORE REMOVED DATA CONFIGURATION-----------------------------------------#
  
  if(restore.removed) {
    
    # info restoration slots
    restore.info <- x@removed.image
    rem.rows <- x@rem.rows
    rem.cols <- x@rem.cols
    
    temp <- matrix(0,nrow = restore.info[[1]][1], ncol = restore.info[[1]][2])
    temp[restore.info[[2]], ] <- NA
    temp[, restore.info[[3]]] <- NA
    l.rr <- length(restore.info[[2]]) != 0 
    l.rc <- length(restore.info[[3]]) != 0
    
    # both cols and rows removed
    if(l.rc && l.rr) {
      temp[-rem.rows, -rem.cols] <- res
      rownames(temp) <- restore.info[[4]]
      colnames(temp) <- restore.info[[5]]
      return(as.data.frame(temp, stringsAsFactors = FALSE))
    } # only columns
    if(l.rc && !l.rr) {
      temp[, -rem.cols] <- res
      rownames(temp) <- restore.info[[4]]
      colnames(temp) <- restore.info[[5]]
      return(as.data.frame(temp, stringsAsFactors = FALSE))
    } # only rows
    if(!l.rc && l.rr) {
      temp[-rem.rows, ] <- res
      rownames(temp) <- restore.info[[4]]
      colnames(temp) <- restore.info[[5]]
      return(temp, stringsAsFactors = FALSE)
    }
  }
  
  # OUTPUT CREATION------------------------------------------------------------#
  
  rownames(res) <- row.names
  colnames(res) <- loc.names
  
  res
}

#########################END INT.GENIND#########################################

