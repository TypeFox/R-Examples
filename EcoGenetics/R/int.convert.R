
#' INTERNAL CONVERSION TOOLS FOR GENETIC DATA
#' @param X Input
#' @param ncod Number of digits coding each allele
#'  (e.g., 1: x, 2: xx, 3: xxx, etc.). If NULL, ncod will we 
#'  obtained from the ploidy and the maximum number of characters
#'  in the data cells.
#' @param sep.in Separator in the input
#' @param sep.out Separator in the output
#' @param ploidy Ploidy of the data.
#' @param chk.names Defalult TRUE. The function make checks of individuals 
#' and loci names during conversion.
#' @param chk.plocod  Defalult TRUE. The function checks coherence 
#' in ploidy and number of digits coding alleles for loci data during conversion.


#--------------------------------------------------------------------#
#-LOCUS MATRIX TO:
# locus matrix to allele matrix
# locus matrix to  locus matrix
# locus matrix to locus list
# locus matrix to allele list

#-ALLELE MATRIX TO:
# allele matrix to locus matrix
# allele matrix to locus list
# allele matrix to allele list

#-LOCUS LIST TO:
# locus list to locus list
# locus list to allele list
# locus list to allele matrix
# locus list to locus matrix

#-ALLELE LIST TO:
# allele list to locus list
# allele list to locus matrix
# allele list to allele matrix


#--------------------------------------------------------------------#
# locus matrix to allele matrix

int.loc2al <- function(X, 
                       ncod = NULL, 
                       ploidy = 2, 
                       sep.in = "", 
                       sep.out = NULL,
                       chk.names = TRUE,
                       chk.plocod = TRUE) {
  
  nind <- nrow(X)
  nloc <- ncol(X)
  
  X <- as.matrix(X, rownames.force = TRUE)
  mode(X) <- "character"
  X <- gsub(sep.in, "", X)
  
  if(chk.names) {
    X <- int.check.colnames(X)
    X <- int.check.rownames(X)
  }
  
  coldata <- int.loc2loc(X, ncod = ncod, ploidy = ploidy, 
                         sep.out = "/", chk.names = FALSE,
                         chk.plocod =  chk.plocod)
  
  # unfolding and creating a matrix with one allele per column
  coldata[is.na(coldata)]<-"NA/NA"
  coldata <- as.vector(t(coldata))
  coldata <- unlist(strsplit(coldata, "/"))
  coldata <- matrix(coldata, ncol = ploidy * nloc, nrow = nind, byrow = TRUE)
  coldata[coldata == "NA"] <- NA
  
  # column names configuration
  if(ploidy != 1) {
    nombres <- lapply(colnames(X), function(x) paste(rep(x, ploidy), 1:ploidy, sep = "."))
    nombres <- unlist(nombres)
  } else {
    nombres <- colnames(X)
  } 
  
  #output mode 
  
  # no.numeric <- grep("[^[:digit:]]", coldata)
  # if(length(no.numeric) == 0) {
  #  suppressWarnings(mode(coldata) <- "numeric")
  # } 
  
  coldata <- aue.rmspaces(coldata)
  coldata <- gsub("(NA)+", NA, coldata)
  
  coldata <- matrix(coldata, nrow = nind)
  
  colnames(coldata) <- nombres
  rownames(coldata) <- rownames(X)
  
  coldata
  
}

#--------------------------------------------------------------------#
# allele matrix to locus matrix

int.al2loc <- function(X, 
                       ncod = NULL, 
                       ploidy, 
                       sep.in = NULL,
                       sep.out = "",
                       chk.names = TRUE,
                       chk.plocod) {
  
  X <- as.matrix(X, rownames.force = TRUE)
  mode(X) <- "character"
  
  if(ploidy != 1) {
    nloc <- ncol(X) / ploidy
  } else {
    nloc <- ncol(X)
  }
  
  nom <- gsub("[.][^.]*$", "", colnames(X))
  nom <- unique(nom)
  
  if(chk.names) {
    X <- int.check.colnames(X)
    X <- int.check.rownames(X)
    nom <- int.check.vnames(nom, len.X = nloc)
  }
  
  xseq <- aue.seqlist(1, ploidy * nloc, by = ploidy)
  
  X.list  <- lapply(1:ncol(xseq), 
                    function(i) {
                      apply(X[, xseq[, i], drop = FALSE], 
                            1, 
                            paste, 
                            sep = "", collapse = sep.out)
                    })
  
  X <- do.call(cbind, X.list)
  
  X <- gsub("(NA)+", NA, X)
  
  colnames(X) <- nom
  
  X
  
}

#--------------------------------------------------------------------#
# locus matrix to locus matrix

int.loc2loc <- function(X, 
                        ncod = NULL, 
                        ploidy = 2, 
                        sep.in = "",
                        sep.out = "",
                        chk.names = TRUE,
                        chk.plocod = TRUE) {
  
  X <- as.matrix(X, rownames.force = TRUE)
  mode(X) <- "character"
  X <- gsub(sep.in, "", X)
  
  # control and configuration
  if(length(sep.out) != 1) {
    stop("sep.out must be a character of length 1")
  }
  
  if(chk.names) {
    X <- int.check.colnames(X)
    X <- int.check.rownames(X)
  }
   
  if(chk.plocod) {
  ncod <- int.check.ncod(X, ploidy = ploidy, ncod = ncod)
  }
  
  # separate alleles with the character "sep.out" 
  X <- gsub(paste("([[:alnum:]]{",ncod,"})",sep = ""), 
            paste("\\1", meta2char(sep.out), sep = ""),  X)
  X <- sub(paste(meta2char(sep.out), "$", sep = ""), "", X)
  X[is.na(X)] <- paste(rep("NA", ploidy), sep.out, collapse = "", sep = "")
  X <- sub(paste(meta2char(sep.out), "$", sep = ""), "", X)
  
  X <- gsub("(NA)+", NA, X)
  
  X
  
}

#--------------------------------------------------------------------#
# locus matrix to locus list

int.loc2list <- function(X,
                         ncod = NULL, 
                         ploidy,
                         sep.in = "",
                         sep.out = "",
                         chk.names = TRUE,
                         chk.plocod = TRUE) {
  
  #matrix contol included loc2loc
  X  <- int.loc2loc(X, ncod = ncod, ploidy = ploidy, 
                    sep.in = sep.in, sep.out = sep.out, 
                    chk.names = chk.names,  chk.plocod =  chk.plocod)
  
  X.list <- list()
  for(i in 1:ncol(X)) {
    X.list[[i]] <- X[, i, drop = FALSE]
  }
  names(X.list) <- colnames(X)
  
  X.list
  
}

#--------------------------------------------------------------------#
# allele matrix to locus list

int.al2list <- function(X, ncod = NULL, ploidy, sep.out = "",  sep.in = NULL) {
  
  X <- int.al2loc(X, ploidy = ploidy, sep.out = sep.out)
  
  X.list <- list()
  for(i in 1:ncol(X)) {
    X.list[[i]] <- X[, i, drop = FALSE]
  }
  names(X.list) <- colnames(X)
  
  X.list
  
}

#--------------------------------------------------------------------#
# locus matrix to allele list

int.loc2listal <- function(X, 
                           ncod = NULL, 
                           ploidy,
                           sep.in = "", 
                           sep.out = NULL,
                           chk.names = TRUE,
                           chk.plocod = TRUE) {
  
  nloc <- ncol(X)
  
  #- to columns
  listnames <- colnames(X)
  if(chk.names) {
    listnames <- int.check.vnames(listnames, len.X = ncol(X))
  }
  
  X <- int.loc2al(X, ncod = ncod, ploidy = ploidy, chk.names = FALSE,
                  chk.plocod = chk.plocod)
  
  X.list <- list()
  
  xseq <- aue.seqlist(1, ploidy * nloc, by = ploidy)
  for(i in 1:ncol(xseq)) {
    X.list[[i]] <- X[, xseq[, i], drop = FALSE]
  }
  
  names(X.list) <- listnames
  X.list
  
}

#--------------------------------------------------------------------#
# allele matrix to allele list

int.al2listal <- function(X, 
                          ncod = NULL, 
                          ploidy, 
                          sep.in = NULL, 
                          sep.out = NULL,
                          chk.names = TRUE,
                          chk.plocod) {
  
  if(ploidy != 1) {
    nloc <- ncol(X) / ploidy
  } else {
    nloc <- ncol(X)
  }
  
  listnames <- gsub("[.][^.]*$", "", colnames(X))
  listnames <- unique(listnames)
  if(chk.names) {
    listnames <- int.check.vnames(listnames, len.X = nloc)
  }
  
  X.list <- list()
  xseq <- aue.seqlist(1, ploidy * nloc, by = ploidy)
  for(i in 1:ncol(xseq)) {
    X.list[[i]] <- X[, xseq[, i], drop = FALSE]
  }
  
  names(X.list) <- listnames
  
  X.list
  
}

#--------------------------------------------------------------------#
# locus list to locus list

int.list2list <- function(X, 
                          ncod = NULL, 
                          ploidy = 2, 
                          sep.in = "",
                          sep.out = "",
                          chk.names = TRUE,
                          chk.plocod = TRUE) {
  
  out <- lapply(X, function(x) int.loc2loc(x, 
                                           ncod = ncod, 
                                           ploidy = ploidy,
                                           sep.in = sep.in,
                                           sep.out = sep.out,
                                           chk.names = chk.names,
                                           chk.plocod =  chk.plocod))
  out
  
}

#--------------------------------------------------------------------#
# locus list to allele list

int.list2listal <- function(X, 
                            ncod = NULL, 
                            ploidy = 2, 
                            sep.in = "",
                            sep.out = NULL,
                            chk.names = TRUE,
                            chk.plocod) {
  
  out <- lapply(X, function(x) int.loc2al(x, 
                                          ncod = ncod, 
                                          ploidy = ploidy,
                                          sep.in = sep.in,
                                          chk.names = chk.names))
  out
}

#--------------------------------------------------------------------#
# allele list to locus list

int.listal2list <- function(X, 
                            ncod = NULL, 
                            ploidy = 2,
                            sep.in = NULL,
                            sep.out = "",
                            chk.names = NULL,
                            chk.plocod) {
  
  out <- lapply(X, 
                function(x) {
                  tmp <- apply(x, 1, paste, sep = sep.out)
                  tmp <- gsub("(NA)+", NA, tmp)
                  tmp
                  })
  out
  
}

#--------------------------------------------------------------------#
# locus list to locus matrix

int.list2loc <- function(X, 
                         ncod = NULL, 
                         ploidy = 2,
                         sep.in = "",
                         sep.out = "",
                         chk.names = TRUE,
                         chk.plocod = TRUE) {
  
  out <- do.call(cbind, data)
  out <- int.loc2loc(data, ncod =  ncod, ploidy = ploidy,
                     sep.in = sep.in, sep.out = sep.out,
                     chk.names = chk.names,  chk.plocod = chk.plocod)
  out
}

#--------------------------------------------------------------------#
# locus list to allele matrix

int.list2al <- function(X, 
                        ncod = NULL, 
                        ploidy = 2,
                        sep.in = "",
                        sep.out = NULL,
                        chk.names = TRUE,
                        chk.plocod = TRUE) {
  
  X <- do.call(cbind, X)
  X<- int.loc2loc(X, ncod =  ncod, ploidy = ploidy,
                  sep.in = sep.in, sep.out = "",
                  chk.names = chk.names,  chk.plocod = chk.plocod)
  X <- int.loc2al(X, ncod = ncod, ploidy = ploidy,
                  chk.names = FALSE)
  X
}

#--------------------------------------------------------------------#
# allele list to locus matrix

int.listal2loc <- function(X, 
                           ncod = NULL,
                           ploidy = 2,
                           sep.in = NULL,
                           sep.out = "",
                           chk.names = TRUE,
                           chk.plocod) {
  
  X <- do.call(cbind, X)
  X <- int.al2loc(data, ploidy = ploidy, sep.out = sep.out, chk.names = chk.names)
  X
}

#--------------------------------------------------------------------#
# allele list to allele matrix

int.listal2al <- function(X, 
                          ncod = NULL, 
                          ploidy = NULL,
                          sep.in = NULL,
                          sep.out = NULL,
                          chk.names = NULL,
                          chk.plocod) {
  
  do.call(cbind, X)
}

#--------------------------------------------------------------------#
