#' \code{projectOntoMGS} 
#' @title projectOntoMGS
#' @description This function takes a list of genes and projects it to a mgs catalogue
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param genebag : vector of gene_ids
#' @param list.mgs : this is the structured MGS information formatted as a list of geneids each corresponding to an MGS
#' @param res.filt.mode : the filtering method, either 'perc' or 'size', default will be size.
#'        perc : genes projected onto mgs are only kept if they exceed a threshold percentage of the mgs size 
#'        size : genes projected onto mgs are only kept if they exceed a threshold number of genes
#' @param res.filt.threshold : the threshold used to retain an mgs according to the filtering mode, either 'perc' or 'size'
#'        if 'perc' : the minimal percentage of genes projected onto mgs (number of genes projected / mgs size) needed to select a mgs 
#'        if 'size' : the minimal number of genes projected onto mgs
#' @param not_projected : whether a last genebag containing not projected genes is to be added to the list of selected mgs
#' @return a list of selected mgs with the geneids according to the projected genes
projectOntoMGS <- function (genebag, list.mgs, res.filt.mode = "size", res.filt.threshold = 50, not_projected = TRUE) {
  res <- list()
  if (res.filt.mode != "perc" & res.filt.mode != "size") {
    stop("Choose 'perc' or 'size' as a filtering mode of the result")
  }
  else {
    for (m in 1:length(list.mgs)) {
      tmp <- genebag[(genebag %in% list.mgs[[m]])]
      if (res.filt.mode == "size") {
        if (length(tmp) >= res.filt.threshold) {
          res[[names(list.mgs)[m]]] <- tmp
        }
      }
      else {
        if (length(tmp) >= (length(list.mgs[[m]]) * res.filt.threshold/100)) {
          res[[names(list.mgs)[m]]] <- tmp
        }
      }
    }
    if (not_projected == TRUE) {
      if(length(unlist(res))==0){
        res[["not_projected"]] <- genebag
      } else {
        res[["not_projected"]] <- genebag[-which(genebag %in% unlist(res))]       
      }
    }
    return(res)
  }
}

#' \code{extractProfiles} 
#' @title extractProfiles
#' @description This function extracts the profiles from a gene profile matrix of a group of genes or a list 
#' of groups of genes. It can also restrict the size of the result
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param genebag : vector or list of gene_ids
#' @param data : raw count or frequency matrix with genes_ids as rawnames
#' @param size.max : default 15000, maximum number of profiles to be extracted
#' @param size.min : default 1, the minimal size threshold above which a group of genes is selected.
#' This is used to extract multiple profiles and filtering the list with a minimal number of genes
#' @param silent : default TRUE, detailling and following computation progress
#' @return a matrix or a list of profile matrixes
extractProfiles <- function(genebag, data, size.max = 15000, size.min = 1, silent=TRUE){
  if(!is.list(genebag)){
    print("Profile extraction from a single gene vector")
    if(length(genebag) < size.max){
      res <- data[match(genebag, rownames(data)),]
    }else{
      res <- data[match(genebag[1:size.max], rownames(data)),]
    }
  }else{
    print("Multiple profile extraction")
    # list of matrixes to be returned
    res <- list()
    for(i in 1:length(genebag)){
      if(i%%10==0 & silent==FALSE){print(i)}
      if(length(genebag[[i]]) < size.max){
        res[[i]] <- data[match(genebag[[i]], rownames(data)),]
      }else{
        res[[i]] <- data[match(genebag[[i]][1:size.max], rownames(data)),]
      }
    }
    names(res) <- names(genebag)
    res <- res[as.numeric(summary(genebag)[,"Length"])>=size.min]
  }
  return(res)
}

#' \code{aggregateProfiles} 
#' @title aggregateProfiles
#' @description This function takes a list of profile matrixes and returns an aggregated big matrix.
#' The individual matrixes can be filtered in size so that the first X rows are used for each of them.
#' This function is used to prepare the data and plot different MGS as barcodes.
#' @author Emmanuelle Le Chatelier
#' @param list.profiles : list of matrix profiles
#' @param max.size : the maximum number of rows to be selected in the final agregated matrix.
#' default   max.size=25.
#' @param min.size : this is the minimum number of rows rows to be selected in the final aggregated matrix. 
#' If a group has less it will be discarded. By default min.size=max.size.
#' @return an aggregated profile matrix.
aggregateProfiles <- function(list.profiles, max.size = 25, min.size = max.size){
  if(!is.list(list.profiles)){
    stop("Error! The object is not a list")
  }
  # list of matrices to be returned
  res <- c()
  for(i in 1:length(list.profiles)){
    if(nrow(list.profiles[[i]]) >= min.size){
      res <- rbind(res, list.profiles[[i]][1:min(max.size, nrow(list.profiles[[i]])),])
    }
  }
  return(res)
}

#' \code{computeFilteredVectors} 
#' @title computeFilteredVectors
#' @description filters and computes verctors based on gene profiles from a single matrix or a 
#' list of matrix profiles
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param profile : list of (or unique) matrix profiles
#' @param type : vectorisation method, vectors are calculated from the mean, median or sum of a given list of genes;
#' default type="mean" otherwise it will be "median" or  "sum"
#' @param filt : filtering threshold in % , rate of positive value in an individual under which matrix is cleaned 
#' and sparse values put to 0 default filt= 0, no filtering
#' @param debug : default is FALSE, when TRUE information on advancement is printed
#' @return a filtered vector or a matrix of filtered vectors
computeFilteredVectors <- function (profile, type = "mean", filt = 0, debug = FALSE) {
  if (is.list(profile)) {
    res <- matrix(data = NA, ncol = ncol(profile[[1]]), nrow = length(profile))
    for (i in 1:length(profile)) {
      if(debug) if (i%%100 == 0) { print(i) }
      if (type == "mean") { # mean
        res[i, ] <- apply(filterMat(as.matrix(profile[[i]]), filt = filt), 2, mean)
      } else {
        if (type == "median") { # median
          res[i, ] <- apply(filterMat(as.matrix(profile[[i]]), filt = filt), 2, median)
        }else { #sum
          res[i, ] <- apply(filterMat(as.matrix(profile[[i]]), filt = filt), 2, sum)
        }
      }
    }
    rownames(res) <- names(profile)
    colnames(res) <- colnames(profile[[1]])
  }
  else {
    if (type == "mean") { # mean
      res <- apply(filterMat(as.matrix(profile), filt = filt), 2, mean)
    }else {
      if(type== "median"){ # median
        res <- apply(filterMat(as.matrix(profile), filt = filt), 2, median)
      }else { # sum
        res <- apply(filterMat(as.matrix(profile), filt = filt), 2, sum)
      }
    }
  }
  return(res)
}

#' \code{buildMgsFinal} 
#' @title buildMgsFinal
#' @@date December 20th 2013
#' @description This function will take a vector of genes (to be transformed into a list of genebags) 
#'      or a list of genebags and will extract the profiles. Next genes well be ordered by connectivity
#'      which is to be computed for each group and the 50 most connected are selected to consitute the 
#'      marker genes. These will be then used to compute the mean vectors. A final object containing all 
#'      this information along with taxonomical annotation will be returned
#' @author Edi Prifti
#' @param genebag : a vector of genes to be projected onto mgs or a list of genebags, default = NULL.
#' @param mgs.cat : MGS catalogue to be used
#' @param mgs.taxo : taxonomy table for the MGS catalogue
#' @param profiles :  the data profile matrix to extract the profiles
#' @param conn : if TRUE the connectivity of a group is to be computed and ordered, default = TRUE.
#' @param silent : print detailled information on progress, default = FALSE.
#' @param filt : filtering based on percentage of prevalence to avoid noise for no signal samples by 
#' computeFilteredVectors.
#' @return a list containing the final elements such as the 50 most connected genes, the mean vectors etc
buildMgsFinal <- function (genebag = NULL, mgs.cat, mgs.taxo, profiles, conn = TRUE, silent = TRUE, filt=20) {
  # if a list of genes than we need to project it onto the mgs.cat
  if (!is.list(genebag)) {
    genebag.list <- projectOntoMGS(genebag = genebag, list.mgs = mgs.cat, res.filt.mode = "size", res.filt.threshold = 50)
    #NOTE: another possibility is to trim the mgs.cat using selectListSize for an upper and lower trim and then projectOntoMGS
  } else { genebag.list <- genebag }
  # extract the profiles for the list of mgs
  genebag.list.dat <- extractProfiles(genebag.list, profiles, silent=silent)
  genebag.list.ordered <- list()
  genebag.list.filtered <- list()
  genebag.list.50 <- list()
  for (i in 1:length(genebag.list)) {
    if (i%%10 == 0 & !silent) {print(i)} # inform on where we are
    dat <- genebag.list.dat[[i]]
    if (conn == TRUE) { # order by connectivity
      con <- connectivity(dat) # connectivity as sum of correlations using default settings.
      dat <- dat[order(con, decreasing = T), ]
    }
    genebag.list.ordered[[i]] <- dat
    genebag.list.filtered[[i]] <- filterMat(as.matrix(dat), filt=filt)
    genebag.list.50[[i]] <- genebag.list.filtered[[i]][1:50, ]
  }
  names(genebag.list.ordered) <- names(genebag.list.dat)
  names(genebag.list.filtered) <- names(genebag.list.dat)
  names(genebag.list.50) <- names(genebag.list.dat)
  # build the resulting object
  res <- list()
  res$genebags <- genebag.list
  res$genebags.all <- mgs.cat[names(genebag.list)]
  res$mgs.profiles <- genebag.list.ordered
  res$mgs.50 <- genebag.list.50
  res$mean_vectors <- computeFilteredVectors(profile = res$mgs.50, type="mean",filt=0) # the filtering was done on the whole MGS
  res$taxonomy <- as.data.frame(mgs.taxo[names(res$genebags), ])
  name <- as.character(res$taxonomy$species)
  name <- paste(names(res$genebags), name)
  name <- gsub(" NA", "", name)
  res$taxonomy$name <- name
  return(res)
}


#' \code{connectivity} 
#' @title connectivity
#' @description This function computes the intra-row correlation and applied a threshold to compute the connections of each row.
#' A connectivity vector is returned.
#' @author Emmanuelle le Chatelier & Edi Prifti
#' @param prof : a profile matrix. This data is used to compute correlations.
#' @param method : the correlation method, default = pearson.
#' @param th : default 0, if >0 than this threshold will be applied to compute a hard threholded connectivity.
#' @param soft : default = FALSE, if TRUE and when th > 0, the connectivity is computed as the soft threholded, sum of correlations 
#' above the threshold
#' @return a vector of connectivity
#' @note this connectivity does not use a hard thresholding but is based on a total correlation score
connectivity <- function (prof, method = "pearson", th = 0, soft=FALSE) {
  if(th < 0 | th >= 1 ) {
    stop("The threhold should be positive and smaller than 1.")
  }else{
    # compute the correlation matrix
    data.cor <- Hmisc::rcorr(t(prof), type = method)$r
    if(soft){ # if real correlations
      if(th==0){
        con <- colSums(data.cor , na.rm = TRUE)
      }else{ # don't take into account values lower than the threshold
        tmp <- data.cor; tmp[tmp<=th] <- 0
        con <- colSums(tmp, na.rm=TRUE)
      }
    }else{ # hard threholding
      con <- colSums((data.cor > th) + 0, na.rm=TRUE)
    }
  }
  return(con)
}


#' \code{selectListSize} 
#' @title selectListSize
#' @description This function will sextract a part of a list based on the length of its components. Typically
#' a genebag list can be used. This function will work for uniclass lists of vectors and data frames.
#' @author Edi Prifti
#' @param l : a list of vectors or matrixes.
#' @param min.size : default = 0 and this has no action. This is the lower threshold for the element's size.
#' @param max.size : default = 0 and this has no action. This is the upper threshold for the element's size.
#' @param names : default = FALSE the function will return a subset of the list, if TRUE the function will 
#' return only the names of the elements of the list as identifiers.
#' @return the trimmed subset of the list.
selectListSize <- function(l, min.size=0, max.size=0, names=FALSE){
  if(min.size==0 & max.size==0){stop("nothing to be done")}
  if(min.size > max.size & max.size!=0){stop("minimum value should be lower than maximum value")}
  if(min.size!=0 | max.size!=0){
    type <- table(unlist(lapply(l, class)))
    if (length(type)==1){
      if(names(type)=="matrix" | names(type)=="data.frame"){
        l.size <- as.numeric(lapply(l, nrow))
      }else {l.size <- as.numeric(lapply(l, length))}
    }else {stop("multi class list case. not handled")}
    
    if(min.size!=0 & max.size!=0) {
      l.ind <- which(l.size >= min.size & l.size <= max.size)
    }else if(max.size!=0){
      l.ind <- which(l.size <= max.size)
    }else if(min.size!=0){
      l.ind <- which(l.size >= min.size)
    }
    if(names){
      return(names(l[l.ind]))
    }else{
      return(l[l.ind])
    }
  }
}


#' \code{profiles2Genebags} 
#' @title profiles2Genebags
#' @description This function will extract from a list of group profiles a list of identifiers of the matrix.
#' @author Edi Prifti
#' @param profiles : a list of dataframes
#' @return a list of genebags
profiles2Genebags <- function(profiles){
    if(!is.list(profiles)){
        stop("Error ! A list of profiles is needed.")
    }else{
      if(dim(profiles[[1]])<2){
        stop("Error ! the elements of the list should be data frames or matrixes.")
      }else{
        return(lapply(profiles, rownames))
      }
    }
}

#' End of section and file