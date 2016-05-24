#' \code{testRelations} 
#' @title testRelations
#' @description This function applies a statistical test either a correlation (spearman, pearson), 
#'      wilcoxon or t.test as a function of a given phenotype. It will return a matrix of probabilities
#'      p and q values along with the correlation coefficient or the enrichment variable when a binary parameter.
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param data : frequency matrix with gene_ids in the rownames
#' @param trait : a vector with the trait to test, binary or numerical variable
#' @param type : a character string indicating the type of test to be applied
#' @param restrict : an optional logical vector to select a subset of the samples to perform the test
#'      default restrict = rep(TRUE, ncol(data)) ie all the samples are selected
#' @param multiple.adjust : type of multiple adjustment default is "BH" i.e. Benjamini & Hochberg method 
#' @param paired : logical with default FALSE wether the test should be paired or not
#' @param debug : default FALSE, when TRUE the progress is printed each 1000 steps
#' @return a matrix with analytical results (correlation tests) indicating rho, rho2, p and q values for each parameter tested
#' @note New addon taking into account a trait for correlation, when it is a two class variable with the same number of elements
#'      a correlation between both groups is performed
testRelations <- function (data, trait, type, restrict = rep(TRUE, ncol(data)), multiple.adjust = "BH", paired = FALSE, debug=FALSE) {
  trait.val <- names(table(trait))
  res <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = 5))
  rownames(res) <- rownames(data)
  colnames(res) <- c("rho", "rho2", "p", "q", "status")
  if (type != "spearman" & type != "pearson") {
    if (length(table(trait)) != 2) {
      stop("Sorry, you can't use this test. The trait should contain only two categories!")
    }
    if (type == "wilcoxon" | type == "t.test") {
      if (type == "wilcoxon") {
        if(debug) print("Executing Wilcoxon test")
        for (i in 1:nrow(data)) {
          if (i%%1000 == 0) {
            print(i)
          }
          tmp <- wilcox.test(data[i, restrict] ~ trait[restrict], 
                             paired = paired)
          res[i, "p"] <- tmp$p.value
          if (mean(data[i, restrict][trait == trait.val[1]]) > 
                mean(data[i, restrict][trait == trait.val[2]])) {
            res[i, "status"] <- trait.val[1]
          }
          else {
            res[i, "status"] <- trait.val[2]
          }
        }
        res[, 4] <- p.adjust(res[, "p"], method = multiple.adjust)
      }
      else {
        if(debug) print("Executing T test")
        for (i in 1:nrow(data)) {
          if (i%%1000 == 0) {
            print(i)
          }
          tmp <- t.test(data[i, restrict] ~ trait[restrict], 
                        paired = paired)
          res[i, "p"] <- tmp$p.value
          if (mean(data[i, restrict][trait == trait.val[1]]) > 
                mean(data[i, restrict][trait == trait.val[2]])) {
            res[i, "status"] <- trait.val[1]
          }
          else {
            res[i, "status"] <- trait.val[2]
          }
        }
        res[, 4] <- p.adjust(res[, "p"], method = multiple.adjust)
      }
    }
    else {
      stop("Sorry, your test does not exist! Available : spearman, pearson, t.test and wilcoxon")
    }
  }
  else {
    if (type == "spearman") {
      if(paired){ # correlate two vectors split by trait in two classes
        if(debug) print("Entering correlation mode between two classes")
        if(length(table(trait))!=2 | table(trait)[1]!=table(trait)[2]){
          stop("Sorry, trait does not seem to be a 2-level categorical variable or with the same prevalence")
        }
        
        # run the correlation
        cl1 <- names(table(trait)[1]); cl1.ind <- trait==cl1; cl2 <- names(table(trait)[2]); cl2.ind <- trait==cl2
        if(debug) print("Executing Spearman correlation")
        for (i in 1:nrow(data)) {
          if (i%%1000 == 0) {
            print(i)
          }
          tmp <- Hmisc::rcorr(data[i, cl1.ind], data[i, cl2.ind], type = "spearman")
          res[i, 1] <- tmp$r[1, 2]
          if (!is.na(res[i, 1])) {
            if (res[i, 1] > 0) {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          }
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$P[1, 2]
        }
        res[, 4] <- p.adjust(res[, "p"], method = multiple.adjust)
      }else{
        if(debug) print("Executing Spearman correlation")
        for (i in 1:nrow(data)) {
          if (i%%1000 == 0) { print(i) }
          tmp <- Hmisc::rcorr(data[i, restrict], trait[restrict], type = "spearman")
          res[i, 1] <- tmp$r[1, 2]
          if (!is.na(res[i, 1])) {
            if (res[i, 1] > 0) {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          }
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$P[1, 2]
        }
        res[, 4] <- p.adjust(res[, "p"], method = multiple.adjust) 
      }
    }
    else {
      if(paired){ # correlate two vectors split by trait in two classes
        if(debug) print("Entering correlation mode between two classes")
        if(length(table(trait))!=2 | table(trait)[1]!=table(trait)[2]){
          stop("Sorry, trait does not seem to be a 2-level categorical variable or with the same prevalence")
        }
        
        # run the correlation
        cl1 <- names(table(trait)[1]); cl1.ind <- trait==cl1; cl2 <- names(table(trait)[2]); cl2.ind <- trait==cl2
        print("Executing Pearson correlation")
        for (i in 1:nrow(data)) {
          if (i%%1000 == 0) {
            print(i)
          }
          tmp <- Hmisc::rcorr(data[i, cl1.ind], data[i, cl2.ind], type = "pearson")
          res[i, 1] <- tmp$r[1, 2]
          if (!is.na(res[i, 1])) {
            if (res[i, 1] > 0) {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          }
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$P[1, 2]
        }
        res[, 4] <- p.adjust(res[, "p"], method = multiple.adjust)
      }else{
        if(debug) print("Executing Pearson correlation")
        for (i in 1:nrow(data)) {
          if (i%%1000 == 0) {
            print(i)
          }
          tmp <- Hmisc::rcorr(data[i, restrict], trait[restrict], type = "pearson")
          res[i, 1] <- tmp$r[1, 2]
          if (!is.na(res[i, 1])) {
            if (res[i, 1] > 0) {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          }
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$P[1, 2]
        }
        res[, 4] <- p.adjust(res[, "p"], method = multiple.adjust)
      }
    }
  }
  return(res)
}


#' \code{hierClust} 
#' @title hierClust
#' @description This function computes the pairwise distance between samples and computes a hierarchical clustering
#' that is further depicted as a heatmap graphic
#' @author Edi Prifti
#' @param data : frequency matrix with gene_ids in the rownames
#' @param side : the distance can be performed on the columns or on the rows
#' @param dist : the type of distance used. By default this is correlation based similarity
#' @param cor.type : when correlation matrix, the default is spearman
#' @param hclust.method : the hierarchical clustering method, by default it is the ward method
#' @param side.col.c : a vector of colors to be applied in the columns, usually depincting a class
#' @param side.col.r : a vector of colors to be applied in the rows, usually depincting a class
#' @param plot : logical default TRUE. It will plot the heatmap of the similarity with the hierchical clustering
#' @return it will return a list of three variables, the correlation matrix, the distance matrix and the hclust object
#' @note updated hierClust functions by elechat april 7th 2015 added options SideColors added + spearman == pearson(rank)
hierClust <- function (data, side = "col", dist = "correlation", cor.type = "spearman", 
                       hclust.method = "ward", side.col.c = NULL, side.col.r = NULL, plot = TRUE) {
  res <- NULL
  if (side == "col") {
    if (sum(colSums(data, na.rm=TRUE) == 0) > 0) {
      warning("Warning there are samples with no signal that need to be discarded")
    }
  }
  else {
    if (sum(rowSums(data, na.rm=TRUE) == 0) > 0) {
      warning("Warning there are genes with no signal that need to be discarded")
    }
  }
  if (dist == "correlation") {
    if (side != "col") { # if rows we need to transpose
      data <- t(data)
    }
    if (cor.type == "spearman") { data <- apply(data, 2, rank) } # compute the rank for spearman
    mat.rho <- Hmisc::rcorr(data, type = "pearson")$r
    diag(mat.rho) <- NA # don't use the diagonal, correlation with themselves
    # compute the distance as 1-correlation
    mat.dist <- as.dist(1 - mat.rho)
    mat.hclust <- hclust(d = mat.dist, method = hclust.method)
    diag(mat.rho) <- 1 # add the diag again
    if (plot) {
      if (is.null(side.col.c) & is.null(side.col.r)) { # if none is provided
        gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                          Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                          margins = c(6, 6), cex.axis = 0.7)
      } else {
        if (is.null(side.col.c)){ # if column class is not provided but the row is
          if (length(side.col.r) != ncol(mat.rho)) {warning("side.col.r must be a character vector of entry length")}
          gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                            Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                            margins = c(6, 6), cex.axis = 0.7, RowSideColors=side.col.r)
        }else{ # if row class is not provided but the culumn is
          if (is.null(side.col.r)){
            if (length(side.col.c) != ncol(mat.rho)) {warning("side.col.c must be a character vector of entry length")}
            gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                              Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                              margins = c(6, 6), cex.axis = 0.7, ColSideColors=side.col.c)
          } else { # if both classes are provided
            if (length(side.col.r) != ncol(mat.rho)) {warning("side.col.r must be a character vector of entry length")}
            if (length(side.col.c) != ncol(mat.rho)) {warning("side.col.c must be a character vector of entry length")}
            gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                              Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                              margins = c(6, 6), cex.axis = 0.7, ColSideColors=side.col.c, RowSideColors=side.col.r)
          }
        }
      }
    }
    res <- list(mat.rho = mat.rho, mat.dist = mat.dist, mat.hclust = mat.hclust)
  } else { stop("Other distances are not yet implemented, only 1-correlation is used here !") }
  return(res)
}


#' \code{filt.hierClust} 
#' @title filt.hierClust
#' @description This function takes as input a square similarity matrix and searches for clusters of samples with strong associations
#' and extracts the sub matrix with the closely related sampless. Only positive correlations are considered here.
#' @author Emmanuelle Le Chatelier & Edi Prifti
#' @param mat.rho : square correlation matrix with ids (can be used for also other than just samples)
#' @param hclust.method : the hierarchical clustering method, by default it is the ward method
#' @param side.col.c : a vector of colors to be applied in the columns, usually depincting a class
#' @param side.col.r : a vector of colors to be applied in the rows, usually depincting a class
#' @param size : the number of samples in the resulting ordered matrix
#' @param plot : logical default TRUE. It will plot the heatmap of the similarity with the hierchical clustering
#' @param filt : default is 0.5 and is the filtering threshold to be applied
#' @return it will return a matrix with samples in rows and their closely related ones on the columns along with the 
#' correlation score.
filt.hierClust <- function (mat.rho, hclust.method = "ward", side.col.c = NULL, side.col.r = NULL, size=10, plot = TRUE, filt = 0.5) {
  rho <- mat.rho # keep it as a backup
  diag(mat.rho) <- 0

  if((nrow(mat.rho)-1) < size){
    size <- nrow(mat.rho)
    warning(paste("There are less than", size, "samples. Setting size to the number of rows"))
    size <- nrow(mat.rho) -1
  }
  
  res <- as.data.frame(matrix(NA, ncol(mat.rho), size*2));  rownames(res) <- colnames(mat.rho)
  colnames(res) <- paste(c("Hit","Hit_rho"),sort(c(1:size,1:size)),sep="_")
  
  # For each row sort the samples and extract the best samples and put in the res the collected data
  for (i in 1:nrow(mat.rho)){
    tmp <- sort(mat.rho[i,], decreasing=T)
    tmp <- tmp[1:size]
    res[i,seq(1,size*2,by=2)] <- names(tmp)
    res[i,seq(2,size*2,by=2)] <- round(tmp,3)
    sel <- colnames(res)[res[i,]==0]
    sel <- gsub("rho_","",sel)
    res[i,sel] <- 0
  }
  
  N <- apply(mat.rho,1,max) # find the maximums
  N <- N>filt # transform it in a logical index
  mat.rho <- rho[N,N] # keep only the sub matrix which is very connected
  if(all(dim(mat.rho)!=0)){
    # sub select the classes if they are not null
    if(!is.null(side.col.c)) side.col.c <- side.col.c[N] 
    if(!is.null(side.col.r)) side.col.r <- side.col.r[N] 
    diag(mat.rho) <- NA
    # commpute distance as 1-correlation
    mat.dist <- as.dist(1 - mat.rho)
    # compute hierarchical clustering
    mat.hclust <- hclust(d = mat.dist, method = hclust.method)
    diag(mat.rho) <- 1
    if (plot) {
      if (is.null(side.col.c) & is.null(side.col.r)) { # if none is provided
        gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                          Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                          margins = c(6, 6), cex.axis = 0.7)
      } else {
        if (is.null(side.col.c)){ # if column class is not provided but the row is
          if (length(side.col.r) != ncol(mat.rho)) {warning("side.col.r must be a character vector of entry length")}
          gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                            Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                            margins = c(6, 6), cex.axis = 0.7, RowSideColors=side.col.r)
        }else{ # if row class is not provided but the culumn is
          if (is.null(side.col.r)){
            if (length(side.col.c) != ncol(mat.rho)) {warning("side.col.c must be a character vector of entry length")}
            gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                              Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                              margins = c(6, 6), cex.axis = 0.7, ColSideColors=side.col.c)
          } else { # if both classes are provided
            if (length(side.col.r) != ncol(mat.rho)) {warning("side.col.r must be a character vector of entry length")}
            if (length(side.col.c) != ncol(mat.rho)) {warning("side.col.c must be a character vector of entry length")}
            gplots::heatmap.2(mat.rho, scale = "none", trace = "none", 
                              Rowv = as.dendrogram(mat.hclust), Colv = as.dendrogram(mat.hclust), 
                              margins = c(6, 6), cex.axis = 0.7, ColSideColors=side.col.c, RowSideColors=side.col.r)
          }
        }
      }
    }
  }else{
    warning(paste("There are no related samples above the threshold",filt))
  }
  return(res)  
}


#' \code{lmp} 
#' @title lmp
#' @description This function will extract the p-value from a linear model object. It is used by phenoPairwiseRelations
#' @author Edi Prifti
#' @param modelobject : a linear model object as produced by lm()
#' @return a p-value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#' \code{phenoPairwiseRelations} 
#' @title phenoPairwiseRelations
#' @description This function will compute all the relations between different variables adapting different statistical tests as a 
#' function of the data type. It will adjust the p-value matrix for multiple testing.
#' @author Edi Prifti
#' @param data : bioclinical data with variables on the rows and samples on the columns
#' @param adjust : method to adjust for multiple testing (default="BH)
#' @param verbose : default=FALSE. If TRUE information will be printed to follow the progression.
#' @return a list of two matrixes containing the p-values and the multiple testing adjustment.
phenoPairwiseRelations <- function(data, adjust="BH", verbose=FALSE){
  #require(nortest)
  #require(Hmisc)
  data.relations <- matrix(NA, ncol(data), ncol(data)); colnames(data.relations) <- colnames(data); rownames(data.relations) <- colnames(data) 
  for(i in 1:ncol(data)){
    if (verbose) print(paste(i,colnames(data)[i]))
    for(j in 1:ncol(data)){
      if (verbose) print(paste(j,colnames(data)[j]))
      x <- data[,i]
      y <- data[,j]
      # test weather the variables are factors
      # a) both factors
      if(is.factor(x)){ 
        if(is.factor(y)){
          # perfom a chisquare
          if(nrow(table(x,y))==1 | ncol(table(x,y))==1 | sum(table(x,y))==0){
            # do nothing
          }else{
            data.relations[i,j] <- chisq.test(table(x,y))$p.value
          }
        }else{
          # b) x is a factor but not y
          # test for normality
          if(nortest::lillie.test(y)$p.value>0.05){ #if normal
            if(length(table(x))>1){
              if(length(table(x))==2){
                if(all(colSums(table(y,x))>1)){
                  data.relations[i,j] <- t.test(y~x)$p.value
                }
              }else{
                if(sum(table(y,x))>1 & length(unique(x[!is.na(x)]))!=1 & sum(colSums(table(y,x))>0)>1){
                  tmp <- lmp(lm(y~x))
                  if(!is.null(tmp)) {
                    data.relations[i,j] <- tmp
                  }
                }
              }
            }
          }else{ # not normal
            # run anova
            if(length(table(x))>1){
              if(length(table(x))==2){
                # should test here the distribution of y if normal than t.test otherwise wilcox
                if(all(colSums(table(y,x))>1)){
                  data.relations[i,j] <- wilcox.test(y~x)$p.value
                }
              }else{
                if(sum(table(y,x))>1 & length(unique(x[!is.na(x)]))!=1 & sum(colSums(table(y,x))>0)>1){
                  tmp <- kruskal.test(y~x)
                  if(!is.null(tmp)) {
                    data.relations[i,j] <- tmp$p.value
                  }
                }
              }
            }
          }
        }
      }else{ # x is not a factor
        if(is.factor(y)){ # but y is a factor
          # test for normality
          if(nortest::lillie.test(x)$p.value>0.05){ #if normal
            if(length(table(y))>1){
              if(length(table(y))==2){
                if(all(colSums(table(x,y))>1)){
                  data.relations[i,j] <- t.test(x~y)$p.value
                }
              }else{
                if(sum(table(x,y))>1 & length(unique(y[!is.na(y)]))!=1 & sum(colSums(table(x,y))>0)>1){
                  tmp <- lmp(lm(x~y))
                  if(!is.null(tmp)) {
                    data.relations[i,j] <- tmp
                  }
                }
              }
            }
          }else{ # non normal
            if(length(table(y))>1){
              if(length(table(y))==2){
                if(all(colSums(table(x,y))>1)){
                  data.relations[i,j] <- wilcox.test(x~y)$p.value
                }
              }else{
                if(sum(table(x,y))>1 & length(unique(y[!is.na(y)]))!=1 & sum(colSums(table(x,y))>0)>1){
                  tmp <- kruskal.test(x~y)
                  if(!is.null(tmp)) {
                    data.relations[i,j] <- tmp$p.value
                  }
                }
              }
            }
          }
        }else{
          # if both are numeric, run a correlation by testing for normality
          if (length(x[!is.na(x)])>4 & length(x[!is.na(y)])>4){
            if(nortest::lillie.test(x)$p.value<0.05 | nortest::lillie.test(y)$p.value<0.05){
              # use spearman
              data.relations[i,j] <- Hmisc::rcorr(x,y,type="spearman")$P[1,2]
            }else{
              # use pearson
              data.relations[i,j] <- Hmisc::rcorr(x,y,type="pearson")$P[1,2]
            }
          }
        }
      }
    }
  }
  # adjust for multiple testing
  data.relations.q <-matrix(p.adjust(data.relations ,method=adjust), nrow=nrow(data.relations)); 
  colnames(data.relations.q) <- colnames(data.relations); rownames(data.relations.q) <- rownames(data.relations)
  return(list(p=data.relations,q=data.relations.q))
}


#' \code{extractSignificant} 
#' @title extractSignificant
#' @description This function will extract a list of vectors p- or q-values from an object produced by phenoPairwiseRelations.
#' @author Edi Prifti
#' @param relation.matrix : a matrix of p- produced by or q-values by phenoPairwiseRelations()
#' @param interest : a vector of variable names of interest.
#' @param threshold : default 0.05 needed to select significant relations
#' @return a list of vectors containing p-values or q-values along with the names of the variables.
extractSignificant <- function(relation.matrix, interest, threshold=0.05){
  res <- list()
  for(i in 1:length(interest)){
    tmp.val <- relation.matrix[,interest[i]]
    tmp <- (tmp.val < threshold)
    tmp[is.na(tmp)] <- FALSE
    res[[i]] <- tmp.val[tmp]
  }
  names(res) <- interest
  return(res)
}


#' End of section and file