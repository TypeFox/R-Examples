#' Detection of outliers of zero-inflated data
#' 
#' detects outliers in compositional zero-inflated data 
#' @param x a data frame 
#' @param impute imputation method internally used
#' @details XXX
#' @return XXX
#' @export
#' @author Matthias Templ
#' @examples 
#' ### Installing and loading required packages
#' data(expenditures)
zeroOut <- function(x, impute="knn"){
  ## @Matthias Templ, TU WIEN, 2012
  rownames(x) <- 1:nrow(x)
  ind <- 1:ncol(x)
  D <- ncol(x)
  ## 1. Imputiere
  if(impute %in% c("impKNNa","knna","KNNa")){
    xi <- impKNNa(x)$xImp
  } else if (impute %in% c("knn", "KNN", "kNN")){
    xi <- kNN(x, imp_var = FALSE)
#  } else if (impute %in% c("fry", "Fry", "FRY")){
#    xi <- rmzero(x, minval=0.01, delta=0.01)
  } else if (impute %in% c("IRMI","irmi","Irmi")){
    xi <- impCoda(x, init="geometricmean")$xImp
  } else {
    stop("wrong method for imputation specified")  
  }	
  x <- cbind(x, ID=1:nrow(x))
  xi <- cbind(xi, ID=1:nrow(x))
  
  ## make sure that xi is a data.frame:
  if(class(xi)=="matrix") xi <- data.frame(xi)
  w <- is.na(x[, ind])
  #    w <- apply(w, 2, as.integer)
  s <- apply(w, 1, paste, collapse=":")
  #    xi <- cbind(xi, id=1:nrow(x)) #new
  #    x <- cbind(x, id=1:nrow(x)) #new    
  xs <- split(xi, s)
  getSortIndex <- function(x, s){
    xs <- split(x, s)
    ## TRUE when zero
    lapply(xs, function(x){
      is.na(x[1,])
    })
  }
  si <- getSortIndex(x[,ind], s)
  zneworder <- xs
  mah <- pval <- mahcorr <- IDlist <- list()
  for(i in 1:length(xs)){
    index <- names(xs[i])
    index <- as.logical(strsplit(index, ":")[[1]])
    sortedxs <- xs[[i]]
    wt <- which(index)
    wf <- which(!index)
    sortedxs <- sortedxs[, c(wt,wf)]
    zneworder <- isomLR(sortedxs)
    zcovs <- robustbase::covMcd(zneworder)
    ## took only last columns of xs
    if(length(wf) == 2){
      p <- ncol(zneworder)
      zscore <- (zneworder[, p] - zcovs$center[p]) / sqrt(zcovs$cov[p,p]) 
      mah[[i]] <- abs(zscore)
      pval[[i]] <- pnorm(mah[[i]])
      mahcorr[[i]] <- mah[[i]] / qnorm(0.975)
      names(mahcorr[[i]]) <- names(pval[[i]]) <- names(mah[[i]]) <- rownames(xs[[i]])
    } else if(length(wf) > 2){
      noneff <- c((length(wt) + 1):(D-1))    
      mah[[i]] <- sqrt(as.numeric(mahalanobis(zneworder[, noneff], center=zcovs$center[noneff], cov=zcovs$cov[noneff, noneff])))
      pval[[i]] <- pchisq((mah[[i]])^2, length(wf))
      mahcorr[[i]] <- mah[[i]] / sqrt(qchisq(0.975, ncol(zneworder[, noneff])))
      names(mahcorr[[i]]) <- names(pval[[i]]) <- names(mah[[i]]) <- rownames(xs[[i]])
    } else{
      mah[[i]] <- NA
      pval[[i]] <- NA
      mahcorr[[i]] <- NA
    }
    IDlist[[i]] <- xs[[i]][ncol(xs[[i]])]
  }
  ## list --> data.frame
  nam <- names(unlist(mahcorr))
  df <- data.frame("mah"=as.numeric(unlist(mah)), "pval"=as.numeric(unlist(pval)), "mahcorr"=as.numeric(unlist(mahcorr)),
                   "ID"=nam)
  df <- merge(x, df, by="ID")
  df <- cbind(df, "outlier"=df$mahcorr > 1)
  return(df)  	
}		
