#' Compares Mahalanobis distances from two approaches
#' 
#' Mahalanobis distances are calculated for each zero pattern.
#' Two approaches are used. The first one estimates Mahalanobis distance for observations belonging to one each zero pattern each.
#' The second method uses a more sophisticated approach described below.
#' @aliases getEstimates print.estimates plot.estimates
#' @param x data frame or matrix
#' @param imp imputation method
#' @param y unused second argument for the plot method
#' @param ... additional arguments for plotting passed through
#' @return \item{df}{a data.frame containing the Mahalanobis distances from the estimation in subgroups, the Mahalanobis distances from the imputation and covariance approach, an indicator specifiying outliers and an indicator specifying the zero pattern} \item{df2}{a groupwise statistics.}
#' @export
#' @import ggplot2
#' @import data.table
#' @author Matthias Templ, Karel Hron
#' @seealso \code{\link{impKNNa}}, \code{\link{isomLR}}
#' @examples
#' 
#' data(arcticLake)
#' # generate some zeros
#' arcticLake[1:10, 1] <- 0
#' arcticLake[11:20, 2] <- 0
#' m <- compareMahal(arcticLake)
#' plot(m)
compareMahal <- function(x, imp="KNNa"){
  if(is.null(dim(x))) stop("x must be a data.frame, data.table or matrix")
  if(imp %in% c("impKNNa","knna","KNNa","kNNa")){
    imp <- "KNNa"
  } else if (imp %in% c("knn", "KNN", "kNN")){
    imp <- "kNN"
  } else if (imp %in% c("IRMI","irmi","Irmi")){
    imp <- "irmi"
  } else {
    stop("wrong method for imputation specified")  
  }	
  
  ## create an ID for latter use
  x <- cbind(x, rn = 1:nrow(x))
  p <- ncol(x)
  ind <- 1:(p-1)
  ## code zeros as NA's
  if(any(is.na(x)) & any(x==0, na.rm=TRUE)) warning("the data includes NA's and zeros. \n Impute the missing values first, otherwise they are treated as zeros")
  x[x == 0] <- NA
  w <- is.na(x[,ind])
  w <- apply(w, 2, as.integer)
  s <- apply(w, 1, paste, collapse="")
  xs <- split(x, s)
  names(xs) <- gsub("0", "x", names(xs))
  names(xs) <- gsub("1", "0", names(xs))
  ## if one group is too small report an error
  check <- as.numeric(unlist(lapply(xs, nrow)))
  w <- which(check < 2*(ncol(x)-1) + 2)
  if(length(w) > 0){
    m <- missPatterns(x[,-ncol(x)])$tabcombPlus
    cat("\n the following subgroups:\n")
    print(m[w,])
    cat("\n are too small for evaluation in each subcomposition")
    stop("subgroups must be larger than 2*ncol(x)+1")
  }
  ## exclude NA columns:
  xs <- lapply(xs, function(x){
    x <- x[,!is.na(x[1,]),drop=FALSE]
    x
  })
  ## mahalanobis distances in each subset:
  calcMahal <- function(x){
    if(ncol(x) > 2){
      xilr <- isomLR(x)
      covs <- robustbase::covMcd(xilr)
      mahcorr <- sqrt(as.numeric(mahalanobis(xilr, center=covs$center, cov=covs$cov)))
      mahcorr <- mahcorr / sqrt(qchisq(0.975, ncol(xilr)))
    } else if(ncol(x) == 2){
      xilr <- isomLR(x)
      covs <- covMcd(xilr)
      zscore <- (as.numeric(xilr) - covs$center) / sqrt(covs$cov) 
      mahcorr <- abs(zscore) / qnorm(0.975)
    } else {
      mahcorr <- NA
    }
    return(mahcorr)
  }
  ## Mahalonibis distances for subsets:
  m1 <- lapply(xs, function(x,...) calcMahal(x[,-ncol(x)]))
  ## Mahalonibis distances via imputation approach:
  resi <- zeroOut(x[,ind,drop=FALSE], imp)
  
  ## prepare data for plot:
  len <- lapply(m1, length)
  id <- lapply(xs, function(x) x[, ncol(x)])
  id <- as.numeric(unlist(id))
  df <- data.frame("sub"=unlist(m1),
                   "group"=factor(rep(names(m1), len)),
                   "cols"=factor(rep(1:length(len), len)),
                   "id"=id
  )
  df <- merge(df, x, by.y="rn",by.x="id")
  df <- cbind(resi, df)
#  df <- merge(df, resi, by.x="id", by.y="ID")
  ## rearange columns:
  ordering <- c("sub","mahcorr","outlier","group")
  df <- df[, ordering]
  colnames(df)[2] <- "imp"
  ## add number of obs to group, easier with data.table 
  df <- data.table::data.table(df)
  ## to avoid CRAN notes:
  obs <- group <- NULL
  ## group:
  df[, obs:=.N, by = list(group)]
  df <- data.frame(df)
  df$group <- apply(df[,c("group","obs")], 1, paste, collapse=",  obs =")
  ## make text:
  outsub <- aggregate(df[,"sub",drop=FALSE], list(df$group), function(x) sum(x > 1, na.rm=TRUE))
  outimp <- aggregate(df[,"imp",drop=FALSE], list(df$group), function(x) sum(x > 1, na.rm=TRUE))
  g1 <- table(df$group)
  df2 <- data.frame("group"=names(table(df$group)),
                    "outsub"=outsub[,2],
                    "outimp"=outimp[,2],
                    "gsize"=as.numeric(g1))
  df2$outrelsub <- df2[, "outsub"] / df2[, "gsize"]
  df2$outrelimp <- df2[, "outimp"] / df2[, "gsize"]  
  df2$x <- rep(1.3, times = nrow(df2))
  df2$y <- rep(5, times = nrow(df2))
  #  ## bring it to right order:
  #  df <- df[order(df$id),]
  result <- list(df=df, df2=df2)
  class(result) <- "mahal"
  return(result)
}

#' @rdname compareMahal
#' @export
#' @method plot mahal
plot.mahal <- function(x, y, ...){
  df <- x$df
  df2 <- x$df2
  ## to avoid CRAN notes:
  imp <- sub <- NULL
  ## ggplot:
  g <- ggplot2::ggplot(aes(x=imp, y=sub), data=df)
  g <- g + ggplot2::geom_point() 
  g <- g + ggplot2::facet_wrap(~group)
  g <- g + ggplot2::geom_hline(yintercept=1, linetype=2, colour="lightgrey")
  g <- g + ggplot2::geom_vline(xintercept=1, linetype=2, colour="lightgrey")
  g <- g + ggplot2::geom_abline(intercept = 0, slope = 1, colour="grey")
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::xlab("imputation subcompositional approach")
  g <- g + ggplot2::ylab("estimation in subcompositions only")
  #  g <- g + geom_text(data=df2[,c("nam","outsub")], 
  #                     aes(x=1, y=5, label=outsub), parse=TRUE)
  #  g <- g + geom_text(aes(x, y, label=outsub, data=df2))
  ## to avoid CRAN notes:
  outsub <- outimp <- NULL
  g <- g + ggplot2::geom_text(data = df2, aes(x=x, y=y, label=outsub), size=4)
  g <- g + ggplot2::geom_text(data = df2, aes(x=y, y=x, label=outimp), size=4)
  #  g <- g + geom_text(data = df2, aes(x=y, y=y, label=paste("obs=",gsize,sep="")))
  print(g)
} 