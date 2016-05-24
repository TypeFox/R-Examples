MLoef <- function(robj, splitcr="median")
{
# performs the Martin-Loef LR-test
# robj... object of class RM
# splitcr... splitting criterion for two groups. "median" (default) and "mean"
#            split items in two groups according to the median/mean or item raw
#            scores.
#            a vector of length k (number of items) containing two different
#            elements signifying group membership of items can be supplied.

  ### no test with missing values rh 19-05-10
  if(any(is.na(robj$X))) stop("Martin-Loef Test with NA currently not available\n")

  wrning <- NULL   # initialize an object for warnings

  if(length(splitcr) == 1){   # generate split-vector if "mean" or "median"
    if(splitcr == "median"){
      raw.scores <- colSums(robj$X,na.rm=T)
      numsplit <- as.numeric(raw.scores > median(raw.scores,na.rm=T))
      ## removed for the time being 2011-09-08 rh
      #if( any(raw.scores == median(raw.scores,na.rm=T)) ){   # Only if one item's raw score == the median, a warning is issued
      #  wrning <- which(raw.scores == median(raw.scores,na.rm=T))   # append a warning-slot to the object for print and summary methods
      #  cat("Item(s)",paste(names(wrning),collapse=", "),"with raw score equal to the median assigned to the lower raw score group!\n")
      #}
    }
    if(splitcr=="mean"){
      raw.scores <- colSums(robj$X,na.rm=T)
      numsplit <- as.numeric(raw.scores > mean(raw.scores,na.rm=T))
      ## removed for the time being 2011-09-08 rh
      #if( any(raw.scores == mean(raw.scores,na.rm=T)) ){   # Only if one item's raw score == the mean, a warning is issued
      #  wrning <- which(raw.scores == mean(raw.scores,na.rm=T))   # append a warning-slot to the object for print and summary methods
      #  cat("Item(s)",paste(names(wrning),collapse=", "),"with raw score equal to the mean assigned to the lower raw score group!\n")
      #}
    }
  } else {   # check if the submitted split-vector is appropriate
    if(length(splitcr) != ncol(robj$X)) stop("Split vector too long/short.")
#    if(length(unique(splitcr)) > 2) stop("Only two groups allowed.")
    if(length(unique(splitcr)) < 2) stop("Split vector must contain at least two groups.")
    numsplit <- splitcr
  }
  sp.groups <- unique(numsplit)
  i.groups <- sapply(sp.groups, function(g){ which(numsplit == g) }, simplify=F)

  # check if any group countains less than 2 items
  if( any(unlist(lapply(i.groups, length)) < 2) ){
    stop("Each group of items must contain at least 2 items.")
  }

  # check if one group contains subject with <=1 valid responses
  if(any(unlist(lapply(i.groups, function(g){ any(rowSums(!is.na(robj$X[,g])) <= 1) })))) stop("Groups contain subjects with less than two valid responses.")

  ### possible missing patterns and classification of persons into groups
#  MV.X <- apply(matrix(as.numeric(is.na(robj$X01)),ncol=ncol(robj$X01)),1,paste,collapse="")
#  MV.p <- sort(unique(MV.X))
#  MV.g <- numeric(length=length(MV.X))
#  g <- 1
#  for(i in MV.p){
#    MV.g[MV.X == i] <- g;
#    g <- g + 1
#  }
#  na.X01 <- list()
#  for(i in 1:length(MV.p)){
#    na.X01[[i]] <- matrix(robj$X01[which(MV.g == i),], ncol=ncol(robj$X01))
#  }

#  res1 <- RM(robj$X01[,i.groups[[1]]])
#  res2 <- RM(robj$X01[,i.groups[[2]]])

  # fitting the submodels
  subModels <- lapply(i.groups, function(g){ PCM(robj$X[,g]) })

  ### calculating the numerator and denominator
  sub.tabs <- as.data.frame(sapply(subModels, function(M){
                rowSums(M$X, na.rm=T)
              }))
  sub.tabs <- table(sub.tabs)

  sub.term <- sub.tabs * (log(sub.tabs) - log(nrow(robj$X)))
  sub.term <- sum(na.omit(as.numeric(sub.term)))

  sub.max <- lapply(i.groups, function(g){ sum(apply(robj$X[,g], 2, max)) })

  full.tab  <- table(rowSums(robj$X, na.rm=T))
  full.term <- sum(na.omit(as.numeric( full.tab * (log(full.tab) - log(nrow(robj$X))) )))

  ML.LR <- 2 * (
             sub.term  + sum(unlist(lapply(subModels, `[[`, "loglik")))
           - full.term - robj$loglik
           )

  df <- prod(unlist(sub.max)+1) - (sum(apply(robj$X, 2, max))+1) - length(sp.groups) + 1

#  ml.num <- ml.den <- df <- numeric()

#  for(i in 1:length(MV.p)){
#    .temp.num <- table(rowSums(na.X01[[i]],na.rm=T))
#    .temp.num <- .temp.num[.temp.num > 0]   ### rh
#    ml.num[i] <- sum( (log(.temp.num)-log(sum(.temp.num)))*.temp.num )
#
#    if(nrow(na.X01[[i]]) > 1){
#      .temp.den <- table(rowSums(na.X01[[i]][,i.groups[[1]]],na.rm=T),
#                         rowSums(na.X01[[i]][,i.groups[[2]]],na.rm=T))
#    }
#    else{
#      .temp.den <- table(sum(na.X01[[i]][,i.groups[[1]]],na.rm=T),
#                         sum(na.X01[[i]][,i.groups[[2]]],na.rm=T))
#    }
#    .temp.den <- .temp.den[.temp.den > 0]   ### rh
#    ml.den[i] <- sum( (log(.temp.den)-log(sum(.temp.den)))*.temp.den )
#
#    k1 <- sum(!is.na(na.X01[[i]][1,i.groups[[1]]])) ### rh
#    k2 <- sum(!is.na(na.X01[[i]][1,i.groups[[2]]])) ### rh
#    df[i] <- k1 * k2 -1                              ### rh
#  }
#
#  a <- sum(ml.num)
#  b <- sum(ml.den)
#  k <- c(length(i.groups[[1]]),length(i.groups[[2]]))
#
#  ML.LR <- -2*( (a + robj$loglik) - (b + res1$loglik + res2$loglik) )
#  DF <- prod(k) - 1
  p.value <- 1 - pchisq(ML.LR, df)

  result <- list(LR=ML.LR, df=df, p.value=p.value,
                 fullModel=robj, subModels=subModels,
                 Lf=robj$loglik,  Ls=lapply(subModels, `[[`, "loglik"),
#                 theta.table.RM=table(rowSums(robj$X01)),                        # both used for the plotting
#                 theta.table.MLoef=table(rowSums(res1$X01),rowSums(res2$X01)),   # routine plot.MLoef
                 i.groups=i.groups,
#                 items1=i.groups[[1]], items2=i.groups[[2]], k=k,
                 splitcr=splitcr, split.vector=numsplit, warning=wrning, call=match.call())
  class(result) <- "MLoef"
  return(result)
}
