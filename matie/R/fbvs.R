# forward backward explanitary variable set selection function
fbvs <- function(dataSet,one=names(d)[length(names(d))],maxv=length(names(d))-1,linear=FALSE){
  
  getRsq <- function(d,one=names(d)[length(names(d))],
                     best=names(d)[1:length(names(d))-1],
                     linear=FALSE){
    # make a new data set and get rid of NaNs
    d <- d[,c(best,c(one))]
    d <- d[complete.cases(d),]
    if(nrow(d)<length(best)+2)return(0.0)
    if(linear) {
      model <- paste(one,"~",best[1])
      if(length(best)>1) for(i in 2:length(best)) model <- paste(model,best[i],sep=" + ")
      fit <- lm(model, data=d)
      return(summary(fit)$adj.r.squared)
    }
    A <- suppressWarnings(ma(cbind(d[,best],d[,one]))$A)
    return(A)
  }
  
  # Remove non-numeric columns from data frames
  d <- dataSet[ , sapply(dataSet, is.numeric)]
  # Remove columns containing NaNs
  # d <- dataSet[ , colSums(is.na(d)) == 0]
  if(! linear) d <- as.data.frame(rwt(d))
  n <- names(d)
  
  # check that the one is still there
  if (!(one %in% n)) stop(paste(one, " is not a name of a numeric variable"))
  
  # prepare print message
  if(linear) mess <- "Rsq" else mess <- "A"
  
  posv <- n[!(n %in% c(one))]
  print(posv)
  
  # find best starting variable
  maxA <- 0
  rmv <- NULL
  for(v in posv){
    A <- getRsq(d,one,c(v),linear)
    if(A>maxA && is.finite(A)) {
      maxA<-A
      rmv <- v
    }
  }
  
  # return if no association found
  if(is.null(rmv))stop("no association found")
  
  # initialize the best group
  bestv <- c(rmv)
  posv <-posv[!(posv %in% c(rmv))]
  print(paste("added... ",rmv," max",mess,"=",round(maxA,2)))
  
  
  while( (length(bestv)<maxv) && (length(posv)>1)){
    # forward step
    rmv <- NULL
    for(v in posv){
      A <- getRsq(d,one,c(bestv,v),linear)
      if(A>maxA && is.finite(A)) {
        maxA<-A
        rmv <- v
      }
    }
    if(is.null(rmv))return(list(one=one,best=bestv,Rsq=maxA))
    bestv <- c(bestv,c(rmv))
    posv <-posv[!(posv %in% c(rmv))]
    print(paste("added... ",rmv," max",mess,"=",round(maxA,2)))
    
    # backward step
    rmv <- NULL
    for(v in bestv){
      A <- getRsq(d,one,bestv[!(bestv %in% c(v))],linear)
      if(A>maxA && is.finite(A)) {
        maxA<-A
        rmv <- v
      }
    }
    if(!(is.null(rmv))){
      bestv <- bestv[!(bestv %in% c(rmv))]
      posv <-c(posv,c(rmv))
      print(paste("removed... ",rmv," max",mess,"=",round(maxA,2)))
    }
  }
  return(list(one=one,best=bestv,Rsq=maxA))
  
}
