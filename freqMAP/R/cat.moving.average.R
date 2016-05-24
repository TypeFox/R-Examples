`cat.moving.average` <-
function(dat,x,hw,cat.names=NULL){

  if(mode(dat[,1])!="numeric")
    stop("First column of dat must be numeric")
  if(mode(dat[,2])!="character")
    stop("Second column of dat must be character")

  if(any(apply(dat,2,function(x){sum(is.na(x))})))
    stop("NA not allowed in dat")

  if( !(mode(cat.names) %in% c("NULL","character")) )
    stop("cat.names must be a vector of strings.")

  if(is.null(cat.names)){
    #The list of unique cat names
    cat.names <- sort(unique(dat[,2]))
  }

  if(any(!(unique(dat[,2] %in% cat.names)))){
    stop("cat.names must contain at least all unique values of dat[,2]")
  }
  
  cat.ma <- data.frame(x=x,n=0)
  
  for(a in cat.names){
    cat.ma[,a] <- NA
    for(i in 1:length(x)){
      in.group <- (dat[,1] >= x[i]-hw) & (dat[,1] <  x[i]+hw)
      n <- sum(in.group)
      if(n>0){
        cat.ma[i,a] <- sum(dat[in.group,2]==a)/n
        cat.ma[i,2] <- n
      }
    }
  }
  
  cat.ma
}

