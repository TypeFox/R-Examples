valid.n<-function(x,na.rm=TRUE) return(ifelse(na.rm,sum(!is.na(x)),length(x)))

propbrk<-function(x,trueval=TRUE,na.rm=TRUE) {
 if(is.na(x) || length(x) == 0) return(0)
 if(na.rm) return(sum(x==trueval,na.rm=TRUE)/valid.n(x))
 else return(sum(x==trueval,na.rm=TRUE)/length(x))
}

sumbrk<-function(x,trueval=TRUE,na.rm=TRUE) {
 return(sum(x==trueval,na.rm=TRUE))
}

binciWu<-function(x,n,alpha=0.05,trueval=TRUE,na.rm=TRUE) {
 if(missing(n)) n<-ifelse(na.rm,valid.n(x),length(x))
 x<-sum(x==trueval,na.rm=TRUE)
 z<-pnorm(1-alpha/2)
 zsq<-z*z
 phat<-ifelse(x<1,x,x/n)
 pest<-phat+zsq/(2*n)
 ci<-(pest+z*sqrt((phat*(1-phat))/n+zsq/(4*n*n)))/(1+zsq/n)
 return(ci)
}

binciWl<-function(x,n,alpha=0.05,trueval=TRUE,na.rm=TRUE) {
 if(missing(n)) n<-ifelse(na.rm,valid.n(x),length(x))
 x<-sum(x==trueval,na.rm=TRUE)
 z<-pnorm(1-alpha/2)
 zsq<-z*z
 phat<-ifelse(x<1,x,x/n)
 pest<-phat+zsq/(2*n)
 ci<-(pest-z*sqrt((phat*(1-phat))/n+zsq/(4*n*n)))/(1+zsq/n)
 return(ci)
}

brkdnNest<-function(formula,data,FUN=c("mean","sd","sd","valid.n"),
 label1="Overall",trueval=TRUE) {

 if(missing(data) || missing(formula)) 
  stop("brkdnNest must be called with a formula for breakdown and a data frame.")
 bn<-as.character(attr(terms(formula),"variables")[-1])
 # number of variables used to break down the leftmost term
 nbn<-length(bn)
 # number of functions to apply at each level of breakdown
 nFUN<-length(FUN)
 brklist<-vector("list",nFUN)
 truevalFUN<-c("propbrk","binciWu","binciWl","sumbrk")
 for(brkfun in 1:nFUN) {
  brklist[[brkfun]]<-vector("list",nbn)
  if(FUN[brkfun] %in% truevalFUN) 
   brklist[[brkfun]][[1]]<-do.call(FUN[brkfun],list(data[[bn[1]]], 
    trueval=trueval,na.rm=TRUE))
  else
   brklist[[brkfun]][[1]]<-do.call(FUN[brkfun],list(data[[bn[1]]],
    na.rm=TRUE))
  names(brklist[[brkfun]][[1]])<-label1
  for(brk in 2:nbn) {
   if(FUN[brkfun] %in% truevalFUN) 
    brklist[[brkfun]][[brk]]<-tapply(data[[bn[1]]],data[bn[2:brk]],
     FUN=match.fun(FUN[brkfun]),trueval=trueval)
   else
    brklist[[brkfun]][[brk]]<-tapply(data[[bn[1]]], 
     data[bn[2:brk]],FUN=match.fun(FUN[brkfun]),na.rm=TRUE)
   names(brklist[[brkfun]][[brk]])<-levels(data[[brkfun[brk]]])
  }
  # get rid of the NAs in valid.n
  if(FUN[brkfun] == "valid.n")
   brklist[[brkfun]]<-
    rapply(brklist[[brkfun]],function(x) ifelse(is.na(x),0,x),how="replace")
 }
 attr(brklist,"class")<-"brklist"
 names(brklist)<-FUN
 return(brklist)
}

sliceArray<-function(x,slice) {
 dimx<-dim(x)
 if(is.null(dimx)) return(x[slice])
 else {
  ndim<-length(dimx)
  slicestring<-
   paste("x[",slice,paste(rep(",",ndim-1),collapse=""),"]",sep="",collapse="")
  newx<-eval(parse(text=slicestring))
  return(newx)
 }
}

print.brklist<-function(x,...) {

 crawlBreakList<-function(x,depth=1) {
  if(length(x)>1) {
   if(depth==1) cat(names(x[[1]]),unlist(x[[1]]),"\n")
   x[[1]]<-NULL
   for(nextbit in 1:length(x[[1]])) {
    newx<-lapply(x,sliceArray,nextbit)
    cat(rep("\t",depth),names(x[[1]][nextbit]),unlist(x[[1]][nextbit]),"\n")
    crawlBreakList(newx,depth=depth+1)
   }
  }
 }

 xnames<-names(x)
 for(func in 1:length(x)) {
  cat(xnames[func],"\n")
  crawlBreakList(x[[func]])
 }
}
