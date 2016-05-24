brkdn<-function(formula,data,maxlevels=10,num.desc=c("mean","var","valid.n"),
 width=10,round.n=2) {

 if(!missing(data) && !missing(formula)) {
  bn<-as.character(attr(terms(formula),"variables")[-1])
  nbn<-length(bn)
  cat("\nBreakdown of",bn[1],"by",bn[nbn],"\n")
  if(!is.numeric(data[[bn[1]]])) stop("\nVariable on left of formula must be numeric")
  if(nbn > 2) {
   by.factor<-as.factor(data[[bn[nbn]]])
   factor.levels<-levels(by.factor)
   factor.labels<-attr(data[,bn[nbn]],"value.labels")
   if(!is.null(names(factor.labels))) {
    factor.levels<-factor.labels
    factor.labels<-names(factor.labels)
   }
   if(is.null(factor.labels)) factor.labels<-factor.levels
   nlevels<-length(factor.levels)
   if(nlevels > maxlevels) {
    nlevels<-maxlevels
    cat("Too many levels - only using first",maxlevels,"\n")
   }
   brkstats<-vector("list",nlevels)
   names(brkstats)<-factor.levels[1:nlevels]
   for(i in 1:nlevels) {
    currentdata<-subset(data,by.factor==factor.levels[i])
    for(j in 1:dim(currentdata)[2])
     attr(currentdata[,j],"value.labels")<-attr(data[,j],"value.labels")
    cat(paste("\nVariable ",bn[1],sep="",collapse=""),"for",bn[nbn],"- level",factor.labels[i],"\n\n")
    cat(formatC(num.desc,width=10),"\n",sep="")
    junk<-describe.numeric(as.matrix(currentdata[bn[1]]),num.desc=num.desc)
    next.formula<-as.formula(paste(paste(bn[1],"~"),paste(bn[2:(nbn-1)],collapse="+")))
    brkstats[[i]]<-brkdn(next.formula,currentdata,maxlevels=maxlevels,num.desc=num.desc)
    class(brkstats[[i]])<-"dstat"
   }
   class(brkstats)<-"dstat"
   invisible(brkstats)
  }
  else {
   by.factor<-as.factor(data[[bn[2]]])
   factor.levels<-levels(by.factor)
   factor.labels<-attr(data[,bn[2]],"value.labels")
   if(!is.null(names(factor.labels))) {
    factor.levels<-factor.labels
    factor.labels<-names(factor.labels)
   }
   if(is.null(factor.labels)) factor.labels<-factor.levels
   nlevels<-length(factor.levels)
   if(nlevels > maxlevels){
    nlevels<-maxlevels
    cat("Too many levels - only using first",maxlevels,"\n")
   }
   gstats<-matrix(NA,ncol=nlevels,nrow=length(num.desc))
   colnames(gstats)<-factor.labels[1:nlevels]
   rownames(gstats)<-num.desc
   if(is.numeric(data[[bn[1]]])) {
    round.ns<-rep(round.n,length(num.desc))
    npos<-match("valid.n",num.desc)
    if(!is.na(npos)) round.ns[npos]<-0
    vwidth<-max(nzchar(factor.labels))
    cat(formatC("Level",width=vwidth))
    cat(formatC(num.desc,width=width),"\n",sep="")
    for(i in 1:nlevels) {
     currentdata<-subset(data[[bn[1]]],by.factor==factor.levels[i])
     if(length(currentdata)){
      gstats[,i]<-describe.numeric(currentdata,num.desc=num.desc)
      cat(formatC(factor.labels[i],width=vwidth),
      formatC(round(gstats[,i],round.ns),width=width),"\n")
     }
    }
    class(gstats)<-"dstat"
    rnames<-rownames(gstats)
   }
   invisible(gstats)
  }
 }
 else cat("Usage: brkdn(formula,data,maxlevels=10,num.desc=c(\"mean\",\"var\",\"valid.n\"))\n")
}
