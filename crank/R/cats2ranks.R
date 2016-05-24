# cats2ranks produces a table of mean ranks for a matrix or data frame of values
# where each value represents one of a number of options as numbers (1:noptions).
# Rows are considered to represent respondents or cases, and columns are
# entry points.
# The item typically contains a small number of entry points labeled:
# "Most important <option>", "Second most important <option>",...
# and the respondent is to enter one option code in each entry point.
# The code entered in the first entry point is considered to have rank 1 and so on.
# The important aspect of this function is that it substitutes the mean of
# unallocated ranks for all options that were not entered.
# That is, if there were ten options and only five entry points, the five options
# not entered by a respondent will receive the mean of the unallocated ranks - 8.
# If the respondent does not fill all entry points, the mean of unallocated ranks
# will be calculated appropriately and NAs will be assigned that mean.
# This procedure corrects for the bias introduced when some options are chosen
# by few respondents, but receive high (numerically low) ranks.
# The user will want to pass the set of option codes if not all appear in x.

cats2ranks<-function(x,cats=NULL) {
 if(!missing(x)) {
  dimx<-dim(x)
  if(is.null(cats)) {
   if(is.null(dimx)) {
    cats<-sort(unique(x))
    nrows<-1
   }
   if(is.data.frame(x)) cats<-sort(unique(unlist(x)))
   if(is.matrix(x)) cats<-sort(unique(as.vector(x)))
  }
 }
 if(is.null(cats)) {
  stopmsg<-paste("Usage: cats2ranks(x,cats=NULL)\n",
   "\twhere x is a vector, matrix or data frame,\n",
   "\tcats is a vector of all possible categories (if not all are in x),\n",
   sep="")
  stop(stopmsg)
 }
 ncats<-length(cats)
 ranks<-1:ncats
 nwranks<-matrix(NA,nrow=dimx[1],ncol=ncats)
 catsum<-catcount<-rep(0,ncats)
 for(row in 1:dimx[1]) {
  n.valid<-sum(!is.na(x[row,]))
  norank<-mean(ranks[-(1:n.valid)])
  nwranks[row,]<-rep(norank,ncats)
  for(col in 1:dimx[2]) {
   if(!is.na(x[row,col])) {
    thisrank<-which(cats==x[row,col])
    nwranks[row,thisrank]<-col
    catcount[thisrank]<-catcount[thisrank]+1
   }
  }
  catsum<-catsum+nwranks[row,thisrank]
 }
 newranks<-list(ranks=nwranks,categories=cats,ranksum=catsum,
  rankcount=catcount)
 class(newranks)<-"cats2ranks"
 return(newranks)
}

print.cats2ranks<-function(x,...) {
 rankdim<-dim(x$ranks)
 if(is.null(rankdim)) nranks<-length(x$ranks)
 else nranks<-rankdim[2]
 meanranks<-unlist(apply(x$ranks,2,mean))
 cat("\n   Mean rank\tN\n")
 rankpos<-order(meanranks)
 for(i in 1:nranks)
  cat(x$categories[rankpos[i]],"\t",round(meanranks[rankpos[i]],2),"\t",
   x$rankcount[rankpos[i]],"\n",sep="")
}
