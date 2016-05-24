getIntersectList<-function(nelem,xnames=NULL,sep="+") {
if(is.null(xnames)) xnames<-LETTERS[1:nelem]
xnamelen<-length(xnames)
if(xnamelen < nelem) {
  extranames<-paste("extra",1:(nelem-xnamelen),sep=sep)
  cat("Not enough names in",xnames,"adding",extranames,"\n")
  xnames<-c(xnames,extranames)
}
intersectList<-vector("list",nelem+2)
total_n<-0
for(comb in 1:nelem) {
  nn<-choose(nelem,comb)
  intersectList[[comb]]<-rep(0,nn)
  currentnames<-names(intersectList[[comb]])<-
   pasteCols(combn(xnames,comb),sep=sep)
  for(intersect in 1:nn) {
   cat("Number of elements in",currentnames[intersect],"- ")
   current_n<-scan(nmax=1,quiet=TRUE)
   intersectList[[comb]][intersect]<-current_n
   total_n<-total_n+current_n
  }
}
cat("Total number of elements (press Enter for ",total_n,") - ",
  sep="")
intersectList[[nelem + 1]]<-scan(nmax=1,quiet=TRUE)
if(length(intersectList[[nelem+1]])==0)
  intersectList[[nelem+1]]<-total_n
intersectList[[nelem+2]]<-xnames
class(intersectList) <- "intersectList"
return(intersectList)
}
