pasteCols<-function(x,sep="") {
 pastestring<-paste("list(",paste("x","[",1:dim(x)[1],",]",
  sep="",collapse=","),")",sep="")
 return(do.call(paste,c(eval(parse(text = pastestring)),sep=sep)))
}

# converts a two column matrix of object identifiers [,1] and attributes [,2]
# into a data frame of TRUE/FALSE values where each row is an object and
# each column is an attribute

categoryReshape<-function(x) {
 dimx<-dim(x)
 if(is.null(dimx) || dimx[2]==1)
  stop("Can only reshape a matrix or data frame with at least two columns")
 row_values<-sort(unique(x[,1]))
 column_values<-sort(unique(x[,2]))
 newx<-
  as.data.frame(matrix(0,nrow=length(row_values),ncol=length(column_values)))
 for(row in 1:dimx[1]) {
  row_index<-which(row_values %in% x[row,1])
  column_index<-which(column_values %in% x[row,2])
  newx[row_index,column_index]<-1
 }
 names(newx)<-column_values
 rownames(newx)<-row_values
 return(newx)
}

# makes an intersectList object from a matrix or data frame of TRUE/FALSE
# values where each row represents an object and each column an attribute
# TRUE indicates that the object has that attribute, FALSE that it does not
# add a "weight" vector that allows a count of objects to be read directly

makeIntersectList<-function(x,xnames=NULL,sep="+") {
 # If any entries in x are not 1/0 OR TRUE/FALSE, assume that x
 # is a two column matrix of object identifiers [,1] and attributes [,2]
 if(any(!(x %in% c(TRUE,FALSE)))) x<-categoryReshape(x)
 if(is.null(xnames)) xnames <- colnames(x)
 dimx<-dim(x)
 if(is.null(xnames)) xnames<-LETTERS[1:dimx[2]]
 intersectList<-vector("list",dimx[2]+2)
 for(intersect in 1:dimx[2])
  intersectList[[1]][intersect]<-sum(rowSums(x)==1 & x[,intersect])
 names(intersectList[[1]])<-xnames
 for(comb in 2:dimx[2]) {
  nn<-choose(dimx[2],comb)
  intersectList[[comb]]<-rep(0,nn)
  currentnames<-
   names(intersectList[[comb]])<-pasteCols(combn(xnames,comb),sep)
  currentcombs<-combn(1:dimx[2],comb,simplify=TRUE)
  for(intersect in 1:nn) {
   combvec<-rep(0,dimx[2])
   combvec[currentcombs[,intersect]]<-1
   intersectList[[comb]][intersect]<-
    sum(colSums(apply(x,1,"==",combvec))==dimx[2])
  }
 }
 intersectList[[dimx[2]+1]]<-dimx[1]
 names(intersectList[[dimx[2] + 1]])<-"Total"
 intersectList[[dimx[2]+2]]<-xnames
 names(intersectList[[dimx[2] + 2]])<-"attributes"
 # drop any empty intersection levels
 for(comb in dimx[2]:1)
  if(sum(intersectList[[comb]])==0) intersectList[[comb]]<-NULL
 class(intersectList)<-"intersectList"
 return(intersectList)
}
