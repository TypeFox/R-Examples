intersectDiagram<-function(x,pct=FALSE,show.nulls=FALSE,xnames=NULL, 
 sep="+",mar=c(0,0,3,0),main="Intersection Diagram",cex=1,col=NULL,
 minspacing=NA,all.intersections=FALSE,include=NULL) {

 matchParts<-function(x,table,ignore.case=TRUE) {
  for(pattern in 1:length(x)) {
   match_index<-grep(x[pattern],table,ignore.case=ignore.case)
   if(length(match_index)) return(match_index)
  }
  return(0)
 }
 if(!match(class(x),"intersectList",0)) {
  if(is.matrix(x) || is.data.frame(x)) {
   if(is.data.frame(x)) 
    x<-as.matrix(x)
    x<-makeIntersectList(x,xnames=xnames,sep=sep)
  }
  if(!match(class(x),"intersectList",0)) 
   stop("x must be a matrix, data frame or intersectList")
 }
 oldmar<-par("mar")
 par(mar=mar)
 # attribute labels
 attributes<-x[[length(x)]]
 # get all the names for the individual attributes
 if(is.null(include)) include<-attributes
 # total number of attributes
 nattributes<-length(attributes)
 # peel off the number of objects and the attributes for display
 x[[length(x)]]<-NULL
 nobjects<-x[[length(x)]]
 x[[length(x)]]<-NULL
 # number of intersection levels with at least one object
 nlevels<-length(x)
 # if no colors specified, use rainbow
 if(is.null(col)) col<-c(rainbow(nattributes),NA)
 else
  if(length(col) < nattributes) col<-rep(col,length.out=nattributes)
 # total number of objects for each intersection level
 objectsums<-sapply(x,sum)
 # index of level with the most objects
 maxlevel<-which.max(objectsums)
 nNonZero<-function(x) return(sum(x > 0))
 # number of intersections with at least one member for each intersection level
 # or all intersections if the somewhat dangerous "show everything" option is TRUE
 if(all.intersections) nintersects<-sapply(x,length)
 else nintersects<-sapply(x,nNonZero)
 # maximum number of intersections in a given level
 maxintersections<-max(nintersects)
 # largest intersection set in x
 maxn<-max(unlist(x))
 # default to a minimum spacing of one tenth of the largest intersection set
 if(is.na(minspacing)) minspacing<-0.1 * maxn
 # x limit that will hold the maximum number of objects and allow
 # spacing for the maximum number of intersections in units of objects 
 maxx<-objectsums[maxlevel] + minspacing * maxintersections
 # have to escape the separator in case it is "+" (default) or
 # some other character that means something to some function
 attsep<-paste("[",sep,"]",sep="")
 # display the empty plot
 plot(0,xlim=c(0,maxx),ylim=c(0,nlevels+show.nulls), 
  main=main,xlab="",ylab="",type="n",axes=FALSE)
 # step through each level of intersections
 for(level in 1:nlevels) {
  # determine the intersect level by the number of elements in the first name
  intersectLevel<-length(unlist(strsplit(names(x[[level]][1]),attsep)))
  # indices of intersections with at least one object in this level
  # or just all of the intersections
  if(all.intersections) intersections<-1:nintersects[[level]]
  else intersections<-which(x[[level]] > 0)
  # get all the names in this level with at least one object
  blocknames<-names(x[[level]])[intersections]
  # spacing between intersection sets in object units
  spacing<-(maxx - objectsums[level])/nintersects[level]
  # left edges of the rectangles in x positions
  leftx<-c(0,cumsum(x[[level]][intersections] + spacing)) + spacing/2
  # now step through the intersections in this level
  for(intersect in 1:length(intersections)) {
   # check if this intersection is to be displayed
   if(matchParts(include,blocknames[intersect])) {
    # make the label for the intersection
    cellqnt<-ifelse(pct,
     paste(round(100*x[[level]][intersections[intersect]]/nobjects,1),"%",sep=""),
     x[[level]][intersections[intersect]])
    # indices of the colors to use for this rectangle
    colindex<-
     which(attributes %in% unlist(strsplit(blocknames[intersect],attsep)))
    # number of colors
    ncol<-length(colindex)
    # width of each color slice
    xinc<-x[[level]][intersections[intersect]]/ncol
    # colors for the slices
    slicecol<-col[colindex]
    # start at the left edge of the sliced rectangle
    offset<-0
    # step through the slices
    for(slice in 1:ncol) {
     # first draw the rectangle with no border
     rect(leftx[intersect]+offset,nlevels-level+show.nulls+0.1,
      leftx[intersect]+offset+xinc,nlevels-level+show.nulls+0.9,
      col=slicecol[slice],border=NA)
     # move to the left edge of the next slice
     offset<-offset+xinc
    }
    # draw a box around the sliced rectangle
    rect(leftx[intersect],nlevels-level+show.nulls+0.1,
     leftx[intersect]+x[[level]][intersections[intersect]], 
     nlevels-level+show.nulls+0.9)
    # display the label for this rectangle
    boxed.labels(leftx[intersect]+x[[level]][intersections[intersect]]/2, 
     nlevels-level+show.nulls+0.5,
     paste(blocknames[intersect],cellqnt,sep="\n"),cex=cex)
   }
  }
 }
 if(show.nulls) {
  # number of objects with no set membership or no attributes
  nonset<-as.numeric(nobjects - sum(objectsums))
  # left edge of the rectangle
  leftnulls<-sum(par("usr")[1:2])/2-nonset/2
  # draw the rectangle
  if(nonset) rect(leftnulls,0.1,leftnulls+nonset,0.9)
  # center of the rectangle
  xpos<-leftnulls+nonset/2
  # display the label
  if(pct) nonset<-paste(round(100*nonset/nobjects,1),"%",sep="")
  boxed.labels(xpos,0.5,paste("Non-members",nonset,sep="\n"),cex=cex)
 }
 # restore the original plot parameters
 par(mar=oldmar)
 # stick the number of objects and attributes back on
 x[[length(x) + 1]]<-nobjects
 x[[length(x) + 1]]<-attributes
 invisible(x)
}
