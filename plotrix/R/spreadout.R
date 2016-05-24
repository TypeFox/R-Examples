spreadout<-function(x,mindist) {
 if(sum(!is.na(x)) < 2) return(x)
 xorder<-order(x)
 goodx<-x[xorder][!is.na(x[xorder])]
 gxlen<-length(goodx)
 start<-end<-gxlen%/%2
 # nicely spread groups of short intervals apart from their mean
 while(start > 0) {
  while(end < gxlen && goodx[end+1] - goodx[end] < mindist) end<-end+1
  while(start > 1 && goodx[start] - goodx[start-1] < mindist) start<-start-1
  if(start < end) {
   nsqueezed<-1+end-start
   newx<-sum(goodx[start:end])/nsqueezed -
    mindist*(nsqueezed%/%2 - (nsqueezed/2 == nsqueezed%/%2) * 0.5)
   for(stretch in start:end) {
    goodx[stretch]<-newx
    newx<-newx+mindist
   }
  }
  start<-end<-start-1
 }
 start<-end<-length(goodx)%/%2+1
 while(start < gxlen) {
  while(start > 1 && goodx[start] - goodx[start-1] < mindist) start<-start-1
  while(end < gxlen && goodx[end+1] - goodx[end] < mindist) end<-end+1
  if(start < end) {
   nsqueezed<-1+end-start
   newx<-sum(goodx[start:end])/nsqueezed -
    mindist*(nsqueezed%/%2 - (nsqueezed/2 == nsqueezed%/%2) * 0.5)
   for(stretch in start:end) {
    goodx[stretch]<-newx
    newx<-newx+mindist
   }
  }
  start<-end<-end+1
 }
 # force any remaining short intervals apart
 if(any(diff(goodx) < mindist)) {
  start<-gxlen%/%2
  while(start > 1) {
   if(goodx[start] - goodx[start-1] < mindist)
    goodx[start-1]<-goodx[start]-mindist
   start<-start-1
  }
  end<-gxlen%/%2
  while(end < gxlen) {
   if(goodx[end+1] - goodx[end] < mindist)
    goodx[end+1]<-goodx[end]+mindist
   end<-end+1
  }
 }
 x[xorder][!is.na(x[xorder])]<-goodx
 return(x)
}
