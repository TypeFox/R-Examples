plotinfo<-function(vecs,info,pos=0,adj=NULL,lift=0,digits=3){
#
nodenum<-length(vecs[,1])
#
#remain<-data$remain
#if (!is.null(remain)){  #if we have cutted, cut also info
#   lenrem<-length(remain)
#   newinfo<-matrix(0,lenrem,1) 
#   for (i in 1:lenrem){
#      point<-remain[i]
#      newinfo[i]<-info[point]
#   }
#   info<-newinfo
##  orinodenum<-length(info)
##  newinfo<-matrix(0,orinodenum,1)
##  ind<-1
##  for (i in 1:orinodenum){
##     if (remain[i]==1){  
##        newinfo[ind]<-info[i]
##        ind<-ind+1
##     }
##  }
##  info<-newinfo[1:nodenum]
#}
##
infolocx<-matrix(nodenum,1)
infolocy<-matrix(nodenum,1)
#
for (i in 1:nodenum){
  infolocx[i]<-vecs[i,3]   #+(vecs[i,3]-vecs[i,1])/2  
  infolocy[i]<-vecs[i,2]+lift
}
info<-format(info,digits=digits)
text(infolocx,infolocy,info,pos,adj)
}


