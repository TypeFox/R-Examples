`sfl2sl` <-
function(sformlist,exclude=NULL,eventlist=NULL){
   sformstats<-list()

   #First two dimensions can be taken from the fixed values of the rows and columns of the
   #sformlist. This is done for optimization purposes. Third dimension is pulled from the actor dimension and is equal to the number
   # of valid events
   
   sfl.dim12<-c(length(sformlist[[1]][1,][[1]]),length(sformlist[[1]][1,]))

   sfl.rownames<-colnames(sformlist[[1]][1,][[1]])
   sfl.colnames<-names(sformlist[[1]][1,])
   for(i in 1:length(sformlist)){
      flm<-list()
      for(k in 1:length(sformlist[[i]][,1])){
      flm[[k]]<-matrix(unlist(sformlist[[i]][k,]),nrow=sfl.dim12[1])}
      flm.ar<-array(unlist(flm), dim = c(dim(flm[[1]]), length(flm)),dimnames=list(sfl.rownames,sfl.colnames,1:length(sformlist[[i]][,1])))
      sformstats[[i]]<-aperm(flm.ar,c(3,1,2))
   }
   if(!is.null(exclude)){
      b<-glapply(eventlist[,1],eventlist[,2],function(x) which(x%in%exclude),regroup=FALSE)
      bv<-lapply(b,function(x) length(x))
      bvd<-which(unlist(bv)>0)
      bvdout<-b[bvd]
      for(i in 1:length(bvd)){
         sformstats[[bvd[i]]]<-sformstats[[bvd[i]]][-bvdout[[i]],,]
         }
      for(i in 1:length(sformstats)){
         dimnames(sformstats[[i]])[[1]]<-1:nrow(sformstats[[i]][,,1])
         }
      }
  sformstats
}

