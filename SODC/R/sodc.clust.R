sodc.clust <-
function(x,centers,l1=-1,l2=-1,cv.num=5,clus=kmeans,boot.num=20,l2.idx=seq(-3, 3, by=6/20),l1.idx=seq(-3, 3, by=6/20)) {
   if(l2==-1){
   
      rlt=odc.optimallambda2(x, centers,cv.num,l2.idx)

      l2=rlt$opt.lambda2
      
      
   }
   if(l1==-1){

    rlt = sodc.optimallambda1.boot.all(x, centers, boot.num,l1.idx)
    l1=rlt$opt.lambda1
   }
   
           
   rlt.sodc = my.lasso.classify( x,centers,l1,l2)
   
   cl=NULL
   clvar=NULL
   if(length(unique(rlt.sodc$Z))>=centers ){
           
      cl=clus(rlt.sodc$Z,centers)
            
      clvar=clus(x[,rlt.sodc$varset],centers) 
      
   }
   
   return(list(cl=cl, clvar=clvar, opt.lambda1=l1, opt.lambda2=l2))
}
