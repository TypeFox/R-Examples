odc.clust <-
function(x,centers,cv.num=5,l2=-1,clus=kmeans,l2.idx=seq(-3, 3, by=6/20)) {  
   
  if(l2 == -1){   
     rlt = odc.optimallambda2(x, centers,cv.num,l2.idx)
   
     rlt.odc = odc.cv(x,centers,rlt$opt.lambda2)
     opt.lambda2 = rlt$opt.lambda2
   }
   else{
       rlt.odc = odc.cv(x,centers,l2)
       
   }
   res = clus(rlt.odc$Z,centers)
   if(l2==-1)
   return(list(res=res, opt.lambda2=opt.lambda2))
   else
   return(res)

}
