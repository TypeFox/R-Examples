recluster.region<-function(mat,tr=50,dist="simpson",method="ward",phylo=NULL, mincl=2,maxcl=3,rettree=FALSE,retmat=FALSE,retmemb=FALSE){
        res<-NULL       
        clusters<-maxcl-mincl+1
        mat2<-as.matrix(mat)    
        rows<-nrow(mat2)    
        tab2<-array(NA,dim=c(rows,tr,clusters))
        rownames(tab2)<-rownames(mat2)        
        rownames(mat2)<-c(1:rows)
        if(data.class(mat)=="dist"){dista<-as.dist(mat2)
                }else{
                dista<-recluster.dist(mat2, phylo=phylo, dist=dist)}
        dista<-as.matrix(dista)                
        for(i in 1:tr){ 
                samp<-sample(as.numeric(as.character(rownames(mat2))))
                dista2<-dista[samp,samp]
                if(method=="pam"){
                                  tree<-NULL}else{
                                  if (method=="diana"){tree<-diana(as.dist(dista2))}else
                                                      {tree<-hclust(as.dist(dista2),method=method)}
                                  }
                for (cut in mincl:maxcl){
                                  if(method=="pam"){
                                  cuttr<-pam(as.dist(dista2),k=cut)$clustering}
                                  else{                
                                  cuttr<-cutree (tree, k=cut)}
                      tab2[,i,(cut-mincl+1)]<-cuttr[order(as.numeric(as.character(rownames(dista2))))]
                      }
                }
        tree<-NULL
        cuttr<-NULL
        matrices<-array(NA,dim=c(rows,rows,clusters))
        for(sel in 1:clusters){
                tabsel<-tab2[,,sel]
                for(cl in 1:rows){
                        for(rw in cl:rows){
                                    vect<-round((tabsel[rw,]-tabsel[cl,])/(tabsel[rw,]-tabsel[cl,]+0.0001),0)
                                     matrices[rw,cl,sel]<-sum(vect,na.rm=T)/tr
                                     }
                        }
                }
        tabsel<-NULL
        if(retmemb){res$memb<-tab2}
        tab2<-NULL
        if(retmat){res$matrices<-matrices}
        pamsol<-matrix(data=NA, nrow=rows,ncol=clusters)
        colnames(pamsol)<-c(mincl:maxcl)
        rownames(pamsol)<-rownames(mat)
        res$solutions<-matrix(data=NA, nrow=clusters,ncol=3)
        colnames(res$solutions)<-c("k","silh","ex.diss")
        res$solutions[,1]<-c(mincl:maxcl)
        for (pamr in 1:clusters){
                  if(method=="pam"){pamsol[,pamr]<-pam(as.dist(matrices[,,pamr]),k=mincl-1+pamr)$clustering}
                          else{
                                   if (method=="diana"){
                                              pami<-diana(as.dist(matrices[,,pamr]))}else{
                                              pami<-hclust(as.dist(matrices[,,pamr]),method=method)}
                                  pamsol[,pamr]<-cutree(pami,k=mincl-1+pamr)
                                  if(rettree){res$tree[[pamr]]<-pami}                                                                    
                 }
                 pami<-NULL
                 res$solutions[pamr,3]<-recluster.expl(dista,pamsol[,pamr])
                 res$solutions[pamr,2]<-mean(silhouette(pamsol[,pamr],dista)[,3])
                 }
        res$grouping<-pamsol
        return(res)
}
