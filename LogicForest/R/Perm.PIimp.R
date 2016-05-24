Perm.PIimp <-
function(resp, fullPIdata, mod, ntrees, wts, nperm)
{
 All.PInms<-colnames(fullPIdata)
 massPI.change<-matrix(0, nrow=nperm, ncol=length(All.PInms))
 colnames(massPI.change)<-All.PInms
 for (t in 1:nperm)
  { 
  n.obs<-nrow(fullPIdata) 
  nCV<-length(mod)
  test.ids<-CV.ids(n=n.obs, kfold=nCV)
  loc.ids<-CV.ids(n=ntrees, kfold=nCV)
  PIchange<-matrix(0, nrow=ntrees, ncol=length(All.PInms))
  colnames(PIchange)<-All.PInms
  for(i in 1:nCV)
   {
   ids<-test.ids[[i]]
   chg.ids<-loc.ids[[i]]
   alphas<-wts[[i]]
   TR.resp<-resp[[i]]
   c.PIdat<-mod[[i]]$PI.TSdata
   for(k in 1:length(c.PIdat))
     {
     spot<-chg.ids[k]
     alpha<-alphas[k]
     tree.PIdat<-c.PIdat[[k]]
     tree.PInms<-colnames(tree.PIdat)
     if (dim(tree.PIdat)[2]==1) {tree.PIdat<-as.vector(tree.PIdat)}
     loc<-which(All.PInms%in%tree.PInms)
     if(is.vector(tree.PIdat))
       {Otree.pred<-tree.PIdat
        O.treemiss<-sum(abs(TR.resp-Otree.pred))
        p.ids<-sample(1:length(tree.PIdat), length(tree.PIdat), replace=F)
        Ptree.pred<-tree.PIdat[p.ids]
        P.treemiss<-sum(abs(TR.resp-Ptree.pred))
        DeltaMiss<-P.treemiss-O.treemiss
        PIchange[spot,loc]<-alpha*DeltaMiss
        }
     if(is.matrix(tree.PIdat))
        {Otree.pred<-ifelse(rowSums(tree.PIdat)>0, 1, 0)
         O.treemiss<-sum(abs(TR.resp-Otree.pred))
         nobs<-length(Otree.pred)
         p.ids<-sample(1:nobs, nobs, replace=F)
         nPIs<-ncol(tree.PIdat)
         for (g in 1:nPIs)
           {
           cloc<-loc[g]
           PPI<-tree.PIdat[p.ids,g]
           if(g==1) {Ptree.PIdat<-cbind(PPI, tree.PIdat[,2:nPIs])}
           if(g>1 & g<nPIs) {Ptree.PIdat<-cbind(tree.PIdat[,1:(g-1)],PPI, tree.PIdat[,(g+1):nPIs])}
           if(g==nPIs) {Ptree.PIdat<-cbind(tree.PIdat[,1:(g-1)],PPI)}
           Ptree.pred<-ifelse(rowSums(Ptree.PIdat)>0, 1, 0)
           P.treemiss<-sum(abs(TR.resp-Ptree.pred))
           DeltaMiss<-P.treemiss-O.treemiss
           PIchange[spot,cloc]<-alpha*DeltaMiss
           }
         }
      }
   }
  PIimp<-colSums(PIchange)/nrow(PIchange)
  massPI.change[t,]<-PIimp
  }
 PI.importance<-colMeans(massPI.change)
 PI.importance
}
