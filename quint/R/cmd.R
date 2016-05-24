cmd <-
function(pmat,dmats2){
#check mean difference per node condition
    selp1<-dmats2==1
    selp2<-dmats2==2
    #vec p1 is sum of nodes belonging to p1 with effect size or mean dif that is not > 0  
    vecp1<- sapply(1:nrow(dmats2),function(kk,pmat,sel){sum(pmat[sel[kk,],3]<=0)},pmat=pmat,sel=selp1)  
    #vecp2 is sum of nodes belonging to p2 with effect size or mean dif that is not > 0 
     vecp2<- sapply(1:nrow(dmats2),function(kk,pmat,sel){sum(-pmat[sel[kk,],3]<=0)},pmat=pmat,sel=selp2)
    condvec<-ifelse(vecp1==0 &vecp2==0,1,0)
    return(condvec)
        }
