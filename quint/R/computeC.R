computeC <-
function(pmat,dmats3,w){
#compute criterion
    # dmats3: designmatrix with admissible assignments of the nodes to the partition classes
    #compute value of partitioning criterion
    selp1<-dmats3==1
    selp2<-dmats3==2
    weight<-pmat[,1]+pmat[,2]
    pmat<-cbind(pmat,weight)
     #weighted average of the effect sizes of the regions belonging to p1 for each possible partition (rows of dmats3)
    dif1<- sapply(1:nrow(dmats3),function(kk,pmat,sel){
    (t(pmat[sel[kk,],4])%*%pmat[sel[kk,],3])/sum(pmat[sel[kk,],4])  },pmat=pmat,sel=selp1)
         #weighted average of the effect sizes of the regions belonging to p1 for each possible partition (rows of dmats3)
    dif2<- sapply(1:nrow(dmats3),function(kk,pmat,sel){
    (t(pmat[sel[kk,],4])%*%-pmat[sel[kk,],3])/sum(pmat[sel[kk,],4]) },pmat=pmat,sel=selp2)
    #compute "Difference in treatment outcomes" component
    crit1<- sapply(1:length(dif1),function(kk,dif1,dif2,w){ w[1]*log(1+dif1[kk])+w[1]*log(1+dif2[kk])},
    dif1=dif1,dif2=dif2,w=w)
    #compute "Difference in cardinality" component
    crit2<-sapply(1:nrow(dmats3),function(kk,pmat,w,sel1,sel2){
         w[2]*log(sum(pmat[sel1[kk,],4]))+w[2]*log(sum(pmat[sel2[kk,],4]))},pmat=pmat,w=w,sel1=selp1,sel2=selp2)
    crittot<-apply(rbind(crit1,crit2),2,sum)
    return(list(crittot=crittot,critdif=crit1,critcard=crit2))}
