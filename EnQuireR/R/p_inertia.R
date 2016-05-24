#Pourcentages d'inertie
p_inertia=function(dataset){
acm.ref=MCA(dataset,graph=F)
p_iner=matrix(0,200,3)
for (k in 1:200){
    tab=dataset
    for (j in 1:ncol(dataset)){
    lev=levels(dataset[,j])
    proba=matrix(summary(dataset[,j]))[,1]/nrow(dataset)
    #tab[,j]=sample(lev,nrow(dataset),replace=T,prob=proba)
    tab[,j]=dataset[sample(1:nrow(dataset),nrow(dataset)),j]
    }
    colnames(tab)=colnames(dataset)
    res=MCA(tab,graph=F)
    p_iner[k,1]=res$eig[[2]][1]
    p_iner[k,2]=res$eig[[2]][2]
    p_iner[k,3]=p_iner[k,1]+p_iner[k,2]
}
colnames(p_iner)=c("Dim.1","Dim.2","Plan.1-2")

#resum=matrix(0,2,2)
#resum[1,1]=mean(p_iner[,1])    
#resum[2,1]=mean(p_iner[,2])
#resum[1,2]=sd(p_iner[,1])
#resum[2,2]=sd(p_iner[,2])
#rownames(resum)=c("Dim.1","Dim.2")
#colnames(resum)=c("Mean","Standard deviation")

resum=matrix(NA,3,2)
rownames(resum)=c("Dim.1","Dim.2","Plan.1-2")
colnames(resum)=c("% of variance","p-value")
resum[1,1]=acm.ref$eig[1,2]
resum[2,1]=acm.ref$eig[2,2]
resum[3,1]=acm.ref$eig[2,3]
resum[1,2]=length(which(p_iner[,1]>acm.ref$eig[1,2]))/200
resum[2,2]=length(which(p_iner[,2]>acm.ref$eig[2,2]))/200
resum[3,2]=length(which(p_iner[,3]>acm.ref$eig[2,3]))/200

res=list()
res$p_iner=p_iner
res$resum=resum           
return(res)
}

#p_iner=p_inertia(tea,c(1:21,23:36))

    
