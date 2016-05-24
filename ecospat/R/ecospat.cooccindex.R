ecospat.co_occurrences <- function (data)
{

#######################
## Co-occurrences 2m ##
#######################

#---------------------------------------------------------------------------------
#This part calculate the matrix P1 of co-occurrence realized on co-occurrences possible
#at a resolution of 2 meters for each pair of species. These results are also put in
#a vector L1, where each comparison is only kept once.
#---------------------------------------------------------------------------------

Nsp1=ncol(data)-4
Nrel1=nrow(data)
P1_n_occ = t(as.matrix(data[,-(1:4)]))%*%as.matrix(data[,-(1:4)])
Tot1_sp = apply(data[,-(1:4)],MARGIN=2,sum)
Min1_n = outer(Tot1_sp, Tot1_sp,FUN=pmin)
P1_n = P1_n_occ / Min1_n

#----------------------------------------------------------------------------------
#This part reunite in columns all the values of co-occurrences relative to each
#species. Hence, one comparison between two pairs of species occur in two columns.
#----------------------------------------------------------------------------------

L1_n = P1_n[lower.tri(P1_n)]
Sp1_n = rep(1:(Nsp1-1),times=((Nsp1-1):1))
Sp2_n = t(outer(1:Nsp1,1:Nsp1,FUN="pmax"))[lower.tri(t(outer(1:Nsp1,1:Nsp1,FUN="pmax")))]
#Sumdata_n = data.frame(cbind(Sp1_n,Sp2_n,L1_n))

#----------------------------------------------------------------------------------
#This part below reunite in columns all the values of co-occurrences relative to each
#species. There is one column by species. Hence, one comparison between two pairs of
#species occur in two columns.
#----------------------------------------------------------------------------------

CoobySp1_n = matrix(data=0, nrow=Nsp1-1, ncol=Nsp1)
CoobySp1_n[lower.tri(CoobySp1_n,diag=T)]=P1_n[lower.tri(P1_n)]
CoobySp1_n[upper.tri(CoobySp1_n)]=P1_n[upper.tri(P1_n)]
CoobySp1_n<-data.frame(cbind(CoobySp1_n))
#colnames(CoobySp1_n)<-Names
boxplot(CoobySp1_n)

return (P1_n)
}
