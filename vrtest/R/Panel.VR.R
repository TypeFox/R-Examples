Panel.VR <-
function(dat,nboot=500){
k=ncol(dat)
vrstat=matrix(NA,nrow=k)

for (i in 1:k) vrstat[i,]=Auto.VR(dat[,i])$stat
vr1 = max(abs(vrstat)); vr2 = sum(vrstat^2);vr3 = sqrt(k)*mean(vrstat);

stats1=matrix(NA,nrow=nboot); stats2=matrix(NA,nrow=nboot); stats3=matrix(NA,nrow=nboot); 
for (ii in 1:nboot){
ys = dat * Mammen(nrow(dat));
vrstats=matrix(NA,nrow=k)
    for (jj in 1:k) vrstats[jj,]=Auto.VR(ys[,jj])$stat
    stats1[ii,]=max(abs(vrstats)); 
    stats2[ii,]=sum(vrstats^2);
    stats3[ii,]=sqrt(k)*mean(vrstats);
}
pboot1 = mean( stats1 > vr1 ); pboot2 = mean( stats2 > vr2 ); pboot3 = mean( abs(stats3) > abs(vr3) );

return(list(MaxAbs.stat=vr1,SumSquare.stat=vr2,Mean.stat=vr3,MaxAbs.pval=pboot1,SunSquare.pval=pboot2,Mean.pval=pboot3))
}
