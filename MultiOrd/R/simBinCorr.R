simBinCorr <-
function(ordPmat, CorrMat, nSim, steps=0.025){

conformity.Check(ordPmat, CorrMat)
p=find.binary.prob(ordPmat) # range of p is checked here
pvec=p$p

Mlocation=p$Mlocation
del.next=CorrMat
change=1
iteration=0
cat("calculating the intermediate binary correlations. \n"); 
while ( sum(change>0.001) >0) {
iteration=iteration+1
ep0 = generate.binary( nSim, pvec, del.next)
Mydata= BinToOrd(pvec, ordPmat, Mlocation, ep0)
if (iteration<15){del.next = del.next + ( CorrMat - Mydata$Corr )*0.9}
else {del.next = del.next + ( CorrMat - Mydata$Corr )*steps}
change = abs(CorrMat- Mydata$Corr)
if(iteration%%40==0){cat("\n")}
cat("."); 
}
cat("\n required ", iteration, " iterations to calculate intermediate binary correlations. \n"); 
return(list(del.next=del.next, Mlocation=Mlocation, pvec=pvec) )
}
