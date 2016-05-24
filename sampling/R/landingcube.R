"landingcube" <-
function(X,pikstar,pik,comment=TRUE) 
# landing phase of the cube method
######################################################
{ 
# extraction of the non-integer values for the landing phase 
EPS=1e-11
p=dim(X)[2]
N=dim(X)[1]
liste=(pikstar>EPS & pikstar<(1-EPS)) 
pikland=pikstar[liste] 
Nland=length(pikland) 
Xland=array(X[liste,] ,c(Nland,p)) 
nland=sum(pikland) 
FLAGI=(abs(nland-round(nland))<EPS)
# construction of the list of possible samples 
if(comment==TRUE){
cat("\n\nBEGINNING OF THE LANDING PHASE\n")
cat("At the end of the flight phase, there remain ",Nland,"non integer probabilities","\n")
cat("The sum of these probabilities is ",nland,"\n")
cat("This sum is ")
if(FLAGI)  cat(" integer\n") else cat(" non-integer\n")
}
if(FLAGI) 
     {
      pikland=round(nland)*pikland/nland
      nland=round(nland);
      SSS= writesample(nland,Nland) 
      }
else 
      SSS= rbind( writesample(trunc(nland),Nland), 
           writesample(trunc(nland)+1,Nland) ) 
lll=nrow(SSS)
if(comment){
cat("The linear program will consider ",lll," possible samples\n")
} 
# computation of the cost 
Asmp=matrix(0,p,lll);
for(i in 1:lll)  {Asmp[,i]=t(Xland/pik[liste]) %*% (SSS[i,]-pikland) }
A=X[pik>EPS,]/pik[pik>EPS]
cost=rep(0,times=lll) 
for(i in 1:lll) 
cost[i]=t(Asmp[,i]) %*% ginv(t(A) %*% A) %*% Asmp[,i] 
# linear programming
V = t(cbind(SSS,rep(1,times=lll)))
b=c(pikland,1)
constdir=rep("==",times=(Nland+1))
x=lp("min",cost,V,constdir,b)$solution
# choice of the sample
u=runif(1,0,1)
i=0
ccc=0
while(ccc<u) {i=i+1;ccc=ccc+x[i]}
if(comment)
    {
     cat("The mean cost is ",mean(cost),"\n")
     cat("The smallest cost is ",min(cost),"\n")
     cat("The largest cost is ",max(cost),"\n")
     cat("The cost of the selected sample is",cost[i])
     }
pikfin=pikstar
pikfin[liste]=SSS[i,];
pikfin
}

