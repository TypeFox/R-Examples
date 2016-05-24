"balancedstratification" <-
function(X,strata,pik,comment=TRUE,method=1)
{
strata=cleanstrata(strata)
H=max(strata)
N=dim(X)[1]
pikstar=rep(0,times=N)
for(h in 1:H) 
{
if(comment==TRUE) cat("\nFLIGHT PHASE OF STRATUM",h)
pikstar[strata==h]=fastflightcube(cbind(X[strata==h,],pik[strata==h]),pik[strata==h],1,comment) 
}
if(comment==TRUE) cat("\nFINAL TREATMENT")
XN=cbind(disjunctive(strata)*pik,X)/pik*pikstar
if(is.null(colnames(X))==FALSE)
    colnames(XN)<-c(paste("Stratum", 1:H, sep = ""),colnames(X))
samplecube(XN,pikstar,1,comment,method) 
}

