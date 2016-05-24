resSim.fun <-
function(i,lambda,covariates, beta,lint,t=NULL,
 tind=TRUE,typeI='Disjoint', typeRes='Pearson',h=NULL)
{
posNH<-simNHP.fun(lambda=lambda)$posNH
if (is.null(t)) t<-c(1:length(lambda))

mod<-fitPP.fun(covariates=covariates, start=beta, posE=posNH,  tind=tind,
tim=t,modCI='FALSE', tit="",modSim=TRUE,dplot='FALSE')


if (typeI=='Disjoint')
{
ResAux<-CalcResD.fun(mlePP=mod, h=h, nint=NULL,lint=lint, 
typeRes=typeRes,modSim='TRUE')

}
else ResAux<-CalcRes.fun(mlePP=mod,  h=h, lint=lint,typeRes=typeRes)

if (typeRes=='Raw') res<-ResAux$RawRes
else res<-ResAux$ScaRes$ScaRes
return(res)
}
