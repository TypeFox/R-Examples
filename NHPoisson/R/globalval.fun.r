globalval.fun <-
function(mlePP,lint=NULL,nint=NULL, 
Xvar=NULL,namXvar=NULL, Xvart=NULL,namXvart=NULL, 
h=NULL, typeRes=NULL, typeResLV='Pearson',typeI='Disjoint', nsim=100,
clevel=0.95,resqqplot=FALSE, nintLP=100, tit='',
flow=0.5, addlow=FALSE,histWgraph=TRUE, plotDisp=c(2,2), indgraph=FALSE,scax=NULL, scay=NULL,
legcex=0.5, cores=1, xlegend='topleft',fixed.seed=NULL)
{
if ((is.null(lint))&(typeI=='Overlapping')) stop('Argument lint must be specified for Overlapping intervals')
if ((is.null(lint))&(is.null(nint))&(typeI=='Overlapping')) stop('one of arguments lint or nint must be specified for Disjoint intervals')
namcovariates<-dimnames(mlePP@covariates)$covnames
covariates<-mlePP@covariates
if(is.null(Xvar)) 
{
Xvar<-covariates
namXvar<-namcovariates
}

if (histWgraph==TRUE)	dev.new(record=TRUE)

t<-mlePP@t
if (is.null(Xvart)==FALSE) Xvart<-as.matrix(Xvart)
if (is.null(Xvar)==FALSE) Xvar<-as.matrix(Xvar)
posE<-mlePP@posE
if (is.null(tit)) tit<-mlePP@tit
inddat<-mlePP@inddat
posEH<-transfH.fun(mlePP)$posEH
res<-unifres.fun(posEH)
graphresU.fun(unires=res$unires, posE=posE, Xvariables=Xvar,
namXv=namXvar,tit=tit, flow=flow, addlow=addlow,histWgraph=FALSE, 
plotDisp=plotDisp,indgraph=indgraph)

if (typeI=='Disjoint')
RRes<-CalcResD.fun(mlePP=mlePP, h=h, nint=nint,lint=lint, typeRes=typeRes)
else
RRes<-CalcRes.fun(mlePP=mlePP, h=h, typeRes=typeRes,lint=lint)

graphrate.fun(objres=RRes, tit=tit,scax=scax,scay=scay, xlegend=xlegend,
histWgraph=FALSE)

graphres.fun(objres=RRes,Xvariables=Xvart,typeRes=typeRes,
namXv=namXvart,histWgraph=FALSE, plotDisp=plotDisp,addlow=addlow,tit=tit,
flow=flow)

graphResCov.fun( Xvar=Xvar,mlePP=mlePP, h=h, nint=nintLP, tit=tit, 
typeRes=typeResLV,namX=namXvar,histWgraph=FALSE, plotDisp=plotDisp)


if (resqqplot==TRUE) resQQplot.fun(nsim=nsim,objres=RRes, 
covariates=covariates, clevel=clevel,tit=tit,cores=cores,fixed.seed=fixed.seed,
histWgraph=FALSE)

return(RRes)
}
