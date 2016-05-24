jointfit <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
PhiFix=NULL,LamFix=NULL,structure="correlated",mord=0,dord=1,convergence=1e-05,Init_Corr=NULL, EstimateCorrelations=TRUE, ZZCorr=NULL) {
    if (structure=="correlated") res<-jointfit_correlated(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations, ZZCorr=ZZCorr)
    if (structure=="shared") res<-jointfit_shared(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       PhiFix=NULL,LamFix=NULL,structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations)
    return(res)
}
