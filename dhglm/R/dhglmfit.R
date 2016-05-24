dhglmfit <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
PhiFix=NULL,LamFix=NULL,mord=0,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06) {
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau))
    res<-dhglmfit_run(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=3)
   return(res)
}
