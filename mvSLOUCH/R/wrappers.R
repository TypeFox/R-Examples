BrownianMotionModel<-function(phyltree,data,predictors=NULL,M.error=NULL,calcCI=FALSE){
    .PhyloSDEestim(phyltree,data,kY=NULL,regimes=NULL,regimes.times=NULL,params=list(EvolModel="bm",EstimParams=list(calcCI=calcCI)),predictors=predictors,M.error=M.error)
}
ouchModel<-function(phyltree,data,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,M.error=NULL,Atype="Invertible",Syytype="UpperTri",calcCI=FALSE,diagA="Positive"){
    if (Atype=="DecomposableNegative"){diagA="Negative"}
    .PhyloSDEestim(phyltree,data,kY=NULL,regimes=regimes,regimes.times=regimes.times,root.regime=root.regime,params=list(EvolModel="ouch",EstimParams=list(Atype=Atype,Syytype=Syytype,diagA=diagA,diagSyy="Positive",calcCI=calcCI)),predictors=predictors,M.error=M.error)
}

mvslouchModel<-function(phyltree,data,kY,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,M.error=NULL,Atype="Invertible",Syytype="UpperTri",calcCI=FALSE,diagA="Positive"){
    if (Atype=="DecomposableNegative"){diagA="Negative"}
    .PhyloSDEestim(phyltree,data,kY=kY,regimes=regimes,regimes.times=regimes.times,root.regime=root.regime,params=list(EvolModel="mvslouch",EstimParams=list(Atype=Atype,Syytype=Syytype,diagA=diagA,diagSyy="Positive",calcCI=calcCI)),predictors=predictors,M.error=M.error)
}

SummarizeBM<-function(phyltree,data, modelParams,t=c(1),dof=NULL,M.error=NULL,predictors=NULL,calcCI=FALSE){
    .SummarizeFullPoint(NULL,dfData=data,PhylTree=phyltree,EvolModel="bm",EstimationParams=list(calcCI=calcCI),regimes.times=NULL,regimes=NULL,modelParams=modelParams,t=t,dof=dof,bShouldPrint=TRUE,LogLik=NULL,maxIter=NULL,tol=NULL,Merror=M.error,predictors=predictors)
}
SummarizeOUCH<-function(phyltree,data,modelParams,regimes=NULL,regimes.times=NULL,t=c(1),dof=NULL,M.error=NULL,predictors=NULL,Atype="Invertible",Syytype="UpperTri",calcCI=FALSE){
    .SummarizeFullPoint(NULL,dfData=data,PhylTree=phyltree,EvolModel="ouch",EstimationParams=list(calcCI=calcCI,Atype=Atype,Syytype=Syytype),regimes.times=regimes.times,regimes=regimes,modelParams=modelParams,t=t,dof=dof,bShouldPrint=TRUE,LogLik=NULL,maxIter=c(10,50),tol=c(0.001,0.0001),Merror=M.error,predictors=predictors)
}
SummarizeMVSLOUCH<-function(phyltree,data,modelParams,regimes=NULL,regimes.times=NULL,t=c(1),dof=NULL,M.error=NULL,predictors=NULL,Atype="Invertible",Syytype="UpperTri",calcCI=FALSE){
    .SummarizeFullPoint(NULL,dfData=data,PhylTree=phyltree,EvolModel="mvslouch",EstimationParams=list(calcCI=calcCI,Atype=Atype,Syytype=Syytype),regimes.times=regimes.times,regimes=regimes,modelParams=modelParams,t=t,dof=dof,bShouldPrint=TRUE,LogLik=NULL,maxIter=10,tol=c(0.001,0.0001),Merror=M.error,predictors=predictors)
}

simulBMProcPhylTree<-function(phyltree,X0,Sigma,dropInternal=TRUE,M.error=NULL,fullTrajectory=FALSE,jumpsetup=NULL){
    evolmodel<-"bm"
    if (fullTrajectory){evolmodel<-"bmStep"}
    Simulparams<-NULL
    if (!is.null(jumpsetup)){Simulparams$jump<-jumpsetup}
    res<-.simulVasicekProcPhylTree(phyltree,evolmodel,modelParams=list(vX0=X0,Sxx=Sigma),EstimationParams=NULL,Simulparams=Simulparams,dropInternal=dropInternal,M.error=M.error,bAllTrajectories=fullTrajectory)
    if (!fullTrajectory){res<-res$ExtantSample}
    res
}
simulOUCHProcPhylTree<-function(phyltree,modelParams,regimes=NULL,regimes.times=NULL,dropInternal=TRUE,M.error=NULL,fullTrajectory=FALSE,jumpsetup=NULL){
    evolmodel<-"ouch"
    if (fullTrajectory){evolmodel<-"ouchStep"}
    Simulparams<-NULL
    if (!is.null(jumpsetup)){Simulparams$jump<-jumpsetup}
    res<-.simulVasicekProcPhylTree(phyltree,evolmodel,modelParams=modelParams,EstimationParams=NULL,regimes=regimes,regimes.times=regimes.times,Simulparams=Simulparams,dropInternal=dropInternal,M.error=M.error,bAllTrajectories=fullTrajectory)
    if (!fullTrajectory){res<-res$ExtantSample}
    res
}

simulMVSLOUCHProcPhylTree<-function(phyltree,modelParams,regimes=NULL,regimes.times=NULL,dropInternal=TRUE, M.error=NULL,fullTrajectory=FALSE,jumpsetup=NULL){
    evolmodel<-"mvslouch"
    if (fullTrajectory){evolmodel<-"mvslouchStep"}
    Simulparams<-NULL
    if (!is.null(jumpsetup)){Simulparams$jump<-jumpsetup}
    res<-.simulVasicekProcPhylTree(phyltree,evolmodel,modelParams=modelParams,EstimationParams=NULL,regimes=regimes,regimes.times=regimes.times,Simulparams=Simulparams,dropInternal=dropInternal,M.error=M.error,bAllTrajectories=fullTrajectory)
    if (!fullTrajectory){res<-res$ExtantSample}
    res
}

.SimulStudyBM<-function(phyltree,X0,Sigma,n=NULL, M.error=NULL){
    .SimulStudy(phyltree,modelParams=list(vX0=X0,Sxx=Sigma),params=list(EvolModel="bm"),n=n,M.error=M.error)
}
.SimulStudyOUCH<-function(phyltree,modelParams,regimes=NULL,regimes.times=NULL,n=NULL, M.error=NULL){
    .SimulStudy(phyltree,modelParams,params=list(EvolModel="ouch"),regimes=regimes,regimes.times=regimes.times,n=n,M.error=M.error)
}
.SimulStudyMVSLOUCH<-function(phyltree,modelParams,regimes=NULL,regimes.times=NULL,n=NULL, M.error=NULL){
    .SimulStudy(phyltree,modelParams,params=list(EvolModel="mvslouch"),regimes=regimes,regimes.times=regimes.times,n=n,M.error=M.error)
}
