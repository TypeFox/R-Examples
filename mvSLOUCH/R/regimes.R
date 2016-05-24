.createRegimes<-function(PhylTree,vRegimes){
    lRegimes<-vector("list",5)
    names(lRegimes)<-c("regimes","regimeTypes","regimeTimes","regimeCoding","regimeVect")
    lRegimes$regimeVect<-vRegimes
    lRegimes$regimeTimes<-sapply(PhylTree@epochs,rev)
    vRegimeTypes<-unique(vRegimes)
    lRegimes$regimeTypes<-as.character(1:length(vRegimeTypes))
    lRegimes$regimeCoding<-cbind(vRegimeTypes,as.character(1:length(vRegimeTypes)))
    lRegs<-sapply(PhylTree@lineages,function(vLin,Regs){vRegs<-rep(NA,length(vLin)-1);if (length(vLin)>1){for(i in 1:length(vRegs)){vRegs[i]<-Regs[vLin[i]]}};rev(vRegs) },Regs=as.character(vRegimes),simplify=FALSE)
    lRegs<-lRegs[(PhylTree@nnodes-PhylTree@nterm+1):PhylTree@nnodes]
    lRegs<-sapply(lRegs,function(vReg,mCode){LinRegs<-sapply(vReg,function(r,mCode){mCode[which(mCode[,1]==r),2]},mCode=mCode,simplify=TRUE);names(LinRegs)<-NULL;LinRegs},mCode=lRegimes$regimeCoding,simplify=FALSE)
    lRegimes$regimes<-lRegs
    lRegimes
}
