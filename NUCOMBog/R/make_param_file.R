make_param_file<-function(WD,par=NULL){
  values<-c(350,0.014,0.014,0.014,250,150,10,50,150,0.2,0.2,0.2,0.99,0.98,0.85,0.7,0.4,0.2,0.2,0.1,37,16,25,35,16,25,35,25,30,35,0.7,0.4,0.4,0.2,0.9,0.24,0.24,0.24,0.13,0.13,0.13,0.00004,0.00004,0.039,0.5,70,0.08,0,0,0,0.06,0,0,0,0.1,0.45,0.4,0.15,0.03,0.02,0.45,0,0.3,0.012,25,2,12,20,1000,-100,0,400,16,100,70,16,100,70,25,100,70,0.7,0.35,0.3,0.35,0.7,0.17,0.17,0.17,0.08,0.06,0.07,0.00004,0.00004,0.00004,0.8,60,0.04,0,0,0,0.01,0.04,0,0,0,0.45,0.25,0.3,0.022,0.014,0.25,0,0.1,0.012,25,5,14,20,600,50,100,300,9,20,40,0.7,1,1.4,0.02,0.02,0.00004,45,0,0.04,0,1,0.032,0.005,0.4,20,5,14,17,500,0,50,200,7,13,25,0.7,1,1.3,0.04,0.04,0.00004,50,0,0.04,0,1,0.036,0.006,0.4,20,3,14,17,350,0,50,150,5,15,30,0.7,1,1.2,0.08,0.08,0.00004,60,0,0.08,0,1,0.042,0.008,0.4,20,0,10,17,200,-100,0,50)
  names<-c("CO2ref","CritNClvms","CritNCacro","CritNCcato","FixDepthtcnp","FixDepthlcnp","MinDepthlvms","MaxDepthlvms","MaxDepthacro","MicrEfflvms","MicrEffacro","MicrEffcato","Mutcnp","Mulcnp","Mulvms","Muacro","Mucato","TDecParlvms","TDecParacro","TDecParcato","ToptDec","gram_BDLvmsleaf","gram_BDLvmsstem","gram_BDLvmsroot","gram_BDAcroleaf","gram_BDAcrostem","gram_BDAcroroot","gram_BDCatoleaf","gram_BDCatostem","gram_BDCatoroot","gram_Beta","gram_CAllocFrleaf","gram_CAllocFrstem","gram_CAllocFrroot","gram_CropFc","gram_DecParLvmsleaf","gram_DecParLvmsstem","gram_DecParLvmsroot","gram_DecParAcroleaf","gram_DecParAcrostem","gram_DecParAcroroot","gram_DecParCatoleaf","gram_DecParCatostem","gram_DecParCatoroot","gram_KExt","gram_MaxGr","gram_MortFrLvmsleaf","gram_MortFrLvmsstem","gram_MortFrLvmsroot","gram_MortFrAcroleaf","gram_MortFrAcrostem","gram_MortFrAcroroot","gram_MortFrCatoleaf","gram_MortFrCatostem","gram_MortFrCatoroot","gram_NAllocFrleaf","gram_NAllocFrstem","gram_NAllocFrroot","gram_NConcMax","gram_NconcMin","gram_ReallFrleaf","gram_ReallFrstem","gram_ReallFrroot","gram_SLA","gram_TMaxGr","gram_TMinGr","gram_TOpt1Gr","gram_TOpt2Gr","gram_WLMax","gram_WLMin","gram_WLOpt1","gram_WLOpt2","eric_BDLvmsleaf","eric_BDLvmsstem","eric_BDLvmsroot","eric_BDAcroleaf","eric_BDAcrostem","eric_BDAcroroot","eric_BDCatoleaf","eric_BDCatostem","eric_BDCatoroot","eric_Beta","eric_CAllocFrleaf","eric_CAllocFrstem","eric_CAllocFrroot","eric_CropFc","eric_DecParLvmsleaf","eric_DecParLvmsstem","eric_DecParLvmsroot","eric_DecParAcroleaf","eric_DecParAcrostem","eric_DecParAcroroot","eric_DecParCatoleaf","eric_DecParCatostem","eric_DecParCatoroot","eric_KExt","eric_MaxGr","eric_MortFrLvmsleaf","eric_MortFrLvmsstem","eric_MortFrLvmsroot","eric_MortFrAcroleaf","eric_MortFrAcrostem","eric_MortFrAcroroot","eric_MortFrCatoleaf","eric_MortFrCatostem","eric_MortFrCatoroot","eric_NAllocFrleaf","eric_NAllocFrstem","eric_NAllocFrroot","eric_NConcMax","eric_NconcMin","eric_ReallFrleaf","eric_ReallFrstem","eric_ReallFrroot","eric_SLA","eric_TMaxGr","eric_TMinGr","eric_TOpt1Gr","eric_TOpt2Gr","eric_WLMax","eric_WLMin","eric_WLOpt1","eric_WLOpt2","humm_BDLvmsshoot","humm_BDAcroshoot","humm_BDCatoshoot","humm_Beta","humm_CAllocFrshoot","humm_CropFc","humm_DecParLvmsshoot","humm_DecParAcroshoot","humm_DecParCatoshoot","humm_MaxGr","humm_MortFrLvmsshoot","humm_MortFrAcroshoot","humm_MortFrCatoshoot","humm_NAllocFrshoot","humm_NConcMax","humm_NConcMin","humm_ReallFrshoot","humm_TMaxGr","humm_TMinGr","humm_TOpt1Gr","humm_TOpt2Gr","humm_WLMax","humm_WLMin","humm_WLOpt1","humm_WLOpt2","lawn_BDLvmsshoot","lawn_BDAcroshoot","lawn_BDCatoshoot","lawn_Beta","lawn_CAllocFrshoot","lawn_CropFc","lawn_DecParLvmsshoot","lawn_DecParAcroshoot","lawn_DecParCatoshoot","lawn_MaxGr","lawn_MortFrLvmsshoot","lawn_MortFrAcroshoot","lawn_MortFrCatoshoot","lawn_NAllocFrshoot","lawn_NConcMax","lawn_NConcMin","lawn_ReallFrshoot","lawn_TMaxGr","lawn_TMinGr","lawn_TOpt1Gr","lawn_TOpt2Gr","lawn_WLMax","lawn_WLMin","lawn_WLOpt1","lawn_WLOpt2","holl_BDLvmsshoot","holl_BDAcroshoot","holl_BDCatoshoot","holl_Beta","holl_CAllocFrshoot","holl_CropFc","holl_DecParLvmsshoot","holl_DecParAcroshoot","holl_DecParCatoshoot","holl_MaxGr","holl_MortFrLvmsshoot","holl_MortFrAcroshoot","holl_MortFrCatoshoot","holl_NAllocFrshoot","holl_NConcMax","holl_NConcMin","holl_ReallFrshoot","holl_TMaxGr","holl_TMinGr","holl_TOpt1Gr","holl_TOpt2Gr","holl_WLMax","holl_WLMin","holl_WLOpt1","holl_WLOpt2")
  allpar<-data.frame(names,values)
  allpar$names<-as.character(allpar$names)

  if(!is.null(par)){
    paruserprovided<-which(par$names %in% allpar$names)
    alluserprovided<-which(allpar$names %in% par$names)
    allpar$values[alluserprovided]<-par$values[paruserprovided]
  }

  file<-paste(WD,"input/param.txt",sep="")
  cat("[bog]",file=file,sep="\n")
  cat("TimeStep = 1",file=file,sep="\n",append=TRUE)

  for(i in 1:nrow(allpar)){
    if(i>=23 && i<=nrow(allpar)){
      if(i == 23){
        cat("[gram]",file=file,sep="\n",append=TRUE)}
      if(i == 73){
        cat("[eric]",file=file,sep="\n",append=TRUE)}
      if(i == 124){
        cat("[humm]",file=file,sep="\n",append=TRUE)}
      if(i == 149){
        cat("[lawn]",file=file,sep="\n",append=TRUE)}
      if(i == 174){
        cat("[holl]",file=file,sep="\n",append=TRUE)}
      cat(paste(substr(allpar[i,1],start=6,stop=nchar(allpar[i,1])),"=",allpar[i,2],sep=""),file=file,append=TRUE,sep="\n")
    }

    if(i<=21){
      cat(paste(allpar[i,1],"=",allpar[i,2],sep=""),file=file,append=TRUE,sep="\n")
    }

  }
}

