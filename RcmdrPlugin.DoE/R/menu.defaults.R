## 2-level factorial menu
.default.design2 <- list(nameVar=gettextRcmdr("Design.1"),
    nrunVar="8",nfacVar="4",nrepVar="1", nblockVar="1", ncenterVar="0",
    cbInitials = as.character(c(0,1,0,1,1,1,0,0,0,0)),
    level1Var="-1",level2Var="1",seedVar=as.character(sample(31999,1)),
    specialrbVariable="none",hardVar="0",designrbVariable="default",
    genVar="NULL", catlgVar="catlg", designVar="NULL",
    resVar="III",qualcritrbVariable="MA",
    comprclassVar="3: all interactions of group 1",
    facnamlist=c("A","B","C","D"),faclev1list=rep(-1,4),faclev2list=rep(1,4),faclablist=c("","","",""),
    estrbVariable="none",comprrbVariable="manual", maxtimeVar="60",est2fislist="",
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))
    
.default.design2pb <- list(nameVar=gettextRcmdr("Design.1"),
    nrunVar="8",nfacVar="4",nrepVar="1", ncenterVar="0",
    cbInitials = as.character(c(0,1,1,1,1,1,0,0,0,0)),
    level1Var="-1",level2Var="1",seedVar=as.character(sample(31999,1)),
    facnamlist=c("A","B","C","D"),faclev1list=rep(-1,4),faclev2list=rep(1,4),faclablist=c("","","",""),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))
    
## orthogonal arrays and numbers of runs are picked from list
## variable is 0-indexed
## colnums is multibox that is only available in case of idVar activated
.default.designoa <- list(nameVar=gettextRcmdr("Design.1"),
    idVar="NULL",nrunVar="NULL",nfacVar="3",nrepVar="1", minrdfVar=0, 
    optimVar="none",
    colnolist=c("","",""),
    cbInitials = as.character(c(0,1,0,1,1,1,0,0,0,0)),
    seedVar=as.character(sample(31999,1)),
    facnamlist=c("A","B","C"),nlevlist=rep("",3), faclevlist=rep("",3),faclablist=c("","",""),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"), fromVar="", toVar="")
    
.default.designfac <- list(nameVar=gettextRcmdr("Design.1"),
    nrunVar="",nfacVar="2",nrepVar="1",nblockVar="1",
    cbInitials = as.character(c(0,1,1,1,1,1,0,0,0,0)),
    seedVar=as.character(sample(31999,1)),
    facnamlist=c("A","B"),nlevlist=rep("",2), faclevlist=rep("",2),faclablist=c("",""),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))
    
.default.designlhs <- list(nameVar=gettextRcmdr("Design.1"),
    nrunVar="20",nfacVar="4",digitsVar="NULL",typerbVariable="optimum",
    cbInitials = as.character(c(0,1,1,1,1,1,0,0,0,0)),
    level1Var="0",level2Var="1",seedVar=as.character(sample(31999,1)),
    facnamlist=c("A","B","C","D"),faclev1list=rep(0,4),faclev2list=rep(1,4),faclablist=c("","","",""),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))

.default.designaugmentlhs <- list(nameVar=gettextRcmdr("Design.1"),
    mVar="1", typerbVariable="optAugment",
    cbInitials = as.character(c(0,1,1,1,1,1,0,0,0,0)),
    seedVar=as.character(sample(31999,1)),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))

.default.designbbd <- list(nameVar=gettextRcmdr("Design.1"),
    nfacVar="4",ncenterVar="4",blockVar="NULL",
    cbInitials = as.character(c(0,1,1,1,1,1,0,0,0,0)),
    level1Var="-1",level2Var="1",seedVar=as.character(sample(31999,1)),
    facnamlist=c("A","B","C","D"),faclev1list=rep(-1,4),faclev2list=rep(1,4),faclablist=c("","","",""),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))
    
.default.designccd <- list(nameVar=gettextRcmdr("Design.1"),
    alphaVar="orthogonal",ncenterVar="4",
    cbInitials = as.character(c(0,1,1,1,1,1,0,0,0,0)),
    seedVar=as.character(sample(31999,1)),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))
    
.default.designDopt <- list(nameVar=gettextRcmdr("Design.1.Dopt"),
    nrunsVar="8",
    rhsVariable="",constraintVar=gettextRcmdr("<all candidate data set rows eligible>"),
    nrepeatVar="5",
    cbInitials = as.character(c(0,1,1,1,1,1,0,0,0,0)),
    seedVar=as.character(sample(31999,1)),
    etyperbVariable="none", decimalrbVariable = "default", 
    dirVar=getwd(), fileVar=gettextRcmdr("Design.1"))
    