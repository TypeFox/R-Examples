FeatureEEG <-
function(data,classes.Id,rec.Id,nselec="default",features="default",
Alpha=0.05, AlphaCorr=0.05,minacc=0.7, Nfea=10, fast = FALSE){

  reps<-rec.Id
  classes<-classes.Id
  
	featype<-c()
	ncomps<-c()
	win<-c()
	stat<-c()
	power<-c()
	abs<-c()
	log<-c()
	mintomax<-c()
  
  mapping<-c("f1","f2","f3","f4","f5","f6","f7")
  mapping<-cbind(mapping,c("windowing", "loadings","windowing of spectrum", "loadings of spectrum",
                           "windowing of principal components", "double windowing", "CWT"))


  if(class(features)=="Features")
  {
  	#### Retrieving parameters of windowing
  	N<-features$feaStatWindowing$N
     if(N>0)
      {
	     featype<-c(featype, rep("f1",N))
	     win<-c(win,features$feaStatWindowing$win)
	     stat<-c(stat,features$feaStatWindowing$stat)
	     power<-c(power,features$feaStatWindowing$power)
	     abs<-c(abs,features$feaStatWindowing$abs)
	     log<-c(log,features$feaStatWindowing$log)
	     mintomax<-c(mintomax,features$feaStatWindowing$mintomax)
      }

	   #### Retrieving parameters of windowing with spec
	    N<-features$feaSpecStatWindowing$N
      if(N>0){
	      featype<-c(featype, rep("f3",N))
	      win<-c(win,features$feaSpecStatWindowing$win)
	      stat<-c(stat,features$feaSpecStatWindowing$stat)
	      power<-c(power,features$feaSpecStatWindowing$power)
	      abs<-c(abs,features$feaSpecStatWindowing$abs)
	      log<-c(log,features$feaSpecStatWindowing$log)
	      mintomax<-c(mintomax,features$feaSpecStatWindowing$mintomax)
      }


	    #### Retrieving parameters of double windowing
	    N<-features$feaDoubleStatWindowing$N
      if(N>0)
      {
  	    featype<-c(featype, rep("f6",N))
	      win<-c(win,features$feaDoubleStatWindowing$win)
  	    stat<-c(stat,features$feaDoubleStatWindowing$stat)
  	    power<-c(power,features$feaDoubleStatWindowing$power)
  	    abs<-c(abs,features$feaDoubleStatWindowing$abs)
  	    log<-c(log,features$feaDoubleStatWindowing$log)
  	    mintomax<-c(mintomax,features$feaDoubleStatWindowing$mintomax)
      }



    	#### Retrieving parameters of cwt
	    N<-features$feaCWT$feaCWT
      if(N!=0) 
      {
          featype<-c(featype, "f7")
	        wavelet<-features$feaCWT$wavelet
	        variance<-features$feaCWT$variance
      }
	    #### Retrieving parameters of PCA
	    ncomps<-features$feaPCA$ncomps
	    if(ncomps>0)
	    {
		    if(features$feaPCA$feaLoadingsPCA)	featype<-c(featype, "f2")
		    if(features$feaPCA$feaSpecLoadingsPCA)	featype<-c(featype, "f4")
    
    
		    if(features$feaPCA$feaSignalsPCA>0)	
      {
		    featype<-c(featype, rep("f5",features$feaPCA$feaSignalsPCA))
		    win<-c(win,features$feaPCA$win)
		    stat<-c(stat,features$feaPCA$stat)
		    power<-c(power,features$feaPCA$power)
		    abs<-c(abs,features$feaPCA$abs)
		    log<-c(log,features$feaPCA$log)
		    mintomax<-c(mintomax,features$feaPCA$mintomax)
		  }  
	    }

  } else if(features[1]=="default"){
    featype=  c("f1",   "f1",   "f1",   "f1",   "f1",   "f1",  "f1",   "f1",   "f1",   "f1",   "f1",   "f1",  "f1",   "f1",   "f1",   "f1",  "f1",   "f2",   "f3",  "f3",   "f3",   "f3",   "f3",   "f3",   "f3",   "f3",   "f3",   "f3",    "f4",   "f5",   "f5",   "f5",   "f5",     "f6"        , "f6"        , "f6",        "f7" )
    stat<-    c("sum",  "sum",  "sum",  "sum",  "sum" , "sum", "sum",  "sum",  "sum",  "sum",  "sum" , "sum", "sum",  "max",  "max",  "max", "max",          "sum", "sum",  "sum",  "sum",  "sum",  "sum" , "max",  "max",  "max",  "max",           "sum",  "sum",  "max",  "max",    "sum","var" , "sum","max" , "max","sum")
    power<-   c(  2,      2,      2,      2,      2,      2,     2,      2,      2,      2,      2,      2,     1,     2,      2,      2,      2,             1,    1,      1,      1,      1,      1,       1,      1,      1,      1,               2,       2,      2,      2,         2  ,1    ,  2, 1       ,  2, 1)
    mintomax<-c(TRUE,   FALSE,  TRUE,   FALSE,  TRUE,   FALSE,  TRUE,   FALSE,  TRUE,   FALSE,  TRUE,   FALSE, TRUE,  FALSE,  TRUE,   FALSE,  TRUE,           TRUE,  FALSE,  TRUE,   FALSE,  TRUE,   FALSE,  FALSE,  TRUE,   FALSE,  TRUE,            TRUE,   TRUE,   TRUE,   TRUE,    FALSE, TRUE ,  FALSE,TRUE ,  FALSE, TRUE)
    win<-     c( 5,       15,     30,     5,      15,     30,   5,       15,     30,     5,      15,     30,    15 ,   15,     15,      15,   15,             5,      15,     30,    5,      15,     30,      15,     15,      15,    15,             15,     15,      15,    15,      15,10       ,   10, 10    ,  10,10)
    abs<-     c(FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE, FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE, FALSE,  FALSE,  FALSE,  FALSE,          FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE, FALSE,  FALSE,  FALSE,  FALSE,           FALSE,  FALSE,  FALSE,  FALSE,   FALSE ,FALSE,  FALSE,FALSE,  FALSE,FALSE)
    log<-     c(FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE, TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   FALSE, FALSE,  TRUE,   TRUE,   FALSE,         FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE,  TRUE,   TRUE,   FALSE,           FALSE,  TRUE,   TRUE,   FALSE,   TRUE , FALSE,  FALSE,FALSE,  FALSE,FALSE)
    wavelet<-"gaussian2"
    variance<-1
    ncomps<-3
  } else if(features[1]=="example"){
    featype=  c("f1",  "f3", "f3")
    stat<-    c("sum",  "sum", "sum")
    power<-   c(  2,      2, 2)
    mintomax<-c(TRUE,   FALSE,  TRUE)
    win<-     c( 5,       10,     10)
    abs<-     c(FALSE,  FALSE,  FALSE)
    log<-     c(FALSE,  FALSE,  FALSE)
    wavelet<-"gaussian2"
    variance<-1
    ncomps<-1   
  } else stop("Invalid parameter format: features")
	
	Alpha<-rep(Alpha,length(featype)*2)
	minacc<-rep(minacc,length(featype)*2)
	
  
  #parametros basicos
  if (is.null(ncol(data))) {nel<-1 ; data<-as.matrix(data)} else nel<-ncol(data)
  which.classes<-unique(classes)
  nclass <- length(which.classes)
  if (nclass!=2) stop("Vector 'classes' must indicate only 2 classes.")
  classes[classes==which.classes[1]]<-1
  classes[classes==which.classes[2]]<-2
  
  nrep <-  c(length(unique(reps[which(classes==1)])),length(unique(reps[which(classes==2)])))
  NR<-sum(nrep)
  L<-nrow(data)/(NR)
  
  if(nselec[1]=="default")
  {
    nselec=c(floor(nrep[1]*0.8),floor(nrep[2]*0.8))
  }
  
  samp1<-sample(1:nrep[1],nselec[1])
  samp2<-sample(1:nrep[2],nselec[2])
  #samp1<-1:6
  #samp2<-1:6

  
  testdata<-rbind(data[which(reps%in%c(1:nrep[1])[-samp1] & classes==1),],
  data[which(reps%in%c(1:nrep[2])[-samp2] & classes==2),])
  testclasses <- c(rep(1,L*(nrep[1]-nselec[1])),rep(2,L*(nrep[2]-nselec[2])))
  testreps <- c(reps[1:(L*(nrep[1]-nselec[1]))],reps[(L*nrep[1]+1):(L*(nrep[1]+nrep[2]-nselec[2]))])
  
  data<-rbind(data[which(reps%in%samp1 & classes==1),],
  data[which(reps%in%samp2 & classes==2),])
  classes <- c(rep(1,L*nselec[1]),rep(2,L*nselec[2]))
  reps <- c(reps[1:(L*nselec[1])],reps[(L*nrep[1]+1):(L*(nrep[1]+nselec[2]))])
  
  nreptest<-c(length(unique(testreps[which(testclasses==1)])),length(unique(testreps[which(testclasses==2)])))
  nrep <-  nselec
  NR<-sum(nrep)
  NRtest <- sum(nreptest)
  
  print("Calculating the features.")
  
  #####SPECTRO
  if (sum(c("f3","f4")%in%featype)>0) {
    r<-.SpecData(data,reps,classes,NR,nrep,nel,L)
    specdata <- r$specdata
    specclasses <- r$specclasses
    specreps <- r$specreps
    specL <- r$L0
    
    r<-.SpecData(testdata,testreps,testclasses,NRtest,nreptest,nel,L)
    testspecdata <- r$specdata
    testspecclasses <- r$specclasses
    testspecreps <- r$specreps
    testspecL <- r$L0
  }
  
  #####WAVELET
  if (c("f7")%in%featype) {
    r<-.CwtData(data,reps,classes,NR,nrep,nel,L,wavelet,variance)
    cwtdata <- r$cwtdata
    cwtclasses <- r$cwtclasses
    cwtreps <- r$cwtreps
    cwtL <- r$L0
    
    r<-.CwtData(testdata,testreps,testclasses,NRtest,nreptest,nel,L,wavelet,variance)
    testcwtdata <- r$cwtdata
    testcwtclasses <- r$cwtclasses
    testcwtreps <- r$cwtreps
    testcwtL <- r$L0
  }



  ####PCA
  if (sum(c("f2","f5")%in%featype)>0) {
    r<-.PcaData(data,reps,classes,ncomps,nrep, nel,NR,L)
    pcadata<-r$pcadata
    loaddata <- r$loaddata
    loadclasses <- r$loadclasses
    loadreps <- r$loadreps
    pcaclasses <- r$pcaclasses
    pcareps <- r$pcareps
    
    
    r<-.PcaData(testdata,testreps,testclasses,ncomps,nreptest, nel,NRtest,L)
    testpcadata<-r$pcadata
    testloaddata <- r$loaddata
    testloadclasses <- r$loadclasses
    testloadreps <- r$loadreps
    testpcaclasses <- r$pcaclasses
    testpcareps <- r$pcareps
  }

  
  ####PCA DO SPEC
  if (c("f4")%in%featype) {
    r<-.PcaData(specdata,specreps,specclasses,ncomps,nrep, nel,NR,specL)
    specpcadata<-r$pcadata
    specloaddata <- r$loaddata
    specloadclasses <- r$loadclasses
    specloadreps <- r$loadreps
    specpcaclasses <- r$pcaclasses
    specpcareps <- r$pcareps
    
    r<-.PcaData(testspecdata,testspecreps,testspecclasses,ncomps,nreptest, nel,NRtest,testspecL)
    testspecpcadata<-r$pcadata
    testspecloaddata <- r$loaddata
    testspecloadclasses <- r$loadclasses
    testspecloadreps <- r$loadreps
    testspecpcaclasses <- r$pcaclasses
    testspecpcareps <- r$pcareps
  }
  
  contWin <-0
  contAlpha <-0
  Wout <- c()
  ParsDiscart <- c()
  FinalFea <-c()
  typeDiscart<-c()

  print("Selecting the features.")
  
  for (FEA in featype){
    
    print(paste("Feature extraction: ",mapping[which(FEA==mapping[,1]),2],", starting.",sep=""))
    
    contAlpha <-contAlpha+1
    
    if (FEA=="f1"){
      contWin<-contWin+1
    
    r<-.WinData(data, reps, classes, win[contWin], stat[contWin], power[contWin], abs[contWin],
     log[contWin], L, nel, nrep,mintomax[contWin])
    usedata <- r$windata
    usereps <- r$winreps
    useclasses <- r$winclasses
    
    r<-.WinData(testdata, testreps, testclasses, win[contWin], stat[contWin], power[contWin], 
    abs[contWin], log[contWin], L, nel, nreptest,mintomax[contWin])
    testusedata <- r$windata
    testusereps <- r$winreps
    testuseclasses <- r$winclasses
    } #FEA=="f1"
    
    if (FEA=="f2"){
      
      usedata <- abs(loaddata)
      usereps <- loadreps
      useclasses <- loadclasses
      
      testusedata <- abs(testloaddata)
      testusereps <- testloadreps
      testuseclasses <- testloadclasses
    }
    if (FEA=="f3"){
      
      contWin<-contWin+1
      
      r<-.WinData(specdata, specreps, specclasses, win[contWin], stat[contWin], power[contWin], abs[contWin],
       log[contWin], specL, nel, nrep,mintomax[contWin])
      usedata <- r$windata
      usereps <- r$winreps
      useclasses <- r$winclasses
      
      r<-.WinData(testspecdata, testspecreps, testspecclasses, win[contWin], stat[contWin], power[contWin], 
      abs[contWin], log[contWin], specL, nel, nreptest,mintomax[contWin])
      testusedata <- r$windata
      testusereps <- r$winreps
      testuseclasses <- r$winclasses
    
    }
    if (FEA=="f4"){
    
      usedata <- abs(specloaddata)
      usereps <- specloadreps
      useclasses <- specloadclasses
      
      testusedata <- abs(testspecloaddata)
      testusereps <- testspecloadreps
      testuseclasses <- testspecloadclasses
    
    }
    if (FEA=="f5"){
    
      contWin<-contWin+1
    
      r<-.WinData(pcadata, pcareps, pcaclasses, win[contWin], stat[contWin], power[contWin], abs[contWin],
       log[contWin], L, ncomps, nrep,mintomax[contWin])
      usedata <- r$windata
      usereps <- r$winreps
      useclasses <- r$winclasses
      
      r<-.WinData(testpcadata, testpcareps, testpcaclasses, win[contWin], stat[contWin], power[contWin], 
      abs[contWin], log[contWin], L, ncomps, nreptest,mintomax[contWin])
      testusedata <- r$windata
      testusereps <- r$winreps
      testuseclasses <- r$winclasses
    
    }
    if (FEA=="f6"){
    
      contWin<-contWin+1
      
      r<-.WinData(data, reps, classes, win[contWin], stat[contWin], power[contWin], abs[contWin],
       log[contWin], L, nel, nrep,mintomax[contWin])
      usedata <- r$windata
      usereps <- r$winreps
      useclasses <- r$winclasses
      L2 <- r$L2
      
      r<-.WinData(testdata, testreps, testclasses, win[contWin], stat[contWin], power[contWin], 
      abs[contWin], log[contWin], L, nel, nreptest,mintomax[contWin])
      testusedata <- r$windata
      testusereps <- r$winreps
      testuseclasses <- r$winclasses
      
      
      contWin<-contWin+1
      
      r<-.WinData(usedata, usereps, useclasses, win[contWin], stat[contWin], power[contWin], abs[contWin],
       log[contWin], L2, nel, nrep,mintomax[contWin])
      usedata <- r$windata
      usereps <- r$winreps
      useclasses <- r$winclasses
      
      r<-.WinData(testusedata, testusereps, testuseclasses, win[contWin], stat[contWin], power[contWin], 
      abs[contWin], log[contWin], L2, nel, nreptest,mintomax[contWin])
      testusedata <- r$windata
      testusereps <- r$winreps
      testuseclasses <- r$winclasses
   
    }
    
    if (FEA=="f7"){
      usedata <- cwtdata
      usereps <- cwtreps
      useclasses <- cwtclasses
      
      testusedata <- testcwtdata
      testusereps <- testcwtreps
      testuseclasses <- testcwtclasses 
    }



    
    P<-.organize(usedata, usereps, useclasses, nrep, NR)
    Ptest<-.organize(testusedata, testusereps, testuseclasses, nreptest, NRtest)
    #print("enter FeaSelec")
    W<-.feaSelect(P,Ptest,nrep,Alpha[contAlpha], AlphaCorr,nreptest,minacc[contAlpha], fast = fast)
    W<-W$Selected

    
    if (!is.null(W)){
      typeDiscart<-c(typeDiscart,contAlpha )
      if (FEA=="f6")  ParsDiscart<-c(ParsDiscart,contWin-1 )
      if (FEA%in%c("f1","f3","f5","f6")) ParsDiscart<-c(ParsDiscart,contWin )
      Wout <- append(Wout,list(W))
      if (length(W)==1) FinalFea <- rbind(FinalFea,append(P[W,],Ptest[W,])) else {
        FinalFea <- rbind(FinalFea,cbind(P[W,],Ptest[W,]))
      }
    } #is.null(W)
    
    
    print(paste("Feature extraction: ",mapping[which(FEA==mapping[,1]),2],", completed.",sep=""))
    cat("\n")



  }##for FEA in featype
  
  label<-factor(LETTERS[c(rep(1,nrep[1]),rep(2,nrep[2]),rep(1,nreptest[1]),rep(2,nreptest[2]))])
  
  if(is.null(FinalFea)) stop("None were selected.")
  P<-cbind(FinalFea[,which(label=='A')],FinalFea[,which(label=='B')])

  nA<-length(which(label=='A'))
  nB<-length(which(label=='B'))
  D.A <- rowSums(P[,1:nA])/nA
  D.B <- rowSums(P[,(nA+1):(nA+nB)])/nB
  Var.A.cwt <- rowSums(P[,1:nA]^2)
  Var.B.cwt <- rowSums(P[,(nA+1):(nA+nB)]^2)
  Var.A.cwt <- Var.A.cwt-nA*D.A^2
  Var.B.cwt <- Var.B.cwt-nB*D.B^2
  sp<-sqrt((Var.A.cwt+Var.B.cwt)/(nA+nB-2))
  t <- sqrt((nA)*(nB)/(nA+nB))*abs(D.A-D.B)/sp    #Teste T.
  
  featype=featype[typeDiscart]
  est=stat[ParsDiscart]
  pow=power[ParsDiscart]
  mm=mintomax[ParsDiscart]
  Wcomp=Wout

  
  T<-c()
  linha=1
  for(i in 1:length(featype)){
    ttam<-length(Wcomp[[i]])
    T<-c(T,max(t[linha:(linha+ttam-1)]))
    linha=linha+ttam
  }
  T<-T+abs(rnorm(length(T))/1000000000)
  
  winners<-c()
  cont=1
  contwin=1
  contfix=1
  contmax=length(featype)


  
  repeat{
    feafix=featype[cont]
    estfix=est[contwin]
    powfix=pow[contwin]
    mmfix=mm[contwin]
    tfix=T[cont]
    
    repeat{
      cont<-cont+1
      if (!(feafix%in%(c("f2","f4")))) contwin=contwin+1
      if (cont>contmax)break;
      feaat=featype[cont]
      estat=est[contwin]
      powat=pow[contwin]
      mmat=mm[contwin]
      tat=T[cont]
      	if(feafix==feaat){
      		if(estfix==estat){
      			if(powfix==powat){
      				if(mmfix==mmat){
      				  if (tat>tfix) tfix=tat
      				} else break;
      			} else break;
      		} else break;
      	}else break;
    }
    winners<-c(winners,which(T==tfix))
    if (cont>contmax)break;
    contfix=cont
    }
    winners<-unique(winners)


    LT<-length(winners)
    if(LT>Nfea){
      winners<-winners[which(T[winners]>=T[winners][order(T[winners],decreasing=TRUE)][Nfea])]
    }


    
    winlinhas<-c()
    cont<-1
    for(i in 1:length(featype)){
      ttam<-length(Wcomp[[i]])
      if(i%in%winners){
        winlinhas<-c(winlinhas,c(cont:(cont+ttam-1)))
      }
      cont<-cont+ttam
    }




    FinalFea=FinalFea[winlinhas,]
    Waux<-c()
    for (i in winners){
    Waux <- append(Waux,list(Wout[[i]]))
  }
  Wout <- Waux

  
  
  winwin<-c()
  cont<-1
  for(i in 1:length(featype)){
    if (featype[i]%in%c("f1","f3","f5")) {
    	if (i%in%winners) {winwin<-c(winwin,cont)}
    	cont<-cont+1
     }
     if (featype[i]%in%c("f6")) {
    	if (i%in%winners) {winwin<-c(winwin,c(cont,cont+1))}
    	cont<-cont+2
     }   
  }

  
  featype=featype[winners]
  stat=stat[ParsDiscart][winwin]
  power=power[ParsDiscart][winwin]
  mintomax=mintomax[ParsDiscart][winwin]
  log=log[ParsDiscart][winwin]
  win=win[ParsDiscart][winwin]
  abs=abs[ParsDiscart][winwin]
  
  
  
  result<-list(FinalFea=FinalFea,ncomps=ncomps, W=Wout, featype=featype,
  win=win,stat=stat,power=power,
  abs=abs,log=log,mintomax=mintomax,label=label,which.classes=which.classes,
  L=L,nch=nel,wavelet=wavelet,variance=variance)
  
  class(result)<-"featureEEG"
  
  return(result)

}

print.featureEEG <-
  function(x,...){
    cat("Number of features selected: ",nrow(x$FinalFea),"\n")
    cat("Use svmEEG to train a classifier and classifyEEG to classify some new data.\n")
  }


summary.featureEEG <-
  function(object,...){
    x<-object
    cat("Number of features selected: ",nrow(x$FinalFea),"\n")
    cat("Use svmEEG to train a classifier and classifyEEG to classify some new data.\n")
  }
