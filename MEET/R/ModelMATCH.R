ModelMATCH<-function(iicc){
    Prob<-iicc$background
    TFBS<-iicc$Transcriptionfactor
    suma<-apply(TFBS,2,function(y){sum(y=="-")})
    TFBS<-TFBS[, suma==0]
    ncolTFBS<-ncol(TFBS)
    logodds <- CalculInformation(TFBS, Prob=Prob)
    vector_informacio <- as.vector(logodds[,1])
    inf <- vector(mode='logical', length=(length(vector_informacio)-4))
      for (j in 1:(length(vector_informacio)-4)){
	    inf[j] <- vector_informacio[j]+vector_informacio[j+1]+vector_informacio[j+2]+vector_informacio[j+3]+vector_informacio[j+4]
	  }
      index <- which.max(inf)
      logodds <- logodds[,2:dim(logodds)[2]]
      print(logodds)
      print(ncol(logodds))
       core <- logodds[(index:(index+4)),]
      print(dim(core))
      minim<-0
      maxim<-0
      minim_core<-0
      maxim_core<-0
      for(j in 1:dim(logodds)[1]){
	  minim <- min(logodds[j,])+minim
	  maxim <- max(logodds[j,])+maxim
	  }
      for (j in 1:dim(core)[1]){
	    minim_core <- minim_core+ min(core[j,])
	    maxim_core <- maxim_core+ max(core[j,])
      }

    parametersModel<-list(posCore=index, minim_core=minim_core, maxim_core=maxim_core, minim=minim, maxim=maxim, Corecut=iicc$parametersIdeal, logodds=logodds, core=core, ncolTFBS=ncolTFBS)

    model=list(model=logodds, parametersModel=parametersModel)

      
}
