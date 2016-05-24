plotEEG <-
function (data, classes.Id, rec.Id, which.classes = "ALL", which.rec="ALL", 
which.channels="ALL", type = 'original', m.a = 1,n.colors=200,
 wavelet="gaussian2", abs=FALSE,Real=TRUE,variance=1) 
{

which.elec = which.channels
  #par(mfrow=c(1,1))
  
  if(length(m.a)==1)
  {
		if(type=="original") ma=m.a
		if(type=="spectrum")
		{
			ma=1
			ma.s=m.a
		}
		if(type=="wavelet" || type=="T.pvalue") 
		{
			mmr=m.a
			mmc=m.a
		}
  }else if(length(m.a)==2)
  {
		if(type=="original") ma=m.a[1]
		if(type=="spectrum")
		{
			ma=m.a[1]
			ma.s=m.a[2]
		}
		if(type=="wavelet" || type=="T.pvalue") 
		{
			mmr=m.a[1]
			mmc=m.a[2]
		}		
  }else
  {
		stop("Parameter m.a: this parameters must be a number or a vector of size 2")
  }
  
  
  if (nrow(data)!=length(rec.Id) ||  length(rec.Id)!=length(classes.Id)) stop("The number
	of rows in 'data' must be equal to the length of the vectors rec.Id and classes.Id.")
  
  if (!(type%in%c("original","spectrum","wavelet","T.pvalue"))) stop("Parameter type: \n'type' must be one of 'original', 'spectrum', 'wavelet' or 'T.pvalue.")
  
	#Caso "ALL" redefinindo os vetores
	if (which.elec[1]=="ALL") which.elec <- 1:ncol(data)
	if (which.classes[1]=="ALL") which.classes <- unique(classes.Id)

	#Numerando as classes.
	cl<-numeric(length(classes.Id))
	lwc <- length(which.classes)
	for (i in 1:lwc)
  {
	  cl[which(classes.Id==which.classes[i])]=i
	}
	which.classes<-1:lwc
	classes.Id<-cl

	if (which.rec[1] == "ALL") which.rec<-lapply(1:lwc,function (g) unique(rec.Id[which(classes.Id==which.classes[g])]))

  if (is.null(ncol(data)))
  {
    L<-length(data[which(rec.Id%in%which.rec[[1]][1] & classes.Id == which.classes[1])])
  }else
  {	
    L<-length(data[which(rec.Id%in%which.rec[[1]][1] & classes.Id == which.classes[1]),1])
  }
  rest<-lapply(1:lwc, function (g) which(classes.Id == g))

  
  #Definindo o banco de dados que deve ser utilizado.
  #Definindo o banco de dados que deve ser utilizado.
  v<-c()
  for (g in 1:lwc)
  {
    v<-c(v,which(rec.Id%in%which.rec[[g]] & classes.Id == which.classes[g]))
  }
  
  if (is.null(ncol(data))){
  dados<-data[v] } else{
  dados<-data[v,which.elec]
  }
  rec.Id<-rec.Id[v]
  classes.Id<-classes.Id[v]
  dados<-cbind(dados,rec.Id)
  
  	
	if (type == "original") 
  {

  	A<-c()
  	a=rep(0,L)
  	for (i in 1:(L-ma+1)){A<-rbind(A,a) }
  	for (i in 1:(L-ma+1))
    {
  	  for (j in 0:(ma-1)){
  	  A[i,i+j]=1
  	  }
  	}

  	P<-c()  #calculando o MA
  	for (el in 1:length(which.elec))
    {
  	  for (class in which.classes) 
      {
  	    for (R in which.rec[[class]])
        {
  	       P<-rbind(P,t(A%*%dados[which(classes.Id==class & rec.Id==R),el]))
  	    }
  	  }
  	}
  	max<-max(P)  #max min para grafico
  	min<-min(P)

  
  	if (lwc > 2) 
    {  #se mais de uma classe colorir por classe
    	plot(1:(L-ma+1),1:(L-ma+1),type='l',col="white",ylim=c(min,max),xlab="",ylab="")
    	cont<-0
    	for (el in 1:length(which.elec)){
    		for (class in which.classes) {
    		  for (R in which.rec[[class]]){
    		   cont <- cont+1
    		   lines(1:(L-ma+1),P[cont,],col=class)
    		  }
    		}
    	}
  	} else 
    { 
    	if (lwc == 1) 
      { # se tiver apenas uma classes colorir por rep.
    	  plot(1:(L-ma+1),1:(L-ma+1),type='l',col="white",ylim=c(min,max),xlab="",ylab="")
    	  cont<-0
    		for (el in 1:length(which.elec)){
    	  	for (class in which.classes) {
    		    for (R in which.rec[[class]]){
    		      cont <- cont+1
    		      lines(1:(L-ma+1),P[cont,],col=R)
    	  	  }
    	  	}
    		}	
    	} else 
      {
    		COL<-list(c(100 ,372, 375, 503 ,506 ,525, 554 ,642, 90 ,93, 501),c(26, 30 ,31,  42,  43, 62, 121, 124, 125, 128, 403))
    		plot(1:(L-ma+1),1:(L-ma+1),type='l',col="white",ylim=c(min,max),xlab="",ylab="")
    		cont<-0
    		cont.R<-1
    		cont.B<-1
    		for (el in 1:length(which.elec)){
      		for (class in which.classes) {
        		for (R in which.rec[[class]]){
          		cont <- cont+1
          		if (class==1){ccol<-cont.R ; cont.R<-cont.R+1 ; if(cont.R>11) cont.R<-1}
          		if (class==2){ccol<-cont.B ; cont.B<-cont.B+1 ; if(cont.B>11) cont.B<-1}
          		lines(1:(L-ma+1),P[cont,],col=colors()[COL[[class]][ccol]])
    	    	}
      		}
    		}
    	}
  	}
	}

	if (type == "spectrum"){
		
		A<-c()   #ma para dados
		a=rep(0,L)
		for (i in 1:(L-ma+1)){A<-rbind(A,a) }
  	for (i in 1:(L-ma+1)){
    	for (j in 0:(ma-1)){
    		A[i,i+j]=1
  	  }
		}

		L2<-length(spectrum(A%*%dados[which(classes.Id==which.classes[1] & rec.Id==which.rec[[1]][1]),1],plot=F)$spec)


		A2<-c() #ma para spec
		a=rep(0,L2)
		for (i in 1:(L2-ma.s+1)){A2<-rbind(A2,a) }
		for (i in 1:(L2-ma.s+1)){
	  	for (j in 0:(ma.s-1)){
	    	A2[i,i+j]=1
		  }
		}


		P<-c()  #calcula ma spec e ma novamente
		for (el in 1:length(which.elec)){
  		for (class in which.classes) {
    		for (R in which.rec[[class]]){
    		  P<-rbind(P,t(A2%*%spectrum(A%*%dados[which(classes.Id==class & rec.Id==R),el],plot=F)$spec))
    	  }
  		}
		}
		max<-max(P)  #max e min para os graficos
		min<-min(P)
		L <- ncol(P)

		if (lwc > 2) {  #se numero de classes maior que 1 colorir por classe
				cont <- 0
				plot(1:L,1:L,type='l',col="white",ylim=c(min,max),xlab="Frequency",ylab="Spectrum")
				for (el in 1:length(which.elec)){
  				for (class in which.classes) {
    				for (R in which.rec[[class]]){
    				  cont <- cont+1
    				  lines(1:L,P[cont,],col=class)
    				}
  				}
				}
			} else {
  			if (lwc == 1) { # se tiver apenas uma classes colorir por rep. #se 1 classe apenas entao colorir por rep
  				cont <- 0
  				plot(1:L,1:L,type='l',col="white",ylim=c(min,max),xlab="",ylab="")
  				for (el in 1:length(which.elec)){
    				for (class in which.classes) {
      				for (R in which.rec[[class]]){
        				cont <- cont+1
        				lines(1:L,P[cont,],col=R)
      				}
    				}
  				}	
			} else {
				COL<-list(c(100 ,372, 375, 503 ,506 ,525, 554 ,642, 90 ,93, 501),c(26, 30 ,31,  42,  43, 62, 121, 124, 125, 128, 403))
				cont <- 0
				cont.R<-1
				cont.B<-1
				plot(1:L,1:L,type='l',col="white",ylim=c(min,max),xlab="",ylab="")
				for (el in 1:length(which.elec)){
  				for (class in which.classes) {
    				for (R in which.rec[[class]]){
      				cont <- cont+1
      				if (class==1){ccol<-cont.R ; cont.R<-cont.R+1 ; if(cont.R>11) cont.R<-1}
      				if (class==2){ccol<-cont.B ; cont.B<-cont.B+1 ; if(cont.B>11) cont.B<-1}
      				lines(1:L,P[cont,],col=colors()[COL[[class]][ccol]])
    				}
  				}
				}	
			}
		}
		
	} #fim do if spectrum

	if (type == "wavelet") {
	
  	tot.rep<-0 #calculando o total de graficos que serao feitos.
  	for (class in which.classes) {tot.rep<-length(which.rec[[class]])+tot.rep}
  	tot.graphs<-length(which.elec)*tot.rep
    
  	if (tot.graphs>5) 	{par(mfrow=c(2,2)); par1=F} else {
  		#par(mfrow=c(1,1))
  		 par1=T}
  	
  	Mcwt<-as.matrix(wavCWT(dados[which(classes.Id==which.classes[1] & rec.Id==which.rec[[1]][1]),1],wavelet=wavelet,variance=variance))
  	ncolcwt<-ncol(Mcwt) #calculando os parametros de linhas colunas e matrizes MA
  	nrowcwt<-nrow(Mcwt)
  	A<-mat.or.vec(nrowcwt-mmr+1, nrowcwt)
  	for (i in 1:(nrowcwt-mmr+1)){
  	  for (j in 0:(mmr-1)){
  	  A[i,i+j]=1
  	  }
  	}
  	A2<-mat.or.vec(ncolcwt-mmc+1, ncolcwt)
  	for (i in 1:(ncolcwt-mmc+1)){
    	for (j in 0:(mmc-1)){
    	  A2[i,i+j]=1
    	}
  	}	
  	
  	P<-array(numeric(0),c(nrowcwt-mmr+1,ncolcwt-mmc+1,tot.graphs)) #define o array para
  	
  
  	# armazenar as matrizes suavizadas de wavelet
  	cont<-0
		for (el in 1:length(which.elec)){ #armazenando as matrizes
  		for (class in which.classes) {
    		for (R in which.rec[[class]]){
    		cont<-cont+1
    		a<-as.matrix(wavCWT(dados[which(classes.Id==class & rec.Id==R),el],wavelet=wavelet,variance=variance))
    		if (abs){a<-abs(a)}
    		b <- sapply(1:ncolcwt,function(g) A%*%a[,g])
    		if (is.null(nrow(b))){b<-rbind(b,b)}
    		
    		P[,,cont] <-  t(sapply(1:(nrowcwt-mmr+1),function(g) A2%*%b[g,]))
    		}
  		}
		}
		if (wavelet == "morlet" & Real) P<-Re(P)  #morlet tem parte imaginaria, utiliza-se
		if (wavelet == "morlet" & !Real) P<-Im(P) #apenas uma dessas
 		max<-max(P)
		min<-min(P)
		brks<-seq(min,max,len=n.colors+1)
		cont<-0
		if(is.null(nrow(P[,,1]))){
			#par(mfrow=c(1,1))
			L<-length(P[,,1])
			plot(1:L,1:L,type='l',col="white",ylim=c(min,max),xlab="",ylab="")
		}
  
    key=1
  
		for (el in 1:length(which.elec)){
			for (class in which.classes) {
			  for (R in which.rec[[class]]){
			  cont<-cont+1
			
			  if(length(P[,,cont])==1){
		    	points(1,P[,,cont],col=class,pch=19)
		  	}else{
			    if(is.null(nrow(P[,,cont]))){
			      lines(1:L,P[,,cont],col=class)
			    }else{        
			      image.plot(P[,,cont], breaks=brks, col=tim.colors(n.colors),main=paste("Class",class,", Channel",which.elec[el],", Rec",R),xlab="",ylab="")
			
			      if(el!=length(which.elec) | class!=tail(which.classes,1) | R!= tail(which.rec[[class]],1))
		      	{
			      	if (par1 | (cont%%4)==0) {
			        	cat("\n","Enter s to stop. Enter any other key to continue.","\n") # prompt
			        	key<-readline() 
			      	}
			      } else {key=="s"}
      
          }
		  	}

			if (key=="s") break}
      if (key=="s") break}
      if (key=="s") break}
		
	}
  
  
  
	if (type == 'T.pvalue') {
  	if (length(which.classes)!=2) stop("In type=T.value the user must
  	specify 2 classes.")
  
  	if (length(which.elec)>5) {par(mfrow=c(2,2));par1=F} else {par1=T}
    
  
  	Mcwt<-as.matrix(wavCWT(dados[which(classes.Id==which.classes[1] & rec.Id==which.rec[[1]][1]),1],wavelet=wavelet,variance=variance))
  	ncolcwt<-ncol(Mcwt) #calculando os parametros de linhas colunas e matrizes MA
  	nrowcwt<-nrow(Mcwt)
  	A<-mat.or.vec(nrowcwt-mmr+1, nrowcwt)
  	for (i in 1:(nrowcwt-mmr+1)){
  	  for (j in 0:(mmr-1)){
  	    A[i,i+j]=1
  	  }
  	}
  	A2<-mat.or.vec(ncolcwt-mmc+1, ncolcwt)
  	for (i in 1:(ncolcwt-mmc+1)){
    	for (j in 0:(mmc-1)){
    	  A2[i,i+j]=1
    	 }
  	}	

  	nA <- length(which.rec[[1]]) #numero de repeticoes grupo A
  	nB <- length(which.rec[[2]]) #numero de repeticoes grupo B
  	if (nA==1 || nB==1)stop("It should have at least 2 repetitions of each class for type=T.pvalue.")
  	EL <- length(which.elec)

		A.cwt <- array(numeric(0),c(nrowcwt-mmr+1,ncolcwt-mmc+1,nA,EL))
		B.cwt <- array(numeric(0),c(nrowcwt-mmr+1,ncolcwt-mmc+1,nB,EL))
	
		for (Rep in 1:nA){  #Faz a CWT para o grupo 1 e armazena no array A.cwt
  		for (el in 1:EL) {
    		a <- as.matrix(wavCWT(dados[which(rec.Id==unique(rec.Id)[Rep] & classes.Id==1),el],wavelet=wavelet,variance=variance))
    		if (abs){a<-abs(a)}		
    		b <- sapply(1:ncolcwt,function(g) A%*%a[,g])
    		A.cwt[,,Rep,el] <- t(sapply(1:(nrowcwt-mmr+1),function(g) A2%*%b[g,]))
  		}
		}

		for (Rep in 1:nB){  #Faz a CWT para o grupo 2 e armazena no array B.cwt
  		for (el in 1:EL) {
    		a <- as.matrix(wavCWT(dados[which(rec.Id==unique(rec.Id)[Rep] & classes.Id==2),el],wavelet=wavelet,variance=variance))
    		if (abs){a<-abs(a)}		
    		b <- sapply(1:ncolcwt,function(g) A%*%a[,g])
    		B.cwt[,,Rep,el] <- t(sapply(1:(nrowcwt-mmr+1),function(g) A2%*%b[g,]))
  		}
		}


		if (nA>1 & nB>1 & EL>1){
  		d.A <- apply(A.cwt[,,,], c(1,2,4), sum)/nA
  		d.B <- apply(B.cwt[,,,], c(1,2,4), sum)/nB
  		var.A.cwt <- apply((A.cwt[,,,])^2, c(1,2,4), sum)
  		var.A.cwt <- (var.A.cwt-nA*d.A^2)/(nA-1)
  		var.B.cwt <- apply((B.cwt[,,,])^2, c(1,2,4), sum)
  		var.B.cwt <- (var.B.cwt-nB*d.B^2)/(nB-1)
		}

		if (nA>1 & nB>1 & EL==1){
  		d.A <- apply(A.cwt[,,,], c(1,2), sum)/nA
  		d.B <- apply(B.cwt[,,,], c(1,2), sum)/nB
  		var.A.cwt <- apply((A.cwt[,,,])^2, c(1,2), sum)
  		var.A.cwt <- (var.A.cwt-nA*d.A^2)/(nA-1)
  		var.B.cwt <- apply((B.cwt[,,,])^2, c(1,2), sum)
  		var.B.cwt <- (var.B.cwt-nB*d.B^2)/(nB-1)
		}
  	d<-d.A-d.B
  	sp <- ((nA-1)*var.A.cwt+(nB-1)*var.B.cwt)/(nA+nB-2)	
  	t <- sqrt(nA*nB/(nA+nB))*d/sqrt(sp)
  	t<-abs(t)
  	pt<-pt(t,nA+nB-2)


		if (EL > 1) {
  		brks<-seq(0.5,1,len=n.colors+1)
  		if (n.colors==2) brks<-c(0.5,0.95,1)
      cont=0
      key=1
  		for (el in 1:EL){
        cont=cont+1
  	  	image.plot(pt[,,el], breaks=brks, col=tim.colors(n.colors),main=paste("Channel",which.elec[el]),xlab="",ylab="")
  		
  			if(el!=EL)
  			{
  				if (par1 | (cont%%4)==0) {
  			  	cat("\n","Enter s to stop. Enter any other key to continue.","\n") # prompt
  			  	key<-readline() 
  				}
  			} else {key=="s"}
  			

        if (key=="s") break  
      }
    } else{
  		brks<-seq(0.5,1,len=n.colors+1)
  		if (n.colors==2) brks<-c(0.5,0.95,1)
  		image.plot(pt, breaks=brks, col=tim.colors(n.colors),main=paste("Channel",which.elec[el]),xlab="",ylab="")
  	}	

		

	} #if T.pvalue

}
