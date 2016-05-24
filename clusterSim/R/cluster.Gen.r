cluster.Gen<-function(numObjects=50,means=NULL,cov=NULL, fixedCov=TRUE, model=1,dataType="m", numCategories=NULL, numNoisyVar=0, numOutliers=0, rangeOutliers=c(1,10) ,inputType="csv2", inputHeader=TRUE, inputRowNames=TRUE, outputCsv="", outputCsv2="", outputColNames=TRUE, outputRowNames=TRUE)
{
	#require(MASS)
	gM<-NULL
	gM1<-NULL
	objects<-numObjects
	covariances<-cov
	fixedCovariance<-fixedCov
	scale<-dataType
	numberOfCategories<-numCategories
	noisyVariables<-numNoisyVar
	outlayers<-numOutliers       
	if(inputRowNames==TRUE)
	{
		inputRowNames<-1
	}
	if(inputRowNames==FALSE)
	{
		inputRowNames<-NULL
	}
	if(inputType!="csv" && inputType!="csv2")
	{
		stop("inputType should be csv or csv2")
	}
	if(dataType=="o" && numOutliers>0)
	{
		stop("numOutliers parameter is for metric ans interval data only")
	}
	if(numOutliers<0)
	{
		stop("numOutliers should be positive integer")
	}
	if(model <=12)
	{
		fixedCovariance <= TRUE
	}
	if(model >=13 && model <=20)
	{
		fixedCovariance = FALSE
	}
	if (model==2 && is.null(dim(means)))
	{    
    means<-as.matrix(means)
	}
	if (model==2 && is.null(dim(covariances)))
	{    
    covariances<-as.matrix(covariances)
	}
	if(model!=1)
	{
		if (model>20)
		{
			if (inputType=="csv")
			{
				try(means<-as.matrix(read.csv(paste("means_",model,".csv",sep=""),header=inputHeader,row.names=inputRowNames)),TRUE)
			}
			else
			{
				try(means<-as.matrix(read.csv2(paste("means_",model,".csv",sep=""),header=inputHeader,row.names=inputRowNames)),TRUE)
			}
			#print(.Last.value)
			if (class(.Last.value)=="try-error")
			{
        model<-3
			}
			else
			{
				if(fixedCovariance && model>20)
				{
					if (inputType=="csv")
					{
						covariances<-as.matrix(read.csv(paste("cov_",model,".csv",sep=""),header=inputHeader,row.names=inputRowNames))      
					}
					else
					{
            #print("czytam")
						covariances<-as.matrix(read.csv2(paste("cov_",model,".csv",sep=""),header=inputHeader,row.names=inputRowNames))      # SPRAWDZIC !!!
						#print(covariances)
						#print(ncol(covariances))
					}
				}
			}
		}
	}
	if(model!=1)
	{
		#print(objects)
		#print(means)
		#print(covariances)
		#start
		############# MODEL 3 #########
		if(model == 3 )
		{
		means<-matrix(c(0,1,0,5),2,2)
		covariances<-matrix(c(1,-0.9,-0.9,1),2,2)
		}
		############# MODEL 4 #########
		if(model == 4 )
		{
		means<-matrix(c(0,1.5,3,0,7,14),3,2)
		covariances<-matrix(c(1,-0.9,-0.9,1),2,2)
		}
		############# MODEL 5 #########
		if(model == 5 )
		{
		means<-matrix(c(1.5,3,4.5,6,12,18,-3,-6,-9),3,3)
		covariances<-matrix(c(1,-0.9,-0.9,-0.9,1,0.9,-0.9,0.9,1),3,3)
		}
		############# MODEL 6 #########
		if(model == 6 )
		{
		means<-matrix(c(5,-3,3,0,-5,5,3,-3,0,-5),5,2)
		covariances<-matrix(c(1,0.9,0.9,1),2,2)
		}
		############# MODEL 7 #########
		if(model == 7 )
		{
		means<-matrix(c(5,-3,3,0,-5,5,3,-3,0,-5,5,-3,3,0,-5),5,3)
		covariances<-matrix(c(1,0.9,0.9,0.9,1,0.9,0.9,0.9,1),3,3)
		}
		############# MODEL 8 #########
		if(model == 8 )
		{
		means<-matrix(c(0,0,5,10,10,0,10,5,0,10),5,2)
		covariances<-matrix(c(1,0,0,1),2,2)
		}
		############# MODEL 9 #########
		if(model == 9 )
		{
		means<-matrix(c(0,10,-10,10,-10,0,10,-10,-10,10,0,10,-10,10,10),5,3)
		covariances<-matrix(c(3,2,2,2,3,2,2,2,3),3,3)
		}
		############# MODEL 10 #########
		if(model == 10 )
		{
		means<-matrix(c(-4,5,14,5,5,14,5,-4),4,2)
		covariances<-matrix(c(1,0,0,1),2,2)
		}
		############# MODEL 11 #########
		if(model == 11 )
		{
		means<-matrix(c(-4,5,14,5,5,14,5,-4,-4,5,14,5),4,3)
		covariances<-matrix(c(1,0,0,0,1,0,0,0,1),3,3)
		}
		############# MODEL 12 #########
		if(model == 12 )
		{
		means<-matrix(c(-2,4,10,16),4,1)
		covariances<-matrix(c(0.5,0.5,0.5,0.5),4,1)
		}
		############# MODEL 13 #########
		if(model == 13 )
		{
		means<-matrix(c(0,1.5,3,0,7,14),3,2)
		covariances1<-matrix(c(1,-0.9,-0.9,1),2,2)
		covariances2<-matrix(c(1.5,0,0,1.5),2,2)
		covariances3<-matrix(c(1,0.5,0.5,1),2,2)
		}
		############# MODEL 14 #########
		if(model == 14 )
		{
		means<-matrix(c(-4,5,14,5,5,14,5,-4,-4,5,14,5),4,3)
		covariances1<-matrix(c(1,0,0,0,1,0,0,0,1),3,3)
		covariances2<-matrix(c(1,-0.9,-0.9,-0.9,1,0.9,-0.9,0.9,1),3,3)
		covariances3<-matrix(c(1,0.9,0.9,0.9,1,0.9,0.9,0.9,1),3,3)
		covariances4<-matrix(c(3,2,2,2,3,2,2,2,3),3,3)
		}
		############# MODEL 15 #########
		if(model == 15 )
		{
		means<-matrix(c(5,-3,3,0,-5,5,3,-3,0,-5,5,-3,3,0,-5),5,3)
		covariances1<-matrix(c(1,-0.9,-0.9,-0.9,1,0.9,-0.9,0.9,1),3,3)
		covariances2<-matrix(c(0.5,0,0,0,1,0,0,0,2),3,3)
		covariances3<-matrix(c(1,0.9,0.9,0.9,1,0.9,0.9,0.9,1),3,3)
		covariances4<-matrix(c(1,0.6,0.6,0.6,1,0.6,0.6,0.6,1),3,3)
		covariances5<-matrix(c(1,0,0,0,1,0,0,0,1),3,3)
		}
		############# MODEL 16 #########
		if(model == 16 )
		{
		means<-matrix(c(0,1,0,5),2,2)
		covariances1<-matrix(c(1,-0.9,-0.9,1),2,2)
		covariances2<-matrix(c(1,0.5,0.5,1),2,2)
		}
		if(length(objects)==1)
		{
			if(model==1)
			{
				objects<-rep(1,objects)
			}
			else
			{
				objects<-rep(objects,nrow(means))
			}
		}
		klasy<-NULL
		for (i in 1: nrow(means))
		{
			if(!fixedCovariance  && model >20)
			{
				if (inputType=="csv")
				{
					covariances<-as.matrix(read.csv(paste("cov_",model,"_",i,".csv",sep=""),header=inputHeader,row.names=inputRowNames))
				}
				else
				{
					covariances<-as.matrix(read.csv2(paste("cov_",model,"_",i,".csv",sep=""),header=inputHeader,row.names=inputRowNames))
				}
			}
			if(!fixedCovariance  && model <=20)
			{
				if(i==1) covariances<-covariances1
				if(i==2) covariances<-covariances2
				if(i==3) covariances<-covariances3
				if(i==4) covariances<-covariances4
				if(i==5) covariances<-covariances5
			}
			if(ncol(covariances)==1)
			{
				gM<-c(gM,rnorm(objects[i],means[i],covariances[i]))
				if (scale=="s")
				{
					gM1<-c(gM1,rnorm(objects[i],means[i],covariances[i]))
				}
			}
			else
			{	
				gM<-rbind(gM,mvrnorm(objects[i],means[i,],covariances))
				if (scale=="s")
				{
					gM1<-rbind(gM1,mvrnorm(objects[i],means[i,],covariances))

				}
			}
			klasy<-c(klasy,rep(i,objects[i]))
		}
    #print("debug 1")
		if(outlayers!=0)
		{
		  if (outlayers>0 && outlayers<1)
		  {
			outlayers=round(sum(objects)*outlayers)
		  }
			outlayersCounter<-0
	
			while (outlayersCounter<outlayers)
			{
			  outlayer<-rep(0,ncol(means))
			  for (i in 1:ncol(means))
			  {
				  if (sample(1:2)[1]==1)
				  {
					outlayer[i]<-runif(1,min=max(gM[1:sum(objects),i])+rangeOutliers[1],max=max(gM[1:sum(objects),i])+rangeOutliers[2])
				  }
				  else
				  {
					outlayer[i]<-runif(1,max=min(gM[1:sum(objects),i])-rangeOutliers[1],min=min(gM[1:sum(objects),i])-rangeOutliers[2])
				  }
				
			  }
			  klasy<-c(klasy,0)
						#print(gM[1,])
						#print(outlayer)
			  gM<-rbind(gM,outlayer);
			  outlayersCounter<-outlayersCounter+1
			}
		}	
	}
	else
	{
		gM<-NULL
		for(i in 1:numNoisyVar)
		{
		 gM<-cbind(gM,runif(objects,min=0,max=1))
		}
		klasy<-1:dim(gM)[1]
	}

	if(model!=1)
	if(ncol(covariances)==1)
	{
		gM<-as.matrix(gM)
		if (scale=="s")
		{
			gM1<-as.matrix(gM1)
		}
	}
	if (model==0 && (is.null(means) || is.null(covariances)))
		stop("Musi byc podany model lub macierze srednich i kowariancji")
	#print(gM)

	if (noisyVariables>0)
	{
		for (i in 1:noisyVariables)		
		{
			if (model!=1)
			gM<-cbind(gM,runif(nrow(gM),min=min(gM),max=max(gM)))
			else
			gM<-cbind(gM,runif(nrow(gM),min=0,max=1))
			if (scale=="s")
			{
				if (model!=1)
				gM1<-cbind(gM1,runif(ncol(gM1),min=min(gM1),max=max(gM1)))
				else
				gM1<-cbind(gM1,runif(ncol(gM1),min=0,max=1))
			}
		}

	}
    #print(gM)
	if(length(numberOfCategories)==1)
	{
		if (model==1)
		{
		numberOfCategories<-rep(numberOfCategories,noisyVariables)
		}
		else
		{
		numberOfCategories<-rep(numberOfCategories,ncol(means)+noisyVariables)
		}
	}
	if (scale=="m")
	{
		resul<-data.frame(gM)
	}
	if (scale=="o")
	{
		#print("debug 10")
		wynik<-gM
		wynik[,]<-0
		numerobiektu<-0
		for (j in 1:dim(gM)[2])
		{
			kolumna<-gM[,j]
			#print(kolumna)
			for (i in 1:dim(gM)[1])
			{
				if (gM[i,j]==max(kolumna))
				{

				wynik[i,j]=numberOfCategories[j]
				}
				else
				{	
	wynik[i,j]<-as.integer((gM[i,j]-min(kolumna))/(max(kolumna)-min(kolumna))*numberOfCategories[j])+1
				}
			}		
		}
		resul<-data.frame(wynik)

	}
	
	if (scale=="s")
	{
    if (nrow(gM)-sum(objects)>0)
    {
      delta<-runif((nrow(gM)-sum(objects))*ncol(gM),0,1)
      dim(delta)<-c(nrow(gM)-sum(objects),ncol(gM))
      gM1<-rbind(gM1,gM[(sum(objects)+1):nrow(gM),]+delta)
		}
		#print(gM1)
		
		wynik<-array(0,c(dim(gM)[1],dim(gM)[2],2))
		dim(wynik)<-c(dim(gM)[1],dim(gM)[2],2)
		for (i in 1:dim(gM)[1])
		for (j in 1:dim(gM)[2])
		{
			if (gM[i,j]<gM1[i,j])
			{
				wynik[i,j,1]<-gM[i,j]
				wynik[i,j,2]<-gM1[i,j]
			}
			else
			{
				wynik[i,j,1]<-gM1[i,j]
				wynik[i,j,2]<-gM[i,j]
			}
		}
	    resul<-wynik
	}
	if(model!=1)
	{
	koncoweNumeryKlas<-sample(1:nrow(means))
	}
	else
	{
	koncoweNumeryKlas<-sample(1:objects)
	}
	#print(koncoweNumeryKlas)
	for(i in 1:(nrow(gM)-outlayers))
	{
		klasy[i]<-koncoweNumeryKlas[klasy[i]]
	}
	#print(klasy)
	if(paste(outputCsv,"",sep="")!="")
	{
		write.table(cbind(1:dim(gM)[1],klasy,resul),file=outputCsv,row.names=outputRowNames,col.names=outputColNames)
	}
	if(paste(outputCsv2,"",sep="")!="")
	{
		write.table(cbind(1:dim(gM)[1],klasy,resul),file=outputCsv2,sep=";",dec=",",row.names=outputRowNames,col.names=outputColNames)
	}
	list(clusters=klasy,data=resul)
}
