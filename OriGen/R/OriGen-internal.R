#This is R code for OriGen made by John Michael O. Ranola ranolaj@uw.edu
#if the function is not to be accessed by users, start with a period(.)

.is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol



ConvertPEDData<-function(PlinkFileName,LocationFileName){
#DataFileName should be the base name of plink ped/map format
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#Location file should be space/tab delimited with columns Label,longitude,latitude
	print("Note: This method assumes each geographical location in the location file has a unique character label in the first column.")
	print('It also assumes the longitude and latitude columns are labeled "longitude" and "latitude".')
	MAPFileName=paste(PlinkFileName,".map",sep="")
	MAPData=read.table(MAPFileName,header=FALSE)
	NumberSNPs=length(MAPData[[1]])
	print(c("NumberSNPs",NumberSNPs))

	LocationData=read.table(LocationFileName,header=TRUE)
	SampleSites=nlevels(LocationData[[1]])
	print(c("SampleSites",SampleSites))

	PEDFileName=paste(PlinkFileName,".ped",sep="")
	PEDData=read.table(PEDFileName,header=FALSE)
	NumberIndividuals=length(PEDData[[1]])

	DataArray=array(0,c(2,SampleSites,NumberSNPs))
	SampleCoordinates=array(0,c(SampleSites,2))
	MembersList=levels(LocationData[[1]])

	#Here we fill in the DataArray
	#set up a vector dividing the sample sites
	SampleSitesLogical=array(FALSE,c(SampleSites,NumberIndividuals))
	for(i in 1:SampleSites){
		SampleSitesLogical[i,]=(as.numeric(LocationData[[1]])==i)
	}

	for(j in 1:NumberSNPs){
		k=2*j-1
		bothlevels=union(levels(PEDData[[6+k]]),levels(PEDData[[7+k]]))
		PEDData[[6+k]]=factor(PEDData[[6+k]],levels=bothlevels)
		PEDData[[7+k]]=factor(PEDData[[7+k]],levels=bothlevels)

		counter=1
		if(length(bothlevels)==3){
			if(bothlevels[1]=="0"){
				counter=2
			}else{
				stop(paste0("3 alleles found at locus ",j))
			}
		}

		tempLogical61=(as.numeric(PEDData[[6+k]])==counter)
		tempLogical62=(as.numeric(PEDData[[6+k]])==(counter+1))
		tempLogical71=(as.numeric(PEDData[[7+k]])==counter)
		tempLogical72=(as.numeric(PEDData[[7+k]])==(counter+1))
		for(i in 1:SampleSites){
			temp1=sum(tempLogical61*(SampleSitesLogical[i,]))
			temp2=sum(tempLogical71*(SampleSitesLogical[i,]))
			DataArray[1,i,j]=temp1+temp2
			temp1=sum(tempLogical62*(SampleSitesLogical[i,]))
			temp2=sum(tempLogical72*(SampleSitesLogical[i,]))
			DataArray[2,i,j]=temp1+temp2
		}
	}

	#here we fill in SampleCoordinates
	foundVector=array(1,SampleSites)

	#first one is free
	temp=as.numeric(LocationData[[1]])[1]
	SampleCoordinates[temp,1]=LocationData$longitude[1]
	SampleCoordinates[temp,2]=LocationData$latitude[1]
	foundVector[temp]=0

	for(i in 2:NumberIndividuals){
		if(sum(foundVector)>0){
			temp=as.numeric(LocationData[[1]])[i]
			if(foundVector[temp]){
				SampleCoordinates[temp,1]=LocationData$longitude[i]
				SampleCoordinates[temp,2]=LocationData$latitude[i]
				foundVector[temp]=0
			}
		}
	}

	ResultsRaw=list(DataArray=DataArray,SampleCoordinates=SampleCoordinates,Membership=as.numeric(LocationData[[1]]),MembersList=MembersList,SampleSites=SampleSites,NumberSNPs=NumberSNPs,NumberIndividuals=NumberIndividuals,PEDFileName=PEDFileName,MAPFileName=MAPFileName,LocationFileName=LocationFileName)

	return(ResultsRaw)
}


FitOriGenModel<-function(DataArray,SampleCoordinates,MaxGridLength=20,RhoParameter=10){
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
#This function takes in the data, fits the model, and returns the allele frequency surfaces
	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}
	NumberSNPs=length(DataArray[1,1,])
	SampleSites=length(DataArray[1,,1])
	GridLength=array(0,2)
	GridCoordinates=array(0.,dim=c(2,MaxGridLength))

	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")

	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=GridAndCoordResults$GridCoordinates

	AlleleFrequencySurfaces=array(0,dim=c(NumberSNPs,GridLength[1],GridLength[2]))
	ResultsRaw=.Fortran("FITORIGENMODEL",AlleleFrequencySurfaces=as.double(AlleleFrequencySurfaces),DataArray=as.integer(DataArray),NumberSNPs=as.integer(NumberSNPs),GridLength=as.integer(GridLength),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),MaxGridLength=as.integer(MaxGridLength),SampleCoordinates=as.double(SampleCoordinates),GridCoordinates=as.double(GridCoordinates),PACKAGE="OriGen")

	ResultsRaw$AlleleFrequencySurfaces=array(ResultsRaw$AlleleFrequencySurfaces,c(NumberSNPs,GridLength[1],GridLength[2]))
	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))

	return(ResultsRaw)
}





ConvertUnknownPEDData<-function(PlinkFileName,LocationFileName,PlinkUnknownFileName){
#DataFileName should be the base name of plink ped/map format
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#Location file should be space/tab delimited with columns Label,longitude,latitude

	print("Note: This method assumes each geographical location in the location file has a unique character label in the first column.")
	print('It also assumes the longitude and latitude columns are labeled "longitude" and "latitude".')
	MAPFileName=paste(PlinkFileName,".map",sep="")
	MAPData=read.table(MAPFileName,header=FALSE)
	NumberSNPs=length(MAPData[[1]])
	print(c("NumberSNPs",NumberSNPs))

	LocationData=read.table(LocationFileName,header=TRUE)
	SampleSites=nlevels(LocationData[[1]])
	print(c("SampleSites",SampleSites))

	PEDFileName=paste(PlinkFileName,".ped",sep="")
	PEDData=read.table(PEDFileName,header=FALSE)
	NumberIndividuals=length(PEDData[[1]])

	UnknownFileName=paste(PlinkUnknownFileName,".ped",sep="")
	UnknownRawData=read.table(UnknownFileName,header=FALSE)
	NumberUnknowns=length(UnknownRawData[[1]])

	DataArray=array(0,c(2,SampleSites,NumberSNPs))
	SampleCoordinates=array(0,c(SampleSites,2))
	MembersList=levels(LocationData[[1]])
	UnknownData=array(0,c(NumberUnknowns,NumberSNPs))

	#performs a check to see whether there is the same number of SNPs in PlinkFileName and PlinkUnknownFileName
	if(length(names(UnknownRawData))!=length(names(PEDData))){
		stop(paste0("Different number of SNPs in ",PlinkFileName, " and ", PlinkUnknownFileName))
	}

	#Here we fill in the DataArray and UnknownData
	#set up a vector dividing the sample sites
	SampleSitesLogical=array(FALSE,c(SampleSites,NumberIndividuals))
	for(i in 1:SampleSites){
		SampleSitesLogical[i,]=(as.numeric(LocationData[[1]])==i)
	}

	for(j in 1:NumberSNPs){
		k=2*j-1
		bothlevels=union(levels(PEDData[[6+k]]),levels(PEDData[[7+k]]))
		PEDData[[6+k]]=factor(PEDData[[6+k]],levels=bothlevels)
		PEDData[[7+k]]=factor(PEDData[[7+k]],levels=bothlevels)

		counter=1
		if(length(bothlevels)==3){
			if(bothlevels[1]=="0"){
				counter=2
			}else{
				stop(paste0("3 alleles found at locus ",j))
			}
		}

		tempLogical61=(as.numeric(PEDData[[6+k]])==counter)
		tempLogical62=(as.numeric(PEDData[[6+k]])==(counter+1))
		tempLogical71=(as.numeric(PEDData[[7+k]])==counter)
		tempLogical72=(as.numeric(PEDData[[7+k]])==(counter+1))
		for(i in 1:SampleSites){
			temp1=sum(tempLogical61*(SampleSitesLogical[i,]))
			temp2=sum(tempLogical71*(SampleSitesLogical[i,]))
			DataArray[1,i,j]=temp1+temp2
			temp1=sum(tempLogical62*(SampleSitesLogical[i,]))
			temp2=sum(tempLogical72*(SampleSitesLogical[i,]))
			DataArray[2,i,j]=temp1+temp2
		}

		#Filling in the UnknownData
		temp1=(UnknownRawData[[6+k]]==bothlevels[counter])
		temp2=(UnknownRawData[[7+k]]==bothlevels[counter])
		UnknownData[,j]=temp1+temp2
	}

	#here we fill in SampleCoordinates
	foundVector=array(1,SampleSites)

	#first one is free
	temp=as.numeric(LocationData[[1]])[1]
	SampleCoordinates[temp,1]=LocationData$longitude[1]
	SampleCoordinates[temp,2]=LocationData$latitude[1]
	foundVector[temp]=0

	for(i in 2:NumberIndividuals){
		if(sum(foundVector)>0){
			temp=as.numeric(LocationData[[1]])[i]
			if(foundVector[temp]){
				SampleCoordinates[temp,1]=LocationData$longitude[i]
				SampleCoordinates[temp,2]=LocationData$latitude[i]
				foundVector[temp]=0
			}
		}
	}

	ResultsRaw=list(DataArray=DataArray,UnknownData=UnknownData,SampleCoordinates=SampleCoordinates,Membership=as.numeric(LocationData[[1]]),MembersList=MembersList,SampleSites=SampleSites,NumberSNPs=NumberSNPs,NumberIndividuals=NumberIndividuals,NumberUnknowns=NumberUnknowns,PEDFileName=PEDFileName,MAPFileName=MAPFileName,LocationFileName=LocationFileName)

	return(ResultsRaw)
}



FitOriGenModelFindUnknowns<-function(DataArray,SampleCoordinates,UnknownData,MaxGridLength=20,RhoParameter=10){
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#by jmor
#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
#UnknownData[NumberUnknowns,NumberSNPs] gives the number of major alleles for the current unknown individual
#This function takes in the data, fits the model, and returns the allele frequency surfaces
	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}
	NumberSNPs=length(DataArray[1,1,])
	SampleSites=length(DataArray[1,,1])
	GridLength=array(0,2)
	GridCoordinates=array(0.,dim=c(2,MaxGridLength))

	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")

	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=GridAndCoordResults$GridCoordinates

	NumberUnknowns=length(UnknownData[,1])
	UnknownGrids=array(0,dim=c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw=.Fortran("FITORIGENMODELFINDUNKNOWNS",UnknownGrids=as.double(UnknownGrids),DataArray=as.integer(DataArray),NumberSNPs=as.integer(NumberSNPs),GridLength=as.integer(GridLength),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),MaxGridLength=as.integer(MaxGridLength),SampleCoordinates=as.double(SampleCoordinates),NumberUnknowns=as.integer(NumberUnknowns),UnknownData=as.integer(UnknownData),GridCoordinates=as.double(GridCoordinates),PACKAGE="OriGen")

	ResultsRaw$UnknownGrids=array(ResultsRaw$UnknownGrids,c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
	ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
	ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))

	ResultsRaw$UnknownData=array(ResultsRaw$UnknownData,c(NumberUnknowns,NumberSNPs))

	return(ResultsRaw)
}



FindRhoParameterCrossValidation<-function(PlinkFileName,LocationFileName,MaxIts=6,MaxGridLength=20){
#DataFileName should be the base name of plink ped/map format
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#Location file should be space/tab delimited with columns ID,Label,AltLabel,Long,Lat
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(MaxIts<3){
		stop("MaxIts must be greater than 3")
	}

	print("Note: This method assumes each geographical location in the location file has a unique character label in the first column.")
	print('It also assumes the longitude and latitude columns are labeled "longitude" and "latitude".')
	MAPFileName=paste(PlinkFileName,".map",sep="")
	MAPData=read.table(MAPFileName,header=FALSE)
	NumberSNPs=length(MAPData[[1]])
	print(c("NumberSNPs",NumberSNPs))

	LocationData=read.table(LocationFileName,header=TRUE)
	SampleSites=nlevels(LocationData[[1]])
	print(c("SampleSites",SampleSites))

	PEDFileName=paste(PlinkFileName,".ped",sep="")
	PEDData=read.table(PEDFileName,header=FALSE)
	NumberIndividuals=length(PEDData[[1]])

	DataArray=array(0,c(2,SampleSites,NumberSNPs))
	SampleCoordinates=array(0,c(SampleSites,2))
	MembersList=levels(LocationData[[1]])

	NumberUnknowns=NumberIndividuals
	UnknownData=array(0,c(NumberUnknowns,NumberSNPs))

	#Here we fill in the DataArray and UnknownData
	#set up a vector dividing the sample sites
	SampleSitesLogical=array(FALSE,c(SampleSites,NumberIndividuals))
	for(i in 1:SampleSites){
		SampleSitesLogical[i,]=(as.numeric(LocationData[[1]])==i)
	}

	for(j in 1:NumberSNPs){
		k=2*j-1
		bothlevels=union(levels(PEDData[[6+k]]),levels(PEDData[[7+k]]))
		PEDData[[6+k]]=factor(PEDData[[6+k]],levels=bothlevels)
		PEDData[[7+k]]=factor(PEDData[[7+k]],levels=bothlevels)

		counter=1
		if(length(bothlevels)==3){
			if(bothlevels[1]=="0"){
				counter=2
			}else{
				stop(paste0("3 alleles found at locus ",j))
			}
		}

		tempLogical61=(as.numeric(PEDData[[6+k]])==counter)
		tempLogical62=(as.numeric(PEDData[[6+k]])==(counter+1))
		tempLogical71=(as.numeric(PEDData[[7+k]])==counter)
		tempLogical72=(as.numeric(PEDData[[7+k]])==(counter+1))
		for(i in 1:SampleSites){
			temp1=sum(tempLogical61*(SampleSitesLogical[i,]))
			temp2=sum(tempLogical71*(SampleSitesLogical[i,]))
			DataArray[1,i,j]=temp1+temp2
			temp1=sum(tempLogical62*(SampleSitesLogical[i,]))
			temp2=sum(tempLogical72*(SampleSitesLogical[i,]))
			DataArray[2,i,j]=temp1+temp2
		}

		#Filling in the UnknownData
		temp1=(PEDData[[6+k]]==bothlevels[counter])
		temp2=(PEDData[[7+k]]==bothlevels[counter])
		UnknownData[,j]=temp1+temp2
	}

	#here we fill in SampleCoordinates
	foundVector=array(1,SampleSites)

	#first one is free
	temp=as.numeric(LocationData[[1]])[1]
	SampleCoordinates[temp,1]=LocationData$longitude[1]
	SampleCoordinates[temp,2]=LocationData$latitude[1]
	foundVector[temp]=0

	for(i in 2:NumberIndividuals){
		if(sum(foundVector)>0){
			temp=as.numeric(LocationData[[1]])[i]
			if(foundVector[temp]){
				SampleCoordinates[temp,1]=LocationData$longitude[i]
				SampleCoordinates[temp,2]=LocationData$latitude[i]
				foundVector[temp]=0
			}
		}
	}

	RhoVector=array(0,c(2,MaxIts))
	RhoParameter=1

#ResultsRaw=.Fortran("LEAVE_ONE_POP_OUT_CROSSVAL_SQUARE",PlinkFileName=as.character(PedFileName),LocationFileName=as.character(LocationFileName),NumberSNPs=as.integer(NumberSNPs),MaxIts=as.integer(MaxIts),MaxGridLength=as.integer(MaxGridLength),RhoVector=as.double(RhoVector),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")

ResultsRaw=.Fortran("LEAVE_ONE_POP_OUT_CROSSVAL_SQUARE2",DataArray=as.integer(DataArray),SampleCoordinates=as.double(SampleCoordinates),UnknownData=as.double(UnknownData),NumberSNPs=as.integer(NumberSNPs),SampleSites=as.integer(SampleSites),NumberUnknowns=as.integer(NumberUnknowns),MaxIts=as.integer(MaxIts),MaxGridLength=as.integer(MaxGridLength),RhoVector=as.double(RhoVector),RhoParameter=as.double(RhoParameter),PACKAGE="OriGen")

ResultsRaw$RhoVector=array(RhoVector,(c(2,MaxIts)))
ResultsRaw2=list(PlinkFileName=PlinkFileName,LocationFileName=LocationFileName,NumberSNPs=as.integer(NumberSNPs),SampleSites=as.integer(SampleSites),NumberUnknowns=as.integer(NumberUnknowns),MaxIts=as.integer(MaxIts),MaxGridLength=as.integer(MaxGridLength),RhoVector=ResultsRaw$RhoVector,RhoParameter=ResultsRaw$RhoParameter)

	return(ResultsRaw2)
}



FitAdmixedModelFindUnknowns<-function(DataArray,SampleCoordinates,UnknownData,MaxGridLength=20,RhoParameter=10,LambdaParameter=100.,MaskWater=TRUE){
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
#UnknownData[NumberUnknowns,NumberSNPs] gives the number of major alleles for the current unknown individual
#This function takes in the data, fits the model, and returns the allele frequency surfaces
	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}
	NumberSNPs=length(DataArray[1,1,])
	SampleSites=length(DataArray[1,,1])
	GridLength=array(0,2)
	GridCoordinates=array(0.,dim=c(2,MaxGridLength))

	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),LambdaParameter=as.double(LambdaParameter),PACKAGE="OriGen")

	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=array(GridAndCoordResults$GridCoordinates,c(2,MaxGridLength))

	IsLand=array(TRUE,dim=c(GridLength[1],GridLength[2]))
	if(MaskWater){
		#change points on water to false here...
		IsLand=.LandArray(GridCoordinates,GridLength)
	}
	NumberUnknowns=length(UnknownData[,1])
	UnknownGrids=array(0,dim=c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw=.Fortran("FITADMIXEDMODELFINDUNKNOWNS",AdmixtureFractions=as.double(UnknownGrids),DataArray=as.integer(DataArray),NumberSNPs=as.integer(NumberSNPs),GridLength=as.integer(GridLength),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),MaxGridLength=as.integer(MaxGridLength),SampleCoordinates=as.double(SampleCoordinates),GridCoordinates=as.double(GridCoordinates),NumberUnknowns=as.integer(NumberUnknowns),UnknownData=as.integer(UnknownData),IsLand=as.logical(IsLand),PACKAGE="OriGen")

	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
	ResultsRaw$AdmixtureFractions=array(ResultsRaw$AdmixtureFractions,c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
	ResultsRaw$UnknownData=array(ResultsRaw$UnknownData,c(NumberUnknowns,NumberSNPs))
	ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))
	ResultsRaw$IsLand=array(ResultsRaw$IsLand,c(GridLength[1],GridLength[2]))

	return(ResultsRaw)
}




RankSNPsLRT<-function(DataArray){
#This function takes in the PED file along with a location file and outputs the Likelihood Ratio ranking
#of each SNP followed by the Likelihood Ratio statistic and the Informativeness for assignment by rosenberg et al..  Note that the statistic is compares the assumption
#that there is just a single global population vs several different sites.

	SampleSites=length(DataArray[1,,1])
	NumberSNPs=length(DataArray[1,1,])

	Rankings=1:NumberSNPs
	LRT=array(0,c(2,NumberSNPs))
	ResultsRaw=.Fortran("CALC_ALL_RANKINGS",DataArray=as.integer(DataArray),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs),Rankings=as.integer(Rankings),LRT=as.double(LRT),PACKAGE="OriGen")

	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
	ResultsRaw$LRT=array(ResultsRaw$LRT,c(2,NumberSNPs))
	return(ResultsRaw)
}




ConvertMicrosatData<-function(DataFileName,LocationFileName){
	print("Note: This method assumes there are two files. DataFileName should have two initial columns, Location Name and Location Number, followed by a single column for each locus.   Each sample is represented by a pair of rows sharing the same Location Name and Location Number giving the two alleles at that locus.  Missing data should be given as -999 and individuals of unknown origin should be given Location Number -1.  The second file, LocationFileName, should have 4 columns: Location Name, Location Number, Latitude, and Longitude.  Both of these files should have headers.")

	#MicrosatData=read.table(DataFileName,header=TRUE)
	MicrosatData=read.table(DataFileName,header=TRUE,colClasses="factor")
	MicrosatData[[2]]=as.numeric(levels(MicrosatData[[2]]))[MicrosatData[[2]]]
	LocationData=read.table(LocationFileName,header=TRUE)
	NumberLoci=length(MicrosatData)-2
	#This counts the number of samples in the training dataset (not including the unknowns)
	NumberUnknowns=sum(MicrosatData[2]==-1)/2
	SampleSites=(length(MicrosatData[[1]]))/2-NumberUnknowns

	#This weeds out the unknown individuals.
	SubMicrosatData=MicrosatData[MicrosatData[[2]]!=-1,]

	AllelesAtLocus=0*1:NumberLoci
	for(i in 1:NumberLoci){
		AllelesAtLocus[i]=length(levels(MicrosatData[[i+2]]))
	}
	MaxAlleles=max(AllelesAtLocus)
	DataArray=array(0,c(MaxAlleles,SampleSites,NumberLoci))
	for(i in 1:NumberLoci){
		for(j in 1:SampleSites){
			DataArray[SubMicrosatData[2*j-1,i+2],j,i]=DataArray[SubMicrosatData[2*j-1,i+2],j,i]+1
			DataArray[SubMicrosatData[2*j,i+2],j,i]=DataArray[SubMicrosatData[2*j,i+2],j,i]+1
		}
	}

	SampleCoordinates=array(0,c(SampleSites,2))
	SampleCoordinates[,1]=LocationData$Longitude
	SampleCoordinates[,2]=LocationData$Latitude

	SubMicrosatData=MicrosatData[MicrosatData[[2]]==-1,]
	UnknownDataArray=array(0,c(NumberUnknowns,2,NumberLoci))
	for(i in 1:NumberLoci){
		for(j in 1:NumberUnknowns){
			UnknownDataArray[j,1,i]=as.integer(SubMicrosatData[2*j-1,i+2])
			UnknownDataArray[j,2,i]=as.integer(SubMicrosatData[2*j,i+2])
		}
	}

	ResultsRaw=list(DataArray=DataArray,SampleCoordinates=SampleCoordinates,AllelesAtLocus=AllelesAtLocus,MaxAlleles=MaxAlleles,SampleSites=SampleSites,NumberLoci=NumberLoci,NumberUnknowns=NumberUnknowns,UnknownDataArray=UnknownDataArray,LocationNames=MicrosatData[[1]],DataFileName=DataFileName,LocationFileName=LocationFileName)

	return(ResultsRaw)
}


FitMultinomialModel<-function(DataArray,SampleCoordinates,MaxGridLength=20,RhoParameter=10){
	#DataArray[MaxAlleles,SampleSites,NumberLoci] Gives the grouped data
	#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
	#This function takes in the data, fits the model, and returns the allele frequency surfaces
	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}
	NumberLoci=length(DataArray[1,1,])
	SampleSites=length(DataArray[1,,1])
	MaxAlleles=length(DataArray[,1,1])

	GridCoordinates=array(0,c(2,MaxGridLength))
	GridLength=array(0,2)

	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")

	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=GridAndCoordResults$GridCoordinates

	AllelesAtLocus=1:NumberLoci
	AllelesAtLocus[]=2
	for(i in 1:NumberLoci){
		for(j in 1:MaxAlleles){
			if(sum(DataArray[j,,i])>0.5){
				AllelesAtLocus[i]=j
			}
		}
	}

	AlleleFrequencySurfaces=array(0.,dim=c(MaxAlleles,GridLength[1],GridLength[2],NumberLoci))

	ResultsRaw=.Fortran("FITMULTINOMIALMODEL",AlleleFrequencySurfaces=as.double(AlleleFrequencySurfaces),DataArray=as.integer(DataArray),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),MaxAlleles=as.integer(MaxAlleles),NumberLoci=as.integer(NumberLoci),SampleCoordinates=as.double(SampleCoordinates),AllelesAtLocus=as.integer(AllelesAtLocus),GridCoordinates=as.double(GridCoordinates),PACKAGE="OriGen")

	ResultsRaw$AlleleFrequencySurfaces=array(ResultsRaw$AlleleFrequencySurfaces,c(MaxAlleles,GridLength[1],GridLength[2],NumberLoci))
	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(MaxAlleles,SampleSites,NumberLoci))
	ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
	ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))

	return(ResultsRaw)
}


.GenerateIsLandMatrix<-function(GridLength,GridCoordinates){
	ans=array(1,dim=c(GridLength[1],GridLength[2]))
	x.vec=rep(GridCoordinates[1,1:GridLength[1]],each=GridLength[2])
	y.vec=rep(GridCoordinates[2,1:GridLength[2]],times=GridLength[1])
	temp.vec=map.where(database="world",x.vec,y.vec)
	results.vec=x.vec*0+1
	for(i in 1:length(temp.vec)){
		if(is.na(temp.vec[i])){
			results.vec[i]=0
		}
	}
	ans=matrix(results.vec,nrow=GridLength[1],ncol=GridLength[2],byrow=TRUE)
	return(ans)
}


GenerateHeatMaps<-function(FitModelOutput,UnknownDataArray,NumberLoci,MaskWater=TRUE){
	NumberUnknowns=length(UnknownDataArray[,1,1])
	#NumberLoci=FitModelOutput$NumberLoci
	UnknownHeatMaps=array(1,dim=c(FitModelOutput$GridLength[1],FitModelOutput$GridLength[2],NumberUnknowns))
	for(i in 1:NumberUnknowns){
		for(j in 1:NumberLoci){
			Allele1=UnknownDataArray[i,1,j]
			Allele2=UnknownDataArray[i,2,j]

			#This takes care of alleles seen in the unknown data that aren't present in the data
			if(Allele1>FitModelOutput$AllelesAtLocus[j]){
				Allele1=0
			}
			if(Allele2>FitModelOutput$AllelesAtLocus[j]){
				Allele2=0
			}
			#Homo vs Hetero
			if(Allele1>0 & Allele2>0){
				if(Allele1==Allele2){
					UnknownHeatMaps[,,i]=UnknownHeatMaps[,,i]*FitModelOutput$AlleleFrequencySurfaces[Allele1,,,j]^2
				}else{
					UnknownHeatMaps[,,i]=UnknownHeatMaps[,,i]*2*FitModelOutput$AlleleFrequencySurfaces[Allele1,,,j]*FitModelOutput$AlleleFrequencySurfaces[Allele2,,,j]
				}
			}else if(Allele1>0){
				UnknownHeatMaps[,,i]=UnknownHeatMaps[,,i]*FitModelOutput$AlleleFrequencySurfaces[Allele1,,,j]
			}else if(Allele2>0){
				UnknownHeatMaps[,,i]=UnknownHeatMaps[,,i]*FitModelOutput$AlleleFrequencySurfaces[Allele2,,,j]
			}

			UnknownHeatMaps[,,i]=UnknownHeatMaps[,,i]/sum(UnknownHeatMaps[,,i])
		}
	}
	if(MaskWater){
		IsLandMatrix=.GenerateIsLandMatrix(FitModelOutput$GridLength,FitModelOutput$GridCoordinates)
		for(i in 1:NumberUnknowns){
			UnknownHeatMaps[,,i]=UnknownHeatMaps[,,i]*IsLandMatrix
		}
	}
	ans=FitModelOutput
	ans$UnknownGrids=UnknownHeatMaps
	ans$UnknownDataArray=UnknownDataArray
	ans$IsLandMatrix=IsLandMatrix
	return(ans)
}






FitMultinomialModelFindUnknowns<-function(DataArray,SampleCoordinates,UnknownDataArray,MaxGridLength=20,RhoParameter=10,MaskWater=TRUE){
	#DataArray[MaxAlleles,SampleSites,NumberLoci] Gives the grouped data
	#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
	#This function takes in the data, fits the model, and returns the allele frequency surfaces
	#UnknownDataArray[NumberUnknowns,2,NumberLoci] gives the unknown data

	NumberLoci=length(DataArray[1,1,])

	Surfaces=FitMultinomialModel(DataArray,SampleCoordinates,MaxGridLength,RhoParameter)
	ResultsRaw=GenerateHeatMaps(Surfaces,UnknownDataArray,NumberLoci,MaskWater)

	# if(!.is.wholenumber(MaxGridLength)){
		# stop("MaxGridLength must be an integer")
	# }
	# if(MaxGridLength<=1){
		# stop("MaxGridLength must be greater than 1")
	# }
	# if(RhoParameter<=0){
		# stop("RhoParameter must be greater than 0")
	# }
	# if(length(SampleCoordinates[1,])!=2){
		# stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	# }
	# NumberLoci=length(DataArray[1,1,])
	# SampleSites=length(DataArray[1,,1])
	# MaxAlleles=length(DataArray[,1,1])

	# GridCoordinates=array(0,c(2,MaxGridLength))
	# GridLength=array(0,2)

	# GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")

	# GridLength=GridAndCoordResults$GridLength
	# GridCoordinates=GridAndCoordResults$GridCoordinates

	# AllelesAtLocus=1:NumberLoci
	# AllelesAtLocus[]=2
	# for(i in 1:NumberLoci){
		# for(j in 1:MaxAlleles){
			# if(sum(DataArray[j,,i])>0.5){
				# AllelesAtLocus[i]=j
			# }
		# }
	# }

	# NumberUnknowns=length(UnknownDataArray[,1,1])
	# UnknownGrids=array(0.,dim=c(GridLength[1],GridLength[2],NumberUnknowns))


	# #AlleleFrequencySurfaces=array(0.,dim=c(MaxAlleles,GridLength[1],GridLength[2],NumberLoci))
	# ResultsRaw=.Fortran("FITMULTINOMIALMODELFIND",UnknownGrids=as.double(UnknownGrids),DataArray=as.integer(DataArray),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),MaxAlleles=as.integer(MaxAlleles),NumberLoci=as.integer(NumberLoci),SampleCoordinates=as.double(SampleCoordinates),GridCoordinates=as.double(GridCoordinates),AllelesAtLocus=as.integer(AllelesAtLocus),NumberUnknowns=as.integer(NumberUnknowns),UnknownDataArray=as.integer(UnknownDataArray),PACKAGE="OriGen")

	# ResultsRaw$UnknownGrids=array(ResultsRaw$UnknownGrids,c(GridLength[1],GridLength[2],NumberUnknowns))
	# ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(MaxAlleles,SampleSites,NumberLoci))
	# ResultsRaw$UnknownDataArray=array(ResultsRaw$UnknownDataArray,c(NumberUnknowns,2,NumberLoci))
	# ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
	# ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))

	return(ResultsRaw)
}


CalcFractionsMultiLoglik<-function(UnknownDataArray,LambdaParameter=100){
	#This function takes the UnknownDataArray which contains allelelic information
	#for individuals WITHIN a single sample site and calculates the resulting
	#fraction loglikelihood for placing all individuals 100% back into their site
	#Note that we force the fortran code to do a single site by pretending it is
	#2 equal sites.
	#UnknownDataArray[NumberUnknowns,2,NumberLoci] lists the two allele numbers of the unknown data
	#Loglikelihoods[NumberUnknowns] lists the loglikelihoods for each unknown individual

	NumberUnknowns=length(UnknownDataArray[,1,1])
	NumberLoci=length(UnknownDataArray[1,1,])
	Loglikelihoods=array(0,dim=c(NumberUnknowns))
	Nodes=2
	MaxAlleles=0
	for(i in 1:NumberLoci){
		TempInt=length(unique(as.vector(UnknownDataArray[,,i])))
		if(TempInt>MaxAlleles){
			MaxAlleles=TempInt
		}
	}

	UnknownDataArrayFact=UnknownDataArray
	AlleleCountsPop=array(0,dim=c(NumberLoci,MaxAlleles,Nodes))
	for(i in 1:NumberLoci){
		TempFact=as.factor(UnknownDataArray[,,i])
		for(j in 1:NumberUnknowns){
			UnknownDataArrayFact[j,1,i]=as.numeric(TempFact[j])
			UnknownDataArrayFact[j,2,i]=as.numeric(TempFact[j+NumberUnknowns])
		}
		for(j in 1:(2*NumberUnknowns)){
			for(k in 1:Nodes){
				AlleleCountsPop[i,TempFact[j],k]=AlleleCountsPop[i,TempFact[j],k]+1
			}
		}
	}
	AlleleFrequencyPop=AlleleCountsPop
	for(i in 1:NumberLoci){
		for(j in 1:MaxAlleles){
			for(k in 1:Nodes){
				AlleleFrequencyPop[i,j,k]=AlleleCountsPop[i,j,k]/sum(AlleleCountsPop[i,,k])
			}
		}
	}

	for(i in 1:NumberUnknowns){
		ResultsRaw=.Fortran("CALC_FRACTIONS_MULTI_LOGLIK_SUB",UnknownIndivFactored=as.integer(t(UnknownDataArrayFact[i,,])),GenomeFractions=as.double(c(1,0)),AlleleFrequencyPop=as.double(AlleleFrequencyPop),Lambda=as.double(LambdaParameter),Nodes=as.integer(Nodes),NumberLoci=as.integer(NumberLoci),MaxAlleles=as.integer(MaxAlleles),Loglik=as.double(Loglikelihoods[i]),PACKAGE="OriGen")

		Loglikelihoods[i]=ResultsRaw$Loglik
	}
	return(Loglikelihoods)
}#end CalcFractionsMultiLoglik



FitMultinomialAdmixedModelFindUnknowns<-function(DataArray,SampleCoordinates,UnknownDataArray,
	MaxGridLength=20,RhoParameter=10,LambdaParameter=100,MaskWater=TRUE,NumberLoci=-1){
	#DataArray[MaxAlleles,SampleSites,NumberLoci] Gives the grouped data
	#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
	#This function takes in the data, fits the model, and returns the allele frequency surfaces
	#UnknownDataArray[NumberUnknowns,2,NumberLoci] lists the two allele numbers of the unknown data

	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}

	if(NumberLoci==-1){
		print("Using all loci")
		NumberLoci=length(DataArray[1,1,])
	}
	SampleSites=length(DataArray[1,,1])
	MaxAlleles=length(DataArray[,1,1])

	GridCoordinates=array(0,c(2,MaxGridLength))
	GridLength=array(0,2)

	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")

	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=array(GridAndCoordResults$GridCoordinates,c(2,MaxGridLength))

	AllelesAtLocus=1:NumberLoci
	AllelesAtLocus[]=2
	for(i in 1:NumberLoci){
		for(j in 1:MaxAlleles){
			if(sum(DataArray[j,,i])>0.5){
				AllelesAtLocus[i]=j
			}
		}
	}

	IsLand=array(TRUE,dim=c(GridLength[1],GridLength[2]))
	if(MaskWater){
		#change points on water to false here...
		#IsLand=.GenerateIsLandMatrix(GridLength,FitModelOutput$Coordinates)
		IsLand=.LandArray(GridCoordinates,GridLength)
	}

	NumberUnknowns=length(UnknownDataArray[,1,1])
	AdmixtureFractions=array(0,dim=c(GridLength[1],GridLength[2],NumberUnknowns))
	Loglikelihoods=array(0,dim=NumberUnknowns)

	#AlleleFrequencySurfaces=array(0.,dim=c(MaxAlleles,GridLength[1],GridLength[2],NumberLoci))
	#ResultsRaw=.Fortran("FITMULTINOMIALMODEL",AlleleFrequencySurfaces=as.double(AlleleFrequencySurfaces),DataArray=as.integer(DataArray),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),MaxAlleles=as.integer(MaxAlleles),NumberLoci=as.integer(NumberLoci),SampleCoordinates=as.double(SampleCoordinates),AllelesAtLocus=as.integer(AllelesAtLocus),GridCoordinates=as.double(GridCoordinates),PACKAGE="OriGen")

	ResultsRaw=.Fortran("FITMULTIADMIXEDMODELFINDUNKNOWNS",AdmixtureFractions=as.double(AdmixtureFractions),DataArray=as.integer(DataArray),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),MaxAlleles=as.integer(MaxAlleles),NumberLoci=as.integer(NumberLoci),SampleCoordinates=as.double(SampleCoordinates),GridCoordinates=as.double(GridCoordinates),AllelesAtLocus=as.integer(AllelesAtLocus),NumberUnknowns=as.integer(NumberUnknowns),UnknownDataArray=as.integer(UnknownDataArray),IsLand=as.logical(IsLand),Loglikelihoods=as.double(Loglikelihoods),PACKAGE="OriGen")

	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(MaxAlleles,SampleSites,NumberLoci))
	ResultsRaw$AdmixtureFractions=array(ResultsRaw$AdmixtureFractions,c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
	ResultsRaw$UnknownDataArray=array(ResultsRaw$UnknownDataArray,c(NumberUnknowns,2,NumberLoci))
	ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))
	ResultsRaw$IsLand=array(ResultsRaw$IsLand,c(GridLength[1],GridLength[2]))

	return(ResultsRaw)
}











#-----------------------------------------------------------------------------------------------------

#The below functions will not be included in V1 of the R package....

#-----------------------------------------------------------------------------------------------------


#This function requires the maps package to work
.IsLand<-function(x.vec,y.vec){
	#require("maps")
	temp.vec=maps::map.where(database="world",x.vec,y.vec)
	result.vec=x.vec*0+1
	for(i in 1:length(temp.vec)){
		if(is.na(temp.vec[i])){
			result.vec[i]=0
		}
	}
	return(result.vec)
}

#this function requires the maps package to work
.MaskWater<-function(GridCoordinates){
	#this short code checks whether the given coordinates are in water and outputs a matrix with 1 meaning land
	#and 0 meaning water...
	#GridCoordinates should be a matrix[x,2] where x is the number of grid points and the first 2 is Long,Lat
	ndiv=length(GridCoordinates[,1])
	latcount=0
	longcount=0
	for(i in 1:ndiv){
		if(GridCoordinates[i,2]>0.001){
			latcount=latcount+1
		}else if(GridCoordinates[i,2]< -0.001){
			latcount=latcount+1
		}
		if(GridCoordinates[i,1]>0.001){
			longcount=longcount+1
		}else if(GridCoordinates[i,1]< -0.001){
			longcount=longcount+1
		}
	}
	temp.mat=mat.or.vec(nc=longcount,nr=latcount)
	temp.mat[,]=1
	for(i in 1:longcount){
		temp.mat[,i]=.IsLand(rep(GridCoordinates[i,1],each=latcount),GridCoordinates[,2])
		}
	#write.table(temp.mat[latcount:1,],file="GridCoordSquare40Water.txt",append=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
	return(temp.mat)
}


#This function requires the maps package to work
.IsLandBool<-function(x.vec,y.vec){
	#require("maps")
	temp.vec=maps::map.where(database="world",x.vec,y.vec)
	result.vec=x.vec
	result.vec[]=TRUE
	for(i in 1:length(temp.vec)){
		if(is.na(temp.vec[i])){
			result.vec[i]=FALSE
		}
	}
	return(result.vec)
}

.LandArray<-function(GridCoordinates,GridLength){
	#this short code checks whether the given coordinates are in water and outputs a matrix with 1 meaning land
	#and 0 meaning water...
	#GridCoordinates should be a matrix[2,x] where x is the number of grid points and the first 2 is Long,Lat
	#ndiv=length(GridCoordinates[1,])

	temp.mat=mat.or.vec(nr=GridLength[1],nc=GridLength[2])
	temp.mat[,]=TRUE
	for(i in 1:GridLength[1]){
		temp.mat[i,]=.IsLandBool(rep(GridCoordinates[1,i],each=GridLength[2]),GridCoordinates[2,1:GridLength[2]])
	}
	return(temp.mat)
}




#this function requires packages ggplot2 and maps to work.  Note that the vectors on the maps package is outdated particularly in europe
#An updated map can be downloaded from http://www.naturalearthdata.com/downloads/50m-cultural-vectors/


PlotAlleleFrequencySurfaceOld<-function(AlleleSurfaceOutput,SNPNumber=1,MaskWater=TRUE){
	#GridCoordinates(2,MaxGridLength)
	print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.")
	#require("maps")
	#require("ggplot2")

	#note this next line is in the function merely to pass R checks.  It serves no other purpose.
	Land=Lat=Long=Frequency=NULL

	TempHM=AlleleSurfaceOutput$AlleleFrequencySurfaces[SNPNumber,,]
	for(i in 1:AlleleSurfaceOutput$GridLength[1]){
		TempHM[i,]=AlleleSurfaceOutput$GridCoordinates[1,i]
	}
	TempOb<-data.frame(Frequency=as.vector(AlleleSurfaceOutput$AlleleFrequencySurfaces[SNPNumber,,]),Long=as.vector(TempHM))
	for(i in 1:AlleleSurfaceOutput$GridLength[2]){
		TempHM[,i]=AlleleSurfaceOutput$GridCoordinates[2,i]
	}
	TempOb$Lat=as.vector(TempHM)
	TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
	subdata=subset(TempOb,Land==1)
	#minp=min(subdata$Frequency)
	minp=0
	#maxp=max(subdata$Frequency)
	maxp=1
	if(MaskWater){
		subdata=subset(TempOb,Land==1)
		#minp=min(subdata$Frequency)
		#minp=0
		#maxp=max(subdata$Frequency)
		p<-ggplot2::ggplot(subset(TempOb,Land==1),aes(Long,Lat))
		} else {
		#minp=min(TempOb$Frequency)
		#minp=0
		#maxp=max(TempOb$Frequency)
		p<-ggplot2::ggplot(TempOb,aes(Long,Lat))
		}
	p+	annotation_map(map_data("world"), fill=NA, colour = "white")+
		geom_tile(aes(fill=Frequency),colour=NA,alpha=1) +
		scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
		annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) +
		ylab("Latitude") + ggtitle(paste0("Allele Frequency Surface SNP:",SNPNumber)) +
		xlab("Longitude")
}



PlotAlleleFrequencySurface<-function(AlleleSurfaceOutput,LocusNumber=1,AlleleNumber=1,MaskWater=TRUE,Scale=FALSE){
#AlleleFrequencySurfaces=array(0.,dim=c(MaxAlleles,GridLength[1],GridLength[2],NumberLoci))
#GridCoordinates(2,MaxGridLength)
print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.")
print("Note: Setting AlleleNumber = 0 when using microsatellites plots all the alleles in a grid.")
#require("maps")
#require("ggplot2")

#note this next line is in the function merely to pass R checks.  It serves no other purpose.
Land=Lat=Long=Frequency=NULL

if(length(dim(AlleleSurfaceOutput$AlleleFrequencySurfaces))==3){
	#We are dealing with the SNP case
	SNPBool=TRUE
	TempHM=AlleleSurfaceOutput$AlleleFrequencySurfaces[LocusNumber,,]
}else{
	#Dealing with markers
	SNPBool=FALSE
	TempHM=AlleleSurfaceOutput$AlleleFrequencySurfaces[AlleleNumber,,,LocusNumber]
}

if(SNPBool==FALSE & AlleleNumber == 0){
	NumberAlleles=AlleleSurfaceOutput$AllelesAtLocus[LocusNumber]
	TempHM=AlleleSurfaceOutput$AlleleFrequencySurfaces[1,,,LocusNumber]

	for(i in 1:AlleleSurfaceOutput$GridLength[1]){
		TempHM[i,]=AlleleSurfaceOutput$GridCoordinates[1,i]
	}

	TempOb2<-data.frame(Frequency=as.vector(AlleleSurfaceOutput$AlleleFrequencySurfaces[1:NumberAlleles,,,LocusNumber]),Long=rep(as.vector(TempHM),each=NumberAlleles),Allele=rep(1:NumberAlleles,times=length(TempHM)))
	for(i in 1:AlleleSurfaceOutput$GridLength[2]){
		TempHM[,i]=AlleleSurfaceOutput$GridCoordinates[2,i]
	}
	TempOb2$Lat=rep(as.vector(TempHM),each=NumberAlleles)
	TempOb2$Land=.IsLand(TempOb2$Long,TempOb2$Lat)

	#subdata=subset(TempOb,Land==1)
	#minp=min(subdata$Frequency)
	minp=0
	#maxp=max(subdata$Frequency)
	maxp=1
	if(Scale){
			maxp = max(TempOb2$Frequency)
		}
	if(MaskWater){
		#subdata=subset(TempOb,Land==1)
		#minp=min(subdata$Frequency)
		#minp=0
		#maxp=max(subdata$Frequency)
		p<-ggplot2::ggplot(subset(TempOb2,Land==1),aes(Long,Lat))
		} else {
		#minp=min(TempOb$Frequency)
		#minp=0
		#maxp=max(TempOb$Frequency)
		p<-ggplot2::ggplot(TempOb2,aes(Long,Lat))
		}
	p+annotation_map(map_data("world"), fill=NA, colour = "white")+
		geom_tile(aes(fill=Frequency),colour=NA,alpha=1) +
		scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
		annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) +
		ylab("Latitude") + ggtitle(paste0("Allele Frequency Surfaces for Locus:",LocusNumber)) +
		xlab("Longitude") + facet_wrap(~Allele)

}else{
	for(i in 1:AlleleSurfaceOutput$GridLength[1]){
		TempHM[i,]=AlleleSurfaceOutput$GridCoordinates[1,i]
	}

	if(SNPBool){
		TempOb<-data.frame(Frequency=as.vector(AlleleSurfaceOutput$AlleleFrequencySurfaces[LocusNumber,,]),Long=as.vector(TempHM))
	}else{
		TempOb<-data.frame(Frequency=as.vector(AlleleSurfaceOutput$AlleleFrequencySurfaces[AlleleNumber,,,LocusNumber]),Long=as.vector(TempHM))
	}

	for(i in 1:AlleleSurfaceOutput$GridLength[2]){
		TempHM[,i]=AlleleSurfaceOutput$GridCoordinates[2,i]
	}
	TempOb$Lat=as.vector(TempHM)
	TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
	subdata=subset(TempOb,Land==1)
	#minp=min(subdata$Frequency)
	minp=0
	#maxp=max(subdata$Frequency)
	maxp=1
	if(Scale){
			maxp = max(TempOb$Frequency)
		}
	if(MaskWater){
		subdata=subset(TempOb,Land==1)
		#minp=min(subdata$Frequency)
		#minp=0
		#maxp=max(subdata$Frequency)
		p<-ggplot2::ggplot(subset(TempOb,Land==1),aes(Long,Lat))
		} else {
		#minp=min(TempOb$Frequency)
		#minp=0
		#maxp=max(TempOb$Frequency)
		p<-ggplot2::ggplot(TempOb,aes(Long,Lat))
		}

		#temp2=data.frame(Frequency=as.vector(AlleleSurfaceOutput$Frequencies[AlleleNumber,,LocusNumber]),Long=as.vector(AlleleSurfaceOutput$SampleCoordinates[,1]),Lat=as.vector(AlleleSurfaceOutput$SampleCoordinates[,2]))
	p+	annotation_map(map_data("world"), fill=NA, colour = "white")+
		geom_tile(aes(fill=Frequency),colour=NA,alpha=1) +
		scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
		#geom_point(data=temp2,aes(Long,Lat,colour=Frequency)) +
		#scale_colour_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
		annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) +
		ylab("Latitude") + ggtitle(paste0("Allele Frequency Surface Locus:",LocusNumber," Allele:",AlleleNumber)) +
		xlab("Longitude")
	}
}



##doing a facet wrap of all the alleles for a given locus
#PlotAlleleFrequencySurfaceAll<-function(AlleleSurfaceOutput,LocusNumber=1,NumberAlleles=0,MaskWater=TRUE,Scale=FALSE){
##Default plots all alleles
##AlleleFrequencySurfaces=array(0.,dim=c(MaxAlleles,GridLength[1],GridLength[2],NumberLoci))
##GridCoordinates(2,MaxGridLength)
#print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.")
##require("maps")
##require("ggplot2")
#if(NumberAlleles==0){
#	NumberAlleles=AlleleSurfaceOutput$AllelesAtLocus[LocusNumber]
#}
#
#TempHM=AlleleSurfaceOutput$AlleleFrequencySurfaces[1,,,LocusNumber]
#
#for(i in 1:AlleleSurfaceOutput$GridLength[1]){
#	TempHM[i,]=AlleleSurfaceOutput$GridCoordinates[1,i]
#}
#
#TempOb2<-data.frame(Frequency=as.vector(AlleleSurfaceOutput$AlleleFrequencySurfaces[1:NumberAlleles,,,LocusNumber]),Long=rep(as.vector(TempHM),each=NumberAlleles),Allele=rep(1:NumberAlleles,times=length(TempHM)))
#for(i in 1:AlleleSurfaceOutput$GridLength[2]){
#	TempHM[,i]=AlleleSurfaceOutput$GridCoordinates[2,i]
#}
#TempOb2$Lat=rep(as.vector(TempHM),each=NumberAlleles)
#TempOb2$Land=.IsLand(TempOb2$Long,TempOb2$Lat)
#
##subdata=subset(TempOb,Land==1)
##minp=min(subdata$Frequency)
#minp=0
##maxp=max(subdata$Frequency)
#maxp=1
#if(Scale){
#    	maxp = max(TempOb2$Frequency)
#    }
#if(MaskWater){
#	#subdata=subset(TempOb,Land==1)
#	#minp=min(subdata$Frequency)
#	#minp=0
#	#maxp=max(subdata$Frequency)
#	p<-ggplot(subset(TempOb2,Land==1),aes(Long,Lat))
#	} else {
#	#minp=min(TempOb$Frequency)
#	#minp=0
#	#maxp=max(TempOb$Frequency)
#	p<-ggplot(TempOb2,aes(Long,Lat))
#	}
#p+annotation_map(map_data("world"), fill=NA, colour = "white")+
#	geom_tile(aes(fill=Frequency),colour=NA,alpha=1) +
#	scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
#	annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) +
#	ylab("Latitude") + ggtitle(paste0("Allele Frequency Surfaces for Locus:",LocusNumber)) +
#	xlab("Longitude") + facet_wrap(~Allele)
#
#}


PlotUnknownHeatMap<-function(HeatMapOutput,UnknownNumber=1,MaskWater=TRUE){
#GridCoordinates(2,MaxGridLength)
print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.")
#require("maps")
#require("ggplot2")

#note this next line is in the function merely to pass R checks.  It serves no other purpose.
Land=Lat=Long=Probability=NULL

TempHM=HeatMapOutput$UnknownGrids[,,UnknownNumber]
for(i in 1:HeatMapOutput$GridLength[1]){
	TempHM[i,]=HeatMapOutput$GridCoordinates[1,i]
}
TempOb<-data.frame(Probability=as.vector(HeatMapOutput$UnknownGrids[,,UnknownNumber]),Long=as.vector(TempHM))
for(i in 1:HeatMapOutput$GridLength[2]){
	TempHM[,i]=HeatMapOutput$GridCoordinates[2,i]
}
TempOb$Lat=as.vector(TempHM)
TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
subdata=subset(TempOb,Land==1)
#minp=min(subdata$Probability)
minp=0
maxp=1
if(MaskWater){
	subdata=subset(TempOb,Land==1)
	#minp=min(subdata$Probability)
	maxp=max(subdata$Probability)
	p<-ggplot2::ggplot(subset(TempOb,Land==1),aes(Long,Lat))
	} else {
	#minp=min(TempOb$Probability)
	maxp=max(TempOb$Probability)
	p<-ggplot2::ggplot(TempOb,aes(Long,Lat))
	}
p+	annotation_map(map_data("world"), fill=NA, colour = "white")+
	geom_tile(aes(fill=Probability),colour=NA,alpha=1) +
	scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
	annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) +
	ylab("Latitude") + ggtitle(paste0("Heat Map Surface Individual:",UnknownNumber)) +
	xlab("Longitude")
}


#.PlotAllUnknowns<-function(HeatMapOutput,NamesList=NULL,MaskWater=TRUE){
##GridCoordinates(2,MaxGridLength)
#print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.")
##require("maps")
##require("ggplot2")

#NumberUnknowns=HeatMapOutput$

#if(is.null(NamesList)){
#	NamesList=1:NumberUnknowns
#}

#TempHM=HeatMapOutput$UnknownGrids[,,UnknownNumber]
#for(i in 1:HeatMapOutput$GridLength[1]){
#	TempHM[i,]=HeatMapOutput$GridCoordinates[1,i]
#}
#TempOb<-data.frame(Probability=as.vector(HeatMapOutput$UnknownGrids[,,UnknownNumber]),Long=as.vector(TempHM))
#for(i in 1:HeatMapOutput$GridLength[2]){
#	TempHM[,i]=HeatMapOutput$GridCoordinates[2,i]
#}
#TempOb$Lat=as.vector(TempHM)
#TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
#subdata=subset(TempOb,Land==1)
##minp=min(subdata$Probability)
#minp=0
#maxp=max(subdata$Probability)
#if(MaskWater){
#	subdata=subset(TempOb,Land==1)
#	#minp=min(subdata$Probability)
#	minp=0
#	maxp=max(subdata$Probability)
#	p<-ggplot(subset(TempOb,Land==1),aes(Long,Lat))
#	} else {
#	#minp=min(TempOb$Probability)
#	minp=0
#	maxp=max(TempOb$Probability)
#	p<-ggplot(TempOb,aes(Long,Lat))
#	}


#BestLocations=data.frame(Labels=NamesList)
#TempLong=1:NumberToPlot
#TempLat=1:NumberToPlot
#for(i in 1:NumberToPlot){
#	#need to find the maximum location of HeatMapOutput$UnknownGrids[,,i]
#	if(MaskWater){
#		#i need a matrix logical showing water...
#		IntLoc=which((HeatMapOutput$UnknownGrids[,,i] == max(HeatMapOutput$UnknownGrids[,,i])&TRUE), arr.ind = TRUE)
#	}else{
#		IntLoc=which(HeatMapOutput$UnknownGrids[,,i]) == max(HeatMapOutput$UnknownGrids[,,i]), arr.ind = TRUE)
#	}
#	TempLong[i]=HeatMapOutput$GridCoordinates[1,IntLoc[1]]
#	TempLat[i]=HeatMapOutput$GridCoordinates[2,IntLoc[2]]
#	}
#BestLocations$Long=TempLong
#BestLocations$Lat=TempLat
#
#p+	annotation_map(map_data("world"), fill=NA, colour = "white")+
#	geom_tile(aes(fill=Probability),colour=NA,alpha=1) +
#	scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
#	annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) +
#	ylab("Latitude") + ggtitle(paste0("Best Locations, Total: ",NumberToPlot)) +
#	xlab("Longitude")
#}




PlotAdmixedSurface<-function(AdmixedOutput,UnknownNumber=1,Percent=FALSE,Title=NULL,MaskWater=TRUE){
#GridCoordinates(2,MaxGridLength)
print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.")
#require("maps")
#require("ggplot2")

#note this next line is in the function merely to pass R checks.  It serves no other purpose.
Rounded=Lat=Long=Fractions=NULL

TempHM=AdmixedOutput$AdmixtureFractions[,,UnknownNumber]
#TempHM=AdmixedOutput$UnknownGrids[,,UnknownNumber]
for(i in 1:AdmixedOutput$GridLength[1]){
	TempHM[i,]=AdmixedOutput$GridCoordinates[1,i]
	}
TempOb<-data.frame(Fractions=as.vector(AdmixedOutput$AdmixtureFractions[,,UnknownNumber]),Long=as.vector(TempHM))
for(i in 1:AdmixedOutput$GridLength[2]){
	TempHM[,i]=AdmixedOutput$GridCoordinates[2,i]
	}
TempOb$Lat=as.vector(TempHM)
TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
if(Percent){
	TempOb$Rounded<-round(100*TempOb$Fractions, digits=1)
}else{
	TempOb$Rounded<-round(TempOb$Fractions, digits=2)
}
if(!is.character(Title)){
	print("Title must be a character vector. Title set to default.")
	Title=NULL
}
if(is.null(Title)){
	if(Percent){
		Title="Admixture Percentages"
	}else{
		Title="Admixture Fractions"
	}
}
subdata=subset(TempOb,Fractions>=0.01)

p<-ggplot(TempOb,aes(Long,Lat))
p+theme(panel.background = element_rect(fill = "lightskyblue1")) +
		annotation_map(map_data("world"), fill="darkolivegreen3", colour = "white")+
		annotation_map(map_data("world"),fill="NA",col="grey10") +
		theme(legend.position="none") +
		geom_text(aes(label=Rounded),alpha=0,size=4) +
		geom_text(data=subdata,aes(Long,Lat,label=Rounded),alpha=1,size=4) +
		theme(legend.position="none") +
		ylab("Latitude")+xlab("Longitude")+ggtitle(paste0(Title))
}
