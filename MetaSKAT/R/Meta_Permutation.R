
ReadPermu_Header<-function(File.MPermu){

	con = file(File.MPermu, "rb")
	header<-readBin(con, integer(), n = 10, size = 8)
	#cat(header, "\n")
	
	# n.permu, n.all, n, nSets, nSNPs, nSNPs.unique
	re = list(con=con, n.permu=header[2], n.all=header[3], n=header[4], nSets=header[5], nSNPs.unique=header[6])
	if(header[1] != 1){
		
		close(con)
		stop("Verion information in File.MPermu is not correct!")
	}
	return(re)

}

# one element of EachInfo.Permu
GetPermu_Score<-function(con, nSNP, nPermu,  StartPosPermu){


	#cat("Start:", StartPosPermu, "\n")
	#StartPosPermu = 8 * 10
	seek(con, where = StartPosPermu, origin = "start")
	out = readBin(con, double(), n = nSNP * (nPermu +1) , size = 8)
	out.m = matrix(out, byrow=TRUE, nrow=nSNP)
	return(out.m)
}


Open_MPermu_File_2Read<-function(File.MPermu.vec, File.MInfo.vec){

	n.cohort<-length(File.MInfo.vec)
	
	if(length(File.MInfo.vec) != length(File.MPermu.vec)){
		stop("Different numbers of Meta Info and Permu files!")
	}
	
	cat("Number of cohorts = ", n.cohort, "\n")
	
	# Check the existence of files 
	for(i in 1:n.cohort){
		File.MPermu.vec[i]<-normalizePath(File.MPermu.vec[i] ,mustWork =FALSE)
		File.MInfo.vec[i]<-normalizePath(File.MInfo.vec[i] ,mustWork =FALSE)	
	
		SKAT:::Check_File_Exists(File.MPermu.vec[i])
		SKAT:::Check_File_Exists(File.MInfo.vec[i])
	}
	
	# Read files
	
	re<-list()
	re.Permu<-list()
	re.SetInfo<-list()
	for(i in 1:n.cohort){
		file.idx<-i-1
		File.MPermu<-File.MPermu.vec[i] 
		File.MetaInfo<-File.MInfo.vec[i] 
		
		data.info<-Read_Info_File(File.MetaInfo)
		re.Permu[[i]]<-ReadPermu_Header(File.MPermu)
		re[[i]]<-data.info
		
	}
	
	# Get unique sets
	Set_unique<-NULL
	for(i in 1:n.cohort){
		Set_unique<-union(Set_unique, re[[i]]$set_unique)
	}	
	
	info<-list(n.cohort=n.cohort, Set_unique = Set_unique, EachInfo=re, EachInfo.Permu=re.Permu)
	
	return(info)

}


GetPermu_obj<-function(Permu.Info, SetID){

	n.cohort = Permu.Info$n.cohort	
	IsExistSNV<-rep(0, n.cohort)
	Info.list<-list()
	Permu.list<-list()
	Score.list<-list()
	N.Permu<-rep(0, n.cohort)
	
	for(i in 1:n.cohort){
		idx<-Permu.Info$EachInfo[[i]]$hash_set[[SetID]]
		
		if(is.null(idx)){
			IsExistSNV[i]<-0
			
		} else {
			IsExistSNV[i]<-1
			nSNP = length(idx)
			N.Permu[i] = Permu.Info$EachInfo.Permu[[i]]$n.permu

			Info.list[[i]]<-Permu.Info$EachInfo[[i]]$Info[idx,]
			StartPosPermu = Info.list[[i]]$StartPOSPermu[1]
			out.m=GetPermu_Score(Permu.Info$EachInfo.Permu[[i]]$con, nSNP , N.Permu[i],  StartPosPermu)
			#out.m1<<-out.m
			#score1<<-Info.list[[i]]$Score
			
			Permu.list[[i]] = out.m[,-1]
			Score.list[[i]] = out.m[,1]
			
			# check score.list
			#cat("Check: ", sum((Score.list[[i]] - Info.list[[i]]$Score)^2), "\n")
			#cat(Score.list[[i]], "\n")
			#cat(Info.list[[i]]$Score, "\n")
		}
	}
	#Info.list1<<-Info.list
	
	obj.oneset = Get_META_Data_OneSet_Align(SMat.list=NULL, Info.list=Info.list, IsExistSNV=IsExistSNV,  n.cohort=n.cohort, Is.SMat=FALSE)
	n.all = obj.oneset$n.all
	
	Permu.list.new<-list()
	Score.list.new<-list()
	
	for(i in 1:n.cohort){
		n1 = N.Permu[i] 
		Permu.list.new[[i]]<-matrix(rep(0, n1*n.all), ncol=n1)
		Score.list.new[[i]]<-rep(0, n.all)
		
		if(IsExistSNV[i] == 1){	
			IDX<-obj.oneset$IDX.list[[i]]
			IDX1<-obj.oneset$IDX1.list[[i]]
			
			Permu.list.new[[i]][IDX,]<-Permu.list[[i]][IDX1,] 
			Permu.list.new[[i]]<-Permu.list.new[[i]]* obj.oneset$Sign.list[[i]]
			
			Score.list.new[[i]][IDX]<-Score.list[[i]][IDX1] 
			Score.list.new[[i]]<-Score.list.new[[i]]* obj.oneset$Sign.list[[i]]
			
		} 
	}
	
	#A2<<-obj.oneset
	re=list(Info.list=obj.oneset$Info.list, Permu.list.new=Permu.list.new, Score.list.new=Score.list.new, n.cohort=n.cohort, N.Permu = N.Permu)
	return(re)	
}


MetaSKAT_MPermu_OneSet<-function(Permu.Info, SetID, n.Resampling=10000, r.corr=0, weights.beta=c(1,25), MAF.cutoff=1){

	n.cohort = Permu.Info$n.cohort	
	n1=rep(1, n.cohort)
	obj=GetPermu_obj(Permu.Info, SetID)
	
	for(i in 1:n.cohort){
		
		idx_miss<-which(is.na(obj$Info.list[[i]]$MAF))
		if(length(idx_miss) > 0){
			obj$Info.list[[i]]$Score[idx_miss] = 0
			obj$Info.list[[i]]$MAF[idx_miss] = 0
		}	
		
		n1[i]<-Permu.Info$EachInfo[[i]]$header$N.Sample
	}
	
	##########################
	# Get Combined MAF
	
	MAF.Combine=0
	weight.list<-list()

	MAF.list<-list()
	for(i in 1:n.cohort){
	
		MAF.list[[i]]<-obj$Info.list[[i]]$MAF
	}
	for(i in 1:n.cohort){
	
		MAF.Combine = MAF.Combine + MAF.list[[i]] * n1[i] / sum(n1)
		
	}
	for(i in 1:n.cohort){
		weight.list[[i]]<-Beta.Weights(MAF.Combine,weights.beta, MAF.cutoff, Is.MAF=TRUE)
	}

	#
	#	Resampling
	#

	#A1<<-obj
	Score1<-NULL
	Score.mat<-NULL
	for(i in 1:n.cohort){
		Score.temp = obj$Score.list.new[[i]]
		
		if(i==1){
			Score1 = Score.temp *weight.list[[i]]
		} else {
			Score1 = Score1 + Score.temp *weight.list[[i]]
		}	
	}
	
	if(r.corr==0){
		TestStat = sum(Score1^2)
	} else {
		TestStat = sum(Score1)^2
	}
	
	nRun<-ceiling(n.Resampling /50000)
	n.Resampling.item = 50000
	TestStat.RA<-NULL
	for(k in 1:nRun){
		for(i in 1:n.cohort){
	
			id<-sample.int(obj$N.Permu[i] ,n.Resampling.item, replace = TRUE)
	
			if(i==1){
				Score.mat<-obj$Permu.list.new[[i]][,id] *weight.list[[i]]
			} else {
				Score.mat = Score.mat + obj$Permu.list.new[[i]][,id] *weight.list[[i]]
			}
	
		}
		if(r.corr==0){
			TestStat.R<-colSums(Score.mat^2)
		} else {
			TestStat.R = colSums(Score.mat)^2
		}
		TestStat.RA = c(TestStat.RA, TestStat.R)
	}
	
	#TestStat1<<-TestStat
	#TestStat.R1<<-TestStat.R
	
	
	pval<-(length(which(TestStat.RA >= TestStat))+1) /(length(TestStat.RA)+1)
	
	return(list(p.value=pval))

}

Permu_Close<-function(Permu.Info){

	n.cohort = Permu.Info$n.cohort
	for(i in 1:n.cohort){
		close(Permu.Info$EachInfo.Permu[[i]]$con)
	}
}

