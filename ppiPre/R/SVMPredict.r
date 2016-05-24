SVMPredict<-function(training_set,predict_set,output="falsePPIs-ppiPre.csv",organism="yeast",drop ="IEA", replaceNA=0)
{
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
	dropcodes <- drop
	trainingsetevidences <- paste("AllEvidences-",training_set)
	predictevidences <- paste("AllEvidences-",predict_set)

	ComputeAllEvidences(input=training_set,output=trainingsetevidences, organism=wh_organism, drop = dropcodes)
	svm.model <- SVMTrain(trainingsetevidences, replaceNA = replaceNA)

	ComputeAllEvidences(input=predict_set,output=predictevidences, organism=wh_organism, drop = dropcodes )
	predictset<- read.csv(predictevidences,header=TRUE,sep=",")
	predictset[is.na(predictset)]<-replaceNA
	predictset$label<-factor(predictset$label)
	predictset_noname<-predictset[,c(-1,-2)] 
	svm.pred<-predict(svm.model, predictset_noname[,-1])
	result <- cbind(predictset,svm.pred) #add the prediction result into the data to predict
	falsePPIs <- result[which (result[3] != result[17]),] # the interactions whose prediction result is different from the original value are the false interactions
	write.csv(falsePPIs,file=output,row.names=FALSE) 
}

SVMTrain<-function(input, replaceNA = 0)
{
	trainingset<-read.csv(file=input,header=TRUE,sep=",")
	trainingset$label<-factor(trainingset$label)
	trainingset<-trainingset[,c(-1,-2)] #remove the name of the proteins
	trainingset [is.na(trainingset)]<-replaceNA
	svm.model<-e1071::svm(label~ .,data=trainingset,cost=100,gamma=1) #train the svm classifier
	return(svm.model)
}

ComputeAllEvidences<-function(input,output="AllEvidences-ppiPre.csv",organism="yeast", drop ="IEA", header=TRUE, sep=",") 
{
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
	dropcodes <- drop
	inputfile<-read.csv(file=input,header=header,sep=sep)
	positivenode<-inputfile[which(inputfile[3]==1),] 
	graph<-graph.data.frame(positivenode,directed=FALSE)   
	nodes<-get.vertex.attribute(graph,"name")
	Sims<-data.frame(protein1=inputfile[1],protein2=inputfile[2],label=inputfile[3],BPWang=0,MFWang=0,CCWang=0,BPTCSS=0,MFTCSS=0,CCTCSS=0,BPIG=0,MFIG=0,CCIG=0,KEGGSim=0,jaccard=0,AA=0,RA=0)
	i<-1
	for(i in 1:length(inputfile[[1]]))
	{
		message(paste("Computing evidences",as.character(Sims[[1]][i]),as.character(Sims[[2]][i])))
		Sims[[4]][i]<-geneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="BP", organism=wh_organism, measure="Wang",drop <- dropcodes )$geneSim
		Sims[[5]][i]<-geneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="MF", organism=wh_organism,measure="Wang",drop <- dropcodes )$geneSim
		Sims[[6]][i]<-geneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="CC", organism=wh_organism,measure="Wang",drop <- dropcodes )$geneSim
		Sims[[7]][i]<-TCSSGeneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="BP", organism=wh_organism,drop <- dropcodes )$geneSim
		Sims[[8]][i]<-TCSSGeneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="MF", organism=wh_organism,drop <- dropcodes )$geneSim
		Sims[[9]][i]<-TCSSGeneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="CC", organism=wh_organism,drop <- dropcodes )$geneSim
		Sims[[10]][i]<-IntelliGOGeneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="BP", organism=wh_organism,drop <- dropcodes )$geneSim
		Sims[[11]][i]<-IntelliGOGeneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="MF", organism=wh_organism,drop <- dropcodes )$geneSim
		Sims[[12]][i]<-IntelliGOGeneSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),ont="CC", organism=wh_organism,drop <- dropcodes )$geneSim
		Sims[[13]][i]<-KEGGSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]))	
		if(Sims[[3]][i]==1)
		{
			Sims[[14]][i]<-JaccardSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),graph)
			Sims[[15]][i]<-AASim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),graph)
			Sims[[16]][i]<-RASim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),graph)
		}
		else{		
			include1<-na.omit(match(as.character(Sims[[1]][i]),nodes))
			include2<-na.omit(match(as.character(Sims[[2]][i]),nodes))
			if(length(include1)!=0&&length(include2)!=0)
			{	
				Sims[[14]][i]<-JaccardSim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),graph)
				Sims[[15]][i]<-AASim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),graph)
				Sims[[16]][i]<-RASim(as.character(Sims[[1]][i]),as.character(Sims[[2]][i]),graph)
			}
			else{
				Sims[[14]][i]<-0
				Sims[[15]][i]<-0
				Sims[[16]][i]<-0
		    	}
		}
		i<-i+1
	}
	write.csv(Sims,file=output,row.names=FALSE) 
}
#ComputeAllEvidences(input="BinaryGS-randomNGS.csv",output="AllEvidences-ppiPre.csv",organism="yeast", drop ="IEA", header=TRUE,sep=",") 
