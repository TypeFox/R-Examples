Rparallel<-function(randfun,distrib,doshuffles,nshuffle,dataIn,returnme,
	boundaries,njobs,qmem){
        if(distrib=="Rparallel"&doshuffles!="NO"){
                ncores<-min(njobs,ifelse(doshuffles=="FROMSCRATCH",
			nshuffle,nshuffle-dataIn$nshuffle))
		cl<-parallel::makeCluster(ncores)
		parallel::clusterEvalQ(cl=cl,expr=library(CORE))
        }
        if(distrib=="Grid"&doshuffles!="NO"){
		ncores<-min(njobs,ifelse(doshuffles=="FROMSCRATCH", 
                         nshuffle,nshuffle-dataIn$nshuffle))
		timeid<-substring(Sys.time(),12,20)
		system(paste("mkdir",timeid))
		setwd(paste(getwd(),"/",timeid,sep=""))
		WrRGrid(ncores)
		nshuffle<-returnme$nshuffle
		save(ncores,doshuffles,dataIn,randfun,boundaries,
			returnme,file=paste(getwd(),"/tempMPI",sep=""))
		sink(paste(getwd(),"/rjob.sh",sep=""))
		cat("R CMD BATCH RGrid.R")
		sink()
		if(is.na(qmem)){
			system(paste("qsub  -cwd -V -sync y -l virtual_free=2G -t 1-",ncores,":1 rjob.sh",sep=""))
			#system(paste("qsub  -cwd -V -sync y -l m_mem_free=2G -t 1-",ncores,":1 rjob.sh",sep=""))
			}
		if(!is.na(qmem)){
			system(paste("qsub  -cwd -V -sync y ",qmem," -t 1-",ncores,":1 rjob.sh",sep=""))
			}
		load(paste(getwd(),"/jobid",sep=""))
		myjobID<-get("job.id")
		load(paste(getwd(),"/mygather.temp.",as.character(myjobID),"1",sep=""))
		myresult<-get("x")
		for(i in 2:ncores){
			load(paste(getwd(),"/mygather.temp.",as.character(myjobID),as.character(i),sep=""))
			myresult<-c(myresult,get("x"))
			}
		if(doshuffles=="FROMSCRATCH")mpidata<-myresult[1:(nshuffle*nrow(returnme$coreTable))]
		if(doshuffles=="ADD")mpidata<-myresult[1:((nshuffle-dataIn$nshuffle)*nrow(returnme$coreTable))]
        }
        if(doshuffles=="FROMSCRATCH"){
                returnme$simscores<-switch(distrib,
                        vanilla=randfun(COREobj=returnme,boundaries=boundaries),
                        Rparallel=matrix(nrow=nrow(returnme$coreTable),
                                data=unlist(parallel::parSapply(cl=cl,X=1:ncores,FUN=randfun,
                                COREobj=returnme,boundaries=boundaries,nprocs=ncores))),
			Grid=matrix(nrow=nrow(returnme$coreTable),data=mpidata)
                        )
        }
        else if(doshuffles=="ADD"){
                returnme$simscores<-
                cbind(dataIn$simscores[1:nrow(returnme$coreTable),,drop=F],switch(distrib,
                        vanilla=randfun(COREobj=returnme,boundaries=boundaries,
                                rngoffset=dataIn$nshuffle),
                        Rparallel=matrix(nrow=nrow(returnme$coreTable),
                                data=unlist(parallel::parSapply(cl=cl,X=1:ncores,
                                FUN=randfun,COREobj=returnme,boundaries=boundaries,
                                rngoffset=dataIn$nshuffle,nprocs=ncores))),
			Grid=matrix(nrow=nrow(returnme$coreTable),data=mpidata)
                ))
        }
        if("simscores"%in%names(returnme))returnme$p<-
                (rowSums(returnme$simscores>returnme$coreTable[,"score"])+1)/
                (ncol(returnme$simscores)+2)
        if(exists("cl"))parallel::stopCluster(cl)
		if(substring(distrib,1,3)=="Gri"){
			returnme$nshuffle<-nshuffle
			setwd("..")
			system(paste("rm -rf",timeid))
			#system(paste("rm mygather.temp.",as.character(myjobID),"*",sep=""))
        		#system("rm tempMPI")
			#system("rm RGrid.*")
			#system("rm jobid")
			#system("rm rjob.sh*")
			}
	return(returnme)
}
