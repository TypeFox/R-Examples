SigTree <-
function(myinput,mystat=c("all","fldc","bldc","fldcc","slb"),mymethod="complete",mymetric="euclidean",
	rand.fun=NA,by.block=NA,distrib=c("vanilla","Rparallel"),Ptail=TRUE,tailmethod=c("ML","MOM"),
	njobs=1,seed=NA,Nperm=ifelse(Ptail,1000,1000*nrow(myinput)),metric.args=list(),rand.args=list()){
	if(!(mymetric %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski",
		"pearson","kendall","spearman"))){
                        define.metric<-get(mymetric)
			save(define.metric,file=paste(getwd(),"/define.metric",sep=""))
                        mymetric<-"define.metric"
                }
	indextable<-TreeStat(myinput,mystat,method=mymethod,metric=mymetric,metric.args=metric.args)
	if(!is.na(rand.fun)){
        #if(class(myinput)!="matrix")stop("Inappropriate input data")
	mystatname<-match.arg(mystat, several.ok = TRUE)
        if(any(mystatname=="all"))mystatname<-c("fldc","bldc","fldcc")
	#size<-ceiling((nrow(indextable)+1)/10)*10
	size<-nrow(indextable)+1
	if(rand.fun!="shuffle.column"&rand.fun!="shuffle.block"){
                        define.function<-get(rand.fun)
                        save(define.function,file=paste(getwd(),"/define.function",sep=""))
                        rand.fun<-"define.function"
        }
	if(!is.na(Nperm))nperm<-Nperm
        if(nperm%%200!=0) batches<-c(rep(200,floor(nperm/200)),nperm%%200)
        if(nperm%%200==0) batches<-c(rep(200,floor(nperm/200)))
	distrib<-match.arg(distrib)
	if(distrib!="Rparallel"){
                if(!(mymetric %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski",
                "pearson","kendall","spearman"))){
                        load(paste(getwd(),"/define.metric",sep=""))
                        mymetric<-"define.metric"
                }
                if(rand.fun!="shuffle.column"&rand.fun!="shuffle.block"){
                        load(paste(getwd(),"/define.function",sep=""))
                        rand.fun<-"define.function"
                }
        }else if(distrib=="Rparallel"){
                ncores<-min(njobs,length(batches),detectCores())
                cl<-parallel::makeCluster(getOption("cl.cores",ncores))
		if(!is.na(seed)) clusterSetRNGStream(cl,iseed=seed)
		#parallel::clusterExport(cl=cl,"myinput")
		parallel::clusterEvalQ(cl=cl,expr={data(myinput,envir=environment())})
		if(!(mymetric %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski",
                "pearson","kendall","spearman"))){
                	parallel::clusterEvalQ(cl=cl,expr={
                        	load(paste(getwd(),"/define.metric",sep=""))})
                	mymetric<-"define.metric"
                }
		if(rand.fun!="shuffle.column"&rand.fun!="shuffle.block"){
			#define.function<-get(rand.fun)
			#save(define.function,file=paste(getwd(),"/define.function",sep=""))
                	parallel::clusterEvalQ(cl=cl,expr={
			load(paste(getwd(),"/define.function",sep=""))})
			rand.fun<-"define.function"
		}
	}
	if(Ptail){
                if(nperm>=10000){
			nperm<-1000
			batches<-10
               		profpack<-vector(mode="list",length=batches)
                	for(pn in 1:batches){
                        	#profpack[[pn]]<-vector(mode="list",length=2)
                        	#names(profpack[[pn]])<-c("myinput","nperm")
                        	#profpack[[pn]]$myinput<-myinput
                        	#profpack[[pn]]$nperm<-100
				profpack[[pn]]<-vector(mode="list",length=1)
                                names(profpack[[pn]])<-c("nperm")
                                profpack[[pn]]$nperm<-100
                	}
		}else if(nperm<10000){
			profpack<-vector(mode="list",length=length(batches))
                	for(pn in 1:length(batches)){
                        	#profpack[[pn]]<-vector(mode="list",length=2)
                        	#names(profpack[[pn]])<-c("myinput","nperm")
                        	#profpack[[pn]]$myinput<-myinput
                        	#profpack[[pn]]$nperm<-batches[pn]
				profpack[[pn]]<-vector(mode="list",length=1)
                                names(profpack[[pn]])<-c("nperm")
                                profpack[[pn]]$nperm<-batches[pn]
                	}
		}
		processed<-switch(distrib,
                        vanilla=lapply(X=profpack,FUN=RandTail,myinput=myinput,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,metric.args=metric.args,rand.args=rand.args),
			Rparallel=parallel::parLapply(cl,X=profpack,fun=RandTail,myinput=myinput,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,metric.args=metric.args,rand.args=rand.args))
                for(pn in 1:length(processed)){
                        if(pn==1)mytailcounts<-processed[[1]]
                        if(pn>1){
			for(i in 1:length(mystatname)){
			#mytailcounts[[i]]<-cbind(mytailcounts[[i]],processed[[pn]][[i]])
			mytailcounts[[i]]<-c(mytailcounts[[i]],processed[[pn]][[i]])
			}
			}
                }
		rm(processed)
		if(any(mystatname!="slb")){
		for(i in 1:length(mystatname)){
			mytailcounts[[i]]<-matrix(nrow=nrow(indextable),data=mytailcounts[[i]])
		}
		names(mytailcounts)<-mystatname
		jointcounts<-matrix(nrow=nrow(indextable),ncol=length(mystatname),
			dimnames=list(c(1:nrow(indextable)),mystatname))
		}else if(mystatname=="slb"){
			mytailcounts[[1]]<-matrix(nrow=1,data=mytailcounts[[i]])
			names(mytailcounts)<-"slb"
			jointcounts<-matrix(ncol=1,nrow=1,dimnames=list(1,"slb"))
		}
		tailmethod<-match.arg(tailmethod)
		if(any(mystatname!="slb")){
		for(statname in mystatname){
			mystat<-vector(mode="list",length=nrow(indextable))
			for(pn in 1:nrow(indextable)){
                        	mystat[[pn]]<-vector(mode="list",length=2)
                        	names(mystat[[pn]])<-c("x","y")
                        	mystat[[pn]]$x<-indextable[pn,statname]
                        	mystat[[pn]]$y<-mytailcounts[[statname]][pn,]
                	}
			jointcounts[,statname]<-unlist(switch(distrib,
				vanilla=lapply(X=mystat,FUN=Pvalue,method=tailmethod,Nexcmax=min(150,nperm/4)),
				Rparallel=parLapply(cl,X=mystat,fun=Pvalue,method=tailmethod,Nexcmax=min(150,nperm/4))))
			jointcounts[jointcounts[,statname]==0,statname]<- .Machine$double.eps
		}
		}else if(mystatname=="slb"){
			mystat<-vector(mode="list",length=1)
			mystat[[1]]<-vector(mode="list",length=2)
			names(mystat[[1]])<-c("x","y")
                        mystat[[1]]$x<-indextable[nrow(indextable),"slb"]
                        mystat[[1]]$y<-mytailcounts[[mystatname]]
			jointcounts[,"slb"]<-unlist(switch(distrib,
                                vanilla=lapply(X=mystat,FUN=Pvalue,method=tailmethod,Nexcmax=min(150,nperm/4)),
                                Rparallel=parLapply(cl,X=mystat,fun=Pvalue,method=tailmethod,Nexcmax=min(150,nperm/4))))
                        jointcounts[jointcounts[,"slb"]==0,"slb"]<- .Machine$double.eps
		}
		if(distrib=="Rparallel")stopCluster(cl)
        }
	if(!Ptail){
		profpack<-vector(mode="list",length=length(batches))
       		for(pn in 1:length(batches)){
                	#profpack[[pn]]<-vector(mode="list",length=2)
                	#names(profpack[[pn]])<-c("myinput","nperm")
                	#profpack[[pn]]$myinput<-myinput
                	#profpack[[pn]]$nperm<-batches[pn]
			profpack[[pn]]<-vector(mode="list",length=1)
                        names(profpack[[pn]])<-c("nperm")
                        profpack[[pn]]$nperm<-batches[pn]
                }
        	processed<-switch(distrib,
                        vanilla=lapply(X=profpack,FUN=RandTree,myinput=myinput,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,metric.args=metric.args,rand.args=rand.args),
                        Rparallel=parallel::parLapply(cl,X=profpack,fun=RandTree,myinput=myinput,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,metric.args=metric.args,rand.args=rand.args))
        	for(pn in 1:length(processed)){
                	if(pn==1)jointcounts<-processed[[1]][[1]]
               		if(pn>1)jointcounts<-jointcounts+processed[[pn]][[1]]
        	}
        	jointcounts<-(jointcounts+1)/(nperm+2)
		if(distrib=="Rparallel")stopCluster(cl)
	}
        if(any(mystatname!="slb")){dimnames(jointcounts)[[2]]<-paste("p",dimnames(jointcounts)[[2]],sep="")
	indextable<-cbind(indextable,jointcounts)}
	else if(mystatname=="slb")indextable<-list(indextable,jointcounts)
	if(rand.fun=="define.function"){system(paste("rm ",getwd(),"/define.function",sep=""))}
	if(mymetric=="define.metric"){system(paste("rm ",getwd(),"/define.metric",sep=""))}
	clusterobj<-vector(mode="list",length=3)
	names(clusterobj)<-c("Call","data","indextable")
	clusterobj$Call<-match.call()
	clusterobj$data<-myinput
        clusterobj$indextable<-indextable
	class(clusterobj)<-"best"
        return(clusterobj)
        }
        if(is.na(rand.fun))return(indextable)
}
