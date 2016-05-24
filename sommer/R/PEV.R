
PEV<-function(PCAs,candidates,Test,ntoselect, npop, nelite, mutprob, niterations, lambda){
  ########### make cross
  
  makeonecross<-function(x1,x2,candidates,mutprob=.5){
    n1<-length(unlist(x1))
    n2<-length(unlist(x2))
    n<-min(c(n1,n2))
    x1x2<-union(unlist(x1),unlist(x2))
    cross<-sample(x1x2,n, replace=FALSE)
    randnum<-runif(1)
    if (randnum<mutprob){
      cross[sample(1:n,1)]<-sample(setdiff(candidates,cross),1)
    }
    return(cross)
  }
  ####### PEV
  PEVMEANAPPROXBYPCA<-function(Train,Test, PCAs, lambda=1e-5){
    PCATrain<-PCAs[rownames(PCAs)%in%Train,]
    PCATest<-PCAs[rownames(PCAs)%in%Test,]
    PEVmean<-mean(diag(PCATest%*%solve(crossprod(PCATrain)+lambda*diag(ncol(PCAs)))%*%t(PCATest)))
    return(PEVmean)
  }
  ######### generate cross
  GenerateCrossesfromElites<-function(Elites, candidates, npop, mutprob){
    newcrosses<-lapply(1:npop, FUN=function(x){
      x1<-Elites[[sample(1:length(Elites),1)]]
      x2<-Elites[[sample(1:length(Elites),1)]]
      return(makeonecross(x1=x1,x2=x2,candidates=candidates,mutprob=mutprob))
    })
    return(newcrosses)
  }
  ########### MAIN FUNCTION
  
	InitPop<-lapply(1:npop, function(x){return(sample(candidates, ntoselect))})
	InitPopFuncValues<-as.numeric(unlist(lapply(InitPop, FUN=function(x){PEVMEANAPPROXBYPCA(Train=x, Test=Test,PCAs=PCAs, lambda=lambda)})))
	orderofInitPop<-order(InitPopFuncValues, decreasing=FALSE)
	ElitePop<-lapply(orderofInitPop[1:nelite], FUN=function(x){return(InitPop[[x]])})
	ElitePopFuncValues<-InitPopFuncValues[orderofInitPop[1:nelite]]
	meanvec<-c()
	for (iters in 1:niterations){
	CurrentPop<-GenerateCrossesfromElites(Elites=ElitePop,candidates=candidates, npop=npop, mutprob=mutprob)	
	
	CurrentPopFuncValues<-as.numeric(unlist(lapply(CurrentPop, FUN=function(x){PEVMEANAPPROXBYPCA(Train=x, Test=Test,PCAs=PCAs, lambda=lambda)})))
	
	orderofCurrentPop<-order(CurrentPopFuncValues, decreasing=FALSE)
	ElitePop<-lapply(orderofCurrentPop[1:nelite], FUN=function(x){return(CurrentPop[[x]])})
	ElitePopFuncValues<-CurrentPopFuncValues[orderofCurrentPop[1:nelite]]
	meanvec<-c(meanvec,min(ElitePopFuncValues))
	##print(sort(ElitePop[[1]]))
	#plot(meanvec)
	}
	return(ElitePop)
}

