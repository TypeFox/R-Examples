"opt.random.par" <-
function(
M,#matrix (network)
k,#number of clusters/groups
n=NULL,#the number of units in each mode (only necessary if mode is larger than 2)
rep,#number of repetitions/different starting partitions to check
approach,
...,
return.all=FALSE,#if 'FALSE', solution for only the best (one or more) partition/s is/are returned
return.err=TRUE,#if 'FALSE', only the resoults of crit.fun are returned (a list of all (best) soulutions including errors), else the resoult is list
maxiter=50,#maximum number of iterations
#m=NULL,#suficient value individual cells
#cut=min(M[M>0]),   #
#BLOCKS=NULL,#array of permissible block types and their ordering for all blocks
trace.iter=FALSE,#save a result of each iteration or only the best (minimal error)
switch.names=NULL,#should partitions that only differ in group names be considert equal (is c(1,1,2)==c(2,2,1))
save.initial.param=TRUE,#should the initial parameters be saved
skip.par=NULL,#the partions that are not allowed or were already checked and should be skiped
save.checked.par=TRUE,#should the checked partitions be saved
merge.save.skip.par=any(!is.null(skip.par),save.checked.par), #should the checked partitions be merged with skiped ones
skip.allready.checked.par=TRUE,#if 'TRUE',the partitions that were already checked when runing 'opt.par' form different statrting points will be skiped
check.skip="iter",#when should the check be preformed:
# "all"  - before every call to 'crit.fun'
# "iter" - at the end of eack iteratiton
# "opt.par"  - before every call to 'opt.par', implemented in opt.these.par and opt.random.par
# "never" - never
#use.for=TRUE, #should fortran rutines be used when possible
print.iter=FALSE, #should the progress of each iteration be printed
max.iden=10, #the maximum number of results that should be saved (in case there are more than max.iden results with minimal error, only the first max.iden will be saved)
seed=NULL,#the seed for random generation of partitions
parGenFun = genRandomPar, #The function that will generate random partitions. It should accept argumetns: k (number of partitions by modes, n (number of units by modes), seed (seed value for random generation of partition), addParam (a list of additional parametres)
mingr=1,	#minimal alowed group size
maxgr=Inf,	#maximal alowed group size
addParam=list(  #list of additional parameters for gerenrating partitions. Here they are specified for dthe default function "genRandomPar"
genPajekPar = TRUE, 	#Should the partitions be generated as in Pajek (the other options is completly random)
probGenMech = NULL),	#Here the probabilities for different mechanizems for specifying the partitions are set. If not set this is determined based on the previous parameter.
maxTriesToFindNewPar=rep*10 	#The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before 
){

  dots<-list(...)
  if(is.null(switch.names)){
  	switch.names<-is.null(dots$BLOCKS)
  }


  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
  optfun<-gen.opt.par(M=M,k=k,maxiter=maxiter, approach=approach,switch.names=switch.names,trace.iter=trace.iter,save.initial.param = save.initial.param,skip.par=skip.par,save.checked.par=save.checked.par,merge.save.skip.par=merge.save.skip.par,check.skip=check.skip,print.iter=print.iter,mingr=mingr,maxgr=maxgr,...)
  eval(optfun)

  nmode<-length(k)

  res<-list(NULL)
  err<-NULL
  nIter<-NULL

  if(nmode==1){
    n<-dim(M)[1]
  } else if(nmode==2){
    n<-dim(M)
  }

  if(!is.null(seed))set.seed(seed)
  
  on.exit({
    res1 <- res[which(err==min(err, na.rm = TRUE))]
    best<-NULL
    best.clu<-NULL
    for(i in 1:length(res1)){
      for(j in 1:length(res1[[i]]$best)){
        if(
          ifelse(is.null(best.clu),
            TRUE,
            if(nmode==1) ifelse(switch.names,
              !any(sapply(best.clu,rand2,clu2=res1[[i]]$best[[j]]$clu)==1),
              !any(sapply(best.clu,function(x)all(x==res1[[i]]$best[[j]]$clu)))
            ) else ifelse(switch.names,
              !any(sapply(best.clu,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(res1[[i]]$best[[j]]$clu))==1),
              !any(sapply(best.clu,function(x)all(unlist(x)==unlist(res1[[i]]$best[[j]]$clu))))
            )
          )
        ){
          best<-c(best,res1[[i]]$best[j])
          best.clu<-c(best.clu,list(res1[[i]]$best[[j]]$clu))
        }
        
        if(length(best)>=max.iden) {
        	warning("Only the first ",max.iden," solutions out of ",length(na.omit(err))," solutions with minimal error will be saved.\n")
        	break
        }
  
      }
    }
  
      names(best)<-paste("best",1:length(best),sep="")
  
    if(any(na.omit(err)==Inf) || ss(na.omit(err))!=0 || length(na.omit(err))==1){
      cat("\n\nOptimization of all partitions completed\n")
      cat(length(best),"solution(s) with minimal error =", min(err,na.rm=TRUE), "found.","\n")
    }else {
      cat("\n\nOptimization of all partitions completed\n")
      cat("All",length(na.omit(err)),"solutions have err",err[1],"\n")
    }
  
    call<-list(call=match.call())
    best<-list(best=best)
    checked.par<-list(checked.par=skip.par)
    if(return.all) res<-list(res=res) else res<-NULL
    if(return.err) err<-list(err=err) else err<-NULL
    if(!exists("initial.param")){
      initial.param<-NULL
    } else initial.param=list(initial.param)
  
    res<-c(list(M=M),res,best,err,list(nIter=nIter),checked.par,call,initial.param=initial.param)
    class(res)<-"opt.more.par"
    return(res)
    })
  

  for(i in 1:rep){
    cat("\n\nStarting optimization of the partiton",i,"of",rep,"partitions.\n")
    find.unique.par<-TRUE
    ununiqueParTested=0
    while(find.unique.par){
      temppar<-parGenFun(n=n,k=k,seed=seed,mingr=mingr,maxgr=maxgr,addParam=addParam)

      find.unique.par<-
      ifelse(is.null(skip.par),
        FALSE,
        if(nmode==1) ifelse(switch.names,
          any(sapply(skip.par,rand2,clu2=temppar)==1),
          any(sapply(skip.par,function(x)all(x==temppar)))
        ) else ifelse(switch.names,
          any(sapply(skip.par,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(temppar))==1),
          any(sapply(skip.par,function(x)all(unlist(x)==unlist(temppar))))
        )
      )
      ununiqueParTested=ununiqueParTested+1
      endFun<-ununiqueParTested>=maxTriesToFindNewPar
      if(endFun) {
      	break
      } else if(ununiqueParTested%%10==0) cat(ununiqueParTested,"partitions tested for unique partition\n")
    }
    
    if(endFun) break
    
    res[[i]]<-opt.par.tmp(
    M=M,
    clu=temppar,
    k=k,
    approach=approach,
    skip.par=skip.par,
    ...
    )

    err[i]<-res[[i]]$best[[1]]$err
    nIter[i]<-res[[i]]$nIter
    if(skip.allready.checked.par && save.checked.par)skip.par<-res[[i]]$checked.par
  }

}

