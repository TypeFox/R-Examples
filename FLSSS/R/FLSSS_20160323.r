FLSSS<-function(len,v,target,ME,sizeNeeded=1L,tlimit="none",catch=NULL,throw=F,LB=1L:len,UB=(length(v)-len+1L):length(v)){
if(len==0||length(v)==0||(is.numeric(sizeNeeded)&sizeNeeded<=0)||(is.numeric(tlimit)&tlimit<=0))return(list())
if(is.numeric(sizeNeeded))
{
    if(is.numeric(tlimit))
    {
        if(is.null(catch))
        {
          if(!throw).Call('FLSSS_FLSSS_SK', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB, sizeNeeded, tlimit)
          else .Call('FLSSS_FLSSS_SK_throw', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB, sizeNeeded, tlimit)
        }
        else
        {
          if(!throw).Call('FLSSS_FLSSS_SK_catch', PACKAGE = 'FLSSS', len, v, ME, catch[[1]],catch[[2]],catch[[3]],catch[[4]], sizeNeeded, tlimit)
          else .Call('FLSSS_FLSSS_SK_catch_throw', PACKAGE = 'FLSSS', len, v, ME,catch[[1]],catch[[2]],catch[[3]],catch[[4]],sizeNeeded, tlimit)
        }
    }
    else
    {
        if(is.null(catch))
        {
          if(!throw).Call('FLSSS_FLSSS_SK', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB, sizeNeeded,-1)
          else .Call('FLSSS_FLSSS_SK_throw', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB, sizeNeeded, -1)
        }
        else
        {
          if(!throw).Call('FLSSS_FLSSS_SK_catch', PACKAGE = 'FLSSS', len, v, ME, catch[[1]],catch[[2]],catch[[3]],catch[[4]], sizeNeeded, -1)
          else .Call('FLSSS_FLSSS_SK_catch_throw', PACKAGE = 'FLSSS', len, v, ME, catch[[1]],catch[[2]],catch[[3]],catch[[4]], sizeNeeded, -1)
        }
    }
}
else
{
    if(is.numeric(tlimit))
    {
        if(is.null(catch))
        {
          if(!throw).Call('FLSSS_FLSSS_SK', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB, 0, tlimit)
          else .Call('FLSSS_FLSSS_SK_throw', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB, 0, tlimit)
        }
        else
        {
          if(!throw).Call('FLSSS_FLSSS_SK_catch', PACKAGE = 'FLSSS', len, v, ME, catch[[1]],catch[[2]],catch[[3]],catch[[4]],0,tlimit)
          else .Call('FLSSS_FLSSS_SK_catch_throw', PACKAGE = 'FLSSS', len, v, ME, catch[[1]],catch[[2]],catch[[3]],catch[[4]], 0, tlimit)
        }
    }
    else
    {
        if(is.null(catch))
        {
          if(!throw).Call('FLSSS_FLSSS_SK', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB,0,-1)
          else .Call('FLSSS_FLSSS_SK_throw', PACKAGE = 'FLSSS', len, v, target, ME, LB, UB, 0, -1)
        }
        else
        {
          if(!throw).Call('FLSSS_FLSSS_SK_catch', PACKAGE = 'FLSSS', len, v, ME, catch[[1]],catch[[2]],catch[[3]],catch[[4]], 0, -1)
          else .Call('FLSSS_FLSSS_SK_catch_throw', PACKAGE = 'FLSSS', len, v, ME, catch[[1]],catch[[2]],catch[[3]],catch[[4]],0,-1)
        }
    }
}
}









# mFLSSS0<-function(len, mV, mTarget, mME, sizeNeeded=1L, tlimit="none", LB=1L:len, UB=(nrow(mV)-len+1L):nrow(mV)){
# # diffKey=diff(mV[[1]])
# bestCol=colSums(cor(mV,method="spearman"))
# bestCol=which(bestCol==max(bestCol))[1]
# mV=data.frame(mV,id=1L:nrow(mV))
# mV=mV[order(mV[[bestCol]]),]
# mVid=mV$id
# mV=data.frame(key=1L:nrow(mV),mV[-ncol(mV)])
#
# d=length(mTarget)
# # if(min(diffKey)<=0){cat("First variable (1st column) is not monotonic increasing\n");return()}
# k=unlist(lapply(mV[-1],function(x)max(0,max(-diff(x)/1))))
# names(k)=NULL
# mV[-1]=mapply(function(x,y)x*mV[[1]]+y,k,mV[-1],SIMPLIFY=F)
#
# keyME=min(mME/k)
#
# keyTarget=seq(sum(1L:len)+as.integer(keyME),sum((nrow(mV)-len+1L):nrow(mV))-as.integer(keyME),by=max(as.integer(2*keyME),1L))
#
# mME=-k*rep(as.integer(keyME),d)+mME
#
# if(as.integer(keyME)<1L)keyME=1e-1
# mME=c(keyME,mME)
#
# # print(mV)
# # print(mME)
# #print(mV)
# rst=list()
# i=1L
# for(keyT in sample(keyTarget,length(keyTarget)))
# #for(keyT in keyTarget[order(abs(keyTarget-(keyTarget[1]+keyTarget[length(keyTarget)])/2))])
# {
#   #cat(keyT," ")
#   mTarget_i=k*rep(keyT,d)+mTarget
#   #mME=k*rep(keyME,d)+mME
#
#   # print(c(keyT,mTarget))
#   # print(c(keyME,mME))
#   if(!is.numeric(tlimit))tlimit=-1
#   realTarget=c(keyT,mTarget_i)
#
#   # print(realTarget)
#
#   rst=c(rst,mFLSSS_SK(len, nrow(mV), mV, realTarget, mME, LB, UB, sizeNeeded, tlimit))
#   if(length(rst)>=sizeNeeded){print(i);break};i=i+1L
# }
#
# lapply(rst,function(x)mVid[x])
# }














mFLSSSpar<-function(len, mV, mTarget, mME, maxCore=8L, totalSolutionNeeded=1L, tlimit=60, singleSolutionNeeded=1L, randomizeTargetOrder=T, LB=1L:len, UB=(nrow(mV)-len+1L):nrow(mV)){

if(ncol(mV)==1L){print("Please call FLSSS for single dimensional set");return(list());}

bestCol=colSums(cor(mV,method="spearman"))

if(abs(bestCol[1]-ncol(mV))<1e-9)
{
  mV=data.frame(mV,id=1L:nrow(mV))
  mV=mV[order(mV[[1L]]),]
  mVid=mV$id
  mV$id=NULL
  rst=.Call('FLSSS_mFLSSS_SK_como', PACKAGE = 'FLSSS', len, nrow(mV), mV, mTarget, mME, LB, UB, totalSolutionNeeded, tlimit)
  return(lapply(rst,function(x)mVid[x]))
}



bestCol=which(bestCol==max(bestCol))[1]
mV=data.frame(mV,id=1L:nrow(mV))
mV=mV[order(mV[[bestCol]]),]
mVid=mV$id
mV=data.frame(key=1L:nrow(mV),mV[-ncol(mV)])

d=length(mTarget)
k=unlist(lapply(mV[-1],function(x)max(0,max(-diff(x)/1))))
names(k)=NULL
mV[-1]=mapply(function(x,y)x*mV[[1]]+y,k,mV[-1],SIMPLIFY=F)

keyME=min(mME/k)

keyTarget=seq(sum(1L:len)+as.integer(keyME),sum((nrow(mV)-len+1L):nrow(mV))-as.integer(keyME),by=max(as.integer(2*keyME),1L))
if(randomizeTargetOrder)keyTarget=sample(keyTarget,length(keyTarget))

mME=-k*rep(as.integer(keyME),d)+mME

if(as.integer(keyME)<1L)keyME=1e-1
mME=c(keyME,mME)

if(maxCore<2)rst=.Call('FLSSS_mFLSSS_SK', PACKAGE = 'FLSSS', len, nrow(mV), mV, keyTarget, mTarget, k, mME, LB, UB, totalSolutionNeeded, singleSolutionNeeded, tlimit)
  else rst=.Call('FLSSS_mFLSSS_SK_par', PACKAGE = 'FLSSS', maxCore, len, nrow(mV), mV, keyTarget, mTarget, k, mME, LB, UB, totalSolutionNeeded, singleSolutionNeeded, tlimit)

lapply(rst,function(x)mVid[x])
}








# mFLSSS<-function(len, mV, mTarget, mME, totalSolutionNeeded=1L, tlimit=60, singleSolutionNeeded=1L, LB=1L:len, UB=(nrow(mV)-len+1L):nrow(mV)){
# bestCol=colSums(cor(mV,method="spearman"))
# bestCol=which(bestCol==max(bestCol))[1]
# mV=data.frame(mV,id=1L:nrow(mV))
# mV=mV[order(mV[[bestCol]]),]
# mVid=mV$id
# mV=data.frame(key=1L:nrow(mV),mV[-ncol(mV)])
#
# d=length(mTarget)
# k=unlist(lapply(mV[-1],function(x)max(0,max(-diff(x)/1))))
# names(k)=NULL
# mV[-1]=mapply(function(x,y)x*mV[[1]]+y,k,mV[-1],SIMPLIFY=F)
#
# keyME=min(mME/k)
#
# keyTarget=seq(sum(1L:len)+as.integer(keyME),sum((nrow(mV)-len+1L):nrow(mV))-as.integer(keyME),by=max(as.integer(2*keyME),1L))
# keyTarget=sample(keyTarget,length(keyTarget))
#
# mME=-k*rep(as.integer(keyME),d)+mME
#
# if(as.integer(keyME)<1L)keyME=1e-1
# mME=c(keyME,mME)
#
# rst=mFLSSS_SK(len, nrow(mV), mV, keyTarget,
#   mTarget, k, mME, LB, UB, totalSolutionNeeded, singleSolutionNeeded, tlimit)
#
# lapply(rst,function(x)mVid[x])
# }















