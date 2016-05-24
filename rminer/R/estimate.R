#-------------------------------------------------------------------------------------------------
# "estimate.R" code by Paulo Cortez 2014@, Department of Information Systems, University of Minho
#
# This file deals with two estimation functions: crossvaldata - k-fold cross validation and holdout - holdout validation
#
#-------------------------------------------------------------------------------------------------

# adapted from the bootstrap library to use x - formula and data!!!
# mode= "stratified", "order" or "random"
crossvaldata<- function(x,data,theta.fit,theta.predict,ngroup=10,mode="stratified",seed=NULL,model,task,feature="none",...)
{
    args=do.call(list,list(...)) # extra arguments (needed for different kernel?)
    #call <- match.call() not needed
    outindex<-output_index(x,names(data))
    y<- data[,outindex]
    n <- length(y)
    ngroup <- trunc(ngroup)

    if( ngroup < 2){
        stop ("ngroup should be greater than or equal to 2")
    }
    if(ngroup > n){
        stop ("ngroup should be less than or equal to the number of observations")
    }
    
    if(ngroup==n) {groups <- 1:n; leave.out <- 1} # leave-one-out
    if(ngroup<n){ leave.out <- trunc(n/ngroup); groups <- crossfolds(y,ngroup,leave.out,mode,seed) } # stratified k-fold or random k-fold

    par=vector("list",length=ngroup)
    if(substr(feature[1],1,4)=="sens" || substr(feature[1],1,4)=="sabs" || substr(feature[1],1,4)=="simp") 
    {
         SEN <- matrix(nrow=ngroup, ncol=ncol(data) )
         if(substr(feature[1],1,4)!="sens")
         { SRESP=TRUE;NCOL=ncol(data);FSRESP=TRUE;}
         else
         { SRESP=FALSE;FSRESP=FALSE;}
    }
    else {SEN<-NULL;SRESP=FALSE}

    if(feature[1]!="none" && substr(feature[1],1,4)!="sens") ATTRIB <- vector("list",ngroup)
    else ATTRIB<-NULL

    #cv.fit <- rep(NA,n) # this line was changed to:
    if(task=="prob") cv.fit <- matrix(ncol=length(levels(y)),nrow=n) 
    else if(task=="class") cv.fit=y
    else cv.fit <- rep(NA,n)

    imethod=switch(feature[1],sabsv=,simpv="sensv",sabsg=,sabs=,simp=,simpg="sensg",sabsr=,simpr="sensr",feature[1])

    for(j in 1:ngroup)
    {
        u=theta.fit(x,data[-groups[[j]], ],task=task,model=model,feature=feature,...)
        if(!is.null(SEN)) 
         {
            #cat("----- j:",j,"\n",sep=" ")
            #aux=(Importance(u,data[-groups[[j]], ],method="sensgrad"))$imp 
            #print(paste(" --- j:",j,"---"))
            #print(aux)
            IMPORTANCE=Importance(u,data[-groups[[j]], ],method=imethod,responses=FSRESP)
#cat(" >> L:",length(IMPORTANCE$sresponses),"NCOL:",NCOL,"\n")
            SEN[j,]=IMPORTANCE$imp 
            if(FSRESP)
            { 
             if(j==1) SRESP=IMPORTANCE$sresponses
             else{ for(l in 1:NCOL) 
                    {
                     if(!is.null(IMPORTANCE$sresponses[[l]]) ) # && !is.null(IMPORTANCE$sresponses[[j]])$x) # $x)) 
                      { if(is.null(SRESP[[l]])) SRESP[[l]]=IMPORTANCE$sresponses[[l]]
                        else{ SRESP[[l]]$x=c(SRESP[[l]]$x,IMPORTANCE$sresponses[[l]]$x);
                              if(task=="prob") SRESP[[l]]$y=rbind(SRESP[[l]]$y,IMPORTANCE$sresponses[[l]]$y,deparse.level=0)
                              else if(task=="class") SRESP[[l]]$y=addfactor(SRESP[[l]]$y,IMPORTANCE$sresponses[[l]]$y)
                              else SRESP[[l]]$y=c(SRESP[[l]]$y,IMPORTANCE$sresponses[[l]]$y)
                            }
                      }
                    }
                }
            }
         }
        if(!is.null(ATTRIB)) ATTRIB[[j]]=u@attributes
        par[[j]]=u@mpar
#print("---")
#cat("cv fit class:",class(cv.fit),"\n")
        if(model=="cubist" && !is.null(args$neighbors)) 
        {
        neighbors=args$neighbors
        if(is.matrix(cv.fit)) cv.fit[ groups[[j]], ] <-  theta.predict(u,data[groups[[j]],],neighbors) # probabilities!
        else cv.fit[ groups[[j]] ] <-  theta.predict(u,data[groups[[j]],],neighbors) # regression or classification, 1 output
        }
        else
        {
        if(is.matrix(cv.fit)) cv.fit[ groups[[j]], ] <-  theta.predict(u,data[groups[[j]],]) # probabilities!
        else cv.fit[ groups[[j]] ] <-  theta.predict(u,data[groups[[j]],]) # regression or classification, 1 output
        }
        #print(cv.fit[groups[[j]]])
    }
    if(leave.out==1) groups <- NULL
    return(list(cv.fit=cv.fit, 
                mpar=par,
                sen=SEN,sresponses=SRESP,
                attributes=ATTRIB,
                ngroup=ngroup, 
                leave.out=leave.out,
                groups=groups))
                #call=call)) 
}

#  auxiliar function, adapted from the bootstap library: only makes the groups, you should not need to use this:
# about mode:
#      > "stratified" - stratified random split if y is factor; else random split
#      > "random" - uses standard random split
crossfolds<-function(y,ngroup,leave.out,mode="stratified",seed=NULL)
{
  n <- length(y)
  groups <- vector("list",ngroup)
  if(!is.null(seed)) set.seed(seed) # use seed to fix groups

  if(is.factor(y) && mode=="stratified") STRATIFIED=TRUE# stratified crossfolds
  else STRATIFIED=FALSE

  if(STRATIFIED) # stratified crossfolds
    groups=factorsample(y,leave.out,ngroup) # groups can be null if stratified is not possible

  if(!STRATIFIED || is.null(groups) ) # pure random crossfolds
  {
   if(mode=="order") o<- 1:n 
   else o <- sample(1:n) # random split

   for(j in 1:(ngroup-1)){
            jj <- (1+(j-1)*leave.out)
            groups[[j]] <- (o[jj:(jj+leave.out-1)])
   }
   groups[[ngroup]] <- o[(1+(ngroup-1)*leave.out):n]
  }
  if(!is.null(seed)) set.seed(NULL) # reset seed
  return(groups)
}


#---------------------------------------------------------
# holdout: create indexes for spliting the data into training and test datasets
#          the holdout is statified if the output is a factor 
# a list is returned with:
#  $tr - indexes of all training examples
#  $ts - indexes of all test examples
#  if internalsplit is TRUE then another holdout if performed on tr:
#  $itr - indexes of tr for internal training 
#  $vtr - indexes of tr for internal testing (validation)
#
# 
# parameters:
# y - is the desired output vector or a vector with a numeric sequence (e.g. c(1,1,2,3,3,3,4,4) )
# ratio is the ratio of training set (in percentage)
# internalsplit if TRUE then another stratified holdout is used within the training set
# mode - the sampling mode used for the holdout
#      > "stratified" - stratified random split if y is factor; else random split
#      > "random" - uses standard random split
#      > "order" - uses the sequential order of y, no random is used. 
#                  the first examples are used as tr while the lattest are used as ts
#      need to check in the future is this makes any sense at all ?
#      > "incremental" - incremental retraining, ratio=batch-size, iter=iterator
# seed : optional, for having the same holdout (note: to reset a seed, use:  set.seed(NULL) )
#---------------------------------------------------------
holdout<-function(y,ratio=2/3,internalsplit=FALSE,mode="stratified",iter=1,seed=NULL, window=10, increment=1)
{ 
  ALLITR=NULL; VAL=NULL;
  NSIZE=length(y)
 
 if(mode=="incremental"||mode=="rolling")
 { 
   aux=window+increment*(iter-1)
   aux=min(aux,NSIZE)
   if(mode=="rolling") iaux=max((aux-window+1),1) else iaux=1
   ALLTR=iaux:aux
   end=aux+ratio
   end=min(end,NSIZE)
   iend=aux+1
   if(iend<end) TS=iend:end else TS=NULL
 }
 else 
 {
  if(ratio>=1) { ratio=(NSIZE-ratio)/NSIZE}

  if(mode=="order")
  { 
    TRS<-round(ratio*NSIZE)
    ALLTR<-1:TRS
    TS<-(TRS+1):NSIZE
    if(internalsplit)
      {
       TRIS<-round(ratio*TRS)
       ALLITR<-1:TRIS
       VAL<-(TRIS+1):TRS
      }
  }
  else # default random holdout
  {
   if(!is.null(seed)) set.seed(seed) # use seed to fix holdout
   ALL=1:NSIZE
    if(is.factor(y) && mode=="stratified") 
     { 
       ALLTR=factorsample(y,ratio*NSIZE)
       NALLTR=length(ALLTR)
       if(internalsplit==TRUE) ALLITR=ALLTR[factorsample(y[ALLTR],ratio*NALLTR)]
     }
    else
     {
       ALLTR<-sample(ALL,(ratio*NSIZE)) 
       if(internalsplit)
         {
          NALLTR<-length(ALLTR)
          ALLITR<-sample(ALLTR,(ratio*NALLTR)) 
         }
     }
   TS<-setdiff(ALL,ALLTR)
   if(internalsplit==TRUE)VAL<-setdiff(ALLTR,ALLITR)       
   if(!is.null(seed)) set.seed(NULL) # reset seed to random
  }
 }
  return(list(tr=ALLTR,itr=ALLITR,val=VAL,ts=TS))
}
#---------------------

# returns indexes of y with size such that frequency of each class is similar (as possible) to the frequency of y 
factorsample=function(y,size,ngroup=1)
{
 L=levels(y)
 NL=length(L)
 YL=length(y)
 I=vector("list",NL)
 if(ngroup==1) S=vector(length=size) 
 else { S=vector("list",ngroup) 
        for(j in 1:(ngroup-1)) S[[j]]=vector(length=size)
        S[[ngroup]]=vector(length=(YL-size*(ngroup-1)))
        rini=1
      }
 ini=1
 for(i in 1:NL)
  {
   I[[i]]=which(y==L[i])
   LS=round( size*length(I[[i]])/YL )
   if(ngroup==1) { S[ini:(ini+LS-1)]=sample( (I[[i]]),size=LS )
                   ini=ini+LS
                }
   else {
         if(length(I[[i]])<LS*ngroup) {S=NULL;break;}
         OS=sample(I[[i]],size=LS*(ngroup-1))
         sini=1
         for(j in 1:(ngroup-1))
            {
             send=sini+(LS-1)
             S[[j]][ini:(ini+LS-1)]=OS[sini:send]
             sini=send+1
            }
         Rest=setdiff(I[[i]],OS)
         LRest=length(Rest)
         S[[ngroup]][rini:(rini+LRest-1)]=Rest
         ini=ini+LS
         rini=rini+LRest
        }
  }
 return(S)
}
#--------------------------------------
