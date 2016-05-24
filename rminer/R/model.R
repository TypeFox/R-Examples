#-------------------------------------------------------------------------------------------------
# "model.R" code by Paulo Cortez 2010-2014@, Department of Information Systems, University of Minho
#
# This file deals will all Data Mining models
#
#  - Several parts of code were gratefully copied/inspired from the kernlab source code package :)
#-------------------------------------------------------------------------------------------------

# personal thoughts of things to improve:
# improved feature selection?
# parallel execution of searches?
# ensembles? 
# use of references: R2.13? R.oo package ? 
# ranking ?
# ordinal classification ?
#
#-------------------------------------------------------------------------------------------------
# libraries that need to be installed:
#library(nnet)   # get nnet: MLP and Multiple Logistic Regression
#library(pls) # get pls methods
#library(MASS) # lda and qda
#library(mda) # mars
#library(rpart)  # get rpart: Decision Trees
#library(randomForest) # get randomForest
#library(adabag) # get bagging and boosting
#library(party) # get ctree
#library(Cubist) # regression, cubist
#library(kknn)   # get kknn: k-nearest neighbours
#library(kernlab)# get the ksvm
#library(e1071) # get naiveBayes

# not used:
##library(neuralnet) # get neuralnet # still experimental stuff
##library(klaR) # get NaiveBayes
##library(RWeka) # get Weka classifiers ?
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
# save/load functions:
#
#-- save a fitted model into a binary (ascii=FALSE) of text (ascii=TRUE) file:
savemodel=function(mm_model,file,ascii=FALSE)
{ #if(mmmodel@model=="wnaivebayes") { .jcache( (mmmodel@object)$classifier ) } # weka objects
  save(mm_model,file=file,ascii=ascii)
}
#-- load a file that was saved using savemodel function
loadmodel=function(file)
{mm_model=FALSE;load(file=file); return(mm_model) }

#-- save a fitted mining into a binary (ascii=FALSE) of text (ascii=TRUE) file:
savemining=function(mmm_mining,file,ascii=TRUE)
{ save(mmm_mining,file=file,ascii=ascii) }

#-- load a fitted mining file that was saved using savemining function
loadmining=function(file)
{mmm_mining=FALSE;load(file=file);return(mmm_mining)}

#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
# set the global object "model", contains the slots:
# @formula - the formula used to create the model
# @outindex - the output index of data
# @attributes - index with all attributes considered 
# @model - the data mining model:
# @object - the data mining object (e.g. nnet, rpart, etc...) 
# @time - the estimation/training time
# @created - string with the date and time (strptime) when the model was build
# @mpar - vector with the hyperparameters of the model (e.g. svm, mlp, knn)
#         other: metric
#         mlp: c(nr,MaxIt,validationmethod,validationpar,metric) 
#         mlp: c(nr,MaxIt,validationmethod,validationpar,metric,PAR) PAR = H if decay is used for the search or DECAY if H is used for the search
#         svm: c(C,Epsilon,validationmethod,validationpar,metric)  if C=NA or Epsilon=NA then heuristics are used
#         knn: c(validationmethod,validationpar,metric)
#              metric = "AUC", "ACC", "MAE", "SSE"
# @task - type of data mining task, one of the following: 
# "default", "class" (or "cla" or "classification"), "prob" (or "probability") - pure class, no probabilities!, "reg" (or "regression")
#     "default" means if is.factor(y) "class" else "reg"
# @scale - if the data needs to be scaled (to use specially with "nnet"). The options are:
#     "none", "inputs", "all", "default"
# @transform - if the output needs to be transformed ("reg"):
#     "log","logpositive","positive","scale"
# @levels - if the output is factor, store the levels. This is needed to recover the level order when predictive results are latter
#           analysed from files.
setClass("model",representation(formula="formula",model="ANY",task="character",mpar="list",attributes="numeric",scale="character",transform="character",created="character",time="numeric",object="ANY",outindex="numeric",levels="character")) 
#--------------------------------------------------------------------------------------
#--- generic and powerful fit method --------------------------------------------------
# most of these parameters are explained in the rminer.R file

# improve this??? 
# todo : search= change to set the parameter and if two or more parameters can be searched and with which function/method?
# todo : specify different kernel for SVM?
# todo: add more models, such as: 
# Relevance Vector Machine (RVM): rvm
# elastic nets: http://cran.r-project.org/web/packages/glmnet/index.html classification and regression
# Conditional inference trees via party
# The party package provides nonparametric regression trees for nominal, ordinal, numeric, censored, and multivariate responses. party: A laboratory for recursive partitioning, provides details: http://www.statmethods.net/advstats/cart.html ctree

# adabag!
# adaboost: http://en.wikibooks.org/wiki/Data_Mining_Algorithms_In_R/Classification/adaboost -> binary classification!
# number of trees, max depth, min split, complexity,xval:
# predict(gdis,iris,type="probs") and predict(gdis,iris,type="vector")
# variable importance plot!
# require(nnet, quietly=TRUE)

# package car ? lm for classification?
# JRIP class: http://en.wikibooks.org/wiki/Data_Mining_Algorithms_In_R/Classification/JRip library(caret)

#setGeneric("fit", function(x, ...) standardGeneric("fit"))
#setMethod("fit",signature(x="formula"),
# change search method to perform: grid, uniform, other function !!!

fit=function(x, data=NULL,model="default",task="default",search="heuristic",mpar=NULL,feature="none",scale="default",transform="none",created=NULL,fdebug=FALSE,...)
{
PTM= proc.time() # start clock
call=match.call() # does not occupy too much memory!
args=do.call(list,list(...)) # extra arguments (needed for different kernel?)
if(is.null(data)) # If x contains the data in formula expressions, then data=NULL 
{ data=model.frame(x) # only formula is used 
     outindex=output_index(x,names(data))
     x=as.formula(paste(names(data)[outindex]," ~ .",sep=""))
}
else outindex=output_index(x,names(data)) 

params=NULL
task=defaultask(task,model,data[1,outindex])
COL=NCOL(data)
attributes=1:COL # set to all
if(is.factor(data[,outindex])) levels=levels(data[1,outindex]) else levels=""

if(!is.list(model)) model=defaultmodel(model,task) 
#if(is.null(mpar)) mpar=defaultmpar(mpar,task)

if(task!="reg" && !is.factor(data[,outindex])) data[,outindex]=factor(data[,outindex])
else if(task=="reg" && transform_needed(transform)) data[,outindex]=xtransform(data[,outindex],transform) #transform output?

feature=defaultfeature(feature) # process feature

search=readsearch(search,model,task,mpar,...)

if(feature_needed(feature[1])) # will treat this better in next big (improved) version: 2/9/2014
  { 
#print(feature)
#print("--")
    smeasure=switch(feature[1],sabsv="v",sabsr="r","g")
    feature[1]=switch(feature[1],sabsv=,sabsr=,sabsg="sabs",feature[1])
    #cat("f:",feature[1],"sm:",smeasure,"\n")
    fstop=as.numeric(feature[2]) # -1 or number of delections
    Runs=as.numeric(feature[3]) 
    vmethod=feature[4] 
    vpar=as.numeric(feature[5]) 
    if(length(feature)>5) { fsearch=suppressWarnings(as.numeric(feature[6])) # primary search parameter
                            if(is.na(fsearch)) fsearch=feature[6]
                          }
    else fsearch=search # this may take a long time...
       # perform both FS and parameter search at the same time...
    RES=bssearch(x,data,algorithm=feature[1],Runs=Runs,method=c(vmethod,vpar),model=model,task=task,search=fsearch,
                 scale=scale,transform="none",fstop=fstop,smeasure=smeasure)
    attributes=RES$attributes
    outindex=output_index(x,names(data)[attributes]) # update the output index if it has changed due to the deletion of variables
    #if(length(fsearch)>1) { search$search=RES$search # check later if this works: 1/10/2014
    #                      }
    #print(search)
    data=data[,attributes] # new
  }
 else attributes=1:NCOL(data); # set to all

if(search$smethod=="grid" || search$smethod=="matrix") # grid or matrix search
{ 
  search=mgrid(search,x,data,model,task,scale,transform,fdebug,args)
}
else if(search$smethod=="2L" || substr(search$smethod,1,2)=="UD") # nested 2-level grid or uniform design
{ 
 if(fdebug) print(" 1st level:")
#print(search)
  if(substr(search$smethod,1,2)=="UD") 
    {
     if(search$smethod=="UD1") points=13 else points=9
     limits=c(search$search$sigma,search$search$C); if(task=="reg") limits=c(limits,search$search$epsilon)
     mud=uniform_design(limits,points=points,center=NULL)
     search$search$sigma=mud[,1];search$search$C=mud[,2]; if(task=="reg") search$search$epsilon=mud[,3]
    }
  saux1=mgrid(search,x,data,model,task,scale,transform,fdebug,args) # 1st level
#print(" > best1 ---:")
#print(saux1)
  if(substr(search$smethod,1,2)=="UD") 
    {
     if(search$smethod=="UD1") points=9 else points=5
     mode="range"
    } 
  else mode="seq"
  search2=midrangesearch(saux1,search,mode=mode) 
  if(substr(search$smethod,1,2)=="UD") 
    {
     limits=c(search2$search$sigma,search2$search$C); if(task=="reg") limits=c(limits,search2$search$epsilon)
     center=c(saux1$search$sigma,saux1$search$C); if(task=="reg") center=c(center,saux1$search$epsilon)
     mud=uniform_design(limits,points=points,center=center)
     search2$search$sigma=mud[,1];search2$search$C=mud[,2]; if(task=="reg") search2$search$epsilon=mud[,3]
    } 
 if(fdebug) print(" 2nd level:")
#print(" > 2nd ---:")
#print(search2)
  saux2=mgrid(search2,x,data,model,task,scale,transform,fdebug,args) # 2nd level
#print(" > best2 ---:")
#print(saux2)
  if(isbest(saux2$bestval,saux1$bestval,search$metric)) search=saux2 else search=saux1
#print(" > final ---:")
#print(search)
}

if(search$smethod=="none") # no search
{ # --- if for all models, assumes one fit per model:
  #args=c(args,search$search) # need to add extra parameters to args? yes but no duplicates should be in args!
  #cat("search2:",search$smethod,"\n")
#if(is.list(search$search) && is.vector(search$search)) # vector list
#    search$search=search$search[[1]] # XXX 
  args=addSearch(args,search$search,model)
  #AAA<<-args
  #print(">> none fit >> args in fit main function:"); print(class(search)); print(args)
  if(is.list(model)) # new way, list(fit=FITFUNCTION,predict=PREDICTFUNCTION,name=NAME)
    M=try( do.call(model$fit,c(list(x,data),args)), silent=TRUE ) # execute model$fit with x, data and extra args (if any!)
  else # model is character
   {
    fargs=args
    # models with different (non do.call) fits:
    if (model=="naive") # naive classification or regression
      { if(substr(task,1,3)=="reg") M=mean(data[,outindex]) # output mean
        else if(is.factor(data[,outindex])) 
         { M=data[1,outindex];M[1]=levels(data[1,outindex])[mostcommon(data[,outindex])] } # most common class
      } 
    else if(model=="kknn") # knn does not have a fit function but the data needs to be stored
      { M=c(list(formula=x,train=data),fargs) } # think about params!
    # -outindex regression models that use x, y ...
    else if(model=="mars"||model=="cubist") # M=mda::mars(data[,-outindex],data[,outindex],...) 
         M=try( do.call(model,c(list(x=data[,-outindex],y=data[,outindex]),fargs)), silent=TRUE) # FIT!!!
    else # equal do.call fits:
    {
     fitfun=model # true for most models
     if(model=="rpart") 
      { if(task=="reg") fargs[["method"]]="anova" else fargs[["method"]]="class" }
     # "ctree", "lda", "qda", "naiveBayes", "bagging", "boosting", "pcr", "plsr", "cppls", "rvm" (sigma but not C or epsilon)
     else if(model=="multinom") fargs[["trace"]]=FALSE 
     # regression fits:
     else if(model=="mr") # multiple/linear regression: uses the nnet function
      { fitfun="mlp.fit"; fargs[["size"]]=0;fargs[["nr"]]=1;fargs[["task"]]="reg";fargs[["scale"]]="none"; }
     else if(model=="randomForest")
      { fitfun="randomForest";fargs[["importance"]]=TRUE;}
     else if(substr(model,1,3)=="mlp")
      {
       if(scale=="default") # removed "svm" from this if due to a "object scal not found" error
         { if(is.factor(data[,outindex])) scale="inputs" else scale="all" }
       if(model=="mlpe") type=2 else type=1
       fitfun="mlp.fit"; 
       fargs[["task"]]=task;
       fargs[["scale"]]=scale;
       fargs[["type"]]=type
      }
     else if(model=="ksvm")
      { fitfun="svm.fit"
        fargs[["task"]]=task;
      }
     #sink("/dev/null")
     #cat(">>fargs:\n")
     #print(fargs)# capture to avoid "maximum number of iterations reached" verbose:

 capture.output({
     M=try( do.call(fitfun,c(list(x,data=data),fargs)), silent=TRUE) 
                }) # do not show any text??? 
     # FIT!!!
     #sink()
     if(class(M)[1]!="try-error") args=modelargs(args,model,M) # update args if needed

    } # end else
   } # end model is character 
} # end no search
if(class(M)[1]=="try-error") { M=M[[1]];msg=paste("fit failed:",M[[1]],capture.output(call));warning(msg);}

TME=(proc.time()-PTM)[3] # computes time elapsed, in seconds
if(is.null(created)) created=format(Sys.time(),"%Y-%m-%d %H:%M:%S")
return(new("model",formula=x,model=model,task=task,mpar=args,attributes=attributes,scale=scale,transform=transform,created=created,time=TME,object=M,outindex=outindex,levels=levels))
}#)

#
addSearch=function(args,search,model)
{
 # decode and treat kpar:
 ##call=match.call() # does not occupy too much memory!
 ##print(call)
 kpar=list()
 if(is.character(model) && (model=="ksvm"||model=="rvm") )
   {
    N=names(search)
    hyper=c("sigma","degree","scale","offset","order","length","lambda","normalized")
    I=NULL;for(i in 1:length(hyper)) I=c(I,which(N==hyper[i]))
    if(length(I)>0)
     {
      for(i in 1:length(I))
        { aname=N[I[i]]
          kpar[[ aname ]]= search [[ aname ]]
          search [[ aname ]] = NULL
        }
      if( is.null(args$kpar)) args$kpar=kpar
      else for(j in 1:length(kpar)) args$kpar[[ names(kpar)[j] ]] = kpar [[j]]
     }
    
   }
 # add search into args:
 L=length(search)
 N=names(search)
 if(L>0) for(i in 1:L) args [[ N[i] ]] = search [[ N[i] ]]
 return(args)
}

modelargs=function(args,model,M) # not sure if needed
{
 if(substr(model,1,3)=="mlp")
  {
   if(!is.null(args$size) && is.na(args$size)) # "NA" -> size
      {
       if(model=="mlp") size=M$mlp$n[2] else size=M$mlp[[1]]$n[2]
       args[["size"]]=size
      }
   if( is.null(args$nr) || (!is.null(args$nr) && is.na(args$nr)) ) # "NA" -> size
      {
       args[["nr"]]=M$nr
      }
  }
 else if(model=="ksvm"||model=="rvm")
  { if(!is.null(args$kpar)) # automatic estimation was used, or exponential scale or normal
      {
       # get kernel:
       #if(model=="ksvm") args[["kpar"]]=M$svm@kernelf@kpar else 
       args[["kpar"]]=M@kernelf@kpar
      }
    if( model=="ksvm" && ( is.null(args$C) || (!is.null(args$C) && is.na(args$C)) ))
       if(!is.null(M@param$C)) args[["C"]]=M@param$C
    if( model=="ksvm" && ( is.null(args$epsilon) || (!is.null(args$epsilon) && is.na(args$epsilon)) ))
       if(!is.null(M@param$epsilon)) args[["epsilon"]]=M@param$epsilon
   }
 else if(model=="randomForest")
  { if(is.null(args$mtry) || is.na(args$mtry)) # automatic estimation was used!
      {
       args[["mtry"]]=M$mtry
      }
  }
 return(args)
}
# 
readsearch=function(search,model,task="reg",mpar=NULL,...) # COL is inputs+output
{
 #S<<-search;M<<-model;TT<<-task;MP<<-mpar
 #mpause()
 #search=S;model=M;task=TT;mpar=MP
 #search="heuristic"; model="ksvm"; task="prob";mpar=NULL;args=list()
 args=do.call(list,list(...)) # extra arguments (needed for different kernel?)
 if(is.list(search) && !is.null(search$smethod) && search$smethod=="normal") search$smethod="grid" # compatibility issue
 if(is.list(search) && !is.null(search$search) && is.character(search$search))
  { search2=readsearch(search$search,model,task,mpar,...)
    search$search=search2$search;search$smethod=search2$smethod
    if(is.null(search$convex)) search$convex=0
    if(is.null(search$method)) search$method=c("holdout",2/3,NA)
    if(is.null(search$metric)) { if(task=="reg") search$metric="SAE" else if(task=="class") search$metric="ACC" else search$metric="AUC" }
  }
 search2=list();smethod="";
 if(is.character(search) && !is.list(model))
 { 
   if(model=="ksvm" && is.null(args$kernel)) # treat the special rbfdot implicit case
      search2[["kernel"]]="rbfdot"
   #if(search=="automatic" && is.null(args$kpar)) { search2$kpar="automatic"; search="heuristic"} # compatibility issue

   if(search=="heuristic"||search=="automatic") 
   { smethod="none"
     if(substr(model,1,3)=="mlp") 
      { if(!is.null(args$size)) search=args$size else search=NA # will use: round((COL-1)/2) in mlp.fit # simple heuristic for hidden nodes
      }
     else if(model=="ksvm") # special case 
      { 
        if( (is.null(args$kernel) || (is.character(args$kernel) && args$kernel=="rbfdot") ) # compatible kernel
            && is.null(args$kpar) && search!="automatic" ) # no kpar definition or no "automatic"
             { if(!is.na(search)) search="automatic" else search=2^-7 # 2^-7 # or "automatic" # use automatic heuristic for the kernel parameter
             }
        else if(is.null(args$kpar)) search="automatic" # other kernels
      }
   }
   # Manuel JMLR paper:
   # args$kernel=="polydot", s = {0.001,0.01,0.1}, offset o = 1, degree d = {1,2,3}, C={0.25,0.5,1}  
   # "rpart", control$cp with 10 values from 0.18 to 0.01: seq(0.01,0.18,length.out=10)
   # "ctree", mincriterion, tuned with the values 0.1:0.11:0.99, control$mincriterion=seq(0.1,0.98,by=0.11)
   else if(search=="heuristic5") # simple heuristic
     { smethod="grid"
       if(model=="kknn") search=seq(1,9,2)
       else if(substr(model,1,3)=="mlp") search=seq(0,8,2)
       else if(model=="ksvm")
         { if(is.null(args$kernel) || is.character(args$kernel) && args$kernel=="rbfdot") search=2^seq(-15,3,4)
           else {cat("Current rminer version does not have an heuristic5 rule for this kernel:");print(args$kernel);}
         }
       else if(model=="randomForest") {search=1:5} 
     }
   else if(search=="heuristic10")
     { smethod="grid"
       if(model=="kknn") search=seq(1,10,1) 
       else if(substr(model,1,3)=="mlp") search=seq(0,9,1) 
       else if(model=="ksvm")
         { if(is.null(args$kernel) || is.character(args$kernel) && args$kernel=="rbfdot") search=2^seq(-15,3,2)
           else {cat("Current rminer version does not have an heuristic10 rule for this kernel:");print(args$kernel);}
         }
       else if(model=="randomForest") {search=1:10}
     }
   else if(model=="ksvm")
     {
      if(substr(search,1,2)=="UD")
        {
         if(is.null(args$kernel) || is.character(args$kernel) && args$kernel=="rbfdot") 
         { 
          smethod=search
          if(task=="reg") { search2$sigma=c(-8,0);search2$C=c(-1,6);search2$epsilon=c(-8,-1) }
          else { search2$sigma=c(-15,3);search2$C=c(-5,15) } # gama, C, epsilon
          search2[["exponentialscale"]]=TRUE
         }
         else {cat("Current rminer version does not have an UD rule for this kernel:");print(args$kernel);}
        }
     }
 }
 if( (is.na(search)||is.character(search)||is.numeric(search)||is.factor(search)) && !is.list(model))
 {
   if(substr(model,1,3)=="mlp") search2[["size"]]=search
   else if(model=="ksvm") 
     { if(is.numeric(search) && (is.null(args$kernel) || (is.character(args$kernel) && args$kernel=="rbfdot")) ) search2$sigma=search
       else if(search=="automatic") search2$kpar=search
       if(is.null(args$kernel)) # treat the special rbfdot implicit case
         search2[["kernel"]]="rbfdot"
     }
   else if(model=="randomForest" && !is.character(search)) search2[["mtry"]]=search
   else if(model=="kknn" && !is.character(search)) search2[["k"]]=search
   if(smethod!="grid" && smethod!="none" && substr(smethod,1,2)!="UD") { if(length(search)<2) smethod="none" else smethod="grid" }
   if(length(search2)==0) search=list(smethod=smethod) else search=list(smethod=smethod,search=search2,convex=0) ###
 }
 else if(is.character(search) && search=="heuristic" && is.list(model)) 
 { # model is a list!
   search=list(smethod="none")
 }
 if(!is.null(mpar))
 {
  Lmpar=length(mpar)
  if(is.null(search$method) && !is.list(model) ) search$method=getmethod(mpar)# fill method
  if(is.null(search$metric) && !is.list(model) ) search$metric=getmetric(mpar,task) # fill metric
  if( !is.list(model) && substr(model,1,3)=="mlp")
   {
     if(is.null(args$nr) && is.null(search$search$nr) && Lmpar>3) 
       { nr=as.numeric(mpar[1])
         if(is.na(nr)) nr=3 # default value 
         search$search$nr=nr
       }
     if(is.null(args$maxit) && is.null(search$search$maxit) && Lmpar>3)
       { maxit=as.numeric(mpar[2])
         if(is.na(maxit)) maxit=100 # default value 
         search$search$maxit=maxit
       }
   }
  if( !is.list(model) && model=="ksvm")
   {
     if(is.null(args$C) && is.null(search$search$C) && Lmpar>0)
       { C=suppressWarnings(as.numeric(mpar[1])) # number or NA
         search$search$C=C
       }
     if(is.null(args$epsilon) && is.null(search$search$epsilon) && Lmpar>1)
       { epsilon=suppressWarnings(as.numeric(mpar[2])) # number or NA
         search$search$epsilon=epsilon
       }
   }
 }
 else 
   { # fill some search components if missing: 
     if(is.list(search))
       { 
        if(is.null(search$smethod) && !is.null(search$search)) 
          { N=length(search$search)
            L=vector(length=N)
            for(i in 1:N) L[i]=length(search$search[[i]])
            LS=prod(L) 
            if(LS>1) search$smethod="grid" else search$smethod="none" 
          }
        if(is.null(search$search)) search$smethod="none"
        if(!is.null(search$smethod) && search$smethod!="none")
          {
           if(is.null(search$convex)) search$convex=0
           if(is.null(search$method)) search$method=c("holdout",2/3,NA)
          }
        if(is.null(search$metric)) { if(task=="reg") search$metric="SAE" else if(task=="class") search$metric="ACC" else search$metric="AUC" }
       }
   }

 return (search)
}

midrangesearch=function(midpoint,search,mode="seq") 
{
 n=length(search$search) # dimension
 for(i in 1:n)
  if(is.numeric(search$search[[i]]) && length(search$search[[i]])>1 ) 
   { Range=diff( range(search$search[[i]]) )/4
     L=length(search$search[[i]])
     Min=max( midpoint$search[[i]]-Range , min(search$search[[i]]) )
     Max=min( midpoint$search[[i]]+Range , max(search$search[[i]]) )
     if(mode=="range") L=2
     search$search[[i]]=seq(Min,Max,length.out=L)
   }
 return(search) 
}

mparheuristic=function(model,n=NA,lower=NA,upper=NA,by=NA,kernel="rbfdot")
{
 if(model=="ksvm" || model=="rvm") model=paste(model,kernel,sep="")

 if(is.na(lower))
 { lower=switch(model,kknn=1,randomForest=1,mlp=,mlpe=0,rvmrbfdot=,ksvmrbfdot=-15,ksvmvanilladot=-2,rvmpolydot=,ksvmpolydot=1,
                      rpart=0.01,ctree=0.1,NA)
 }
 if(is.na(upper))
 { upper=switch(model,kknn=10,randomForest=10,mlp=,mlpe=9,rvmrbfdot=,ksvmrbfdot=3,ksvmvanilladot=7,rvmpolydot=,ksvmpolydot=3,
                      rpart=0.18,ctree=0.98,NA)
 } 
 #if(is.na(by) && is.na(n) && !is.na(upper) && !is.na(lower))
 #{ by=switch(model,kknn=1,randomForest=1,mlp=,mlpe=1,ksvmrbfdot=2,ksvmvanilladot=1,ksvmpolydot=1,
 #                     rpart=0.017,ctree=0.11,NA)
 #} 

 # not na(n) 
 if(!is.na(n) && n>1) by=(upper-lower)/(n-1)

 if(is.na(n) && is.na(by)) n=1

 if(!is.na(by)) 
  { s=seq(lower,upper,by=by) 
    if(model=="kknn"||model=="mlp"||model=="mlpe"||model=="randomForest") s=round(s) # integer for these models
    s=unique(s) # remove duplicates if any
    n=length(s) 
  }

 if(n>1)
  {
   if(model=="kknn"||model=="rvm") l=list(k=s)
   else if(model=="randomForest") l=list(mtry=s)
   else if(model=="mlp" || model=="mlpe") l=list(size=s)
   else if(model=="ksvmrbfdot" || model=="rvmrbfdot") l=list(kernel=kernel,sigma=2^s) # SVMLIB authors suggestion
   else if(model=="ksvmvanilladot") l=list(kernel=kernel,C=2^s) # Fernandez-Delgado et al, 2014 
   else if(model=="ksvmpolydot" || model=="rvmpolydot") l=list(kernel=kernel,scale=(1/10)^s,offset=1,degree=s) # Fernandez-Delgado et al, 2014
   else if(model=="rpart")
    {
     vl=vector("list",n) 
     names(vl)=rep("cp",n) # same cp name 
     for(i in 1:n) vl[[i]]=s[i] # cycle needed due to [[]] notation
     l=list(control=vl)
    }
   else if(model=="ctree")
    {
     vl=vector("list",n) 
     for(i in 1:n) vl[[i]]=party::ctree_control(mincriterion=s[i]) 
     l=list(controls=vl)
    }
   else l=NULL # not available
  }
 else  # special case: n=1
  {
   l=NULL
   if(substr(model,1,4)=="ksvm" || substr(model,1,3)=="rvm") l=list(kernel=kernel,kpar="automatic") # default heuristic
   if(substr(model,1,3)=="mlp") l=list(size=NA) # default heuristic
  }
 return(l)
}

# -------------------------------------------------------------------------
#defaultmpar=function(mpar,task="reg")
#{ if(task=="reg") metric="SAE" else if (task=="prob") metric="AUC" else metric="ACC"
#  mpar=c("holdout",2/3,metric)
#  return(mpar)
#}

defaultfeature=function(feature)
{ if(!is.na(feature[1]))
   { if(feature[1]=="s") return("simp")
     else if(substr(feature[1],1,4)=="sens" || substr(feature[1],1,4)=="simp" || feature[1]=="none") return(feature[1]) 
   }
  else if(is.na(feature[1]) && length(feature)==1) return("none")
  LF=length(feature);DLF=5
  if(LF<DLF || sum(is.na(feature))!=0) 
  {
   dfeature=c("sabs",-1,1,"holdout",2/3)
   for(i in 1:DLF) { if(is.na(feature[i])) feature[i]=dfeature[i]
                     if(i>1 && !is.na(feature[i-1]) && feature[i-1]=="kfold" && !is.na(feature[i]) && feature[i]==2/3) feature[i]=3
                   }
  }
  return(feature)
}

# ====== Feature Selection Methods (beta development) =======================
# future: put some control in feature, in order to change default smeasure = gradient???
bssearch=function(x,data,algorithm="sabs",Runs,method,model,task,search,scale,transform,smeasure="g",fstop=-1,debug=FALSE)
{ 
 #debug=TRUE
 #cat("alg:",algorithm,"smeasure:",smeasure,"\n")
 #debug=TRUE; cat("debug:",debug,"\n")
 metric=search$metric
 OUT=output_index(x,names(data)) # output
 Z=1:NCOL(data);BZ=Z;
 LZ=length(Z); if(fstop==-1)fstop=(LZ-2);LZSTOP=min(LZ-2,fstop)
 #cat("LZSTOP:",LZSTOP,"\n")
 t=0 # iteration
 NILIM=LZ #convex 0 # LZ # 2 # I need to program this better, transform this parameter into a input user parameter!!!
 JBest=worst(metric);BK=search; # worst error
 notimprove=0 # iterations with improvement
 if(algorithm=="sabs") imethod=switch(smeasure,g="sensg",v="sensv",r="sensr",a="sensa")
#cat("smeasure:",smeasure,"imethod:",imethod,"\n") 
#print(BK)
#debug=TRUE ###
 stop=FALSE;t=0;
 while(!stop)
 {
   if(algorithm=="sabs") 
   { 
    M=mining(x,data[,Z],model=model,task=task,method=method,feature=c(imethod),search=search,scale=scale,transform=transform,Runs=Runs)
    J=mean(M$error); K=centralpar(M$mpar,medianfirst=TRUE)
    LZ=length(Z);Imp=vector(length=LZ);for(i in 1:LZ) Imp[i]=mean(M$sen[,i]);
   }
   else if(algorithm=="sbs") # initial step with all Z values
   { if(t==0) { M=mining(x,data[,Z],Runs=Runs,model=model,method=method,task=task,feature="none",search=search,scale=scale,transform=transform)
               J=mean(M$error); K=centralpar(M,medianfirst=TRUE)
              }
   }
   if(isbest(J,JBest,metric)) {JBest=J;BZ=Z;notimprove=0;BK=K;BT=t;} else notimprove=notimprove+1;
   if(debug){ cat(">> t:",t,"cerr:",J,"Best:",JBest,"BT:",BT,"Z:",Z,"S:\n",sep=" ");print(K)}
   if((t==LZSTOP || notimprove==NILIM)) stop=TRUE
   else # select next attribute to be deleted 
    { 
#cat("Z:",Z,"\n")
     if(algorithm=="sabs"){IOUT=which(Z==OUT);Zimp=setdiff(1:length(Z),IOUT);IDel=which.min(Imp[Zimp]);Del=Z[Zimp[IDel]]; Z=setdiff(Z,Del) 
#cat("imp:",round(Imp[Zimp],digits=3),"\n")
#cat("iout:",IOUT,"Zi:",Zimp,"Idel:",IDel,"del:",Del,"Z:",Z,"\n")
                          }
     else if(algorithm=="sbs"){
                               J=worst(metric);Zb=Z;
                               for(i in setdiff(Z,OUT))
                               {
                                Zi=setdiff(Z,i)
                                M=mining(x,data[,Zi],Runs=Runs,model=model,method=method,task=task,feature="none",search=search,scale=scale,transform=transform)
                                Ji=mean(M$error); Ki=centralpar(M,medianfirst=TRUE); 
                                if(isbest(Ji,J,metric)) {J=Ji;K=Ki;Zb=Zi}
                                #cat("   >> i:",i,"Zi",Zi,"Ki:",Ki,"err:",Ji,"best:",J,"\n")
                               }
                               Z=Zb
                              }
    }
   t=t+1
 }
 #debug=TRUE
 if(debug) {cat(" | t:",BT,"Best:",JBest,"BZ:",BZ,"BS:\n",sep=" "); print(BK) }
 return(list(attributes=BZ,search=BK)) # attributes and search parameter
}
# ====== End of Feature Selection Methods ==================================

# fast uniform design: only works for some uniform design setups: 2, 3 factors and 5,9,13 points
uniform_design=function(limits,factors=length(limits)/2,points,center=NULL)
{
 ud=matrix(nrow=points,ncol=factors) # gamma, C (and epsilon?)
 
 SEQ=vector("list",length=factors)
 for(i in 1:factors)
  { ini=2*(i-1)+1;SEQ[[i]]=seq(limits[ini],limits[ini+1],length=points) }

# 2, 5 
if(factors==2)
{ if(points==5) m=t(matrix(c(1,2,2,5,4,1,5,4,3,3),ncol=5)) 
  else if(points==9) m=t(matrix(c(5,5,1,4,7,8,2,7,3,2,9,6,8,3,6,1,4,9),ncol=9))
  else if(points==13) m=t(matrix(c(5,4,12,3,2,11,9,10,7,7,6,13,3,2,11,12,13,8,10,5,1,6,4,9,8,1),ncol=13))
}
else if(factors==3)
{ if(points==5) m=t(matrix(c(5,3,3,4,4,5,3,1,1,2,5,2,1,2,4),ncol=5))
  else if(points==9) m=t(matrix(c(3,9,6,9,4,7,7,1,4,2,2,8,1,6,3,8,8,2,4,3,1,5,5,5,6,7,9),ncol=9))
  else if(points==13) m=t(matrix(c(9,9,1,7,4,5,8,13,10,2,3,11,4,6,2,13,7,9,6,1,8,10,5,12,12,11,6,3,12,4,1,8,7,5,10,13,11,2,3),ncol=13))
}
for(i in 1:nrow(m))
for(j in 1:ncol(m))
 {
  val=SEQ[[j]][m[i,j]]
  ud[i,j]=val
 }
if(length(center)>0) # delete center point
{ ALL=1:NROW(ud)
  I=which(ud[ALL,1]==center[1])
  ud=ud[setdiff(ALL,I),]
} 
return (ud)
}

# --- the best fit internal function ---------------------------------------------------
mgrid=function(search,x,data,model,task,scale,transform,fdebug,...)
{
 #fdebug=1
 metric=search$metric
 bestval=worst(metric)
 stop=FALSE; notimprove=0

 # check type of search$search parameter: vector or object
 if(!is.vector(search$search[[1]][[1]])) object=TRUE else object=FALSE
 ###cat("object:",object,"\n")
 N=length(search$search); si=list(smethod="none",convex=search$convex,metric=metric); nsi=names(search$search)
 saux=list()
 # set saux names:
 for(i in 1:N) saux[[ nsi[i] ]]= search$search [[ nsi[i] ]] [1]
 L=vector(length=N)
 for(i in 1:N) L[i]=length(search$search[[i]])  

 if(search$smethod=="grid"||search$smethod=="2L") # explore combinations of combinations
 {
  LS=prod(L)
  s=matrix(ncol=N,nrow=LS) # set the search space 
  for(i in 1:N) #
    {
     if(i==1) E=1 else E=E*L[i-1]
     s[,i]=rep(1:length(search$search[[i]]),length.out=LS,each=E) #
    }
 }
 else if(search$smethod=="matrix" ||substr(search$smethod,1,2)=="UD") # generate s
 {
  LS=max(L)
  for(i in 1:N) if(length(search$search[[i]])==1) search$search[[i]]=rep(search$search[[i]],LS)

  s=matrix(ncol=N,nrow=LS) # set the search space 
  for(i in 1:N) #
    {
     s[,i]=1:LS
    }
 }
 if(fdebug) { cat(search$smethod,"with:",LS,"searches") 
              if(is.character(metric)) cat(" (",metric," values)\n",sep="") else cat("\n")}
 i=1 
 while(!stop)
   {
    for(j in 1:N) 
      { 
        if(object) #{ cat("object IF:\n"); 
                     saux[[j]]=search$search[[j]][[s[i,j]]] #}
        else #{ cat("object ELSE:\n"); 
           saux[[ j ]]= search$search[[ j ]] [s[i,j]] # }
      }
    si$search=saux
 
    eval=try( mining(x,data,Runs=1,method=search$method,model=model,task=task,scale=scale,search=si,transform=transform,...)$error, silent= TRUE) 
    if(class(eval)[1]=="try-error") eval=worst(metric)
    if(isbest(eval,bestval,metric)) {bestval=eval;notimprove=0;best=si} 
    else notimprove=notimprove+1
    if(fdebug) { cat("i:",i,"eval:",eval,"best:",bestval,"\n") }
    i=i+1
    if(search$convex>0 && notimprove==search$convex) {stop=TRUE; if(fdebug) cat (" -> stop since not improve searches =",notimprove,"\n") }
    else if(i>LS) stop=TRUE
    } 
 best$bestval=bestval
 return (best)
}
# -------- end of best fit -----------

output_index=function(x,namesdata)
{ # check also: y = model.extract(m, response) for multi targets ? => future work?
 #x=as.character(x)[2]
 #T=terms(x,data=data); Y=names(attr(T,"factors")[,1])[1] 
 #return(which(names(data)==Y))
 return(which(namesdata==(as.character(x)[2])))
}

transform_needed=function(transform) { return (switch(transform,log=,logpositive=,positive=,scale=TRUE,FALSE)) }
feature_needed=function(feature) { return (switch(feature,sbs=,sabsg=,sabsv=,sabsr=,sabs=TRUE,FALSE)) }

defaultask=function(task="default",model="default",output=1)
{ if(task=="default") 
    { 
      if(is.character(model)) task=switch(model,bagging=,boosting=,lda=,multinom=,naiveBayes=,qda="prob",cubist=,mr=,mars=,pcr=,plsr=,cppls=,rvm="reg","default")
      if(task=="default") { if (is.factor(output)) task="prob" else task="reg" }
    }
  else if(substr(task,1,1)=="c") task="class" else if(substr(task,1,1)=="p") task="prob"
  return (task)
}

defaultmodel=function(model,task="reg")
{ 
 if(model=="lr" || model=="logistic") model="multinom" 
 else if(model=="dt") model="rpart" 
 else if(model=="knn") model="kknn" 
 else if(model=="svm") model="ksvm" 
 else if(model=="randomforest") model="randomForest" 
 else if(model=="naivebayes") model="naiveBayes" 
 else if(model=="default")
  { if(substr(task,1,3)=="reg") model="mr" else model="rpart" 
  }
 return(model)
}

#--- end of generic fit ------------------------------------------------------

#--- generic predict method --------------------------------------------------
# @@@ search pattern for this function
# object - a model created with fit
# newdata - data.frame or matrix or vector or factor
#           if formula, please use: predict(...,model.frame(data))

#if(!isGeneric("predict")){
#  if (is.function("predict"))
#    fun = predict 
#  else fun = function(object) standardGeneric("predict")
#  setGeneric("predict", fun)
#}
setGeneric("predict")
setMethod("predict",signature(object="model"),
function(object,newdata,...){

if(NCOL(newdata)>length(object@attributes)) newdata=newdata[,object@attributes] # due to feature selection
# --- if code for the models ---
if(is.list(object@model)) P=object@model$predict(object@object,newdata)
else # test all character models
{
 object@model=defaultmodel(object@model,object@task) # due to compatibility issues (e.g. "svm" -> "ksvm")

 # if scale is needed:
 if(substr(object@model,1,3)=="mlp"){ if(object@scale=="inputs" || object@scale=="all") newdata=scaleinputs2(newdata,object@object$cx,object@object$sx) }

 if(object@model=="naive")
  { if(object@task=="reg") P=rep(object@object,length=nrow(newdata)) # numeric
    else { L=levels(object@object[1])
           P=matrix(nrow=nrow(newdata),ncol=length(L)); P[,]=0
           P[,which(L==object@object)]=1 # matrix
           if(object@task=="class") P=majorClass(P,L) # factor
         }
  }
 else if(object@model=="ctree" || object@model=="rpart" || object@model=="randomForest")
  { if(object@model=="rpart") type=switch(object@task,class="class",prob="prob","vector") # rpart
    else type=switch(object@task,reg=,class="response",prob="prob") # ctree or randomForest
    P=predict(object@object,newdata,type=type)
    if(object@model=="randomForest" && object@task=="prob") attr(P,"class")=NULL
  }
 else if(object@model=="multinom") 
  {
    type=switch(object@task,class="class","probs")
    P=predict(object@object,newdata,type=type)
    if(length(object@levels)==2 && object@task=="prob") {MM=matrix(ncol=2,nrow=length(P));MM[,2]=P;MM[,1]=1-P;P=MM;}
  }
 else if(object@model=="bagging" || object@model=="boosting")
  { P=predict(object@object,newdata)
    if(object@task=="class") P=P$class else P=P$prob
  }
 # pure regression:
 else if(object@model=="mr")
    { P=predict(object@object$mlp,newdata)[,1] }
 else if(object@model=="mars")
    { P=predict(object@object,newdata[,-object@outindex])[,1] }
 else if(object@model=="cubist")
    { P=predict(object@object,newdata[,-object@outindex],...)}
 else if(object@model=="pcr"||object@model=="plsr"||object@model=="cppls")
    { P=predict(object@object,newdata,comp=1:(object@object$ncomp)) }
 # pure classification:
 else if(object@model=="lda" || object@model=="qda")
   { P=predict(object@object,newdata)
     if(object@task=="class") P=P$class else P=P$posterior 
   }
 else if(object@model=="naiveBayes") # || object@model=="wnaivebayes")
   { if(is.null(object@object)) P=newdata[,-object@outindex] # only 1 input
     else { type=switch(object@task,class="class","raw")
            P=suppressWarnings(predict(object@object,newdata[,-object@outindex],type=type))
          }
   }
 # reg or cla:
 else if(object@model=="kknn") # k-nearest neighbour, uses kknn implementation
   { 
     P=do.call(kknn::kknn,c(object@object,list(test=newdata)))
     if(object@task=="reg") P=P$fitted.values
     else P=switch(object@task,prob=P$prob,P$fitted.values)
   }
 else if(object@model=="mlp") 
   { if(object@task=="reg") 
       { P=predict(object@object$mlp,newdata)[,1]
         if(object@scale=="all") P=invtransform(P,"scale",object@object$cy,object@object$sy)
       }
     else
      { type=switch(object@task,class="class","raw")
        P=predict(object@object$mlp,newdata,type=type)
        if(length(object@levels)==2 && object@task=="prob") {MM=matrix(ncol=2,nrow=length(P));MM[,2]=P;MM[,1]=1-P;P=MM;}
      }
   }
 else if(object@model=="mlpe") # 
   {  MR=object@mpar$nr
      if(object@task!="reg"){L=length(object@levels); P=matrix(0,ncol=L,nrow=NROW(newdata))} 
      else P=rep(0,NROW(newdata))
      for(i in 1:MR)
      {
       if(object@task=="reg") 
        { P1=predict(object@object$mlp[[i]],newdata)[,1]
          if(object@scale=="all") P1=invtransform(P1,"scale",object@object$cy,object@object$sy)
          P=P+P1
        }
       else
       { 
        P1=predict(object@object$mlp[[i]],newdata)
        if(L==2) {MM=matrix(ncol=2,nrow=length(P1));MM[,2]=P1;MM[,1]=1-P1;P1=MM;}
        P=P+P1
       }
      }
      P=P/MR
      if(object@task=="class") P=majorClass(P,L=object@levels)
   }
 else if(object@model=="ksvm") 
   {
      if(object@task=="prob") P=predict(object@object,newdata,type="probabilities")
      else if(object@task=="class") P=predict(object@object,newdata) 
      else P=predict(object@object,newdata,type="response")[,1]
   }
 else if(object@model=="rvm") P=predict(object@object,newdata,type="response")[,1]

 # transform P into the same standard format: 
 if(object@task=="reg")
        { if(is.matrix(P)) P=P[,1] 
          if( !is.null(names(P)) ) attr(P,"names")=NULL
        } 
 else if(object@task=="prob")
        {
         if(is.list(P)) P=t(matrix(data=unlist(P),nrow=length(object@levels))) 
         if(is.matrix(P)) { V=vector("list",2);V[[1]]=NULL;V[[2]]=object@levels;attr(P,"dimnames")=V } 
        }
 else # cla
        { if(is.character(P)) P=factor(P,levels=object@levels)
          if(!is.null(names(P))) names(P)=NULL
        }
} # end test all model character types
 if(substr(object@task,1,3)=="reg" && transform_needed(object@transform))P=invtransform(P,object@transform)
 return(P)
})

svm.fit = function(x,data,task,exponentialscale=FALSE,...)
{ 
 args=do.call(list,list(...)) # extra arguments (needed for different kernel?)
 # to avoid this error: object "scal" not found, I will adopt SCALED=TRUE
 SCALED=TRUE
 outindex=output_index(x,names(data))

 # think about kernel change? 
 if(!is.null(args$C) && is.na(args$C) ) # set the C value using heuristic of CheMa04
  {
   if(substr(task,1,3)=="reg") # scale the output
   { CY=mean(data[,outindex])
     SY=sd(data[,outindex])
     Y=xtransform(data[,outindex],"scale",A=CY,B=SY)
   }
   else{CY=0; SY=0; Y=data[,outindex]}
   if(is.factor(Y)) {MEANY=0;SDY=1;}
   else { MEANY=mean(Y); SDY=sd(Y) }
   C=max(abs(MEANY+3*SDY),abs(MEANY-3*SDY))
   args[["C"]]=C
  }
 
 if(!is.null(args$epsilon) && is.na(args$epsilon) && substr(task,1,3)=="reg") # set the epsilon, in regression
 {
  data2=data; data2[,outindex]=Y
  N=nrow(data)
  K3=(kknn::kknn(x,train=data2,test=data2,k=3))$fitted.values # 3NN, as argued by CheMa04
  SUM=sum((Y-K3)^2)
  SD=sqrt(1.5/N * SUM)
  if(N>100) epsilon=3*SD*sqrt(log(N)/N) # eq.14 of CheMa04: use for large N 
  else epsilon=SD/sqrt(N) # eq.13 of CheMa04 : ok if N small 
  args[["epsilon"]]=epsilon
 }

 if(!is.null(args$sigma))
  { if(args$sigma=="automatic") args[["kpar"]]="automatic"
    else if(is.null(args$kpar)) args[["kpar"]]=list(sigma=args$sigma)
    else if(args$kpar!="automatic") args[["kpar"]]=c(args$kpar,list(sigma=args$sigma))
    args[["sigma"]]=NULL
  }

 if(exponentialscale) 
  { 
   if(!is.null(args$C)) args$C=2^args$C
   if(!is.null(args$epsilon)) args$epsilon=2^args$epsilon
   if(!is.null(args$nu)) args$nu=2^args$nu
   if(!is.null(args$kpar) && is.list(args$kpar))
     { for(i in 1:length(args$kpar))
         if(is.numeric(args$kpar[[i]])) args$kpar[[i]]=2^args$kpar[[i]]
     }
  }
 if(task=="prob") args[["prob.model"]]=TRUE
 args[["scaled"]]=SCALED;
 #sink("/dev/null")
 # capture to avoid "maximum number of iterations reached" verbose:
 capture.output( { M = do.call("ksvm",c(list(x,data),args)) } )
 #sink()
 return(M)
}

#------------------------------------------------------------------------------------------
# Create and fit a MLP
# current version: performs NetRuns trainings and selects the MLP 
#                  with lowest minimum criterion (SSE+decay?) 
# size - number of hidden nodes 
                # decay - decay value (in 0 ... 1)
# nr - NetRuns - number of MLP trainings
# maxit   - maximum number of training epochs
# type = 1 mlp, 2 mlpe, 3 mlp neuralnet (rprop)
mlp.fit= function(x,data,size=2,nr=3,task,scale,type=1,...) # metric="SSE")
{
 if(is.na(size)) { size=round( (ncol(data)-1) /2) # simple heuristic for hidden nodes
                   if(task!="reg" && size>10) size=10 # nnet does not allow huge size for classification tasks
                 }

 if(type==2) { MinNet=vector("list",nr) } # ensemble

 outindex=output_index(x,names(data))
 if(scale=="inputs" || scale=="all") # scale the inputs if needed 
   { S=scaleinputs(data,outindex)
     data=S$data; cx=S$cx; sx=S$sx
   }
 else{cx=NULL; sx=NULL}
 if(substr(task,1,3)=="reg" && scale=="all") # scale the output if needed
   { CY=mean(data[,outindex])
     SY=sd(data[,outindex])
     data[,outindex]=xtransform(data[,outindex],"scale",A=CY,B=SY)
   }
 else{CY=0; SY=0}

 MinError=Inf # maximum real number
 if(substr(task,1,3)=="reg") LINOUT=TRUE else LINOUT=FALSE

 goon=TRUE;i=1;
 if(size>0) skip=FALSE else skip=TRUE
 size=round(size) # because of "2L" !
 while(goon)#for(i in 1:nr)
    {
      if(type==1)
      { Net=nnet::nnet(x,data=data,size=size,trace=FALSE,skip=skip,linout=LINOUT,...)
        err=Net$value #---- this code used minimum penalized error (SSE+decay...) #err=Sse(data[,outindex],Net$fitted.values)
      }
      else if(type==2)
        MinNet[[i]]=nnet::nnet(x,data=data,size=size,trace=FALSE,skip=skip,linout=LINOUT,...)
      #else # type==3
      #{
      # #if(metric=="SSE") ERR="sse" 
      # #else ERR="sad" # sse2 #sad
      # NM=names(data); ALL=setdiff(1:ncol(data),outindex) 
      # fx=as.formula(paste(NM[outindex], paste(NM[ALL],collapse='+'), sep='~'))
      # Net=neuralnet(fx,data=data,hidden=H,linear.output=LINOUT,rep=1) #,err.fct=ERR)
      # #cat(data[1,outindex],Net$response[1],"\n")
      # err=Sse(data[,outindex],Net$response)
      #}
      if(type!=2 && MinError>err){ MinError=err; MinNet=Net;}
      i=i+1
      if(i>nr || MinError==0) goon=FALSE
    }
 NNmodel=list(mlp=MinNet,cx=cx,sx=sx,cy=CY,sy=SY,nr=nr)
 return (NNmodel)
}
#---------------------------------------------------------------------------------
# Performs the Data Mining stage.
# The dataset (x - inputs, y - output) is used to fit the model.
# It returns the performance of a model using all data, holdout or kfold.
# PARAMETERS:
# x - a formula (default) or input dataframe, matrix, vector or factor
# data - NULL or the dataframe (default) or the output vector or factor 
# Runs - number of runs (e.g. 1, 10, 20, 30, 100, 200)
# method - vector with 2 arguments (see point A) or list (point B)
#      A) 1st: estimation method: "all", "holdout", "holdoutorder", "holdoutinc", "kfold", "kfoldo"
#          - holdoutorder is similar to holdout, except that a sequencial ordered split is used
#          - holdoutfor is equal to holdoutorder but is used for time series forecasting, special holdout that does not use out-samples
#            instead of the random split. This is quite useful for temporal data.
#          - holdoutinc - incremental retraining method (vpar=batchsize)
#         2nd: the estimation method parameter (e.g. 2/3, 3/4 - holdout or 3, 5 or 10 - kfold)
# model - the type of model used. Currently: "lm", "svm","mlp","logistic","naivebayes", "knn", "rpart","randomforest"
# task - "default", "class" (or "classification", "cla"), "reg" (or "regression")
# search - NULL or a vector of seaching hyperparameters (used for "mlp", "svm", "knn")
# mpar - vector with the internal searching parameters (for "mlp", "svm", "knn"):
# mpar[1] - knn: internal estimation method ("all", "kfold", "holdout", "holdoutorder") 
# mpar[2] - knn: internal split ratio or number of kfolds
# mpar[1] - svm or mlp: 1st internal parameter (svm - C, mlp - nr)
# mpar[2] - svm or mlp: 2nd internal parameter (svm - Epsilon, mlp - maxit) 
# mpar[3] - svm or mlp: internal estimation method ("all", "kfold", "holdout", "holdoutorder")
# mpar[4] - svm or mlp: internal split ratio or number of kfolds 
# scale   - "none", "inputs", "all" - if a 0 mean and 1 std scale is used
#           "default" - the best generic method is automatically applied
# transform - if "log" and task="reg" then a logistic transformation is used in the output  
#             if "logpositive" is equal to "log" except that negative output values are set to zero.
# feature - "default" - "none" or "sens" if model="mlp","svm",...
#           "sens" - sensitivity needs to be stored for each fit... ("sensg" - sensitivity with gradient)
#           "simpv" "simpg" - similar to above but stores more info (sresponses) (v -variance, g- gradient)
#           "none" - if no feature selection is used
#            when using a feature selection the following scheme should be adopted:
#            -  feature=c(FM,Runs,Method,MethodPar,Search), where
#                  FM - "sabs", "sfs", "sffs", "sbs", "sbfs"
#                       "sabsv" - sabs with variance, "sabsg" - sabs with gradient
#                  Runs/Fixed - number of feature selection runs OR "fN", where N is the number of fstop deletions
#                  Method - "holdout" or "kfold"
#                  MethodPar - splitratio or kfold (vpar)
#                  Search - primary search parameter (mlp, svm, knn)
#---------------------------------------------------------------------------------
mining=function(x,data=NULL,Runs=1,method=NULL,model="default",task="default",search="heuristic",mpar=NULL,feature="none",scale="default",transform="none",debug=FALSE,...)
{
 args=do.call(list,list(...)) # extra arguments (needed for neighbors)
 if(is.null(data)) 
   { data=model.frame(x) # only formula is used 
     outindex=output_index(x,names(data))
     x=as.formula(paste(names(data)[outindex]," ~ .",sep=""))
   }
 #else { if(feature[1]=="aries" && length(feature)>1) outindex=as.numeric(feature[2]) 
 else outindex=output_index(x,names(data)) 

 firstm = proc.time()
 previoustm=firstm

 if(!is.list(model)) model=defaultmodel(model,task) # due to compatibility issues

 time=vector(length=Runs)
 error=vector(length=Runs)
 PRED=vector("list",Runs); TEST=vector("list",Runs)

 feature=defaultfeature(feature)
 task=defaultask(task,model,data[1,outindex])
 if(task!="reg" && !is.factor(data[,outindex])) data[,outindex]=factor(data[,outindex])

 vmethod=readmethod(method,Runs) 
#cat("vmethod:\n")
#print(vmethod)
 if(vmethod$m=="holdout" && (vmethod$mode=="incremental" || vmethod$mode=="rolling")) 
    Runs=floor( (nrow(data)-vmethod$w+vmethod$p)/vmethod$i ) 
#cat("runs:",Runs,"\n")

 par=vector("list",length=Runs)
 #if(is.null(mpar)) mpar=defaultmpar(mpar,task)
 search=readsearch(search,model,task,mpar,...)
 metric=search$metric

 if(substr(feature[1],1,4)=="sens" || substr(feature[1],1,4)=="sabs" || substr(feature[1],1,4)=="simp") 
    { 
       if(vmethod$m=="kfold") MULT=vmethod$p
       else MULT=1
       SEN=matrix(ncol=NCOL(data),nrow=(Runs*MULT) )
       if(substr(feature[1],1,4)!="sens") { SRESP=TRUE;FSRESP=TRUE;}
       else {SRESP=NULL; FSRESP=FALSE;}
       imethod=switch(feature[1],sabsv=,simpv="sensv",sabsg=,sabs=,simp=,simpg="sensg",sabsr=,simpr="sensr",simpa="sensa",feature[1])
    }
 else {SEN=NULL; SRESP=NULL;FSRESP=FALSE;imethod="none"}

 #cat("is.null sen?",is.null(SEN),"\n")
 if(feature[1]!="none" && substr(feature[1],1,4)!="sens") # store features?
  { if(vmethod$m=="kfold") attrib=vector("list",Runs*vmethod$p)
    else attrib=vector("list",Runs)
  }
 else attrib=NULL 
 HOLD=FALSE

 for(i in 1:Runs)
 {
  if(vmethod$m=="all") {V=list();V$tr=1:nrow(data);V$ts=V$tr;HOLD=TRUE}
  else if(vmethod$m=="holdout") #all holdout types  
     { if(length(vmethod$seed)==1) si=1 else si=i
       V=holdout(data[,outindex],ratio=vmethod$p,mode=vmethod$mode,iter=i,seed=vmethod$s[si]);HOLD=TRUE
     }
  if(HOLD)     
     { 
       L=fit(x=x,data=data[V$tr,],model=model,task=task,search=search,scale=scale,transform=transform,feature=feature,...)
       if(model=="cubist" && !is.null(args$neighbors)) P=predict(L,data[V$ts,],neighbors=args$neighbors)
       else P=predict(L,data[V$ts,])
       TS=data[V$ts,outindex]
       if(!is.null(SEN)) IMPORTANCE=Importance(L,data,method=imethod,responses=FSRESP) # store sen
     }
  else # KFOLD 
     { if(length(vmethod$seed)==1) si=1 else si=i
       L=crossvaldata(x,data,fit,predict,ngroup=vmethod$p,mode=vmethod$mode,seed=vmethod$s[si],model=model,task=task,feature=feature,
                      scale=scale,search=search,transform=transform,...)
       P=L$cv.fit
       TS=data[,outindex]
     }

  if(!is.null(SEN))
     { 
       if(vmethod$m=="kfold")
         {  Begi=(i-1)*MULT+1; Endi=Begi+vmethod$p-1;
            SEN[Begi:Endi,]=L$sen # store sen
            if(!is.null(SRESP)) 
                 { 
                  if(i==1) SRESP=L$sresponses # store sen
                  else{ for(j in 1:length(SRESP))
                           {  
                            if(!is.null(L$sresponses[[j]]) ) # && !is.null(IMPORTANCE$sresponses[[j]])$x) # $x)) 
                              { if(is.null(SRESP[[j]])) SRESP[[j]]=L$sresponses[[j]]
                                else{ SRESP[[j]]$x=c(SRESP[[j]]$x,L$sresponses[[j]]$x);
                                      if(task=="prob") SRESP[[j]]$y=rbind(SRESP[[j]]$y,L$sresponses[[j]]$y,deparse.level=0)
                                      else if(task=="class") SRESP[[j]]$y=addfactor(SRESP[[j]]$y,L$sresponses[[j]]$y)
                                      else SRESP[[j]]$y=c(SRESP[[j]]$y,L$sresponses[[j]]$y)
                                    }
                              }
                           }
                      }
                 }
         }
       else
         { 
           SEN[i,]=IMPORTANCE$imp # store sen
           if(!is.null(SRESP)) 
            { if(i==1) SRESP=IMPORTANCE$sresponses # store sen
              else{ for(j in 1:length(SRESP)) 
                       { 
                         if(!is.null(IMPORTANCE$sresponses[[j]]) ) # && !is.null(IMPORTANCE$sresponses[[j]])$x) # $x)) 
                           { if(is.null(SRESP[[j]])) SRESP[[j]]=IMPORTANCE$sresponses[[j]]
                             else{ SRESP[[j]]$x=c(SRESP[[j]]$x,IMPORTANCE$sresponses[[j]]$x);
                                   if(task=="prob") SRESP[[j]]$y=rbind(SRESP[[j]]$y,IMPORTANCE$sresponses[[j]]$y,deparse.level=0)
                                   else if(task=="class") SRESP[[j]]$y=addfactor(SRESP[[j]]$y,IMPORTANCE$sresponses[[j]]$y)
                                   else SRESP[[j]]$y=c(SRESP[[j]]$y,IMPORTANCE$sresponses[[j]]$y)
                                 }
                           }
                       }
                  }
            }
         }
     }

  if(!is.null(attrib))
  { if(vmethod$m=="kfold") for(j in 1:vmethod$p) attrib[[ (i-1)*vmethod$p+j ]]=L$attributes[[j]]
    else attrib[[i]]=L@attributes
  }
  #if(npar>0) { if(vmethod$m=="kfold"){ ini=(i-1)*vmethod$p+1;end=ini+vmethod$p-1;par[ini:end,]=L$mpar} else par[i,]=as.numeric(L@mpar[1,]) }
  if(is.list(L)) par[[i]]=L$mpar else par[[i]]=L@mpar
  # store P and TS:
  PRED[[i]]=P
  TEST[[i]]=TS
  error[i]=do.call("mmetric",c(list(y=TS,x=P),metric))

  #print(PRED[[i]]); #print(TEST[[i]]) #cat(" e:",error[i],"\n")
  curtm= proc.time()
  time[i]=(curtm-previoustm)[3] # execution time for run i
  previoustm=curtm
  if(debug) { 
              cat("Run: ",i,", ",error[i],sep="")
              if(is.character(metric)) cat(" (",metric," values)",sep="") 
              cat(" time:",round((curtm-firstm)[3],digits=1),"\n",sep="") 
            }
 }
 if(!is.null(SRESP)) { for(i in 1:length(SRESP)) if( !is.null(SRESP[[i]]) && is.factor(data[,i])) SRESP[[i]]$x=factor(SRESP[[i]]$x) }
 #if(npar>0) { par=data.frame(par); names(par)=modelnamespar(model);}
 return(list(time=time,test=TEST,pred=PRED,error=error,mpar=par,model=model,task=task,method=vmethod,sen=SEN,sresponses=SRESP,runs=Runs,attributes=attrib,feature=feature))
}

readmethod=function(method,Runs=1)
{
 seed=NULL;VFLAG=FALSE;mode=NULL;increment=1;window=10
 if(is.null(method)) method="holdout"
 # read vmethod and assign mode:
 vmethod=method[1]
 if(substr(vmethod,1,6)=="kfoldo") {vmethod="kfold";mode="order";vpar=10}
 else if(substr(vmethod,1,6)=="kfoldr") {vmethod="kfold";mode="random";vpar=10}
 else if(substr(vmethod,1,5)=="kfold") {vmethod="kfold";mode="stratified";vpar=10}
 else if(substr(vmethod,1,8)=="holdouto") {vmethod="holdout";mode="order";vpar=2/3}
 else if(substr(vmethod,1,9)=="holdoutro") {vmethod="holdout";mode="rolling";vpar=1}
 else if(substr(vmethod,1,8)=="holdoutr") {vmethod="holdout";mode="random";vpar=2/3}
 else if(substr(vmethod,1,8)=="holdouti") {vmethod="holdout";mode="incremental";vpar=1}
 else if(vmethod=="holdout") {mode="stratified";vpar=2/3}
 else if(vmethod=="all") {mode="order";vpar=1}
 else # none of above, then use default: 
 { vmethod="holdout";mode="stratified";vpar=2/3;VFLAG=TRUE}

 if(length(method)>1 && !VFLAG) vpar=as.numeric(method[2])
 if(mode=="incremental"||mode=="rolling") 
  { 
   if(length(method)>2 && !is.na(method[3])) window=as.numeric(method[3])
   if(length(method)>3 && !is.na(method[4])) increment=as.numeric(method[4])
  }
 else if(length(method)>2 && !is.na(method[3])) seed=as.numeric(method[3:length(method)])

 if(!is.null(seed) && length(seed)<Runs) # need to build a distinct seed for each run:
 { set.seed(seed)
   seed=sample(1:10000,Runs) 
 } 
 return(list(m=vmethod,p=vpar,s=seed,mode=mode,i=increment,w=window))
}

# mpar is a mining object component
# medianfirst is only used for multiple runs of feature selection!
centralpar=function(mpar,medianfirst=FALSE)
{
 Runs=length(mpar) # runs
 # has folds? 
 if( is.null( names(mpar[[1]]) ) ) {par=mpar[[1]][[1]];Folds=length(mpar[[1]])} else {par=mpar[[1]];Folds=1}
 Args=length(par)
 A=sapply(mpar,unlist) # 
 # AA<<-A
 RA=nrow(A)
 DIST=RA/Folds
 # compute now the median value:

 k=1 
 for(a in 1:Args)
     if(is.list(par[[a]])) # kpar issue!
       { 
        L=length(par[[a]])
        for(i in 1:L) 
          {
             if( is.numeric(par[[a]][[i]]) || is.logical(par[[a]][[i]]) || is.character(par[[a]][[i]]) ) 
               par[[a]][[i]]=centralaux(par[[a]][[i]],(k+i-1),Folds,RA,DIST,A,medianfirst)
             k=k+1
          }        
       } 
     else if(is.numeric(par[[a]]) || is.logical(par[[a]]) || is.character(par[[a]]) ) 
       {
        par[[a]]=centralaux(par[[a]],k,Folds,RA,DIST,A,medianfirst)
        k=k+1
       }
     else k=k+1
 return(par)
}

centralaux=function(elem,k,Folds,Rows,BY,A,medianfirst=FALSE)
{
 if(Folds==0) ROWS=k else ROWS=seq(k,Rows,by= BY )
 aux=A[ROWS,]
 if(is.numeric(elem)) 
  { aux=as.numeric(aux) 
    if(medianfirst) res=medianfirst(aux)$val
    else res=median(as.numeric(aux)) # central value ?
  }
 else
  { T=table(aux) # mode: most common element
    if(is.logical(elem)) res=as.logical(names(which.max(T)))
    else res=names(which.max(T)) # character or factor
  }
 return(res) 
}

# 
getmetric=function(mpar,task="reg") 
{ 
  if(is.null(mpar)) defaultmetric=TRUE
  else { metric=mpar[length(mpar)]
         if(is.numeric(metric)||is.na(metric)) defaultmetric=TRUE
         else if(!is.mmetric(metric)) defaultmetric=TRUE
         else defaultmetric=FALSE
       }
  if(defaultmetric) 
   { if(task=="reg") metric="SAE" else if(task=="prob") metric="AUC" else metric="ACC" }
  return(metric)
}

getmethod=function(mpar)
{
 L=length(mpar)
 if(L==1) {vmethod=mpar[1];vmethodp=NA}
 else if(L==2 || L==3) {vmethod=mpar[1];vmethodp=mpar[2]}
 else if(L>3) {vmethod=mpar[3];vmethodp=mpar[4]}
 else {vmethod="holdout";vmethodp=2/3}

 if(vmethod=="all") vmethodp=NA   
 else if(substr(vmethod,1,7)=="holdout" && is.na(vmethodp)) vmethodp=2/3
 else if(substr(vmethod,1,5)=="kfold" && is.na(vmethodp)) vmethodp=3
 method=readmethod(c(vmethod,vmethodp))
 return (c(method$m,method$p,NA))
}

# TO DO: remove sampling argument, not used!!!
#-- example code of variant to allow direct use of models built with nnet, ksvm, randomForest, etc...
#-- needs more coding... still in beta phase
#-- type should be explicitly defined, as it depends on the predict function.
#-- for instance: nnet - type="raw" or "class", kvsm - type ="response", "votes, etc...

# method="randomforest" - Leo Breiman method; "sens","sensv", "sensg", "sensv", Embrechts method ; "DSA" - uses all data.frame
# sampling - "regular" - grid spread from min to max or "quantile" - based on the distribution percentiles
# responses - if TRUE then the full response (y) values are stored
# interactions - default NULL else is a vector of one or two variables to explore 2 variable interactions (VECC and VEC3)...
# baseline - only used by "sens" method: "mean" (pure kewley method), "median" or row of base values

# To be used only if model is not of "model" type:
# outindex - the output column index 
# task - NULL or "class" or "prob" or "reg" 
# PRED - new prediction function, needs to return a vector for "reg" or matrix of probabilities for "prob"

# to do: think about I-D sensitivity, where I is the number of features? 
# measure - "variance", "range", "gradient" (new: "lenght" !!!)

# to do: check model.matrix(~ a + b, dd)
# to do: new method="shist" that uses h$mids => average(y) => sensitivity?
# to do: input something that distinguishes between L = max levels for factor and all levels for factor are used...
# aggregation: number of B boxplot statistics, if 1 then average, else if 3 then min, average, max!
# LRandom: number of M montecarlo samples, -1 means all
#
# to do: new parameter Lfactor: FALSE - up to RealL, TRUE - up to Levels
# method="SA" or "sens";
#       ="DSA" 
#       ="MSA" 
#       ="CSA"
# responses =TRUE, FALSE or "ALL"
Importance=function(M,data,RealL=7,method="1D-SA",measure="AAD",sampling="regular",baseline="mean",responses=TRUE,
                    outindex=NULL,task="default",PRED=NULL,interactions=NULL,Aggregation=-1,LRandom=-1,MRandom="discrete",Lfactor=FALSE)
{
#model=M;data=d;RealL=6;method="sens";measure="variance";sampling="regular";responses=TRUE;outindex=NULL;task=NULL;
#model<=M;DDD<=data;RealL<=RealL;method<=method;measure<=measure;sampling<=sampling;responses<=responses; 
 SDD=NULL
 if(class(M)=="model"){outindex=M@outindex; task=M@task; Attributes=M@attributes} # rminer model
 else # another R supervised learning model
 { 
  task=defaultask(task,model="default",data[1,outindex])
  if(is.null(PRED)) 
  { nclass=class(M);nclass=nclass[length(nclass)] 
    if(task=="prob")
     {PRED=switch(nclass,
                       lda=,qda=function(M,D){predict(M,D)$posterior},
                       randomForest=function(M,D){predict(M,D,type="prob")},
                       # and so on, to be completed...
                       NULL)
     }
    else if(task=="class")
     {PRED=switch(nclass,
                       lda=,qda=function(M,D){predict(M,D)$class},
                       randomForest=function(M,D){predict(M,D,type="response")},
                       # and so on, to be completed...
                       NULL)
     }
    else{ PRED=switch(nclass,
                       lm==function(M,D){predict(M,D)},
                       randomForest=function(M,D){predict(M,D)},
                       # and so on, to be completed...
                       NULL)}
  }
  Attributes=1:ncol(data)
 }

 # use the method proposed by Leo Breiman 2001:
 if(method=="randomForest" || method=="randomforest") # user should be sure this is a randomforest!
 { if(class(M)=="model") Sv=randomForest::importance(M@object,type=1)
   else Sv=randomForest::importance(M,type=1) # true randomforest
   # need to change this code!!! it seems that negative values mean less important inputs when compared with positive ones
   # the problem is how to scale and get a percentage ? 0% for the most negative input?
   ASv=abs(Sv)
   imp=100*abs(ASv)/sum(ASv) 
   RESP=NULL
   measure=""
 }
 #else if(substring(method,1,4)=="sens" || method=="DSA" || method=="MSA") # set SPREAD
 else # all other methods
 {
  measure=switch(method,sensv="variance",sensg="gradient",sensr="range",sensa="AAD",measure)
  method=switch(method,sensa=,sensv=,sensg=,sensr="sens",method)

  if(method=="SA") method="sens"
  if(method=="sens" && is.null(interactions)) method="1D-SA" else if(method=="sens") method="GSA"

  if(responses) YSTORE=TRUE else YSTORE=FALSE

  if(method=="DSA"||method=="MSA"||method=="CSA"||method=="GSA") { if(Aggregation<1 && substr(task,1,3)=="reg") Aggregation=3 else Aggregation=1 }

  #if(method=="sensi") INTERACT=TRUE else INTERACT=FALSE
  INTERACT=FALSE

  if(method=="DSA"||method=="MSA"){if(LRandom<1) Lsize=NROW(data) else Lsize=LRandom}
  Dsize=NCOL(data); 
  # set ML and NLY if needed
  if(is.factor(data[1,outindex])){Lev=levels(data[1,outindex]);NLY=length(Lev);ML=table(data[,outindex])[];ML=ML/sum(ML)}
  else{ML=NULL;NLY=0;}
  IRANGE=setdiff(Attributes,outindex)

  if(method!="CSA") # set SPREAD and other initializations
  { SPREAD=vector("list",Dsize); MR=0; # set SPREAD
    for(i in IRANGE) # 
    { L=length(levels(data[1,i])); if(L>0){if(!Lfactor) L=min(RealL,L)} else L=RealL
#cat(">>> Factor:",Lfactor,"L:",L,"Real:",RealL,"\n")
      SPREAD[[i]]=rmsample(data[,i],L=L,index=FALSE,sampling=sampling)$x
      MR=max(MR,L)
    } 
    if(!is.null(interactions)) # set LINT and JDOMAIN
    { LINT=length(interactions)
      i1=interactions[1]; L=length(SPREAD[[i1]]);
      if(LINT==2) {Dsize=1;JDOMAIN=interactions[2];}
      else if(LINT>2){Dsize=1; INTERACT=TRUE;}
      else JDOMAIN=setdiff(IRANGE,i1)
    }else LINT=0
    if(method=="1D-SA"||method=="GSA") # set the baseline if needed
    { 
      if(class(baseline)=="data.frame") v=baseline
      else #--- set the average/median input
      { v=data[1,]
        for(i in IRANGE) 
        {  
         if(is.factor(data[1,i]))  
         { if(is.ordered(data[1,i])) ModeClass=middleclass(data[,i],method=baseline) # "mean" or "median"
           else ModeClass=mostcommon(data[,i])
           v[1,i]=levels(data[1,i])[ModeClass]
         }
         else if(baseline=="mean") v[1,i]=mean(data[,i])
         else if(baseline=="median") v[1,i]=median(data[,i])
        } #end of for cycle 
      }
      # new data frame with MR rows, for speeding up the Importance computation:
      data=v[1,]
      if(LINT>0) 
      { L2=2;
         if(!INTERACT) {for(j in JDOMAIN) L2=max(L2,length(SPREAD[[j]]));MR=L*L2}
         else{ MR=1;for(j in 1:LINT) MR=MR*length(SPREAD[[interactions[j]]])
             }
      }
      data=v[rep(1,MR),]
    }
  } # ----------------------------------------------------------------------------
 Llevels=rep(NA,length(Attributes))
 if(is.null(interactions)) # method=="sens", "SA", "1D-SA" || method=="sensgrad" || method=="DSA" || method=="CSA" || method="MSA"
 {
  Sv=rep(0,Dsize)
  if(responses) RESP=vector("list",length=Dsize) else RESP=NULL
  YY=NULL
  LRDATA=nrow(data)
  DD=NULL
  if(method=="DSA"){if(Lsize<LRDATA){S=sample(1:LRDATA,Lsize);DD=data[S,]} else DD=data}
  else if(method=="MSA"){DD=data[1:Lsize,];DD=MCrandom(DD,IRANGE,SPREAD,mode=MRandom)}
  else if(method=="CSA"){ if(class(M)=="model") y=predict(M,data) else y=PRED(M,data)
                          if(is.factor(y[1]))   y=one_of_c(y)
                        }
  if(YSTORE){SDD=DD}
  if(method=="DSA"||method=="MSA"||method=="CSA") value=vector(length=Aggregation)

  for(i in IRANGE)
   { 
    if(method=="DSA"||method=="MSA")
      { L=length(SPREAD[[i]])
        if(responses) { xinp=SPREAD[[i]];xname=names(DD)[i];}
        MEM=DD[,i]
        if(NLY==0) Y=matrix(ncol=L,nrow=Aggregation) else Y=matrix(nrow=L,ncol=Aggregation*NLY)
        if(YSTORE) { z=1;if(NLY>0) YY=matrix(nrow=nrow(DD),ncol=(L*NLY)) else YY=matrix(nrow=nrow(DD),ncol=L)} # store all y values
        for(k in 1:L)
        {
           if(is.factor(DD[1,i])) DD[,i]=factor(rep(SPREAD[[i]][k],Lsize),levels=levels(DD[1,i]))
           else DD[,i]=rep(SPREAD[[i]][k],Lsize)
           if(class(M)=="model") y=predict(M,DD) else y=PRED(M,DD)
           if(is.factor(y[1])) y=one_of_c(y)
           if(YSTORE){if(NLY>0) {for(j in 1:NLY) {YY[,z]=y[,j];z=z+1;} } else {YY[,z]=y;z=z+1;}}
           if(NLY==0) Y[,k]=yaggregate(y,Aggregation)
           else  
           {  for(j in 1:NLY)
              { sss=seq(j,NLY*Aggregation,NLY)
                Y[k,sss]=yaggregate(y[,j],Aggregation)
              }
           }
        }
        DD[,i]=MEM # restore previous values
        if(NLY==0) { for(k in 1:Aggregation)  value[k]=s_measure(Y[k,],measure,x=SPREAD[[i]],Levels=ML) }
        else { for(k in 1:Aggregation)
                  { ini=(k-1)*NLY+1;end=ini+NLY-1
                    value[k]=s_measure(Y[,ini:end],measure,x=SPREAD[[i]],Levels=ML) 
                  }
             }
        Sv[i]=mean(value)
      }
      else if(method=="CSA")
      { L=length(levels(data[1,i])); if(L>0){if(!Lfactor)L=min(RealL,L)} else L=RealL
        IR=rmsample(data[,i],L=L,MAX=Inf,sampling=sampling) 
        if(responses) { xinp=IR$x; xname=names(data)[i];}
        IR=IR$index; K=1:length(IR); LK=length(K)
        if(NLY==0) Y=matrix(ncol=LK,nrow=Aggregation) else Y=matrix(nrow=LK,ncol=Aggregation*NLY)
        k1=1 
        for(k in K)
        {
          if(NLY==0) Y[,k1]=yaggregate(y[IR[[k]]],Aggregation)
          else
          {  
             for(j in 1:NLY)
             { sss=seq(j,NLY*Aggregation,NLY)
               Y[k1,sss]=yaggregate(y[IR[[k]],j],Aggregation)
             }
          }
          k1=k1+1
        } 
        if(NLY==0) { for(k in 1:Aggregation)  value[k]=s_measure(Y[k,],measure,x=xinp,Levels=ML) }
        else { for(k in 1:Aggregation)
                  { ini=(k-1)*NLY+1;end=ini+NLY-1
                    value[k]=s_measure(Y[,ini:end],measure,x=xinp,Levels=ML) 
                  }
             }
        Sv[i]=mean(value)
      }
      else # "sens"
      {   L=length(SPREAD[[i]])
          if(responses) {xinp=SPREAD[[i]];xname=names(data)[i];}
          MEM=data[(1:L),i] 
          data[(1:L),i]=SPREAD[[i]]
          if(class(M)=="model") Y=predict(M,data[(1:L),]) else Y=PRED(M,data[(1:L),])
          data[(1:L),i]=MEM # restore previous values
          Sv[i]=s_measure(Y,measure,x=SPREAD[[i]],Levels=ML)
      }
      if(responses) RESP[[i]]=list(n=xname,l=L,x=xinp,y=Y,yy=YY)
      Llevels[i]=L

   } # cycle for?
  Sum=sum(Sv)
  if(Sum==0) # no change in the model, equal importances to all attributes
  {imp=rep(1/(Dsize-1),length=Dsize);imp[outindex]=0;}
  else imp=Sv/Sum
  } #end of if null interactions --------------------------------------------------------------
  else if(!INTERACT) # LINT<3) # if(!is.null(interactions))
  {
   Sv=rep(0,Dsize)
   if(responses) RESP=vector("list",length=Dsize) else RESP=NULL
   if(method=="DSA") MR=NROW(data) 
   for(j in JDOMAIN)
    {
     if(method=="1D-SA"||method=="GSA")
     { MEM=data[,j] 
       L2=length(SPREAD[[j]]); MR=L*L2;
       for(k in 1:L2) {ini=1+(k-1)*L;end=ini+L-1;data[ini:end,i1]=SPREAD[[i1]];data[ini:end,j]=SPREAD[[j]][k];}
       if(responses) {xinp=data[1:MR,c(i1,j)];xname=c(names(data)[i1],names(data)[j]);}
       if(class(M)=="model") Y=predict(M,data[(1:MR),]) else Y=PRED(M,data[(1:MR),])
       data[,j]=MEM # restore previous values
     }
     else if(method=="DSA") 
     { L2=length(SPREAD[[j]]);Y=NULL
       MEM=data[,c(i1,j)];
       if(responses) {xname=c(names(data)[i1],names(data)[j]);xinp=data[1:(L*L2),c(i1,j)];o=1;}
       for(k in 1:L2) 
       for(l in 1:L) 
       { 
         if(responses){xinp[o,1]=SPREAD[[i1]][l];xinp[o,2]=SPREAD[[j]][k];o=o+1;}
         if(is.factor(data[1,i1])) data[,i1]=factor(SPREAD[[i1]][l],levels=levels(data[1,i1])) else data[,i1]=SPREAD[[i1]][l]
         if(is.factor(data[1,j])) data[,j]=factor(SPREAD[[j]][k],levels=levels(data[1,j])) else data[,j]=SPREAD[[j]][k]
         if(class(M)=="model") y=predict(M,data) else y=PRED(M,data)
         Y=c(Y,mean(y)) 
         #data[,i1]=MEM[,1];data[,j]=MEM[,2] # restore previous values
       }
       data[,i1]=MEM[,1];data[,j]=MEM[,2] # restore previous values
     }
     if(length(interactions)==2) k=1 else k=j
     Sv[k]=s_measure(Y,measure,x=SPREAD[[k]])
     if(responses) RESP[[k]]=list(n=xname,l=L,x=xinp,y=Y) 
     Llevels[k]=L
    }
#cat("SV:",Sv,"\n")
    Sum=sum(Sv)
    if(Sum==0) # no change in the model, equal importances to all attributes
    {imp=rep(1/(Dsize-1),length=Dsize);imp[outindex]=0;}
    else imp=Sv/Sum
  }
  else # LINT>=3, INTERACT
  {
   Sv=rep(0,1)
   if(responses) {RESP=vector("list",length=1);L=vector(length=LINT);} else RESP=NULL
   # interactions / SPREAD / data
   E=1;
#cat("LINT:",LINT,"MR:",MR,"\n")
   for(i in 1:LINT)
   { k=interactions[[i]];LS=length(SPREAD[[k]]);T=MR/(E*LS)
#cat("i:",i,"LS:",LS,"T:",T,"\n")
     if(is.factor(data[,k][1])) data[,k]=rep(factor(SPREAD[[k]],levels=levels(data[,k])),each=E,times=T) else data[,k]=rep(SPREAD[[k]],each=E,times=T)
     E=E*LS;
     if(responses) L[i]=LS
     Llevels[i]=L[i]
   }
   if(class(M)=="model") y=predict(M,data) else y=PRED(M,data)
#rmboxplot(y,MEAN=TRUE,LINE=2,MIN=TRUE,MAX=TRUE,ALL=TRUE)
#mypause()
   Sv=s_measure(y,measure)
   if(responses) RESP[[1]]=list(n=names(data)[interactions],l=L,x=data[,interactions],y=y) 
   imp=1
  }
 }
 return (list(value=Sv,imp=imp,sresponses=RESP,data=SDD,method=method,measure=measure,agg=Aggregation,nclasses=NLY,inputs=setdiff(Attributes,outindex),Llevels=Llevels,
              interactions=interactions))
}

# -- 
avg_imp_1D=function(I,measure="variance")
{
#II<=I
 x=I$sresponses[[1]]$x 
 NC=NCOL(x)
 
 R=list(responses=vector("list",length=NC),value=vector(length=NC),imp=rep(1/NC,NC))
 for(i in 1:ncol(x))
 {
  I1=avg_imp(I,i,measure=measure)
  R$sresponses[[i]]$x=I1$sresponses[[1]]$x
  y=I1$sresponses[[1]]$y
  R$sresponses[[i]]$y=y
#print(y)
  R$value[i]=s_measure(y,measure)
 }
#print("RVALUE:")
#print(R$value)
 SUM=sum(R$value)
 if(SUM>0) R$imp=R$value/SUM 
 return(R)
}

# I, AT=attribute
avg_imp=function(I,AT,measure="variance")
{
 x=I$sresponses[[1]]$x
 y=I$sresponses[[1]]$y
 X1=unique(x[,AT[1]])

 if(length(AT)>1)
 { 
   X2=unique(x[,AT[2]])
   LX=length(X2)*length(X1)
   if(is.matrix(y)) my=matrix(ncol=NCOL(y),nrow=LX)
   else if(is.factor(y)) my=factor(rep(levels(y)[1],LX),levels=levels(y))
   else my=vector(length=LX)
   Im=vector(length=LX)
   k=1;
   for(i in X1)
     for(j in X2)
     {
      W=which(x[,AT[1]]==i & x[,AT[2]]==j)
      Im[k]=W[1]
      if(is.matrix(y)) my[k,]=mean_resp(y[W,])
      else my[k]=mean_resp(y[W])
      k=k+1
     }
 }
 else # length(AT)==1
 { LX=length(X1)
   if(is.matrix(y)) my=matrix(ncol=NCOL(y),nrow=LX)
   else if(is.factor(y)) my=factor(rep(levels(y)[1],LX),levels=levels(y))
   else my=vector(length=LX)
   Im=vector(length=LX); k=1;
   for(i in X1)
     {
      W=which(x[,AT[1]]==i)
      Im[k]=W[1]
      if(is.matrix(y)) my[k,]=mean_resp(y[W,])
      else my[k]=mean_resp(y[W])
      k=k+1
     }
 }
 I$sresponses[[1]]$x=x[Im,AT]
 I$sresponses[[1]]$y=my
 I$sresponses[[1]]$n=I$sresponses[[1]]$n[AT]
 I$sresponses[[1]]$l=I$sresponses[[1]]$l[AT]
 I$value=s_measure(my,measure)
 return(I)
}

# experimental:
# need to check if using the output variable at the first column produces errors in other rminer functions:
agg_matrix_imp=function(I,INP=I$inputs,measure=I$measure,Aggregation=I$agg,method=I$method,outdata=NULL,L=I$Llevels,ameth="xy",Tolerance=0.1)
{
 N=length(INP)
 m1=matrix(0,ncol=N,nrow=N)
 m2=m1
 for(i in 1:N)
  for(j in 1:N)
   { if(INP[i]!=INP[j]) { 
                #cat("i:",i,"j:",j,"Li:",L[i],"Lj:",L[j],"\n")  
                II=aggregate_imp(I,AT=c(INP[i],INP[j]),measure=measure,Aggregation=Aggregation,method=method,outdata=outdata,L=L,ameth=ameth,Tolerance=Tolerance)
                m1[i,j]=mean(II$value); m2[i,j]=mean(II$value2)
              }
   }
 return(list(m1=m1,m2=m2))
}


# need to improve this function:
# I, AT=attribute
# -----
# remove outdata? use ML table in I$ML?
# remove avg_imp?
aggregate_imp=function(I,AT,measure="variance",Aggregation=1,method="sens",outdata=NULL,L,ameth="xy",Tolerance=0.1)
{
 if(length(AT)==1){
 x=I$sresponses[[ AT[1] ]]$x
 y=I$sresponses[[ AT[1] ]]$y
                  }
 else{
      if(substr(method,1,4)=="sens"||method=="GSA")
      {
       x=I$sresponses[[ 1 ]]$x
       y=I$sresponses[[ 1 ]]$y
      }
      else{
            x=I$sresponses[[ AT[1] ]]$x
            y=I$sresponses[[ AT[1] ]]$y # check if this is right???
          }
     }

 if(is.factor(outdata[1])){ Lev=levels(outdata[1]);NLY=length(Lev);ML=table(outdata)[];ML=ML/sum(ML)}
 else {ML=NULL;NLY=0;}

 if(length(AT)==2)
 { 
   if(substr(method,1,4)=="sens"||method=="GSA"||method=="1D-SA"||method=="SA")
   {  
      X1=unique(x[,AT[1]])
      LX=length(X1)
      xlab=I$sresponses[[1]]$n[AT[1]]
      ylab=I$sresponses[[1]]$n[AT[2]]
      X2=unique(x[,AT[2]])
      LY=length(X2)
      # nao sei o que fazer aqui, testar mais tarde...
      if(NLY==0) J=Aggregation else J=NLY
      mm=matrix(nrow=LX,ncol=LY*J)
      for(i in 1:LX)
       for(j in 1:LY)
         {
           W=which(x[,AT[1]]==X1[i] & x[,AT[2]]==X2[j])
           sss=seq(j,(LY*(J-1)+j),LY)
           mm[i,sss]=yaggregate(y[W],Aggregation)
         }
   }
   else # check LY e LX instead of AT 1 and 2? "MSA" / "DSA"
   { 
      X1=unique(x)
      LX=length(X1)
      xlab=I$sresponses[[AT[1]]]$n
      ylab=I$sresponses[[AT[2]]]$n
#cat("kmsample> xlab:",xlab,"ylab:",ylab,":\n")
      K=rmsample(I$data[,AT[2]],L=L[AT[2]],index=TRUE,sampling="regular",Tolerance=0.1) # think about zero slots?
      X2=K$x
      LY=length(X2)
#cat("X2:",X2,"\n")
     
      yy=I$sresponses[[AT[1]]]$yy
      if(NLY==0) J=Aggregation else J=NLY
      mm=matrix(nrow=LX,ncol=LY*J)
     for(i in 1:LX)
      for(j in 1:LY)
       {
         sss=seq(j,(LY*(J-1)+j),LY)
#yy<=yy
#K<=K
         mm[i,sss]=yaggregate(yy[K$index[[j]],i],Aggregation)
       } 
   }  
#cat("done kmsample\n")
#mm<=mm
#NLY<=NLY
#ML<=ML
#ameth<=ameth;measure<=measure
   # s_measure global OU em x + em y ?
   value=vector(length=Aggregation)
   value2=vector(length=Aggregation)
value3=vector(length=Aggregation)
    if(NLY==0) { for(k in 1:Aggregation) { ini=(k-1)*LY+1;end=ini+LY-1;
# switch value
#cat("ini:",ini,"end:",end,"measure:",measure,"xy:",ameth,"\n")
                                  val=s_measure(mm[,ini:end],measure,Levels=ML,xy=ameth) 
value3[k]=s_measure(mm[,ini:end],measure,Levels=ML,xy="global") 
                                  value[k]=val[2];value2[k]=val[1]
                                }
#cat("done\n")
               }
   # think better about the class case:
   else{ for(k in 1:Aggregation)
                  { 
#cat(">> k:",k,"\n")
                    ini=(k-1)*NLY+1;end=ini+NLY*(L[AT[2]])-1;
                    val=s_measure(mm[,ini:end],measure,Levels=ML,xy=ameth) 
                    value[k]=val[1];value2[k]=val[2]
                  }
       }
   #I$sresponses[[1]]$x=x[Im,AT] # ???
   #I$sresponses[[1]]$y=mm # ????
   #I$sresponses[[1]]$n=I$sresponses[[1]]$n[AT]
   #I$sresponses[[1]]$l=I$sresponses[[1]]$l[AT]

   # new importance object:
   I=list()
   I$value=value
   I$value2=value2
#for(i in 1:3)cat(round(value2[i],digits=2),",",round(value[i],digits=2),"\n")
#for(i in 1:3)cat(round(value3[i],digits=2),"\n")
#cat("val: ",round(mean(I$value),digits=2),",",round(mean(I$value2),digits=2),"\n")
   I$Y=mm
   I$x1=X1
   I$x2=X2
   I$agg=Aggregation
   I$NLY=NLY
   I$xy=ameth
   I$xlab=xlab
   I$ylab=ylab
#cat("L:",L,"\n")
   I$L=c(L[AT[1]],LY)
   return(I)
 }
 else # length(AT)==1 ## there are huge problems with this if code, which does not work, later i need to rethink and recode this!!!
 { 
   if(is.numeric(x)) X1=unique(x) else X1=unique(x[,AT[1]])
   #X1=unique(x)
   LX=length(X1)
   Im=vector(length=LX); k=1;
   Min=vector(length=LX)
   Max=vector(length=LX)
   z=1

  if(I$method=="GSA") 
  {
   LY=LX
   mm=vector(length=LX)
   for(i in X1)
     {
      W=which(x[,AT[1]]==i)
      Im[k]=W[1]
      if(NLY==0) mm[k]=mean(y[W])
      else{ for(j in 1:NLY) {mm[z]=mean(y[W,j]);z=z+1;}
          }
      k=k+1
     }
    value=1
    I$method="1D-SA"
  }
  else
  {
   LY=round(length(y)/LX)
   if(NLY==0) mm=matrix(ncol=LX,nrow=LY)
   else mm=matrix(ncol=LX*NLY,nrow=LY)
   for(i in X1)
     {
      W=which(x[,AT[1]]==i)
      Im[k]=W[1]
      if(NLY==0) mm[,k]=y[W]
      else{ for(j in 1:NLY) {mm[,z]=y[W,j];z=z+1;}
          }
      k=k+1
     }
   if(NLY==0) Y=matrix(ncol=LX,nrow=Aggregation) else Y=matrix(nrow=LX,ncol=Aggregation*NLY)
   z=1
   for(i in 1:LX)
   {
    if(NLY==0) Y[,i]=yaggregate(mm[,i],Aggregation)
    else
    {
      for(j in 1:NLY)
               {
                sss=seq(j,NLY*Aggregation,NLY)
                Y[i,sss]=yaggregate(mm[,z],Aggregation)
                z=z+1
               }
    }
   }
  
   value=vector(length=Aggregation)
   if(NLY==0) { for(k in 1:Aggregation)  value[k]=s_measure(Y[k,],measure,Levels=ML) }
   else{ for(k in 1:Aggregation)
                  { 
                    ini=(k-1)*NLY+1;end=ini+NLY-1
                    value[k]=s_measure(Y[,ini:end],measure,Levels=ML) 
                  }
       }
   }
}
 I$sresponses[[1]]$x=x[Im,AT]
 I$sresponses[[1]]$y=mm
 I$sresponses[[1]]$n=I$sresponses[[1]]$n[AT]
 I$sresponses[[1]]$l=I$sresponses[[1]]$l[AT]
# avalue=AAD_responses(value)
#avalue=switch(measure,AAD=AAD_responses(value),variance=variance_responses(value),gradient=gradient_responses(value),range=range_responses(value))
#cat("min:", Min,"\n")
#cat("max:", Max,"\n")
#cat("area:",curvearea(cbind(Min,Max)),"minarea:",((min(Max)-min(Min))*LX),"\n")
 #I$value=area; #-mean(Min)
 #I$value=avalue;
 I$value=mean(value)
 return(I)
}




#----------------------------------------------------------------------------------------------------
#-- Auxiliary private functions, in principle you should not need to use these:
#----------------------------------------------------------------------------------------------------
mean_resp=function(y)
{
 if(is.matrix(y)) return (colMeans(y))
 else if(is.ordered(y)) return (levels(y)[middleclass(y,method="mean")])
 else if(is.factor(y)) return (levels(y)[mostcommon(y)])
 else return (mean(y))
}

# xy = if there is a sensitivity pair and y is a matrix related with regression
s_measure=function(y,measure,x=1:length(y),ABS=TRUE,center="median",Levels=NULL,xy=FALSE)
{
 if(xy=="xy")  { rx=0;N=NROW(y);for(i in 1:N) rx=rx+s_measure(y[i,],measure=measure,x=x,ABS=ABS,Levels=NULL,center=center);rx=rx/N;
                 ry=0;N=NCOL(y);for(i in 1:N) ry=ry+s_measure(y[,i],measure=measure,x=x,ABS=ABS,Levels=NULL,center=center);ry=ry/N;
                 return(c(rx,ry))
               }
 else{ 
      if(xy=="global"){dim(y)=NULL} # convert to vector
      return(switch(measure,range=range_responses(y,Levels=Levels),gradient=gradient_responses(y,ABS=ABS,Levels=Levels),variance=variance_responses(y,Levels),
                            AAD=AAD_responses(y,center=center,Levels=Levels),balanced=balanced_responses(y))) 
     }
}

range_responses=function(y,Levels=NULL)
{
  if(is.ordered(y)) return(range_responses(as.numeric(y)))
  else if(is.character(Levels)) return(sum(y>0)) # y is table of frequencies
  else if(is.numeric(Levels) && is.factor(y)) return(range_responses(one_of_c(y),Levels=Levels)) 
  else if(is.factor(y)) return((length(unique(y))-1)) # future: test in clustering?:)
  else # matrix of probabilities
  { LV=NCOL(y) 
    if(LV==1) return (abs(diff(range(y))))
    else{
         res=0; if(is.null(Levels)) Levels=rep(1/LV,LV)
         for(i in 1:LV) res=res+Levels[i]*abs(diff(range(y[,i]))); 
         return(res) 
        }
  }
}

# deprecated....
#-- compute the gradient of response y: vector, matrix of numeric data (probabilities) or factor
#-- ordered is only used if y is factor, true if it is an ordered factor or else false
# XXX think here!
# Euclidean distance: http://en.wikipedia.org/wiki/Magnitude_(mathematics)
pathlength_responses=function(y,x=1:length(y))
{
x=1:length(y)
 if(is.factor(x)) { x=1:length(y) } # change this later to include quantile analysis? 
 if(is.ordered(y)) return(pathlength_responses(as.numeric(y),x))
 else if(is.factor(y)) return(pathlength_responses(one_of_c(y),x))
 else{ dist=0; LY=NROW(y); CY=NCOL(y)
       for(i in 2:LY)
       { dy=(x[i]-x[i-1])^2 
         if(CY>1) { for(j in 1:CY) dy=dy+(y[i,j]-y[i-1,j])^2 }
         else dy=dy+(y[i]-y[i-1])^2 
         dist=dist+sqrt(dy)
       }
     }
 return (dist-(x[LY]-x[1]))
}


gradient_responses=function(y,ABS=TRUE,Levels=NULL)
{
 if(is.ordered(y)) return(gradient_responses(as.numeric(y)))
 else if(is.character(Levels)) return(sum(abs(diff(y)))/(length(y)-1)) # y table of frequencies, not adequately defined?
 else if(is.numeric(Levels) && is.factor(y)) return(gradient_responses(one_of_c(y),ABS=ABS,Levels=Levels)) 
 else if(is.factor(y)) 
 { G=0;LY=length(y)
   for(i in 2:LY)
   { if(y[i]!=y[i-1]) G=G+1
   } 
   return (G/(LY-1))
 }
 else{ if(ABS) FUNC=abs else FUNC=identity
       LV=NCOL(y) 
       if(LV==1) return (mean(FUNC(diff(y))))
       else{
            if(is.null(Levels)) Levels=rep(1/LV,LV)
            res=0
            for(i in 1:LV) res=res+Levels[i]*abs(diff(range(y[,i]))); 
            return (res)
           }
     }
}
 
# compute the variance of response y (vector or matrix or data.frame)
variance_responses=function(y,Levels=NULL)
{ 
 if(is.ordered(y)) return(variance_responses(as.numeric(y)))
 else if(is.factor(y)) return(variance_responses(one_of_c(y),Levels=Levels))
 else
  { LV=NCOL(y)
    if(LV==1) return (var(y))
    else # LV>2 
    {
      if(is.null(Levels)) Levels=rep(1/LV,LV)
      res=0; for(i in 1:LV) res=res+Levels[i]*var(y[,i]) 
      return(res) 
    }
  }
}

AAD_responses=function(y,x,center="median",Levels=NULL)
{ 
 if(is.ordered(y)) return(AAD_responses(as.numeric(y)))
 else if(is.character(Levels)) {t=y;L=sum(t);C=length(t);N=c(L,rep(0,C-1));return( mean(abs(t-N)) ) } # y table of frequencies, not adequately defined?
 else if(is.numeric(Levels) && is.factor(y)) return(AAD_responses(one_of_c(y),center=center,Levels=Levels)) 
 else if(is.factor(y)) 
 { t=sort(table(y)[],decreasing=TRUE); L=length(y);C=length(levels(y[1]))
   N=c(L,rep(0,C-1))
   return( mean(abs(t-N)) )
 }
 else{ # matrix of probabilities
      LV=NCOL(y)
      if(LV==1) {
                 center=switch(center,mode=mostcommon(y),mean=mean(y),median=,median(y))
                 return(mean(abs(y-center)))
                } 
      else{
#print(Levels)
#print("---")
#print(y)
#print("===")
      if(is.null(Levels)) Levels=rep(1/LV,LV)
      res=0
      for(i in 1:LV)
        {
         center2=switch(center,mode=mostcommon(y[,i]),mean=mean(y[,i]),median=,median(y[,i]))
         res=res+Levels[i]*mean(abs(y[,i]-center2))
#cat("i:",i,"lev:",Levels[i],"aad:",mean(abs(y[,i]-center2)),"AAD:",Levels[i]*mean(abs(y[,i]-center2)),"\n")
        }
      return(res)
     }
    }
}

# experimental, only for factor data!
balanced_responses=function(y)
{
 t=sort(table(y)[],decreasing=TRUE)
 C=length(levels(y[1]))
 L=length(y);
 S=sum(abs(diff(t)))
 R=1*(L%%C>0);D=L-R
 return((L-S)/D)
 #return(L-mean(abs(diff(t))))
}

# --- auxiliary functions, do not use directly:
# x is a $sresponses object
# problem if task="class", y is factor and not a matrix, need to adapt rminer for dealing with factors,
# which means changing several functions, including in plots.R, metrics.R (vaveraging), etc...
resp_to_list=function(x,TC=-1)
{
 LX=length(x$x);lines=LX/x$l;
 RES=vector("list",lines)
 for(i in 1:lines)
    { 
      ini=(i-1)*x$l+1;end=ini+x$l-1
#      if(is.factor(x$y)
#      {
#       if(TC==-1) TC=length(
#       M=cbind(x$x[ini:end],x$y[ini:end])
#      }
#      else
#      { }
       #if(TC==-1) M=data.frame(x=x$x[ini:end],y=x$y[ini:end])
       if(TC==-1) M=cbind(x=x$x[ini:end],y=x$y[ini:end])
       else M=cbind(x$x[ini:end],x$y[ini:end,TC])
#      }
      RES[[i]]=M   
    }
 return (RES)
}
# ---
mids=function(x,L=2,R=range(x))
{
 LR=R[2]-R[1];Step=LR/L;Start=R[1]+LR/(2*L)
 return (Start+(0:(L-1)*Step))
}
ranges=function(x,L=2)
{R=range(x)
 return (seq(R[1],R[2],length.out=L+1))
}

# auxiliary function
# x = factor, L= levels
rmtable=function(x,Levels)
{
 t=table(x)
 L=length(levels)
 res=vector(length=L)
 for(i in 1:length(Levels))
 {
  res[i]=t[][which(Levels[i]==names(t))] 
 }
 return(res)
}

# adapt to factor data?
# iregular and equal may not work for ordered data?
# need to correct this function for: X=c(1,2.1,3.3,4.0,7.0);rmsample(X,L=5,Tolerance=0.0)
rmsample=function(x,L=2,MAX=Inf,Tolerance=1,sampling="regular",index=TRUE)
{
 #MAX=Inf;Tolerance=0.1;sampling="regular";index=TRUE;L=4;x=c(1,2.1,3.3,4.0,7.0)
 IL=NULL
 # LU=length(unique(x)); if(LU<L) L=LU ; # add this line?
 #if(is.ordered(x)) { return (rmsample(as.numeric(x),L=L,MAX=MAX,sampling=sampling,index=index)) } # add this line?
 if(!is.ordered(x) && is.factor(x))
                 { t=table(x); sx=sort.int(t[],decreasing=TRUE,index.return=TRUE);
                   MID=names(t);LI=length(MID)
                   if(LI>L) MID=MID[ - sx$ix[(L+1):LI] ] 
                   if(index){
                   IL=vector("list",L)
                   for(i in 1:L)
                   {
                    range=which(x==MID[i])
                    if(length(range)>MAX) {range=range[1:MAX]}
                    IL[[i]]=range
                   }        }
                 }
 else
 {
 #if(is.factor(x)){OL=levels(x); x=as.numeric(x);} else OL=NULL
 if(sampling=="equal") { 
                            #L=5; MAX=1;
                            sx=sort.int(x,index.return=TRUE)
                            LX=length(sx$ix)
                            R=trunc(LX/L)
                            if(index) IL=vector("list",L)
                            MID=rep(NA,L)
                            if(R>MAX) RL=MAX else RL=R
                            for(i in 1:(L-1)) 
                              { 
                                # compute mid first, choose neighbours around mid!
                                ii=1+(i-1)*R;ie=ii+R-1
                                range=sx$ix[ii:ie]
                                if(R%%2==1) MID[i]=x[ range[(R-1)/2+1] ]
                                else {R1=R/2;MID[i]=(x[range[R1]]+x[range[R1+1]])/2}
                                if(index){ if(R>MAX){ x1=sort(abs(x[range]-MID[i]),index.return=TRUE);range=range[x1$ix[1:RL]] }
                                           IL[[i]]=range
                                         }
                              }
                            ii=1+(L-1)*R;ie=LX;R=LX-(L-1)*R
                            range=sx$ix[ii:ie]
                            if(R%%2==1) MID[L]=x[range[(R-1)/2+1] ]
                            else {R1=R/2;MID[L]=(x[range[R1]]+x[range[R1+1]])/2}
                            if(index){ if(R>MAX){ x1=sort(abs(x[range]-MID[L]),index.return=TRUE); range=range[x1$ix[1:RL]] }
                                       IL[[L]]=range
                                     }
                        }
 else if(sampling=="iregular" || sampling=="regular" || sampling=="quantile"){

 #x=fo; L=3; sampling="regular"

 if(is.ordered(x))
 { LEV=levels(x[1]);LLEV=length(LEV)
   if(L>LLEV) {MID=LEV;L=LLEV}
   else if(sampling=="regular") { NAUX=1:LLEV; MID=LEV[unique(round(seq(NAUX[1],NAUX[LLEV],length.out=L)))] }
   else if(sampling=="quantile"){ MID=LEV[unique(quantile(as.numeric(x),seq(0,1,length.out=L)))] }
   if(index) 
   { IL=vector("list",L); for(i in 1:L) {range=which(x==MID[i]);LR=length(range);if(LR>MAX) range=range[1:MAX];IL[[i]]=range} }
 }
 else{
 R=range(x)
 if(sampling=="iregular"){ MID=mids(x,L,R=R)}
 else if(sampling=="quantile"){ MID=unique(quantile(x,seq(0,1,length.out=L)));attr(MID,"names")=NULL}
 else if(sampling=="regular"){ MID=seq(R[1],R[2],length.out=L)}

#print("MID:")
#print(MID)
#xxx<=x
 #index=TRUE;sampling="regular";L=3;
 if(index){
 S=seq(R[1],R[2],length.out=(L+1))
#print("S:")
#print(S)
 IL=vector("list",L);Iaux=1:length(x)
 ki=1
 if(Tolerance<1) TTolerance=Tolerance*(MID[2]-MID[1])
 for(i in 1:L)
 { 
  I=Iaux[which(x[Iaux]<=S[i+1])]
  Iaux=setdiff(Iaux,I)
  LI=length(I)
  if(LI>MAX || Tolerance<1) 
             { x1=abs(x[I]-MID[ki]) 
               if(LI>MAX) {sx1=sort.int(x1,index.return=TRUE);
                           I=I[sx1$ix[1:MAX]] # check?
                          }
               else { sx1=which(x1<=TTolerance)
                      if(length(sx1)>0)I=I[sx1]else{I=NULL;LI=length(I);}
                    }
               #cat("i:",i,"tol:",TTolerance,"lsx1:",length(sx1),"sx1:",sx1,"x1:",x1,"\n")
               #mypause()
             }
  if(LI>0){IL[[ki]]=I;ki=ki+1;}else{IL[[ki]]=NULL;MID=MID[-ki]}
 }
          }
 }
 }
 else if(sampling=="kmeans"){ # || sampling=="pam"){

 if(is.ordered(x))
 { LEV=levels(x[1]);LLEV=length(LEV)
   if(L>LLEV) {MID=LEV;L=LLEV}
   else { 
         if(sampling=="kmeans") {K=kmeans(as.numeric(x),centers=L);SK=sort.int(round(K$centers[,1]),index.return=TRUE)}
         #else if(sampling=="pam") {K=pam(as.numeric(x),k=L);SK=sort.int(K$medoids[,1],index.return=TRUE)}
         MID=LEV[as.numeric(SK$x)] 
        }
   if(index) 
   { IL=vector("list",L); for(i in 1:L) {range=which(K$cluster==SK$ix[i]);LR=length(range);if(LR>MAX) range=range[1:MAX];IL[[i]]=range} }
 }
 else{
 if(sampling=="kmeans"){K=kmeans(x,centers=L);SK=sort.int(K$centers[,1],index.return=TRUE); MID=as.numeric(SK$x) }
 #else { K=pam(x,k=L); SK=sort.int(K$medoids[,1],index.return=TRUE); MID=as.numeric(SK$x) } # pam
 if(index){
 IL=vector("list",L)
 for(i in 1:L)
 { 
  IL[[i]]= which(K$cluster==SK$ix[i])
  LI=length(IL[[i]])
  if(LI>MAX) { x1=abs(x[ IL[[i]] ]-MID[i]) 
               sx1=sort.int(x1,index.return=TRUE)
               IL[[i]]=IL[[i]][sx1$ix[1:MAX]] # check?
             }
 }
          }
 }
 }
 }
 #if(!is.null(OL)) { cat("MID:",MID,"\n")
 #                   MID=OL[MID] }
 return (list(index=IL,x=MID))
}

# mode=="discrete" or "continuous"
#MCrandom=function(d,IRANGE,SPREAD,mode="continous") # mode="discrete")
MCrandom=function(d,IRANGE,SPREAD,mode="discrete")
{
 LR=nrow(d)
 {
  for(i in IRANGE)
  {
   if(is.factor(d[1,i])) d[,i]=sample(factor(SPREAD[[i]],levels=levels(d[1,i])),LR,replace=TRUE) 
   else{
        if(mode=="discrete") d[,i]=sample(SPREAD[[i]],LR,replace=TRUE) 
        else d[,i]=runif(LR,SPREAD[[i]][1],SPREAD[[i]][length(SPREAD[[i]])])
       }
  }
 }
 return(d)
}
# --- end of auxiliary functions --------------

