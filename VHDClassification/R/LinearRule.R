
# Author: robin
###############################################################################


setClass(
		Class="LinearRule",
		representation=representation(
				normalVector="numeric",
                normalIndex="integer",
				centerVector="numeric",
                proportions="numeric",
                prior="logical",
				labels="factor"
		)
)

.learnLinearRule<-function(x1,x2,labels=factor(c(1,0)),procedure='noThresh',covariance='diagonal',ql=1*(procedure=='Fisher')+0.05,prior=FALSE)
{
    # dim(xi) is p x ni  
    # Procedure can be 'Fisher', covariance can be 'full' or 'diagonal'
    #Testing if the data are ok

    procedure=as.character(procedure)
    #declaration of variables
    p=length(x1[1,]); n1=length(x1[,1]); n2=length(x2[,1]);  n=n1+n2;
    center=array(0,p);
    hatmu=array(0,c(2,p));
    hatC=array(0,p)
    s=array(0,p)
    F=array(0,p) 
    # start by learning covariance and group mean vectors
    hatmu[1,]= colMeans(x1);  
    hatmu[2,]= colMeans(x2);
    hatC=1/(n1+n2-1)*(colSums((x1-hatmu[1,])^2)+colSums((x2-hatmu[2,])^2))
    indices=integer(0);
    F=(hatmu[1,]-hatmu[2,])/(hatC)^(1/2);
    s=(hatmu[1,]+hatmu[2,])/2;
    # Comput the normal vectors depending on the used procedure.
    switch(procedure,
            'noThresh'={indices=1:p},
            'FANThresh'={
                indices=.FAN_Selection(abs(F),hatC^(-1),n);
                F[!1:p%in%indices]=0;
            },  
            'UnivThresh'={
                lambda=array(as.numeric(sqrt(2*as.numeric(ql)*log(p)/n)),p);
                indices=1:p
                indices=indices[abs(F)>lambda]
                F[!1:p%in%indices]=0;
            },
            'FDRThresh'={
                indices=.FDR_Selection(sqrt(n)*abs(F),q=as.numeric(ql))
                F[!1:p%in%indices]=0;
            },
            'FDRstudent'={
                indices=.FDRS_Selection(sqrt(n)*abs(F),q=as.numeric(ql),n)
                F[!1:p%in%indices]=0;
            } ,
            'Otherwise'={stop('Trying to use a non-implemented procedure')}
    )
    F=F/(hatC)^(1/2);
    return(new(Class="LinearRule",normalVector=F,centerVector=s,
                    normalIndex=as.integer(indices),
                    proportions=as.numeric(n1/(n1+n2)),
                    prior=prior,labels=labels))
}

.learnLinearRulefortune<-function(x,y,...){
    y=ordered(y)
    classes=levels(y)
    return(.learnLinearRule(x[y==classes[1],],x[y==classes[2],],labels=factor(levels(y)),...))
}


.LDA<-function(x,F,s,p=1/2)
	return((sum(F*(x-s))>log((1-p)/p)))


setMethod(
        f="plotClassifRule",
        signature=c("LinearRule"),
        definition=function(object,...){
            Index=object@normalIndex
            NormalVector=object@normalVector[Index]
            Center=object@centerVector[Index]
            Mydata=data.frame(NormalVector=NormalVector,Center=Center,Index=Index)
            A=TRUE;
            return(xyplot(NormalVector+Center~Index,data=Mydata,auto.key = A,...))
        }
)

setMethod(
        f="show",
        signature="LinearRule",
        definition=function(object){
                    cat('Normal Vector: \n')
            print(object@normalVector)
                    cat('Center Vector: \n')
            print(object@centerVector)
        }
)
setMethod(
        f="getLogLikeRatio",
        signature="LinearRule",
        definition=function(object){
            return(list(NormalVector=object@normalVector,VecterVector=object@centerVector))
        }
)
setMethod(
        f="predict",
        signature="LinearRule",
        definition=function(object,newdata){
            if (!(is(newdata, "matrix")||is(newdata, "data.frame")))  
                stop("newdata must be a matrix or a data frame")
            p=length(object@normalVector);
            n=nrow(newdata);
            if (object@prior)
                res=apply(newdata,1,.LDA,object@normalVector,object@centerVector,object@proportions)
            else
                res=apply(newdata,1,.LDA,object@normalVector,object@centerVector)
            
            y=factor(object@labels)
            resultat=res
            resultat[res==1]=levels(y)[1]
            resultat[res==0]=levels(y)[2]
            return(resultat)
        }
)

#setMethod(f="getNormalVector",signature='LinearRule',
#		definition=function(object){return(object@normalVector)}
#)
#
#setMethod(f="getCenterVector",signature='LinearRule',
#		definition=function(object){return(object@centerVector)}
#)

setMethod(f=".EvaluateLogLikeRatio",signature=c(x='numeric',object='LinearRule'),
        definition=function(x,object){
            return(sum(object@normalVector*(x-object@centerVector)))        
        }
)

setMethod(f=".minus",signature=c(object='LinearRule'),
        definition=function(object){
                    return(new(Class="LinearRule",normalVector=-object@normalVector,
                                    centerVector=object@centerVector,
                                   normalIndex=object@normalIndex,
                                   proportions=object@proportions,
                                   prior=object@prior))
        }
)

.tune.LDA<-function(x,y,procedure='FDRThresh',ql=10^(-1:-7),prior=prior,...)
{ #### for internal use only
  #### y is an array of 2 factors
    y=ordered(y)
    call<-match.call()
    ranges<-list(ql=ql,prior=prior,procedure=procedure)
	ranges[sapply(ranges,is.null)]<-NULL
	tunecontrol=tune.control(sampling='cross',best.model=FALSE)
    if(length(ranges)<1) ranges=NULL
	modeltmp<-tune('.learnLinearRulefortune',train.x=x,train.y=y,ranges=ranges,
                 tunecontrol=tunecontrol,predict.func=predict,...)
    besti=length(ranges$ql)+1-which.min(rev(modeltmp$performances[['error']]))
    return(.learnLinearRulefortune(x,y,procedure=procedure,ql=modeltmp$performances$ql[besti],prior=prior))
}