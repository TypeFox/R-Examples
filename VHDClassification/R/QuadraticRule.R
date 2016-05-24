# TODO: Add comment
# 
# Author: robingirard
###############################################################################



setClass(
		Class="QuadraticRule",
		representation=representation(
				formVector="numeric",
                formIndex="integer",
				constant="numeric"
		),
		contains="LinearRule"
)

.learnQuadraticRule<-function(x1,x2,labels=factor(c(1,0)),procedure='Fisher',covariance='diagonal',ql=1*(procedure=='Fisher')+0.05,qq=1*(procedure=='Fisher')+0.05,prior=FALSE)
{
    # dim(x) is p x n 
	# Procedure can be 'Fisher', covariance can be 'full' or 'diagonal'
	#Testing if the data are ok
	procedure=as.character(procedure)
	#declaration of variables
	p=length(x1[1,]); n1=length(x1[,1]); n2=length(x2[,1]);  
    center=array(0,p);
    hatmu=array(0,c(2,p));
    hatC=array(0,c(2,p))
    s=array(0,p)
    F=array(0,p) 
    A=array(0,p)
    # start by learning covariance and group mean vectors
    
    hatmu[1,]= colMeans(x1);  
    hatmu[2,]= colMeans(x2);
    hatC[1,]=1/(n1-1)*colSums((x1-hatmu[1,])^2)
    hatC[2,]=1/(n2-1)*colSums((x2-hatmu[2,])^2)
    indices=integer(0);
    nindices=integer(0);
    F=(hatmu[1,]-hatmu[2,])*(1/hatC[1,]+1/hatC[2,])^(1/2);
    s=(hatmu[1,]+hatmu[2,])/2;
	# start by learning covariance and group mean vectors

	TF=F/((1+hatC[1,]/hatC[2,])/n1+(1+hatC[2,]/hatC[1,])/n2)^(1/2)
	# Comput the normal vectors depending on the used procedure.

	switch(procedure,
            'noThresh'={nindices=1:p},
			'UnivThresh'={
				lambda=as.numeric(sqrt(2*as.numeric(ql)*log(p)));
				tmp=1:p
				nindices=tmp[(abs(TF)>lambda)]
				F[!1:p%in%nindices]=0;
			},
			'FDRThresh'={
				nindices=.FDR_Selection(abs(TF),q=as.numeric(ql))
				F[!1:p%in%nindices]=0;
			},
            'Otherwise'={stop('Trying to use a non-implemented procedure')}
	)
	F=F/2*(1/hatC[1,]+1/hatC[2,])^(1/2);
	A=1/hatC[1,]-1/hatC[2,]
	Tw=(hatC[1,]-hatC[2,])/(2*hatC[1,]^2/(n1-1)+2*hatC[2,]^2/(n2-1))^(1/2)
	# Comput the form vectors depending on the used procedure.

	switch(procedure,
            'noThresh'={indices=1:p},
			'UnivThresh'={
				lambda=as.numeric(sqrt(2*qq*log(p)));
				tmp=1:p
				indices=tmp[(abs(Tw)>lambda)]
				A[!1:p%in%indices]=0;
			},
			'FDRThresh'={
				indices=.FDR_Selection(abs(Tw),q=as.numeric(qq))
				A[!1:p%in%indices]=0;
			},
            'Otherwise'={stop('Trying to use a non-implemented procedure')}
    )
    lambda=as.numeric(sqrt(2*log(p)));
    if (length(indices)==1)
    {
        if ((abs(Tw[indices])<lambda))
        {
            A[indices]=0;
            indices=integer(0);
        }
    }
	if (length(indices)>0)
    {
         c=sum(1/8*A[indices]*(hatmu[1,indices]-hatmu[2,indices])^2)+1/2*sum(log(hatC[1,indices]/hatC[2,indices]))
	}
    else{c=0;}
         return(new(Class="QuadraticRule",normalVector=F,
                         centerVector=s,
                         constant=c,
                         formVector=A,
                         normalIndex=as.integer(nindices),
                         formIndex=as.integer(indices),
                         proportions=as.numeric(n1/(n1+n2)),
                         prior=prior,
						 labels=labels))
}
.learnQuadraticRulefortune<-function(x,y,...){
	y=ordered(y)
    classes=levels(y)
    return(.learnQuadraticRule(x[y==classes[1],],x[y==classes[2],],labels=factor(levels(y)),...))
}

setMethod(
		f="predict",
		signature="QuadraticRule",
		definition=function(object,newdata){
			
			if (!(is(newdata, "matrix")||is(newdata, "data.frame")))  
				stop("newdata must be a matrix or a data frame")
            p=length(object@normalVector);
            n=nrow(newdata);
            if (object@prior)
                res=apply(newdata,1,.QDA,object@normalVector,object@centerVector,object@formVector,object@constant,object@proportions)
            else
                res=apply(newdata,1,.QDA,object@normalVector,object@centerVector,object@formVector,object@constant,1/2)
			y=object@labels;
            resultat=res
            resultat[res==1]=levels(y)[1]
            resultat[res==0]=levels(y)[2]
            return(resultat)
		}
)

setMethod(
        f="plotClassifRule",
        signature=c("QuadraticRule"),
        definition=function(object,...){
            Index=union(object@normalIndex,object@formIndex)
            NormalVector=object@normalVector[Index]
            Center=object@centerVector[Index]
            FormVector=object@formVector[Index]
            Mydata=data.frame(NormalVector=NormalVector,FormVector=FormVector,
                    Center=Center,Index=Index)
            return(xyplot(NormalVector+Center+FormVector~Index,data=Mydata,auto.key = TRUE,...))
        }
)

.QDA<-function(x,G,s,A,c,p=1/2)
{
    res=sum(-1/2*A*(x-s)^2+G*(x-s))-c-log((1-p)/p)
    #cat('s=',sum(-1/2*A*(x-s)^2+G*(x-s)),' c=',c,' s-c=',res, ' \n')
    return((res>0))
}


setMethod(f=".EvaluateLogLikeRatio",signature=c(x='numeric',object='QuadraticRule'),
        definition=function(x,object){
            return(-object@constant+sum(-1/2*object@formVector*(x-object@centerVector)^2+object@normalVector*(x-object@centerVector)))        
        }
)
setMethod(f=".minus",signature=c(object='QuadraticRule'),
        definition=function(object){
        return(new(Class="QuadraticRule",normalVector=-object@normalVector,
                         centerVector=object@centerVector,
                         constant=-object@constant,
                         formVector=-object@formVector,
                         normalIndex=object@normalIndex,
                         formIndex=object@formIndex,
                         proportions=object@proportions,
                         prior=object@prior))
        }
)
setMethod(
        f="show",
        signature="QuadraticRule",
        definition=function(object){
        cat('Normal Vector: \n')
            print(object@normalVector)
        cat('Center Vector: \n')
            print(object@centerVector)
        cat('Form Vector: \n')
            print(object@formVector)
        cat('constant: \n')
        print(object@constant)
        }
)


setMethod(
        f="getLogLikeRatio",
        signature="QuadraticRule",
        definition=function(object){
            return(list(	NormalVector=object@normalVector,
                            VecterVector=object@centerVector,
                            FormVector=object@formVector,
                            Constant=object@constant))
        }
)
.tune.QDA<-function(x,y,procedure='FDRThresh',ql=10^(-1:-7),qq=10^(-1:-2),prior=prior,...)
{
	y=ordered(y)
	call<-match.call()
	ranges<-list(ql=ql,qq=qq,prior=prior,procedure=procedure)
	ranges[sapply(ranges,is.null)]<-NULL
	tunecontrol=tune.control(sampling='cross',best.model=FALSE)
	if(length(ranges)<1) ranges=NULL
	modeltmp<-tune('.learnQuadraticRulefortune',train.x=x,train.y=y,ranges=ranges,
                 tunecontrol=tunecontrol,predict.func=predict,...)
	
	besti=length(ranges$ql)*length(ranges$qq)+1-which.min(rev(modeltmp$performances[['error']]))
    return(.learnQuadraticRulefortune(x,y,
                    procedure=procedure,ql=modeltmp$performances$ql[besti],
					qq=modeltmp$performances$qq[besti],prior=prior))  
}	