# TODO: Add comment
# 
# Author: robin
###############################################################################



setClass(
        Class="PartitionWithLLR",
        representation=representation(
                LogLikeRatio="list",
                labels="ordered",
                ThreshProc='character',
                ql='numeric',
                qq='numeric'                
        )
)

setMethod(
        f="predict",
        signature="PartitionWithLLR",
        definition=function(object,newdata){
            if (!(is(newdata, "matrix")||is(newdata, "data.frame")))  
                stop("newdata must be a matrix or a data frame")
            p=ncol(newdata);
            n=nrow(newdata);
            y=object@labels;
            coherence=array(0,c(n,nlevels(y),nlevels(y)),
                                dimnames=list(1:n,levels(y),levels(y)))	
            for (k in levels(y))
            {
                for (l in levels(y))
                {
                    if (l!=k)
                    {
                        LLR=getBinaryRule(object,k,l)
                        if(LLR@prior)
                            p=LLR@proportions
                        else
                            p=1/2
                        priorthresh=log((1-p)/p)
                        coherence[,k,l]=(apply(newdata,1,.EvaluateLogLikeRatio,LLR)>priorthresh)
                    }
                }
            }
            return(apply(coherence,1,.coherence2group))
        }
)


setMethod(
        f="show",
        signature="PartitionWithLLR",
        definition=function(object){
            labels=object@labels;
            cat('Partition of input space into ', nlevels(labels),'areas \n')
            if (class(object@LogLikeRatio[[1]][[2]])=='LinearRule')
            {
                cat('Obtained with Linear separations\n')
                cat('used thresolding procedure: ',object@ThreshProc,'\n')
                if (length(object@ql)>1)
                     cat('with 10 fold cross validation on parameters ql with grid: ',object@ql,'\n')
                else cat('with parameter ql= ',object@ql,'\n')
            }
            if (class(object@LogLikeRatio[[1]][[2]])=='QuadraticRule')
            {
                cat('Obtained with Quadratic separations\n')
                cat('used thresolding procedure: ',object@ThreshProc,'\n')
                if (length(object@ql)>1)
                    cat('with 10 fold cross validation on parameters ql with grid: ',object@ql,'\n')
                else cat('with parameter ql= ',object@ql,'\n')
                if (length(object@qq)>1)
                    cat('with 10 fold cross validation on parameters qq with grid: ',object@qq,'\n')
                else cat('with parameter qq= ',object@qq,'\n')
            }
            

        }
)

setMethod(
        f="plotClassifRule",
        signature=c("PartitionWithLLR"),
        definition=function(object,...){
            labels=object@labels;
            Mydata=data.frame()
            for (k in labels)
            {
                for (j in labels)
                {
                   if (k<j) 
                   {
                       Rule=object@LogLikeRatio[[k]][[j]]
                   
                       if (class(Rule)=='LinearRule')
                       {
    	                   Index=Rule@normalIndex
                           NormalVector=Rule@normalVector[Index]
                           Center=Rule@centerVector[Index]
    	                   CrossLabels=factor(array(paste(k,' Vs ',j,sep=""),length(Index)))
                           Mydata=rbind(Mydata,
                                   data.frame(NormalVector=NormalVector,Center=Center,Index=factor(Index),
                                   CrossLabels= CrossLabels
                                            )
                                    )
                       }
                       if  (class(Rule)=='QuadraticRule')
                       {
                           Index=union(Rule@normalIndex,Rule@formIndex)
                           NormalVector=Rule@normalVector[Index]
                           Center=Rule@centerVector[Index]
                           FormVector=Rule@formVector[Index]
                           CrossLabels=factor(array(paste(k,' Vs ',j,sep=""),length(Index)))
                           Mydata=rbind(Mydata,
                                  data.frame(NormalVector=NormalVector,FormVector=FormVector,
                                           Center=Center,Index=factor(Index),
                                           CrossLabels=CrossLabels
                                            )
                                        )                   
                        }
                    }               
                }
            }
           
            mykey=list(space="right",
                    text=list(levels(Mydata$CrossLabels)),
                    pch=22,
                    title="labels",
                    cex=1.5,
                    columns=1)
            return(xyplot(NormalVector+Center~Index,data=Mydata,
                                                groups=CrossLabels,
                                                #panel=my.Partition.panel(),
                                                auto.key=list(points=FALSE,rectangles=TRUE,space="right",title="labels",columns=1),
                                                xlab='index in feature space',
                                                ylab="",
                                                pch=22,
                                                cex=1.5,
                                                col='black',
                                                aspect='xy',...))
        }
)

learnPartitionWithLLR<-function(x,y,type='linear',procedure='FDRThresh',
        ql=NULL,qq=NULL,BinaryLearningProcedure=NULL,prior=FALSE)
{# type can be linear or quadratic
 # procedure can be FDRThresh, Fisher or 'UnivThresh' in quadratic mode
    
############ Default ql and qr
    if (is.null(ql))
    {
        if (procedure=='UnivThresh'){ql=1:5/4}
        if (!(procedure%in%c('noThresh','UnivThresh'))){ql=10^(-1:-7)}
    }
    if (is.null(qq))
    {
        if (procedure=='UnivThresh'){qq=1:2/2}
        if (!(procedure%in%c('noThresh','UnivThresh'))){qq=10^(-1:-7)}    
    }

############ Setting BinaryLearningProcedure
    if (is.null(BinaryLearningProcedure))
    {
        if (type=='linear')
            {
                if (length(ql)>1){ BinaryLearningProcedure=.tune.LDA;}
                else { BinaryLearningProcedure=.learnLinearRulefortune;}
            } 
        if (type=='quadratic')
        {
            if ((length(ql)>1)|(length(qq)>1)){ BinaryLearningProcedure=.tune.QDA;}
            else { BinaryLearningProcedure=.learnQuadraticRulefortune;}
        }
        if ((type!='linear')&(type!='quadratic'))
        {
            stop('The type of procedure has to be linear or quadratic')
        }
    }
######### Learning

    y=ordered(y)
    f=list()
    if (type=='quadratic')
    {
        for (k in levels(y))
            {
                    f[[k]]=list()
                    for (l in levels(y))
                    {
                            if (l!=k)
                            {
                                    x.train=x[(y==k)|(y==l),]
                                    ytmp=y[(y==k)|(y==l)]
                                    y.train=factor(array(c(1,0),length(ytmp)))
                                    y.train[ytmp==l]=1
                                    y.train[ytmp==k]=0
                                    f[[k]][[l]]=.tune.QDA(x.train,y.train,procedure,ql,qq,prior=prior)
                            }
                
                    }
            }
    }
    if (type=='linear')
    {
        for (k in levels(y))
            {
                    f[[k]]=list()
                    for (l in levels(y))
                    {
                            if (l!=k)
                            {
                                    x.train=x[(y==k)|(y==l),]
                                    ytmp=y[(y==k)|(y==l)]
                                    y.train=factor(array(c(1,0),length(ytmp)))
                                    y.train[ytmp==l]=1
                                    y.train[ytmp==k]=0
                                    f[[k]][[l]]=.tune.LDA(x.train,y.train,procedure,ql,prior=prior)
                            }
                
                    }
            }
    }
    return(new(Class="PartitionWithLLR",LogLikeRatio=f,labels=y,ThreshProc=procedure,ql=ql,qq=qq)) 
}

setMethod(f="getBinaryRule",signature='PartitionWithLLR',
        definition=function(object,k,l){
            #i and j should be ordered factors
            if (k!=l) return(object@LogLikeRatio[[k]][[l]])
            #if (k>l) return(.minus(object@LogLikeRatio[[l]][[k]]))
            if (k==l) return(NULL)
        }
)



.coherence2group<-function(coherence)
{    
    X=rowMeans(coherence)  
    indexes=which(X==max(X))
    if (length(indexes)==1){return(dimnames(coherence)[[2]][indexes])}
    if (length(indexes)==2)
    {
        if (coherence[indexes[1],indexes[2]]==1)
            return(dimnames(coherence)[[2]][indexes[1]])
        else
            return(dimnames(coherence)[[2]][indexes[2]])
    }
    if (length(indexes)>2)
        return(dimnames(coherence)[[2]][sample(indexes,1)])
}