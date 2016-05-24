#setGeneric("predict", function(object, ...) standardGeneric("predict"))
setGeneric(
        name="getLogLikeRatio",
        def=function(object){standardGeneric('getLogLikeRatio')}
)
setGeneric(
        name=".EvaluateLogLikeRatio",
        def=function(x,object){standardGeneric('.EvaluateLogLikeRatio')}
)
setGeneric(
        name=".minus",
        def=function(object){standardGeneric('.minus')}
)
setGeneric(
        name="plotClassifRule",
        def=function(object,...){standardGeneric('plotClassifRule')}
)
setGeneric(
        name="getBinaryRule",
        def=function(object,k,l){standardGeneric('getBinaryRule')}
)

learnBinaryRule<-function(x,y,type='linear',procedure='FDRThresh',covariance='diagonal',ql=NULL,qq=NULL,prior=FALSE)
{
    
    if (is.null(ql))
    {
        if (procedure=='UnivThresh'){ ql=1:5/4}
        if (!(procedure%in%c('noThresh','UnivThresh'))){ql=10^(-1:-7)}
    }
    if (is.null(qq))
    {
        if (procedure=='UnivThresh'){qq=1:2/2}
        if (!(procedure%in%c('noThresh','UnivThresh'))){qq=10^(-1:-7)}    
    }
    
    
    if (type=='linear')
    {
        if (length(ql)>1){ BinaryLearningProcedure=.tune.LDA;}
        else { BinaryLearningProcedure=.learnLinearRulefortune;}
        return(BinaryLearningProcedure(x,y,procedure=procedure,covariance=covariance,ql=ql,prior=prior))
    } 
    if (type=='quadratic')
    {
        if ((length(ql)>1)|(length(qq)>1)){ BinaryLearningProcedure=.tune.QDA;}
        else { BinaryLearningProcedure=.learnQuadraticRulefortune;}
        return(BinaryLearningProcedure(x,y,procedure=procedure,covariance=covariance,ql=ql,qq=qq,prior=prior))
    }
    if ((type!='linear')&(type!='quadratic'))
    {
        stop('The type of procedure has to be linear or quadratic')
    }
    
    
}
    
