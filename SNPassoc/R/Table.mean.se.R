`Table.mean.se` <-
function(var, dep, subset = !is.na(var))
{
    var <- as.factor(var)
    n <- ifelse(is.na(tapply(dep[subset],var[subset],FUN=length)),0,tapply(dep[subset],var[subset],FUN=length))
    me <- ifelse(is.na(tapply(dep[subset],var[subset],FUN=mean)),0,tapply(dep[subset],var[subset],FUN=mean))
    se <- ifelse(is.na(tapply(dep[subset],var[subset],FUN=function(x){sd(x)/sqrt(length(x))})),0,tapply(dep[subset],var[subset],FUN=function(x){sd(x)/sqrt(length(x))}))
    ta <- cbind(n=n,me=me,se=se)
    rownames(ta) <- names(n)
    list(tp = ta)
}

