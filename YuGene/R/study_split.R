# Author : F.Rohart
# created 07-07-2014
# last modified 07-07-2014
#
# split the data according to the study factor, results are in a list


#=============================================================
study_split = function(X, study)
{
    
    X = as.matrix(X)
    
    M=length(levels(study))
    P=ncol(X)
    
    #---------------------- split data
    X.list.study = split(X,study)
    study.name = split(rownames(X),study)
    
    for(m in 1:M)
    {
        
        X.list.study[[m]]=matrix(X.list.study[[m]], ncol=P)
        colnames(X.list.study[[m]]) = colnames(X)
        rownames(X.list.study[[m]]) = study.name[[m]]

    }
    
    
    result = list(X.list.study = X.list.study)
    
    return(invisible(result))
    
}