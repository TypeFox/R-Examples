discrete.bayes=
function (df, prior, y, ...) 
{
    param = as.numeric(names(prior))
    lk=function(j)
       prod(df(y,param[j],...))
    likelihood=sapply(1:length(param),lk)
    pred = sum(prior * likelihood)
    prob = prior * likelihood/pred
    obj = list(prob = prob, pred = pred)
    class(obj) <- "bayes"
    obj
}
