`roc.area` <-
function (obs, pred) 
{
   id <- is.finite(obs)&is.finite(pred)
   obs<- obs[id]
   pred <- pred[id]
        n1 <- sum(obs)
        n<- length(obs)
   A.tilda <-        (mean(rank(pred)[obs == 1]) - (n1 + 1)/2)/(n - n1)
    stats <- wilcox.test(pred[obs == 1], pred[obs == 0], alternative = "great")
    return(list( A = A.tilda, n.total = n, n.events = n1, n.noevents = sum(obs == 0), 
        p.value = stats$p.value))
}
