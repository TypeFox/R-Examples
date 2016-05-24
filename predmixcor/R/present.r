present <- function(prob.pred,x)
{
    category=character(10)
    meanprob=prop1=no=rep(0,10)
    for( i in 1:10 )
    {
        category[i]=paste((i-1)/10,"-",i/10,sep="")
        if(i == 1)
            label <- (prob.pred <=i/10 )
        else
        	label=(prob.pred<=i/10) * (prob.pred>(i-1)/10)
        no[i]=sum(label)
        meanprob[i]=sum(prob.pred*label) / no[i]
        prop1[i]=sum(x*label)/no[i]
    }
    data.frame(Category=category,no,Pred=meanprob,Actual=prop1)    
}
