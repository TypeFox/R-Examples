`jstats` <-
function(d)
  {
    s1 = boxplot.stats(d, do.conf=FALSE, do.out=FALSE);
    m1 = mean(d, na.rm = TRUE)
    m2 = sqrt(var(d, na.rm = TRUE))
    
    return(list(mean=m1, std=m2, med=s1$stats[3], qdist=s1$stats[4]-s1$stats[2] ,
                bstats=s1$stats, mstats=c(m1,m2), N=length(d)))
    
    
  }

