kids.g <-
function(n.kids, p.gene){ 
ptmp <- sum((p.gene-c(1,0))*c(3,1))
return(sample(1:3, n.kids, replace=TRUE, prob=Pgene(1:3, ptmp)))
}
